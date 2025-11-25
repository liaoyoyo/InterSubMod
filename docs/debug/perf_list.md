# Performance and Logic Issues List

## 1. Performance Issues

### 1.1 Repeated BAM File Opening and Index Loading

- **Location**: `RegionProcessor::process_single_region` (src/core/RegionProcessor.cpp)
- **Issue**: `BamReader` is instantiated inside `process_single_region`, which is called for every SNV region. The constructor calls `sam_index_load`, which reads the `.bai` file.
- **Impact**: High I/O overhead and latency, especially for large numbers of regions. The index is re-loaded thousands of times.
- **Improvement**:
  - Refactor `RegionProcessor` to manage `BamReader` instances at the thread level (e.g., in `process_all_regions` loop).
  - Pass an open `BamReader` (or a pool) to `process_single_region`.

### 1.2 Inefficient Data Structure in MatrixBuilder

- **Location**: `MatrixBuilder` (include/core/MatrixBuilder.hpp, src/core/MatrixBuilder.cpp)
- **Issue**: Uses `std::map<int, std::map<int32_t, float>> read_methyl_map_` to store temporary methylation calls.
- **Impact**: High memory overhead and slow insertion/traversal due to many small allocations and tree balancing.
- **Improvement**:
  - Use `std::vector<std::vector<std::pair<int32_t, float>>>` (read_id -> list of calls).
  - Since we iterate sequentially to build the final matrix, a vector is sufficient and much faster.

### 1.3 Redundant CIGAR Traversal

- **Location**: `ReadParser::determine_alt_support` and `MethylationParser::parse_read`
- **Issue**: Both methods traverse the CIGAR string to map sequence positions to reference positions. `ReadParser` does it to find the SNV, `MethylationParser` does it to find CpGs.
- **Impact**: Double work for every read.
- **Improvement**:
  - Compute the `seq_to_ref` mapping once (e.g., in `ReadParser` or a shared utility) and pass it to both.
  - Or, since `ReadParser` only needs one position, maybe keep it simple, but `MethylationParser` definitely needs the full map.

### 1.4 Vector Allocation in Loop

- **Location**: `MethylationParser::build_seq_to_ref_map`
- **Issue**: Allocates `std::vector<int32_t> seq_to_ref(b->core.l_qseq, -1)` for every read.
- **Impact**: Frequent heap allocations.
- **Improvement**:
  - Use a thread-local vector and `clear()`/`resize()` it to reuse memory.

## 2. Logic and Correctness Issues

### 2.1 Reverse Strand Methylation Ignored

- **Location**: `MethylationParser::parse_read` (src/core/MethylationParser.cpp)
- **Issue**: The code iterates the read sequence and checks `if (base == 'C')`. For reads mapped to the reverse strand, the BAM SEQ is reverse complemented (so original 'C' becomes 'G'). The current logic skips these 'G's, effectively ignoring methylation on the reverse strand.
- **Impact**: Loss of ~50% of methylation data (assuming 50/50 strand balance).
- **Improvement**:
  - Check `bam_is_rev(b)`.
  - If reverse strand, look for 'G' (which corresponds to 'C' on the original strand) and check if it has modification.
  - OR, understand how `MM` tag behaves for reverse strand. Standard says `MM` refers to SEQ. If SEQ has 'G', `MM` should index 'G's?
  - Need to verify if the upstream caller (Dorado/Guppy) produces `C+m?` or `G-m?` for reverse reads, or if it normalizes everything.
  - **Action**: Add logic to handle reverse strand 'G's if they carry the methylation signal.

### 2.2 Potential Type Mismatch

- **Location**: Various headers.
- **Issue**: Mixing `int32_t` and `uint32_t` for positions.
  - `SomaticSnv::pos` is `uint32_t`.
  - `RegionOfInterest::win_start_pos` is `int32_t`.
  - `ReadInfo::align_start` is `int32_t`.
- **Impact**: Potential overflow or comparison warnings.
- **Improvement**: Standardize on `int32_t` (as HTSlib uses `int32_t` for positions) or `int64_t` to be safe, but `int32_t` is standard for BAM.

## 3. Code Quality and Readability

### 3.1 Missing Comments

- **Location**: General
- **Issue**: Some complex logic (like CIGAR parsing) lacks detailed comments explaining *why* it's done that way.
- **Improvement**: Add English comments explaining the logic, especially for the coordinate transformations.

### 3.2 Error Handling

- **Location**: `BamReader::fetch_reads`
- **Issue**: Returns empty vector on error, but doesn't distinguish between "no reads" and "error".
- **Improvement**: Log a warning or throw an exception on error.
