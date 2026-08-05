<!--
建立時間: 2026-07-23 00:00 +08:00
目標: 逐步記錄 exact-PS strict endpoint read-linkage production 修正的設計決定、偏離、折衷與未決問題
處理範圍: Python/C++ grouping、runner、adapter、tests、HCC1395與7-dataset validation
關聯檔案:
  - InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/pre-decision-audit.md
-->

# Production exact-PS strict read-linkage — Implementation Notes

- status: `validated_strict_region_report_with_completion_boundary_audit`（strict directed topology 仍列 V001）
- spec_id: `production_exact_ps_strict_read_linkage`
- task_type: `C_production_deployment`
- evidence_scope: `7 technical datasets / 6 biological cell lines / chr1–22；154/154 extraction + strict receipts；正式 W=85,621；L1=7/7、production strict topology=0/7、clone/parent/fusion=0/7；HTML READY`

## 🔵 Design Decisions

### D001 — Production scientific region contract

- Status: Accepted（user-specified）
- Tier: L1
- Decision: We will define a scientific region as a connected component of sSNV nodes inside one chromosome × exact nonmissing PS × HP container, where each retained endpoint edge has at least the declared number of distinct canonical molecules with fixed R/A calls at both endpoints.
- Rationale: Different PS blocks have unresolved relative HP orientation; cut-span threshold can include a site without a threshold-supported endpoint edge.
- Revisit if: a signed cross-PS SAME/FLIP orientation model and independent validation are implemented.

<!-- BEGIN USER-SPECIFIED -->
**DO NOT change**: exact PS and HP are hard primary-container boundaries; distance alone cannot define a scientific region.
<!-- END USER-SPECIFIED -->

### D002 — Region and bounded block remain separate units

- Status: Accepted
- Tier: L1
- Decision: We will retain `source_component_id` across every k≤12 bounded block; block counts will never be reported as biological region counts.
- Revisit if: an exact solver removes the k bound without segmentation.

### D003 — Read crossing definition

- Status: Accepted
- Tier: L1
- Decision: We will count a molecule as supporting edge i–j only when both endpoints have fixed R/A calls after canonical alignment deduplication and MAPQ/BQ gates. Coordinate span or X/O/D/S/L at either endpoint is not edge support.
- Revisit if: an explicit missingness/error emission replaces the hard endpoint gate.

### D004 — Singleton component is not a tree region

- Status: Accepted（user-corrected）
- Tier: L1
- Decision: A k=1 graph component is retained only as `ABSTAIN_SINGLETON_UNLINKED` for locus/mass conservation. It must never enter k≤12 partitioning, Steiner enumeration, topology denominators, or clone/subclone interpretation. Only k>1 components with at least one threshold-qualified direct endpoint edge are `PRIMARY_PS_AWARE` tree-eligible regions.
- Rationale: A singleton is a valid graph-theory component but contains no co-read relation and therefore cannot identify mutation order or branching.
- Revisit if: never; a different single-locus analysis may consume it, but that is not regional tree reconstruction.

<!-- BEGIN USER-SPECIFIED -->
**DO NOT change**: isolated single sSNVs are not tree-building regions; preserve them only as an explicit abstention branch.
<!-- END USER-SPECIFIED -->

### D005 — 50 kb is a diagnostic threshold, not a region boundary

- Status: Accepted
- Tier: L1 implementation + L2 descriptive validation
- Decision: No coordinate-distance cutoff is allowed in primary edge creation or connected-component construction. Report `W span>50 kb`, `max adjacent gap>50 kb`, and direct-edge distance separately as QC.
- Rationale: Fixed endpoint calls on the same canonical molecule are direct linkage evidence. A global distance cutoff can discard such evidence; conversely, a long coordinate envelope does not imply one molecule spans the full W.
- HCC1395 evidence: 1,064/11,462 W have total span>50 kb, only 4 have an adjacent gap>50 kb, and 47/76,202 retained direct edges exceed 50 kb. The long edges occur in 22 W; removing them changes 4 W partitions and loses 4 memberships to singleton status even though the W count remains 11,462.
- Revisit if: a held-out error model supports a calibrated distance-dependent edge likelihood; it must enter as uncertainty/weight rather than an undocumented hard split.

### D006 — Report funnels preserve units and denominators

- Status: Accepted
- Tier: L1
- Decision: The all-dataset report uses three separate conserved tracks: physical loci, PS×HP memberships, and graph components. It must not draw `S loci → W regions` as one funnel.
- Verification: `build_all7_region_report_data.py` independently checks component partition, membership mass, physical-locus roles, 7×22 scope, and threshold grid before emitting a ready package.

### D007 — HTML data package and graph rows fail closed

- Status: Accepted
- Tier: L1
- Decision: The report builder must independently validate every primary component, membership and retained edge, including exact dataset/chromosome/PS/HP scope, component-ID uniqueness, `membership count=k`, component geometry, edge-coordinate equality, `support_total=RR+RA+AR+AA`, within-component endpoints and graph connectivity. The artifact builder separately verifies the canonical 7×22 scope plus the signed report-data receipt before any chart/table staging.
- Verification: `32 passed` across strict graph, strict builder, all-7 report and release-hardening tests; C++ strict graph suite `PASS`. SQLite sources are executed queries and are atomically published; failed staging is preserved under an audit archive rather than deleted.

### D008 — Portable HTML is checksum-pinned and responsive

- Status: Accepted
- Tier: L1
- Decision: The delivery wrapper is pinned to the audited `datascience-mcp-widgets@0.2.8` entrypoints, embedded reader bundle parts, bundled validator server, browser helper/CLI modules, package metadata, and Chromium headless-shell SHA-256 values; caller-supplied plugin/browser overrides are rejected. The artifact JSON is the final atomic commit marker, and failed DB/query/artifact staging is moved to an audit archive.
- Responsive fix: the official 390 px sandbox exposed a Recharts legend with an intrinsic 570 px width. The wrapper constrains `.chart-legend` to its card and permits wrapping; it also retains the required top-bar width normalization. No chart data or scales are modified.
- Verification: synthetic full-contract preflight passes canonical artifact validation, 11 charts, 7 tables, 37 blocks, 22 light/dark inline SVG variants, official 1440/390 browser verification, source-dialog interaction, JavaScript on/off fallback checks, and custom Playwright QA with zero horizontal overflow. The fixture duplicates HCC1395 only for structural QA and must never be cited as cross-dataset science.

### D009 — Missing PS tokens are normalized before container construction

- Status: Accepted and regression-tested
- Tier: L1
- Decision: Python and C++ both trim the PS field, compare it case-insensitively against `""`, `.`, `NA`, `N/A`, `NAN`, `NONE`, and `NULL`, and exclude matching rows before any exact-PS container key is created. Nonmissing tokens are retained after trimming.
- Rationale: Checking only Python truthiness allowed standard text missing tokens such as `.` or `NULL` to become false exact-PS containers. Current extraction serializes missing PS as empty, and the report reader already rejected the broader missing-token set, so no completed formal dataset was shown to contain this leakage; nevertheless the production contract must be correct at its first enforcement point.
- Verification: Python missing-token regression and strict graph suite `15 passed`; C++ strict graph suite `PASS`. Post-fix builder SHA-256 is `912721f934bae6a58ccfe66b872706f035b658d1fc1b53db4025a998916e4b4d`.
- Compatibility probe: Rebuilding pre-hotfix HCC1937 chr21 with the post-fix builder produced byte-identical decompressed component (811 rows), membership (1,313), edge (285), and container (309) TSVs plus identical counts/threshold summaries/checks. The receipt bytes differ only because time and builder/output identities are expected to change.
- Provenance: The exact pre-hotfix builder bytes used by earlier strict receipts are preserved at `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/source_snapshots/7260a7631b30cbb5e4878159583b8a5b27a153de07e8d001303417dd2f29aedd_build_strict_ps_hp_regions.py`; its recomputed SHA-256 is exactly the filename prefix.

### D010 — Mixed builder provenance is accepted only as data-specific no-trigger equivalence

- Status: Accepted with explicit disclosure
- Tier: L1
- Decision: The formal all-7 report may use the 132 pre-hotfix chromosome artifacts for HCC1395, HCC1395_DORADO, COLO829, H1437, H2009 and HCC1937 together with the 22 post-hotfix HCC1954 artifacts because the patched missing-PS branch is provably not exercised by the earlier data.
- Verification: All 156,316 rows in the first six datasets' production container TSVs were scanned against the normalized missing-token set and yielded zero matches. A post-hotfix rebuild of HCC1937 chr21 was byte-identical after decompression for component, membership, edge and container TSVs. The all-7 report independently validates 154/154 extraction receipts, 154/154 strict receipts and every primary graph row before publication.
- Disclosure boundary: This supports numerical equivalence for the current data, not same-builder provenance. Any release that requires all 154 chromosome artifacts to share one builder SHA must create new roots and rerun the earlier 132 chromosomes; existing artifacts must not be overwritten.

### D011 — Completion layers are reported separately and fail closed

- Status: Accepted and independently audited
- Tier: L1 status audit
- Decision: The report must display four separate completion layers: L1 strict read-linkage regions, earlier/legacy candidate-tree references, current production strict directed topology, and cellular clone/parent/fusion inference. A PASS at one layer must never be promoted to a later layer.
- Latest result: L1 is complete for 7/7 datasets with W=85,621. No full-scope eligible v4 strict topology receipt exists, so current production strict directed topology is 0/7. Cellular clone count, exact cellular parent-child edges, and cross-HP/cross-technical fused trees are 0/7 validated.
- HCC1395 boundary: the earlier exact-PS topology observation is a technical pilot whose upstream receipt is `PARTIAL` and `validation_evidence_eligible=false`; it is not bound to the current 11,462-W production strict root.
- Legacy boundary: the older all-7 candidate-tree census contains mixed-PS/50-kb-era regions and remains historical reference only.
- Verification: machine-readable audit, signed receipt and TSV are in `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/topology_status_audit/`.

## 🟠 Deviations

### V001 — Existing exact-PS pilot is not the production runner

- Status: Open
- Tier: L1
- Deviation: The v4 strict production entrypoint now consumes the exact-PS strict contract, but the latest executed v4 receipt is only an HCC1395 chr1 partition smoke with `topology=null`. The older HCC1395 exact-PS topology pilot validates a side route and is not eligible production evidence.
- Required action: execute the v4 strict topology stage over the preregistered full dataset×chr1–22 scope, publish per-sample topology receipts, then separately validate ranking and L4 cellular fusion claims.

### V002 — First strict artifact version mislabeled singleton components as primary

- Status: Corrected; v2 recomputation and all-7 validation complete
- Tier: L1
- Deviation: `hcc1395_strict_regions_v1` correctly found graph components but marked k=1 rows `PRIMARY_PS_AWARE`.
- Correction: schema 1.1.0 labels k=1 `ABSTAIN_SINGLETON_UNLINKED`; v1 is superseded and must not be used as partition input.

## 🟡 Trade-offs

### T001 — Fixed pairwise threshold versus distributed partial evidence

- Status: Accepted for primary; sensitivity required
- Tier: L2
- Trade-off: `threshold=3` is conservative and auditable but can split evidence distributed across several low-count endpoint pairs. We will retain edge weights and report threshold 1/2/3/5 sensitivity rather than hide the distinction.
- Provenance boundary: threshold=3 is independently recomputed from component, membership and edge rows; threshold 1/2/5 values are reaggregated from the current 22 chromosome receipts and reconciled field-by-field to the dataset summary, but are not a second read of every non-primary-threshold TSV row.
- Revisit if: held-out molecule likelihood demonstrates a calibrated edge-confidence alternative.

## 🔴 Open Questions

### Q001 — Full seven-dataset input readiness

- Status: Resolved for the strict-region report
- Tier: L1
- Question: Do all seven canonical datasets have source-bound exact HP/PS molecule rows required for a no-BAM-rescan production rerun?
- Resolution: 7/7 datasets × chr1–22 are complete. The report independently verifies 154/154 extraction receipts and 154/154 strict receipts, reconciles all primary graph rows, and publishes only after artifact validation plus JavaScript on/off desktop/mobile browser QA. This resolves input readiness for the current report; it does not close V001 or establish a unique topology.

### Q002 — Production runtime and memory

- Status: Open
- Tier: L5
- Question: Can pair support be accumulated per exact container without materializing pathological O(k²) pairs for high-k reads, while preserving exact molecule counts?

## 📚 Lore

- HP numeric labels are local to a phase set unless an external orientation bridge proves SAME/FLIP.
- A read spanning genomic coordinates is not an allele-informative bridge.
- Current cut-span support is O(calls+sites) but has different scientific semantics from endpoint co-call support.
- Same-read marginal VAF must not be added again as an independent likelihood term.
- Raising the endpoint-support threshold guarantees that retained edges do not increase and the graph partition only refines. It does **not** guarantee that the count of k>1 W decreases, because one larger W can split into multiple k>1 W.
- HCC1395 v1 preliminary threshold=3 census (component geometry is valid, role label superseded): 39,846 all components = 28,384 unlinked singletons + 11,462 read-linked multisite regions; S=79,687, 62,651 exact-container node memberships, 76,202 qualifying endpoint edges.
- HCC1395 v2 singleton leakage audit: strict memberships conserve exactly as 62,651 = 28,384 abstain + 34,267 tree-eligible; partition has 11,462 units with min(k)=2, k=1=0, singleton component-ID overlap=0, and eligible/partition component-ID sets are identical. At physical-locus level, 22,689/36,384 active loci enter at least one linked membership and 13,695 are singleton-only; 5,717 occur in both roles across different PS/HP memberships, so exclusion is membership-scoped rather than global-locus deletion.
- Independent child-block audit: partition produces 11,712 bounded blocks, including 12 k=1 child blocks derived exclusively from k=13–153 source components; none intersect source singleton IDs, and all 12 have zero retained weight/pattern support. They are computational fragments rather than regions and are explicitly skipped by the adapter before mutation-state tree construction.
- All-7 HTML report implementation uses the canonical Data Analytics artifact contract, native SVG charts, semantic tables, and a repo-local strict delivery wrapper. The wrapper imports checksum-pinned official build/extract/verify modules, adds deterministic top-bar and narrow-screen legend containment without modifying plugin cache, requires light/dark SVG extraction for every chart, and runs official 1440/390 browser verification before atomic publication.
- Formal all-7 result: dataset-level totals are S=469,849, active loci=342,374, all components=255,752, k=1 abstain=170,131, W=85,621, W memberships=443,349 and retained direct edges=1,197,530. Removing all >50 kb direct edges changes 1,172 W partitions, fully loses 228 W, returns 961 memberships to k=1 and changes the summed W count to 85,872 because some W split into multiple k≥2 components. The READY bundle is `InterSubMod/research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/READY.json`.
- Latest completion census: strict read-linkage L1=7/7; current production strict directed topology=0/7; cellular clone count=0/7; exact cellular parent-child=0/7; cross-HP/cross-technical fusion tree=0/7. Endpoint edges remain undirected linkage evidence.
