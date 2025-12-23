# InterSubMod Documentation

Welcome to the InterSubMod documentation hub. This directory contains all technical, architectural, and operational documents for the project.

## Directory Structure

We organize documentation by **purpose** and **lifecycle**:

### Core
*   **[`architecture/`](architecture/)**: High-level system design, module interactions, and trade-off analysis.

### Development Process (`development/`)
*   **[`plans/`](development/plans/)**: Concrete implementation plans, roadmaps, and strategies.
*   **[`analysis/`](development/analysis/)**: Design evaluations, algorithm research, and feasibility studies.
*   **[`guides/`](development/guides/)**: Implementation guides, "how-to" docs, and checklists.
*   **[`notes/`](development/notes/)**: Implementation notes, specific algorithm details, and working drafts.

### Reporting (`reports/`)
*   **[`investigation/`](reports/investigation/)**: Deep dives into specific bugs, anomalies, or logic inquiries.
*   **[`verification/`](reports/verification/)**: Validation reports ensuring features work as intended.
*   **[`reviews/`](reports/reviews/)**: User feedback, code reviews, and design critiques.
*   **[`progress/`](reports/progress/)**: General progress updates and organizational reports.

### Other
*   **[`issues/`](issues/)**: Bug analyis, known issues, and fix logs.
*   **[`manual/`](manual/)**: User manuals and operational guides.
*   **[`api/`](api/)**: API specifications.
*   **[`data/`](data/)**: Data format descriptions and dictionaries.
*   **[`archive/`](archive/)**: Superseded or obsolete documents.
*   **[`garbage/`](garbage/)**: Temporary holding area for files pending deletion.

## Naming Conventions

To ensure chronologically sortable and searchable filenames, we enforce the following rules:

1.  **Date Prefix**: `YYYYMMDD_filename.md`
    *   *Required* for: Reports, roadmaps, implementation plans, and time-sensitive analysis.
    *   *Optional* for: Long-living architectural diagrams or strictly evergreen manuals (though versioning is preferred).
2.  **Snake Case**: Use `lower_snake_case` for all filenames and directories.
    *   Example: `20251223_significance_analysis_strategy.md`
3.  **Format**: Markdown (`.md`) is the standard format.

## Guidelines

*   **New Documents**: Place new design docs in `development/plans` or `development/analysis` initially.
*   **Updates**: Do not overwrite old reports. Create a new report with a new date if the findings change significantly.
*   **Deletion**: Do not delete files. Move them to `archive/` or `garbage/` if they are no longer needed.

For more details, see [Documentation Guide](documentation_guide.md).
