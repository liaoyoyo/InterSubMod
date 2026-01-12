# Documentation Guide

This guide outlines the structure, naming conventions, and organization rules for the `InterSubMod` documentation.

## Directory Structure

The documentation is organized by **role** and **purpose**:

* **`architecture/`**: High-level system design, module interactions, and data flow. For architects and lead developers.
* **`development/`**: Implementation details, setup guides, roadmaps, and developer-specific instructions.
* **`reports/`**: One-off reports including progress updates, test results, design reviews, and verification analysis.
* **`issues/`**: Issue tracking, bug analysis, and known issues.
* **`manual/`**: User manuals and operational guides.
* **`requirements/`**: (Optional) System and user requirements.
* **`api/`**: (Optional) API documentation.
* **`data/`**: (Optional) Data formats, dictionaries, and source descriptions.
* **`operations/`**: (Optional) Deployment guides, runbooks, and monitoring.
* **`decisions/`**: (Optional) Architecture Decision Records (ADR).
* **`archive/`**: Deprecated or superseded documents.

## Naming Conventions

* **Files**: Use `snake_case.md` (e.g., `system_overview.md`, `phase1_report.md`).
* **Directories**: Use `snake_case` (e.g., `design_reviews`, `known_issues`).
* **Dates**: When including dates in filenames, use `YYYY_MM_DD` format (e.g., `2025_11_30_verification_analysis.md`).

## Organization Rules

1. **Architecture vs. Development**:
    * Place long-term, high-level design documents in `architecture/`.
    * Place implementation details, "how-to" guides, and specific plans in `development/`.

2. **Reports**:
    * Group reports by type (`progress`, `tests`, `verification`, `design_reviews`).
    * Use `misc/` for reports that don't fit into specific categories.

3. **Issues**:
    * Maintain `known_issues.md` and `resolved_fixes.md` as living documents.
    * Place detailed analysis of specific bugs in `issues/analysis/`.

4. **Archive**:
    * Do not delete old documents. Move them to `archive/` if they are no longer relevant but historically interesting.

## Verification Documents

Verification reports should be placed in `reports/verification/` and named with the date and topic, e.g., `YYYY_MM_DD_topic_verification.md`.
