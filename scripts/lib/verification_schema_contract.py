#!/usr/bin/env python3
"""Typed, fail-closed helpers for Verification/RegionStratification consumers.

This module is intentionally independent of InterSubMod data paths.  Consumers
provide already-loaded pandas objects; the helpers validate schema identity and
return an explicit view instead of guessing biological or historical meaning.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Sequence, Tuple
import warnings

import pandas as pd


VERIFICATION_SCHEMA_VERSION = 2
REGION_STRATIFICATION_SCHEMA_VERSION = 1

CURRENT_CLASSES_V2: Tuple[str, ...] = (
    "Strong_Bidirectional",
    "ClusterFirstOnly",
    "LOH-Structure",
    "MultiGroupNoLabel",
    "LabelShift",
    "PermanovaLocation",
    "StructureNoLabel",
    "DispersionStructure",
    "Noise_Uniform",
    "Noise_Chaotic",
    "Noise_Uncorrelated",
)
UNKNOWN_CURRENT_CLASS = "UnknownCurrentClass"

V1_CURRENT_CLASSES: Tuple[str, ...] = (
    "Strong",
    "LOH-Structure",
    "MultiGroupNoLabel",
    "LabelShift",
    "PermanovaLocation",
    "StructureNoLabel",
    "DispersionStructure",
    "Noise_Uniform",
    "Noise_Chaotic",
    "Noise_Uncorrelated",
)

CURRENT_TO_EVIDENCE_PATH = dict(
    zip(
        CURRENT_CLASSES_V2,
        (
            "BIDIRECTIONAL",
            "CLUSTER_FIRST_ONLY",
            "LOH_STRUCTURE",
            "WITHIN_HP_MULTIGROUP",
            "LABEL_SHIFT",
            "PERMANOVA_LOCATION",
            "HP_AUC_STRUCTURE_NO_LABEL",
            "DISPERSION_STRUCTURE",
            "NOISE_UNIFORM",
            "NOISE_CHAOTIC",
            "NOISE_UNCORRELATED",
        ),
    )
)

LEGACY_CLASSES: Tuple[str, ...] = ("Strong", "Subclone", "Weak", "Noise")
UNKNOWN_LEGACY_CLASS = "UnknownLegacyClass"

EVIDENCE_PATHS: Tuple[str, ...] = (
    "BIDIRECTIONAL",
    "CLUSTER_FIRST_ONLY",
    "LOH_STRUCTURE",
    "WITHIN_HP_MULTIGROUP",
    "LABEL_SHIFT",
    "PERMANOVA_LOCATION",
    "HP_AUC_STRUCTURE_NO_LABEL",
    "DISPERSION_STRUCTURE",
    "NOISE_UNIFORM",
    "NOISE_CHAOTIC",
    "NOISE_UNCORRELATED",
)
EVIDENCE_DERIVATIONS: Tuple[str, ...] = ("LIVE", "LEGACY_CLASS")

LOH_LEGACY_CLASSES: Tuple[str, ...] = (
    "None",
    "LOH_Noise",
    "LOH_Weak",
    "LOH_Strong",
    "LOH_Subclone",
)

REGION_STATUS_VALUES: Tuple[str, ...] = (
    "VALID",
    "INSUFFICIENT_REGIONS",
    "NOT_APPLICABLE_TUMOR_ONLY",
    "FAILED",
)
REGION_ASSIGNMENTS = {
    0: ("BaselineLowAsm", "BASELINE_LOW_ASM"),
    1: ("HighHpAsm", "HIGH_HP_ASM"),
    2: ("LohFlagged", "LOH_FLAGGED"),
    3: ("HighSampleAsm", "HIGH_SAMPLE_ASM"),
}
REGION_UNASSIGNED_REASON_BY_STATUS = {
    "VALID": "INELIGIBLE_REGION",
    "INSUFFICIENT_REGIONS": "INSUFFICIENT_REGIONS",
    "NOT_APPLICABLE_TUMOR_ONLY": "NOT_APPLICABLE_TUMOR_ONLY",
    "FAILED": "FAILED",
}

VERIFICATION_PROVENANCE_COLUMNS: Tuple[str, ...] = (
    "VerificationSchemaVersion",
    "VerificationClass",
    "VerificationClass_V1_Deprecated",
    "VerificationClass_Legacy",
    "LabelFirstSupport",
    "ClusterFirstSupport",
    "WithinHPSupport",
    "DispersionWarning",
    "EvidencePath",
    "EvidenceDerivation",
    "LOH_Subtype_LegacyVC",
    "LOH_Subtype",
)


class SchemaContractError(ValueError):
    """Raised when a consumer cannot safely determine the requested schema."""


@dataclass
class ClassView:
    values: pd.Series
    field: str
    schema_status: str
    categories: Tuple[str, ...]
    unknown_counts: Dict[str, int]
    warning_messages: Tuple[str, ...] = ()

    def metadata(self) -> Dict[str, object]:
        return {
            "selection_field": self.field,
            "schema_status": self.schema_status,
            "categories": list(self.categories),
            "unknown_counts": dict(self.unknown_counts),
            "warnings": list(self.warning_messages),
        }


@dataclass
class RegionStratumView:
    frame: pd.DataFrame
    status: Dict[str, str]
    assignment_count: int
    n_occupied_region_strata: int


def _require_columns(df: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise SchemaContractError(f"{context}: missing required columns: {', '.join(missing)}")


def _clean_string_series(series: pd.Series) -> pd.Series:
    def clean(value: object) -> object:
        if pd.isna(value):
            return pd.NA
        # Enum identity is exact.  In particular, do not silently accept
        # surrounding whitespace or normalize case.
        text = str(value)
        return pd.NA if text == "" else text

    return series.map(clean)


def _unknown_counts(values: pd.Series, allowed: Sequence[str]) -> Dict[str, int]:
    allowed_set = set(allowed)
    counts = values.dropna().astype(str).value_counts()
    return {str(key): int(value) for key, value in sorted(counts.items()) if str(key) not in allowed_set}


def _require_single_integer_version(df: pd.DataFrame, field: str, expected: int, context: str) -> None:
    _require_columns(df, [field], context)
    if any(isinstance(value, bool) for value in df[field].tolist()):
        raise SchemaContractError(f"{context}: {field} must be an integer, not boolean")
    numeric = pd.to_numeric(df[field], errors="coerce")
    if numeric.isna().any():
        raise SchemaContractError(f"{context}: {field} contains missing or non-numeric values")
    versions = sorted(set(float(value) for value in numeric.tolist()))
    if versions != [float(expected)]:
        raise SchemaContractError(f"{context}: expected {field}={expected}, observed {versions}")


def select_current_view(df: pd.DataFrame, allow_unversioned_raw: bool = False) -> ClassView:
    """Select canonical VerificationClass v2 or an explicitly allowed raw view."""

    _require_columns(df, ["VerificationClass"], "current verification view")
    raw = _clean_string_series(df["VerificationClass"])

    if "VerificationSchemaVersion" not in df.columns:
        if not allow_unversioned_raw:
            raise SchemaContractError(
                "current verification view: VerificationSchemaVersion is missing; "
                "explicit unversioned raw mode is required"
            )
        categories = tuple(sorted(set(str(value) for value in raw.dropna().tolist())))
        message = "UNVERSIONED: plotting raw VerificationClass values without v2 interpretation"
        return ClassView(
            values=raw,
            field="VerificationClass",
            schema_status="UNVERSIONED",
            categories=categories,
            unknown_counts={},
            warning_messages=(message,),
        )

    _require_single_integer_version(
        df, "VerificationSchemaVersion", VERIFICATION_SCHEMA_VERSION, "current verification view"
    )
    unknown = _unknown_counts(raw, CURRENT_CLASSES_V2)
    normalized = raw.where(raw.isin(CURRENT_CLASSES_V2), UNKNOWN_CURRENT_CLASS)
    normalized = normalized.where(raw.notna(), UNKNOWN_CURRENT_CLASS)
    if raw.isna().any():
        unknown["<MISSING>"] = int(raw.isna().sum())
    return ClassView(
        values=normalized,
        field="VerificationClass",
        schema_status="V2",
        categories=CURRENT_CLASSES_V2 + (UNKNOWN_CURRENT_CLASS,),
        unknown_counts=unknown,
    )


def select_legacy_view(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
    unversioned_unknown_policy: str = "fail",
) -> ClassView:
    """Select the explicit legacy four-state field, with opt-in H1 fallback."""

    if unversioned_unknown_policy not in {"fail", "exclude"}:
        raise ValueError("unversioned_unknown_policy must be 'fail' or 'exclude'")

    if "VerificationClass_Legacy" in df.columns:
        raw = _clean_string_series(df["VerificationClass_Legacy"])
        unknown = _unknown_counts(raw, LEGACY_CLASSES)
        if raw.isna().any():
            unknown["<MISSING>"] = int(raw.isna().sum())
        if unknown:
            raise SchemaContractError(f"legacy verification view: invalid values: {unknown}")
        return ClassView(
            values=raw,
            field="VerificationClass_Legacy",
            schema_status="LEGACY_EXPLICIT",
            categories=LEGACY_CLASSES,
            unknown_counts={},
        )

    if not allow_unversioned_v1:
        raise SchemaContractError(
            "legacy verification view: VerificationClass_Legacy is missing; "
            "cannot infer legacy classes from current VerificationClass"
        )
    if "VerificationSchemaVersion" in df.columns:
        raise SchemaContractError(
            "legacy verification view: versioned input lacks VerificationClass_Legacy; fallback is forbidden"
        )

    _require_columns(df, ["VerificationClass"], "unversioned v1 fallback")
    raw = _clean_string_series(df["VerificationClass"])
    unknown = _unknown_counts(raw, LEGACY_CLASSES)
    if raw.isna().any():
        unknown["<MISSING>"] = int(raw.isna().sum())
    if unknown and unversioned_unknown_policy == "fail":
        raise SchemaContractError(f"unversioned v1 fallback: invalid legacy values: {unknown}")

    message = (
        "UNVERSIONED: using VerificationClass as legacy four-state input under explicit "
        "--allow-unversioned-v1 authorization"
    )
    warnings.warn(message, UserWarning, stacklevel=2)
    selected = raw.where(raw.isin(LEGACY_CLASSES), pd.NA)
    return ClassView(
        values=selected,
        field="VerificationClass",
        schema_status="UNVERSIONED_V1",
        categories=LEGACY_CLASSES,
        unknown_counts=unknown,
        warning_messages=(message,),
    )


def ordered_class_crosstab(
    labels: pd.Series,
    view: ClassView,
    normalize: bool = False,
    margins: bool = False,
) -> pd.DataFrame:
    """Build a stable crosstab that retains every class, including zero-frequency classes."""

    if not labels.index.equals(view.values.index):
        raise SchemaContractError("class crosstab: label and verification indexes differ")

    valid = labels.notna() & view.values.notna()
    row_order = list(dict.fromkeys(labels[valid].tolist()))
    normalize_arg = "index" if normalize else False
    table = pd.crosstab(labels[valid], view.values[valid], normalize=normalize_arg)
    table = table.reindex(index=row_order, columns=list(view.categories), fill_value=0)
    if margins:
        if normalize:
            raise ValueError("margins are not supported for normalized crosstabs")
        table["All"] = table.sum(axis=1)
        all_row = table.sum(axis=0)
        all_row.name = "All"
        table = pd.concat([table, all_row.to_frame().T])
    return table


def _parse_typed_bool(series: pd.Series, field: str, allow_na: bool) -> pd.Series:
    parsed = []
    invalid: Dict[str, int] = {}
    for value in series.tolist():
        if pd.isna(value) or str(value) == "NA":
            if allow_na:
                parsed.append(pd.NA)
                continue
            invalid["<MISSING>"] = invalid.get("<MISSING>", 0) + 1
            parsed.append(pd.NA)
            continue
        if isinstance(value, bool):
            parsed.append(value)
            continue
        text = str(value)
        if text == "true":
            parsed.append(True)
        elif text == "false":
            parsed.append(False)
        else:
            invalid[str(value)] = invalid.get(str(value), 0) + 1
            parsed.append(pd.NA)
    if invalid:
        raise SchemaContractError(f"evidence field {field}: invalid boolean values: {invalid}")
    return pd.Series(parsed, index=series.index, dtype="boolean", name=field)


def read_evidence(df: pd.DataFrame) -> pd.DataFrame:
    """Return typed v2 evidence fields without deriving any value from class names."""

    fields = (
        "LabelFirstSupport",
        "ClusterFirstSupport",
        "WithinHPSupport",
        "DispersionWarning",
        "EvidencePath",
        "EvidenceDerivation",
    )
    _require_columns(df, fields, "verification evidence")
    result = pd.DataFrame(index=df.index)
    result["LabelFirstSupport"] = _parse_typed_bool(df["LabelFirstSupport"], "LabelFirstSupport", False)
    result["ClusterFirstSupport"] = _parse_typed_bool(
        df["ClusterFirstSupport"], "ClusterFirstSupport", False
    )
    result["WithinHPSupport"] = _parse_typed_bool(df["WithinHPSupport"], "WithinHPSupport", True)
    result["DispersionWarning"] = _parse_typed_bool(df["DispersionWarning"], "DispersionWarning", True)

    paths = _clean_string_series(df["EvidencePath"])
    path_unknown = _unknown_counts(paths, EVIDENCE_PATHS)
    if paths.isna().any():
        path_unknown["<MISSING>"] = int(paths.isna().sum())
    if path_unknown:
        raise SchemaContractError(f"verification evidence: invalid EvidencePath values: {path_unknown}")
    result["EvidencePath"] = paths

    derivations = _clean_string_series(df["EvidenceDerivation"])
    derivation_unknown = _unknown_counts(derivations, EVIDENCE_DERIVATIONS)
    if derivations.isna().any():
        derivation_unknown["<MISSING>"] = int(derivations.isna().sum())
    if derivation_unknown:
        raise SchemaContractError(
            f"verification evidence: invalid EvidenceDerivation values: {derivation_unknown}"
        )
    result["EvidenceDerivation"] = derivations
    return result


def select_loh_legacy(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> ClassView:
    """Select legacy-derived LOH subtype and enforce the deprecated alias contract."""

    canonical = "LOH_Subtype_LegacyVC"
    deprecated = "LOH_Subtype"
    if canonical in df.columns:
        _require_columns(df, [deprecated], "legacy LOH view")
        raw = _clean_string_series(df[canonical])
        unknown = _unknown_counts(raw, LOH_LEGACY_CLASSES)
        if raw.isna().any():
            unknown["<MISSING>"] = int(raw.isna().sum())
        if unknown:
            raise SchemaContractError(f"legacy LOH view: invalid values: {unknown}")
        alias = _clean_string_series(df[deprecated])
        mismatch = ~((raw == alias) | (raw.isna() & alias.isna()))
        if mismatch.any():
            indices = [str(index) for index in df.index[mismatch][:10].tolist()]
            raise SchemaContractError(
                "legacy LOH view: deprecated LOH_Subtype is not an exact alias at rows "
                + ", ".join(indices)
            )
        return ClassView(raw, canonical, "LOH_LEGACY_EXPLICIT", LOH_LEGACY_CLASSES, {})

    if not allow_unversioned_v1:
        raise SchemaContractError(
            "legacy LOH view: LOH_Subtype_LegacyVC is missing; explicit unversioned fallback is required"
        )
    if "VerificationSchemaVersion" in df.columns:
        raise SchemaContractError("legacy LOH view: versioned input cannot fall back to deprecated LOH_Subtype")
    _require_columns(df, [deprecated], "unversioned legacy LOH fallback")
    raw = _clean_string_series(df[deprecated])
    unknown = _unknown_counts(raw, LOH_LEGACY_CLASSES)
    if raw.isna().any():
        unknown["<MISSING>"] = int(raw.isna().sum())
    if unknown:
        raise SchemaContractError(f"unversioned legacy LOH fallback: invalid values: {unknown}")
    message = "UNVERSIONED: using deprecated LOH_Subtype under explicit fallback authorization"
    warnings.warn(message, UserWarning, stacklevel=2)
    return ClassView(raw, deprecated, "UNVERSIONED_V1", LOH_LEGACY_CLASSES, {}, (message,))


def validate_region_strata(df: pd.DataFrame, status: Mapping[str, object]) -> RegionStratumView:
    """Validate RegionStratification v1 status, sentinels, labels and deprecated alias."""

    if not df.index.is_unique:
        raise SchemaContractError("region stratification: dataframe index must be unique")

    required_status = (
        "status",
        "reason",
        "schema_version",
        "eligible_region_count",
        "min_regions_required",
        "assignment_count",
        "n_occupied_region_strata",
        "warning_count",
        "generated_at",
    )
    missing_status = [field for field in required_status if field not in status]
    if missing_status:
        raise SchemaContractError("region status: missing fields: " + ", ".join(missing_status))
    status_value = str(status["status"])
    if status_value not in REGION_STATUS_VALUES:
        raise SchemaContractError(f"region status: invalid status {status_value!r}")

    def nonnegative_integer(field: str) -> int:
        value = status[field]
        if isinstance(value, bool):
            raise SchemaContractError(f"region status: {field} must be a non-negative integer")
        try:
            numeric = float(value)
            integer = int(numeric)
        except (TypeError, ValueError, OverflowError) as exc:
            raise SchemaContractError(f"region status: {field} must be a non-negative integer") from exc
        if numeric != integer or integer < 0:
            raise SchemaContractError(f"region status: {field} must be a non-negative integer")
        return integer

    schema_version = nonnegative_integer("schema_version")
    eligible_region_count = nonnegative_integer("eligible_region_count")
    min_regions_required = nonnegative_integer("min_regions_required")
    expected_assignment_count = nonnegative_integer("assignment_count")
    expected_occupied = nonnegative_integer("n_occupied_region_strata")
    nonnegative_integer("warning_count")
    if schema_version != REGION_STRATIFICATION_SCHEMA_VERSION:
        raise SchemaContractError(
            f"region status: expected schema_version=1, observed {schema_version}"
        )
    if pd.isna(status["reason"]) or str(status["reason"]) == "":
        raise SchemaContractError("region status: reason must not be empty")
    if pd.isna(status["generated_at"]) or str(status["generated_at"]) == "":
        raise SchemaContractError("region status: generated_at must not be empty")
    if expected_occupied > expected_assignment_count:
        raise SchemaContractError(
            "region status: n_occupied_region_strata cannot exceed assignment_count"
        )
    if status_value == "VALID" and expected_assignment_count != eligible_region_count:
        raise SchemaContractError(
            "region status: VALID assignment_count must equal eligible_region_count"
        )
    if status_value == "VALID" and eligible_region_count == 0:
        raise SchemaContractError("region status: VALID cannot have zero eligible regions")
    if status_value == "VALID" and eligible_region_count < min_regions_required:
        raise SchemaContractError(
            "region status: VALID eligible_region_count is below min_regions_required"
        )
    if (
        status_value == "INSUFFICIENT_REGIONS"
        and eligible_region_count > 0
        and eligible_region_count >= min_regions_required
    ):
        raise SchemaContractError(
            "region status: INSUFFICIENT_REGIONS conflicts with eligible_region_count"
        )

    columns = (
        "RegionStratificationSchemaVersion",
        "RegionStratum_ID",
        "RegionStratum_Label",
        "RegionStratum_Reason",
    )
    _require_columns(df, columns, "region stratification")
    if not df.empty:
        _require_single_integer_version(
            df,
            "RegionStratificationSchemaVersion",
            REGION_STRATIFICATION_SCHEMA_VERSION,
            "region stratification",
        )

    ids = pd.to_numeric(df["RegionStratum_ID"], errors="coerce")
    if ids.isna().any() or (~ids.isin([-1, 0, 1, 2, 3])).any():
        raise SchemaContractError("region stratification: RegionStratum_ID contains invalid values")
    ids = ids.astype(int)
    labels = _clean_string_series(df["RegionStratum_Label"])
    reasons = _clean_string_series(df["RegionStratum_Reason"])
    errors = []
    for index in df.index:
        stratum_id = int(ids.loc[index])
        label = labels.loc[index]
        reason = reasons.loc[index]
        if stratum_id == -1:
            expected_label = "Unassigned"
            expected_reason = REGION_UNASSIGNED_REASON_BY_STATUS[status_value]
        else:
            if status_value != "VALID":
                errors.append(f"row {index}: assigned stratum under status {status_value}")
                continue
            expected_label, expected_reason = REGION_ASSIGNMENTS[stratum_id]
        if label != expected_label or reason != expected_reason:
            errors.append(
                f"row {index}: expected ({expected_label},{expected_reason}), observed ({label},{reason})"
            )
    if errors:
        raise SchemaContractError("region stratification: " + "; ".join(errors[:10]))

    if "Subclone_ID" in df.columns:
        deprecated_ids = pd.to_numeric(df["Subclone_ID"], errors="coerce")
        mismatch = deprecated_ids.isna() | (deprecated_ids.astype("Int64") != ids.astype("Int64"))
        if mismatch.any():
            raise SchemaContractError("region stratification: Subclone_ID is not an exact deprecated alias")

    assignment_count = int((ids >= 0).sum())
    occupied = int(ids[ids >= 0].nunique())
    if assignment_count != expected_assignment_count:
        raise SchemaContractError(
            f"region stratification: assignment_count expected {expected_assignment_count}, observed {assignment_count}"
        )
    if occupied != expected_occupied:
        raise SchemaContractError(
            f"region stratification: occupied count expected {expected_occupied}, observed {occupied}"
        )
    if status_value != "VALID" and (assignment_count != 0 or occupied != 0):
        raise SchemaContractError("region stratification: non-VALID status must have zero assignments")

    normalized = df.copy()
    normalized["RegionStratum_ID"] = ids
    return RegionStratumView(
        frame=normalized,
        status={key: str(value) for key, value in status.items()},
        assignment_count=assignment_count,
        n_occupied_region_strata=occupied,
    )


def extract_provenance_frame(
    df: pd.DataFrame,
    allow_unversioned_v1: bool = False,
) -> pd.DataFrame:
    """Validate and copy the authoritative provenance columns for derived tables."""

    if not df.index.is_unique:
        raise SchemaContractError("provenance export: dataframe index must be unique")

    if "VerificationSchemaVersion" not in df.columns:
        if not allow_unversioned_v1:
            raise SchemaContractError(
                "provenance export: VerificationSchemaVersion missing; explicit unversioned fallback required"
            )
        legacy = select_legacy_view(df, allow_unversioned_v1=True, unversioned_unknown_policy="fail")
        columns = [legacy.field]
        if "LOH_Subtype" in df.columns:
            select_loh_legacy(df, allow_unversioned_v1=True)
            columns.append("LOH_Subtype")
        result = df[columns].copy()
        result.attrs.update(legacy.metadata())
        return result

    current = select_current_view(df)
    legacy = select_legacy_view(df)
    evidence = read_evidence(df)
    select_loh_legacy(df)
    _require_columns(df, VERIFICATION_PROVENANCE_COLUMNS, "provenance export")

    if current.unknown_counts:
        raise SchemaContractError(
            f"provenance export: canonical VerificationClass has unknown values: {current.unknown_counts}"
        )

    deprecated = _clean_string_series(df["VerificationClass_V1_Deprecated"])
    deprecated_unknown = _unknown_counts(deprecated, V1_CURRENT_CLASSES)
    if deprecated.isna().any():
        deprecated_unknown["<MISSING>"] = int(deprecated.isna().sum())
    if deprecated_unknown:
        raise SchemaContractError(
            "provenance export: invalid VerificationClass_V1_Deprecated values: "
            f"{deprecated_unknown}"
        )

    relation_errors = []
    for index in df.index:
        current_value = str(current.values.loc[index])
        legacy_value = str(legacy.values.loc[index])
        deprecated_value = str(deprecated.loc[index])
        evidence_path = str(evidence.loc[index, "EvidencePath"])

        expected_deprecated = (
            "Strong"
            if current_value in {"Strong_Bidirectional", "ClusterFirstOnly"}
            else current_value
        )
        if deprecated_value != expected_deprecated:
            relation_errors.append(
                f"row {index}: deprecated={deprecated_value}, expected {expected_deprecated}"
            )

        expected_path = CURRENT_TO_EVIDENCE_PATH[current_value]
        if evidence_path != expected_path:
            relation_errors.append(
                f"row {index}: EvidencePath={evidence_path}, expected {expected_path}"
            )

        if current_value == "Strong_Bidirectional" and legacy_value != "Strong":
            relation_errors.append(
                f"row {index}: Strong_Bidirectional requires legacy Strong"
            )
        elif current_value == "ClusterFirstOnly" and legacy_value != "Subclone":
            relation_errors.append(
                f"row {index}: ClusterFirstOnly requires legacy Subclone"
            )
        elif current_value not in {"Strong_Bidirectional", "ClusterFirstOnly"} and legacy_value not in {
            "Weak",
            "Noise",
        }:
            relation_errors.append(
                f"row {index}: {current_value} requires legacy Weak or Noise"
            )

        expected_label_first = legacy_value in {"Strong", "Weak"}
        expected_cluster_first = legacy_value in {"Strong", "Subclone"}
        if bool(evidence.loc[index, "LabelFirstSupport"]) != expected_label_first:
            relation_errors.append(
                f"row {index}: LabelFirstSupport conflicts with legacy {legacy_value}"
            )
        if bool(evidence.loc[index, "ClusterFirstSupport"]) != expected_cluster_first:
            relation_errors.append(
                f"row {index}: ClusterFirstSupport conflicts with legacy {legacy_value}"
            )
        if evidence.loc[index, "EvidenceDerivation"] == "LEGACY_CLASS":
            if pd.notna(evidence.loc[index, "WithinHPSupport"]):
                relation_errors.append(
                    f"row {index}: LEGACY_CLASS requires WithinHPSupport=NA"
                )
            if pd.notna(evidence.loc[index, "DispersionWarning"]):
                relation_errors.append(
                    f"row {index}: LEGACY_CLASS requires DispersionWarning=NA"
                )

    if relation_errors:
        raise SchemaContractError("provenance export: " + "; ".join(relation_errors[:10]))

    result = df[list(VERIFICATION_PROVENANCE_COLUMNS)].copy()
    result.attrs["schema_status"] = "V2"
    return result
