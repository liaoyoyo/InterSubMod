import pandas as pd
import os
import argparse
import glob
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.verification_schema_contract import (  # noqa: E402
    SchemaContractError,
    select_legacy_view,
)


LEGACY_SELECTION_COLUMN = "_VerificationClass_Legacy_Selected"

def setup_args():
    parser = argparse.ArgumentParser(description="Find verification candidates from significance summary.")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory (containing significance_summary.csv)")
    parser.add_argument("--top_n", type=int, default=5, help="Number of candidates to list per category")
    return parser.parse_args()

def find_heatmap(base_site_path):
    # Try to find heatmap in distance/METRIC/*.png
    # or just recursively *.png in the site dir
    if not os.path.exists(base_site_path):
        return None
        
    # Search for any heatmap_distance.png
    files = glob.glob(os.path.join(base_site_path, "**", "heatmap_distance.png"), recursive=True)
    if files:
        return files[0]
        
    # Search for any png in distance/ folder
    files = glob.glob(os.path.join(base_site_path, "distance", "**", "*.png"), recursive=True)
    if files:
        return files[0]
        
    return None


def prepare_legacy_candidates(df):
    """Validate and attach the explicit four-state legacy taxonomy."""
    view = select_legacy_view(df, allow_unversioned_v1=False)
    selected = df.copy()
    selected[LEGACY_SELECTION_COLUMN] = view.values
    return selected, view

def analyze_candidates(df, args):
    categories = {
        "Weak": {"filter": df[LEGACY_SELECTION_COLUMN] == "Weak", "sort": [("LabelP", True), ("LabelDelta", False)]},
        "Subclone": {"filter": df[LEGACY_SELECTION_COLUMN] == "Subclone", "sort": [("GlobalP", True), ("LabelP", False)]},
        "Strong": {"filter": df[LEGACY_SELECTION_COLUMN] == "Strong", "sort": [("GlobalP", True)]},
        "Noise": {"filter": df[LEGACY_SELECTION_COLUMN] == "Noise", "sort": [("GlobalP", False)]}
    }

    print("# Legacy verification candidates")
    print("Selection field: `VerificationClass_Legacy`")
    print(f"Input Directory: {args.output_dir}\n")

    for cat_name, logic in categories.items():
        subset = df[logic["filter"]].copy()
        print(f"## Category: {cat_name} (Total: {len(subset)})")
        
        if subset.empty:
            print("No candidates found.\n")
            continue

        # Sorting
        sort_cols = [x[0] for x in logic["sort"]]
        sort_asc = [x[1] for x in logic["sort"]]
        subset = subset.sort_values(by=sort_cols, ascending=sort_asc)
        
        top_candidates = subset.head(args.top_n)
        
        for _, row in top_candidates.iterrows():
            chrom = row["Chr"]
            pos = row["Pos"]
            site_id = f"{chrom}_{pos}"
            
            # Construct base path
            # Search for the directory starting with site_id prefix in the chrom directory
            # Because exact folder might be chr1_123_window...
            chrom_dir = os.path.join(args.output_dir, "filtered_snv_tp", chrom)
            
            site_dir = None
            region_dir = None
            
            if os.path.exists(chrom_dir):
                # pattern: chr1_123* (to catch chr1_123_456_789)
                # Matches the SITE folder (chr_pos)
                site_folder_candidates = glob.glob(os.path.join(chrom_dir, f"{site_id}")) # Exact match usually for site folder
                if site_folder_candidates:
                    site_dir = site_folder_candidates[0]
                    # Format: .../chr1_POS/chr1_START_END/
                    # Find inner subdirectory
                    inner_candidates = glob.glob(os.path.join(site_dir, "*"))
                    # Filter for directories
                    inner_dirs = [d for d in inner_candidates if os.path.isdir(d)]
                    if inner_dirs:
                        region_dir = inner_dirs[0]

            heatmap_path = "File not found"
            gen_cmd = "Region directory not found"
            
            if region_dir:
                gen_cmd = f"python3 tools/plot_distance_heatmap.py --region-dir {region_dir} --metric BERNOULLI"
                found_file = find_heatmap(region_dir)
                if found_file:
                    heatmap_path = found_file
                    gen_cmd = "" # Found, no need to gen
            
            print(f"- **Site**: {site_id}")
            print(f"  - Stats: GlobalP={row['GlobalP']:.3e}, LabelP={row['LabelP']:.3e}, LabelDelta={row['LabelDelta']:.3f}, LegacyClass={row[LEGACY_SELECTION_COLUMN]}")
            if heatmap_path != "File not found":
                print(f"  - Heatmap: `{heatmap_path}`")
            else:
                print(f"  - Heatmap: MISSING. Generate command: `{gen_cmd}`")
        print("\n")

def main():
    args = setup_args()
    
    csv_path = os.path.join(args.output_dir, "significance_summary.csv")
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        sys.exit(1)
        
    try:
        df = pd.read_csv(csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        sys.exit(1)
        
    # Ensure columns exist
    required_cols = ["Chr", "Pos", "VerificationClass_Legacy", "GlobalP", "LabelP", "LabelDelta"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Error: Missing columns in CSV: {missing}")
        sys.exit(1)
        
    try:
        selected, _ = prepare_legacy_candidates(df)
    except SchemaContractError as exc:
        print(f"Error: legacy verification schema contract failed: {exc}", file=sys.stderr)
        sys.exit(2)

    analyze_candidates(selected, args)

if __name__ == "__main__":
    main()
