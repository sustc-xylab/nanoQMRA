#!/usr/bin/env python3
import argparse
import pandas as pd


def safe_num(x, default=1):
    if pd.isna(x) or str(x).strip() == "":
        return default
    try:
        return float(x)
    except Exception:
        return default


def is_plasmid_row(row):
    source = str(row.get("source", "")).strip()
    plasmid_name = str(row.get("plasmid.name", "")).strip()
    return source == "plasmid.like" or plasmid_name != ""


def get_raw_host_label(row):
    # 这里读取原始宿主分类信息，不能读取改写后的 species
    raw_species = str(row.get("species", "")).strip()
    raw_genus = str(row.get("genus", "")).strip()

    if raw_species:
        return raw_species
    if raw_genus:
        return raw_genus
    return "unclassified"


def get_display_species_label(row):
    if is_plasmid_row(row):
        return "plasmid"

    raw_species = str(row.get("species", "")).strip()
    raw_genus = str(row.get("genus", "")).strip()

    if raw_species:
        return raw_species
    if raw_genus:
        return raw_genus
    return "unclassified"


def build_risk_score_list(df):
    out = df.copy()

    out["risk_score"] = out.apply(
        lambda r: (
            safe_num(r.get("ARG_rank_score", 1), 1)
            * safe_num(r.get("WHO_priority_score", 1), 1)
            * safe_num(r.get("mobile_factor", 1), 1)
            * safe_num(r.get("conjugative_plasmid_factor", 1), 1)
        ),
        axis=1
    )

    # host 必须先从原始分类列取
    out["host"] = out.apply(
        lambda r: get_raw_host_label(r) if is_plasmid_row(r) else "",
        axis=1
    )

    # species 是展示列；若为质粒则统一写 plasmid
    out["species_display"] = out.apply(get_display_species_label, axis=1)

    out = out[["species_display", "ARG_type", "risk_score", "host"]].copy()
    out.columns = ["species", "ARG_type", "risk_score", "host"]

    out = out.sort_values(
        by=["risk_score", "species", "ARG_type", "host"],
        ascending=[False, True, True, True]
    ).reset_index(drop=True)

    out.insert(0, "index", range(1, len(out) + 1))
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Generate risk_score_list.tsv from mother_table.scored.tsv"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input mother_table.scored.tsv"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output risk_score_list.tsv"
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", dtype=str).fillna("")
    out = build_risk_score_list(df)
    out.to_csv(args.output, sep="\t", index=False)

    print(f"[OK] wrote: {args.output}")


if __name__ == "__main__":
    main()