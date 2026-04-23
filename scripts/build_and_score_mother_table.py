#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import pandas as pd


WHO_SCORE = {
    "Critical Priority": 10,
    "High Priority": 5,
    "Medium Priority": 3,
    "Not listed": 1,
    "": 1,
}

ARG_RANK_SCORE = {
    "Rank I": 5,
    "Rank II": 4,
    "Rank III": 3,
    "Rank IV": 2,
    "": 1,
}


def load_tsv(path, header=True, colnames=None):
    if header:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return pd.read_csv(path, sep="\t", dtype=str, header=None, names=colnames).fillna("")


def load_table_auto(path):
    return pd.read_csv(path, sep=None, engine="python", dtype=str).fillna("")


def normalize_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalize_species_name(x):
    x = normalize_text(x)
    alias_map = {
        "Salmonella enterica serovar Typhi": "Salmonella Typhi",
        "Salmonella enterica subsp. enterica serovar Typhi": "Salmonella Typhi",
        "Salmonella enterica": "Salmonella enterica",
        "Klebsiella_pneumoniae": "Klebsiella pneumoniae",
        "Serratia_marcescens": "Serratia marcescens",
        "Sinorhizobium_meliloti": "Sinorhizobium meliloti",
        "Yersinia_pseudotuberculosis": "Yersinia pseudotuberculosis",
        "Citrobacter_freundii": "Citrobacter freundii",
        "Enterobacter_cloacae": "Enterobacter cloacae",
        "Enterobacter_sp_MGH_14": "Enterobacter sp. MGH_14",
        "Enterobacter_sp_MGH_26": "Enterobacter sp. MGH_26",
        "Gallibacterium_anatis": "Gallibacterium anatis",
        "Escherichia coli": "Escherichia coli",
        "Klebsiella pneumoniae": "Klebsiella pneumoniae",
        "Haemophilus influenzae": "Haemophilus influenzae",
        "Shigella boydii": "Shigella boydii",
        "Enterococcus faecium": "Enterococcus faecium",
        "Staphylococcus aureus": "Staphylococcus aureus",
        "Neisseria gonorrhoeae": "Neisseria gonorrhoeae",
        "Streptococcus pyogenes": "Streptococcus pyogenes",
        "Streptococcus pneumoniae": "Streptococcus pneumoniae",
        "Streptococcus agalactiae": "Streptococcus agalactiae",
        "Pseudomonas aeruginosa": "Pseudomonas aeruginosa",
        "Mycobacterium tuberculosis": "Mycobacterium tuberculosis",
    }
    return alias_map.get(x, x.replace("_", " "))


def normalize_genus_name(x):
    return normalize_text(x).replace("_", " ")


def normalize_order_name(x):
    x = normalize_text(x)
    alias_map = {
        "Enterobacteriales": "Enterobacterales",
    }
    return alias_map.get(x, x)


def load_rules(rule_dir):
    who_rules = load_tsv(rule_dir / "who_rules.tsv", header=True)
    arg_map = load_tsv(rule_dir / "arg_to_phenotype.tsv", header=True)
    priority_order = load_tsv(rule_dir / "priority_order.tsv", header=True)

    for df in [who_rules, arg_map, priority_order]:
        for col in df.columns:
            df[col] = df[col].apply(normalize_text)

    priority_score = {
        row["priority"]: int(row["score"])
        for _, row in priority_order.iterrows()
    }
    return who_rules, arg_map, priority_score


def match_arg_to_phenotypes(arg_subtype, arg_map_df):
    arg_subtype = normalize_text(arg_subtype)
    if not arg_subtype:
        return []

    hits = []

    for _, row in arg_map_df.iterrows():
        pattern = normalize_text(row["arg_pattern"])
        match_type = normalize_text(row["match_type"])

        matched = False
        if match_type == "exact":
            matched = (arg_subtype == pattern)
        elif match_type == "prefix":
            matched = arg_subtype.startswith(pattern)

        if matched:
            hits.append({
                "phenotype": normalize_text(row["phenotype"]),
                "confidence": normalize_text(row["confidence"]),
                "map_id": normalize_text(row["map_id"]),
                "note": normalize_text(row["note"]),
            })

    # 兼容 grouped 注释
    if arg_subtype in {"CTX-M-group1", "CTX-M-group2", "CTXM-group1", "CTXM-group2"}:
        hits.append({
            "phenotype": "third-generation-cephalosporin-resistant",
            "confidence": "high",
            "map_id": "LOCAL_CTXM_GROUP",
            "note": "Grouped CTX-M annotation mapped locally",
        })

    if arg_subtype in {"TEM-group1", "TEM_group1"}:
        hits.append({
            "phenotype": "ampicillin-resistant",
            "confidence": "medium",
            "map_id": "LOCAL_TEM_GROUP1",
            "note": "Grouped TEM annotation mapped locally",
        })

    dedup = {}
    for h in hits:
        key = (h["phenotype"], h["confidence"], h["map_id"])
        dedup[key] = h

    return list(dedup.values())


def record_matches_rule(record, phenotype, who_rule):
    taxon_level = normalize_text(who_rule["taxon_level"])
    taxon_name = normalize_text(who_rule["taxon_name"])
    rule_phenotype = normalize_text(who_rule["phenotype"])
    rule_id = normalize_text(who_rule["rule_id"])

    record_species = normalize_species_name(record.get("species", ""))
    record_genus = normalize_genus_name(record.get("genus", ""))
    record_order = normalize_order_name(record.get("order", ""))

    if phenotype != rule_phenotype:
        return False

    if taxon_level == "species":
        return record_species == taxon_name

    elif taxon_level == "genus":
        if rule_id == "WHO009":
            return (record_genus == "Salmonella") and (record_species != "Salmonella Typhi")
        return record_genus == taxon_name

    elif taxon_level == "order":
        return record_order == taxon_name

    return False


def annotate_row(record, who_rules_df, arg_map_df, priority_score):
    arg_subtype = normalize_text(record.get("ARG_subtype", ""))
    phenotype_hits = match_arg_to_phenotypes(arg_subtype, arg_map_df)

    if not phenotype_hits:
        return {
            "phenotype": "",
            "phenotype_confidence": "",
            "phenotype_map_ids": "",
            "phenotype_notes": "",
            "WHO_rule_id": "",
            "WHO_all_rule_ids": "",
            "WHO_pathogen_label": "",
            "WHO_all_pathogen_labels": "",
            "WHO_priority": "Not listed",
            "WHO_all_priorities": "",
            "WHO_match_note": "",
            "WHO_implemented": "",
        }

    matched_rules = []
    all_phenotypes = []
    all_confidences = []
    all_map_ids = []
    all_notes = []

    for ph in phenotype_hits:
        phenotype = ph["phenotype"]
        all_phenotypes.append(phenotype)
        all_confidences.append(ph["confidence"])
        all_map_ids.append(ph["map_id"])
        all_notes.append(ph["note"])

        for _, rule in who_rules_df.iterrows():
            if record_matches_rule(record, phenotype, rule):
                matched_rules.append({
                    "rule_id": normalize_text(rule["rule_id"]),
                    "priority": normalize_text(rule["priority"]),
                    "pathogen_label": normalize_text(rule["pathogen_label"]),
                    "match_note": normalize_text(rule["match_note"]),
                    "implemented": normalize_text(rule["implemented"]),
                })

    dedup_rules = {}
    for r in matched_rules:
        dedup_rules[r["rule_id"]] = r
    matched_rules = list(dedup_rules.values())

    if not matched_rules:
        return {
            "phenotype": ";".join(sorted(set(all_phenotypes))),
            "phenotype_confidence": ";".join(sorted(set(all_confidences))),
            "phenotype_map_ids": ";".join(sorted(set(all_map_ids))),
            "phenotype_notes": "; ".join(sorted(set(all_notes))),
            "WHO_rule_id": "",
            "WHO_all_rule_ids": "",
            "WHO_pathogen_label": "",
            "WHO_all_pathogen_labels": "",
            "WHO_priority": "Not listed",
            "WHO_all_priorities": "",
            "WHO_match_note": "",
            "WHO_implemented": "",
        }

    def rule_score(r):
        return priority_score.get(r["priority"], -999)

    best_rule = sorted(matched_rules, key=rule_score, reverse=True)[0]

    return {
        "phenotype": ";".join(sorted(set(all_phenotypes))),
        "phenotype_confidence": ";".join(sorted(set(all_confidences))),
        "phenotype_map_ids": ";".join(sorted(set(all_map_ids))),
        "phenotype_notes": "; ".join(sorted(set(all_notes))),
        "WHO_rule_id": best_rule["rule_id"],
        "WHO_all_rule_ids": ";".join(sorted(r["rule_id"] for r in matched_rules)),
        "WHO_pathogen_label": best_rule["pathogen_label"],
        "WHO_all_pathogen_labels": "; ".join(sorted(set(r["pathogen_label"] for r in matched_rules))),
        "WHO_priority": best_rule["priority"],
        "WHO_all_priorities": ";".join(sorted(set(r["priority"] for r in matched_rules))),
        "WHO_match_note": best_rule["match_note"],
        "WHO_implemented": best_rule["implemented"],
    }


def canonicalize_arg_name(x):
    """
    用于 ARG subtype 标准化匹配：
    - 去空格
    - 统一下划线/连字符
    - 去掉常见 bla 前缀
    - 将 CTX-M 统一成 CTXM
    """
    x = normalize_text(x)
    if not x:
        return ""

    x = x.replace("_", "-")
    x = x.replace("–", "-")
    x = x.replace("—", "-")
    x = x.replace(" ", "")

    xl = x.lower()
    if xl.startswith("bla"):
        x = x[3:]

    x = x.replace("CTX-M", "CTXM")
    x = x.replace("ctx-m", "CTXM")
    x = x.replace("Ctx-M", "CTXM")

    return x


def load_arg_rank_map(arg_csv_path):
    arg_df = load_table_auto(arg_csv_path)

    for c in arg_df.columns:
        arg_df[c] = arg_df[c].apply(normalize_text)

    if "ARGsubtype" not in arg_df.columns or "Level" not in arg_df.columns:
        raise ValueError("ARG.csv must contain columns: ARGsubtype and Level")

    best_map = {}

    for _, row in arg_df.iterrows():
        raw_name = normalize_text(row["ARGsubtype"])
        level = normalize_text(row["Level"])
        canon = canonicalize_arg_name(raw_name)

        old = best_map.get(canon, "")
        if ARG_RANK_SCORE.get(level, 1) > ARG_RANK_SCORE.get(old, 1):
            best_map[canon] = level

    return best_map


def match_arg_rank(arg_subtype, arg_rank_map):
    raw = normalize_text(arg_subtype)
    if not raw:
        return "", "none"

    canon = canonicalize_arg_name(raw)

    # 1) exact canonical match
    if canon in arg_rank_map:
        return arg_rank_map[canon], "exact"

    # 2) prefix/family fallback, choose highest score
    candidates = []
    for ref_name, level in arg_rank_map.items():
        if canon.startswith(ref_name) or ref_name.startswith(canon):
            candidates.append((ref_name, level))

    if candidates:
        candidates = sorted(
            candidates,
            key=lambda x: (ARG_RANK_SCORE.get(x[1], 1), len(x[0])),
            reverse=True
        )
        return candidates[0][1], "prefix"

    return "", "none"


def calc_mobile_factor(plasmid_class):
    plasmid_class = normalize_text(plasmid_class)
    if plasmid_class in {"Conj", "mob_unconj"}:
        return 2
    return 1


def calc_conj_factor(plasmid_class):
    plasmid_class = normalize_text(plasmid_class)
    if plasmid_class == "Conj":
        return 5
    return 1


def build_mother_table(arg_taxa_path, plasmid_sum_path, rules_dir):
    arg_taxa = load_tsv(arg_taxa_path, header=True)
    plasmid_sum = load_tsv(
        plasmid_sum_path,
        header=False,
        colnames=["query", "plasmid_class", "plasmid_arg_summary"]
    )

    for col in arg_taxa.columns:
        arg_taxa[col] = arg_taxa[col].apply(normalize_text)

    for col in plasmid_sum.columns:
        plasmid_sum[col] = plasmid_sum[col].apply(normalize_text)

    if "subtype" not in arg_taxa.columns:
        raise ValueError("Input ARG+taxa table must contain column: subtype")
    if "type" not in arg_taxa.columns:
        raise ValueError("Input ARG+taxa table must contain column: type")
    if "query" not in arg_taxa.columns:
        raise ValueError("Input ARG+taxa table must contain column: query")

    for col in ["species", "genus", "order"]:
        if col not in arg_taxa.columns:
            arg_taxa[col] = ""

    arg_taxa["species"] = arg_taxa["species"].apply(normalize_species_name)
    arg_taxa["genus"] = arg_taxa["genus"].apply(normalize_genus_name)
    arg_taxa["order"] = arg_taxa["order"].apply(normalize_order_name)

    arg_taxa["ARG_subtype"] = arg_taxa["subtype"]
    arg_taxa["ARG_type"] = arg_taxa["type"]

    mother = arg_taxa.merge(plasmid_sum, on="query", how="left").fillna("")

    who_rules_df, arg_map_df, priority_score = load_rules(Path(rules_dir))

    annotations = mother.apply(
        lambda row: pd.Series(annotate_row(row, who_rules_df, arg_map_df, priority_score)),
        axis=1
    )

    mother = pd.concat([mother, annotations], axis=1)
    return mother


def score_mother_table(mother_df, arg_csv_path):
    arg_rank_map = load_arg_rank_map(arg_csv_path)

    arg_rank_list = []
    arg_rank_score_list = []
    arg_rank_match_method_list = []
    who_score_list = []
    mobile_factor_list = []
    conj_factor_list = []
    risk_score_list = []

    for _, row in mother_df.iterrows():
        arg_subtype = normalize_text(row.get("ARG_subtype", ""))
        who_priority = normalize_text(row.get("WHO_priority", ""))
        plasmid_class = normalize_text(row.get("plasmid_class", ""))

        arg_rank, match_method = match_arg_rank(arg_subtype, arg_rank_map)

        who_score = WHO_SCORE.get(who_priority, 1)
        rank_score = ARG_RANK_SCORE.get(arg_rank, 1)
        mobile_factor = calc_mobile_factor(plasmid_class)
        conj_factor = calc_conj_factor(plasmid_class)

        risk_score = who_score * rank_score * mobile_factor * conj_factor

        arg_rank_list.append(arg_rank)
        arg_rank_score_list.append(rank_score)
        arg_rank_match_method_list.append(match_method)
        who_score_list.append(who_score)
        mobile_factor_list.append(mobile_factor)
        conj_factor_list.append(conj_factor)
        risk_score_list.append(risk_score)

    df = mother_df.copy()
    df["ARG_rank"] = arg_rank_list
    df["ARG_rank_score"] = arg_rank_score_list
    df["ARG_rank_match_method"] = arg_rank_match_method_list
    df["WHO_priority_score"] = who_score_list
    df["mobile_factor"] = mobile_factor_list
    df["conjugative_plasmid_factor"] = conj_factor_list
    df["risk_score"] = risk_score_list

    preferred_cols = [
        "query",
        "ARG_subtype", "ARG_type",
        "ARG_rank", "ARG_rank_score", "ARG_rank_match_method",
        "WHO_priority", "WHO_priority_score",
        "plasmid_class", "mobile_factor", "conjugative_plasmid_factor",
        "risk_score",
        "subtype", "type",
        "q.start", "q.end",
        "kingdom", "phylum", "class", "order", "family", "genus", "species",
        "plasmid.name", "source",
        "plasmid_arg_summary",
        "phenotype", "phenotype_confidence", "phenotype_map_ids", "phenotype_notes",
        "WHO_rule_id", "WHO_all_rule_ids",
        "WHO_pathogen_label", "WHO_all_pathogen_labels",
        "WHO_all_priorities",
        "WHO_match_note", "WHO_implemented"
    ]

    existing = [c for c in preferred_cols if c in df.columns]
    remaining = [c for c in df.columns if c not in existing]
    df = df[existing + remaining]

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Build and score mother_table.tsv from ARG+taxa, plasmid classification summary, WHO rules, and ARG rank CSV."
    )
    parser.add_argument("--arg_taxa", required=True, help="Path to test.fa.uniq_arg.w.taxa.tab")
    parser.add_argument("--plasmid_sum", required=True, help="Path to test_plasmids_classification_sum.txt")
    parser.add_argument("--rules", required=True, help="Rules directory")
    parser.add_argument("--arg_csv", required=True, help="Path to ARG.csv")
    parser.add_argument("-o", "--output", required=True, help="Output final scored mother_table.tsv")
    parser.add_argument("--mother_only_output", default="", help="Optional output path for intermediate mother_table.tsv")
    args = parser.parse_args()

    mother = build_mother_table(args.arg_taxa, args.plasmid_sum, args.rules)

    if args.mother_only_output:
        Path(args.mother_only_output).parent.mkdir(parents=True, exist_ok=True)
        mother.to_csv(args.mother_only_output, sep="\t", index=False)
        print(f"[OK] Intermediate mother table written to: {args.mother_only_output}")

    scored = score_mother_table(mother, args.arg_csv)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.output, sep="\t", index=False)

    print(f"[OK] Final scored mother table written to: {args.output}")
    print(f"[OK] Total rows: {len(scored)}")


if __name__ == "__main__":
    main()