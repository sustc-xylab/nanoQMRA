# Implementation notes for WHO priority matching

## Rule design
WHO priority is assigned only when both conditions are satisfied:

1. taxon/pathogen group matches
2. resistance phenotype matches

This is a conservative implementation based on pathogen + resistance phenotype rather than pathogen name alone.

## Important limitations

### 1. Carbapenem resistance
Only strong carbapenemase families are included by default:
- KPC
- NDM
- VIM
- IMP

OXA-family genes are not globally mapped to carbapenem-resistant because many OXA variants do not reliably indicate carbapenem resistance across taxa.

### 2. Third-generation cephalosporin resistance
Only conservative markers are included by default:
- CTX-M family
- CMY family
- VEB family

TEM and SHV families are not globally mapped because many variants do not confer third-generation cephalosporin resistance.

### 3. Fluoroquinolone resistance
Current implementation only uses plasmid-mediated quinolone resistance genes:
- qnrA
- qnrB
- qnrS

Chromosomal mutations such as gyrA/parC are not included unless a mutation-calling pipeline is available.

### 4. Rifampicin-resistant Mycobacterium tuberculosis
Not implemented in the current ARG-only pipeline because rifampicin resistance usually depends on rpoB mutations rather than mobile ARG markers.

### 5. Penicillin-resistant Group B Streptococcus
Not implemented in the current ARG-only pipeline because penicillin resistance is usually associated with target-site alterations rather than standard ARG markers.

### 6. Ampicillin-resistant Haemophilus influenzae
ROB-1 is treated as high-confidence.
TEM-1 is treated as medium-confidence and should be interpreted cautiously.

### 7. Macrolide resistance
erm genes are treated as high-confidence markers.
mph/msr/lnu genes are included but may require cautious interpretation depending on species and phenotype definition.

## Recommended output columns
Recommended fields in the merged output table:

- ARG_subtype
- ARG_type
- ARG_rank
- phenotype
- phenotype_confidence
- WHO_rule_id
- WHO_pathogen_label
- WHO_priority
- WHO_match_note
- WHO_implemented

## Priority assignment behavior
If multiple WHO rules match the same record, keep:
1. the highest priority level
2. all matched rule IDs in a semicolon-separated field

Priority order:
- Critical Priority
- High Priority
- Medium Priority
- Not listed

## Recommended future extensions
- Add OXA subtype whitelist for Acinetobacter baumannii and selected taxa
- Add ESBL subtype whitelist for SHV/TEM where variant-level annotation is available
- Add rpoB mutation support for Mycobacterium tuberculosis
- Add gyrA/parC mutation support for fluoroquinolone resistance
- Add pbp mutation support for penicillin-resistant Streptococcus agalactiae