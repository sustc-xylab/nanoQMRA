# nanoQMRA

`nanoQMRA` is a bioinformatics tool for quantitative microbial risk assessment (QMRA) based on nanopore sequencing data. It is designed to calculate QMRA-related risk values for both individual reads in species level.

For each read, `nanoQMRA` performs independent scoring and indicates:

- whether the read is associated with a pathogenic host
- whether the corresponding genetic element has mobility
- whether the read carries antibiotic resistance genes (ARGs)

By integrating these features, `nanoQMRA` generates:

- a per-read scored result table
- a ranked per-read risk list
- a final sample-level QMRA risk score

Therefore, `nanoQMRA` can be used as a nanopore-based bioinformatics workflow for microbial risk evaluation by combining **pathogen association**, **mobility potential**, and **ARG carriage**.

## Download and environment setup

```bash
git clone https://github.com/your-username/nanoQMRA.git
cd nanoQMRA

nanoQMRA is designed to run in a Linux environment.

Basic tools required:

bash
python3
awk
grep
find
sort
cp
mv
External tools required:

ARGpore2
Plascad
Please also make sure that:

all required downstream scripts are present
all required reference files are available in the input/ directory
the conda environment for Plascad is correctly installed if needed
Example of preparing a basic conda environment:
conda create -n nanoqmra python=3.10 -y
conda activate nanoqmra
