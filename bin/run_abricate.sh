#!/bin/bash

# Usage: ./run_abricate.sh <SampleName> <ConsensusFASTA> <SeroDB> <VFDB>

set -euo pipefail

# Input arguments
SampleName="$1"
Consensus="$2"         # e.g., path/to/sample_assembly.fasta
SeroDB="$3"            # e.g., path/to/abricate_db/sero
VFDB="$4"              # e.g., path/to/abricate_db/vf

# Define default values
DefaultLine="${SampleName}\t${SampleName}_contig_1\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone"
HeaderOnly=1  # Expected line count if only the header is present

# Output filenames
SeroOut="${SampleName}_sero.csv"
VFOut="${SampleName}_vf.csv"
AMROut="${SampleName}_AMR.csv"

# Run serotype search
abricate --datadir "$SeroDB" --db Ssuis_serotype --minid 60 --mincov 60 --quiet "$Consensus" > "$SeroOut"
sed -i 's,_assembly.fasta,,g' "$SeroOut"
if [ "$(wc -l < "$SeroOut")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$SeroOut"
fi

# Run virulence factor search
abricate --datadir "$VFDB" --db Ssuis_vfdb "$Consensus" > "$VFOut"
sed -i 's,_assembly.fasta,,g' "$VFOut"
if [ "$(wc -l < "$VFOut")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$VFOut"
fi

# Run AMR search using CARD
abricate --db card "$Consensus" > "$AMROut"
sed -i 's,_assembly.fasta,,g' "$AMROut"
if [ "$(wc -l < "$AMROut")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$AMROut"
fi

