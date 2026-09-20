"""
GSE303448 - Step 01: inspect the raw files before loading anything.


Two questions this answers:
  1. What does the alignment metadata say?
  2. Do the IFNa and PBS samples use the same gene list? The features files
     differ in size (161688 vs 148633 bytes), which suggests they do not.
     If they differ, genes missing from one reference would silently become
     zeros on merge - indistinguishable from "not expressed".
"""

import pandas as pd
from pathlib import Path

RAW = Path("data/raw")
SAMPLES = ["IFNa1", "IFNa2", "PBS1", "PBS2"]

print("\n=== samples_alignment.csv ===")
print(pd.read_csv(RAW / "GSM9126968_samples_alignment.csv.gz").to_string())

features = {}
for s in SAMPLES:
    features[s] = pd.read_csv(
        RAW / f"GSM9126968_{s}_features.tsv.gz",
        sep="\t", header=None,
    )

print("\n=== features per sample ===")
for s in SAMPLES:
    print(f"{s}: {features[s].shape[0]} features, {features[s].shape[1]} columns")

print("\n=== first 3 rows, IFNa1 ===")
print(features["IFNa1"].head(3).to_string())
print("\n=== first 3 rows, PBS1 ===")
print(features["PBS1"].head(3).to_string())

# column 1 = ensembl id, column 2 = gene symbol (0-indexed: 0 and 1)
ifna = set(features["IFNa1"][1])
pbs = set(features["PBS1"][1])

print("\n=== feature set comparison (gene symbols) ===")
print(f"In IFNa only: {len(ifna - pbs)}")
print(f"In PBS only : {len(pbs - ifna)}")
print(f"Shared      : {len(ifna & pbs)}")

print("\nExamples in IFNa only:")
print(sorted(ifna - pbs)[:25])
print("\nExamples in PBS only:")
print(sorted(pbs - ifna)[:25])

# Also check the third column - 10x uses it for feature type
# (Gene Expression, Antibody Capture, CRISPR Guide Capture...)
if features["IFNa1"].shape[1] >= 3:
    print("\n=== feature types ===")
    for s in SAMPLES:
        print(f"{s}: {features[s][2].value_counts().to_dict()}")

print("\nSTOP. Read the above before running 02_load.py.")
print("If the IFNa-only features are a transgene (GFP, PyMT, a viral")
print("sequence), that is a deliberate custom reference - decide")
print("consciously whether to keep or drop it.")