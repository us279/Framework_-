import os
import pandas as pd

# --- User parameters ---
folder1 = "/Users/uday03/Desktop/Framework_repository/Framework_-/Building Framework/Kahan Summation/errors"
folder2 = "/Users/uday03/Desktop/Framework_repository/Framework_-/Building Framework/No preconditioning/errors"
output_folder = "/Users/uday03/Desktop/Framework_repository/Framework_-"
os.makedirs(output_folder, exist_ok=True)

# Base names of your files (without “.csv”)
file_bases = ["advection_errors", "heat_errors", "burgers_errors"]

# Columns to compare
error_cols = ["L1 Error", "L2 Error", "Linf Error"]

for base in file_bases:
    f2 = os.path.join(folder1,  base + ".csv")
    f1 = os.path.join(folder2,  base + ".csv")

    # Load
    df1 = pd.read_csv(f1)
    df2 = pd.read_csv(f2)

    # Sanity check: same precisions in same order
    if not df1["precision"].equals(df2["precision"]):
        raise ValueError(f"Precision mismatch in {base}")

    # Compute percent difference: (new - old) / old * 100
    diffs = (df2[error_cols] - df1[error_cols]) / df1[error_cols] * 100

    # Build output DF
    out = pd.concat([df1[["precision"]], diffs], axis=1)
    out.columns = ["precision"] + [c + " (%)" for c in error_cols]

    # Write
    out_path = os.path.join(output_folder, base + "_percent_diff.csv")
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
