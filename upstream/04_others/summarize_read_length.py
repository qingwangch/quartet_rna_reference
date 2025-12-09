#!/usr/bin/env python3
"""
Summarize average read lengths from *_read_length.txt files
(each file has sample_id on line 1, and average_length on line 2).
"""

import os
import glob
import sys
import pandas as pd

def collect_read_lengths(directory, pattern="*_read_length.txt"):
    file_paths = glob.glob(os.path.join(directory, "**", pattern), recursive=True)
    if not file_paths:
        return pd.DataFrame()  # 空
    records = []
    for fp in file_paths:
        with open(fp) as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        if len(lines) < 2:
            continue
        sample_id = lines[0]
        try:
            avg_len = float(lines[1])
        except ValueError:
            avg_len = lines[1]
        records.append({"sample_id": sample_id,
                        "average_length": avg_len})
    return pd.DataFrame.from_records(records)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Summarize average read lengths from *_read_length.txt files"
    )
    parser.add_argument("-i","--indir", required=True,
                        help="Directory root to search for *_read_length.txt")
    parser.add_argument("-o","--out", default="read_length_summary.csv",
                        help="Output CSV file")
    args = parser.parse_args()

    df = collect_read_lengths(args.indir)
    if df.empty:
        print(f"No '*_read_length.txt' files found under {args.indir}", file=sys.stderr)
        sys.exit(1)

    # 检查必要列
    expected = {"sample_id", "average_length"}
    missing = expected - set(df.columns)
    if missing:
        print("Error: missing columns:", missing, file=sys.stderr)
        print("Found columns:", list(df.columns), file=sys.stderr)
        sys.exit(1)

    # 排序并写出
    df = df.sort_values("sample_id").reset_index(drop=True)
    df.to_csv(args.out, index=False)
    print(f"Wrote summary of {len(df)} samples to {args.out}")
    