#!/usr/bin/env python3
"""
Merge all *_read_distribution.txt files into one table keyed by sample_id.

Each input file looks like:

Total Reads        702622
Total Tags         11173967
Total Assigned Tags11062524
-------------------------------------------------------------------------
Group               Total_bases  Tag_count  Tags/Kb
CDS_Exons           40127211     713208     177.74
5'UTR_Exons         44767799     856687     19.14
...

Output: a CSV with one row per sample_id and columns like:
   sample_id,
   CDS_Exons_Total_bases,
   CDS_Exons_Tag_count,
   CDS_Exons_Tags/Kb,
   5'UTR_Exons_Total_bases,
   ...
"""

import os
import glob
import pandas as pd

def merge_read_distributions(root_dir, pattern="**/*_read_distribution.txt"):
    # 1) 找到所有文件
    paths = glob.glob(os.path.join(root_dir, pattern), recursive=True)
    all_tables = []
    for fp in paths:
        sample_id = os.path.basename(fp).replace("_read_distribution.txt","")
        
        # 2) 读取文件内容
        with open(fp) as f:
            lines = f.readlines()
        
        # 3) 找到以 "Group" 开头的表头行号
        header_idx = None
        for i, line in enumerate(lines):
            if line.strip().startswith("Group"):
                header_idx = i
                break
        if header_idx is None:
            continue
        
        # 4) 用 pandas 读取从 header_idx 行开始的表格
        df = pd.read_csv(
            fp,
            delim_whitespace=True,
            header=header_idx
        )
        df["sample_id"] = sample_id
        all_tables.append(df)
    
    # 5) 合并所有小表（长格式）
    long = pd.concat(all_tables, ignore_index=True)
    
    # 6) 转宽格式：以 sample_id 为行，以 Group×Metric 为列
    wide = long.pivot_table(
        index="sample_id",
        columns="Group",
        values=["Total_bases","Tag_count","Tags/Kb"]
    )
    # 7) flatten MultiIndex
    wide.columns = [
        f"{grp}_{metric}" for metric, grp in wide.columns
    ]
    wide = wide.reset_index()
    return wide

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Merge RSeQC read_distribution files into one table"
    )
    parser.add_argument(
        "-i","--indir", required=True,
        help="Root directory to search for *_read_distribution.txt"
    )
    parser.add_argument(
        "-o","--out", default="read_distribution_summary.csv",
        help="Output CSV filename"
    )
    args = parser.parse_args()

    summary_df = merge_read_distributions(args.indir)
    summary_df.to_csv(args.out, index=False)
    print(f"Wrote merged table ({summary_df.shape[0]} samples × {summary_df.shape[1]-1} metrics) to {args.out}")
    