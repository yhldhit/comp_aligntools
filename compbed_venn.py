#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare two BED files and visualize shared/different sites using a Venn diagram.

Example:
    python compare_bed_venn.py A.bed B.bed output.png
"""

import sys
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import gzip
import numpy as np

# bed1 = "/home/lijia/yanhongliang/dataset/humb_meth/bed/SRR13832659.dedup.bed"
# bed2 = "/home/lijia/yanhongliang/dataset/humb_meth2/bed/SRR13832659.dedup.bed"
bed1 = "/home/lijia/yanhongliang/proj/comp_aligntools/parabrick_bed/bed/SRR15064000.dedup.bed"
bed2 = "/home/lijia/yanhongliang/dataset/humbv9/bed_bak/SRR15064000.final.bed"
out_png = "/home/lijia/yanhongliang/proj/comp_aligntools/ven/SRR15064000_comp_bed.png"

# ======================
# 2. 读取 BED 文件
# ======================
def read_bed(path):
    df = pd.read_csv(path, sep='\t', header=None, usecols=[0,1,2,3,4],
                     names=['chr','start','end','meth','count'])
    # 创建唯一位点ID
    df['pos'] = df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str)
    df['state'] = np.where(df['meth'] / df['count'] > 0.5, 1, 0)
    return df

df1 = read_bed(bed1)
df2 = read_bed(bed2)

# ======================
# 3. 构建集合
# ======================
set1 = set(df1['pos'])
set2 = set(df2['pos'])

# 交集
common = set1 & set2
only1 = set1 - set2
only2 = set2 - set1

# ======================
# 4. 判断状态是否相同
# ======================
df_common1 = df1[df1['pos'].isin(common)][['pos','state']]
df_common2 = df2[df2['pos'].isin(common)][['pos','state']]
merged = df_common1.merge(df_common2, on='pos', suffixes=('_A','_B'))

same_state = set(merged.loc[merged['state_A'] == merged['state_B'], 'pos'])
diff_state = set(merged.loc[merged['state_A'] != merged['state_B'], 'pos'])

# ======================
# 5. 打印统计结果
# ======================
print("=== Comparison summary ===")
print(f"Total in A: {len(set1)}")
print(f"Total in B: {len(set2)}")
print(f"Common sites: {len(common)}")
print(f"  ├─ Same state: {len(same_state)}")
print(f"  └─ Different state: {len(diff_state)}")
print(f"Unique to A: {len(only1)}")
print(f"Unique to B: {len(only2)}")

# ======================
# 6. 绘制韦恩图
# ======================
plt.figure(figsize=(6,6))
v = venn2(subsets=(len(only1), len(only2), len(common)), set_labels=('parabrick', 'bismark'))

# 在交集区域额外标注状态差异
v.get_label_by_id('11').set_text(f"Common:\nSame={len(same_state)}\nDiff={len(diff_state)}")

plt.title("BED site comparison")
plt.tight_layout()
plt.savefig(out_png, dpi=300)
print(f"[Saved] Venn diagram → {out_png}")


