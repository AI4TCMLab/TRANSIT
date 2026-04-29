import pandas as pd
import gseapy as gp
import sys
import os
import numpy as np

Lfile = sys.argv[1]   # lincs drug profile
Sfile = sys.argv[2]   # disease signature

# ========== 1. 读取疾病 signature ==========
S_df = pd.read_csv(Sfile, sep='\t')

# 检查必要列
required_cols = ['Gene_Symbol', 'logFC', 'pv']
for col in required_cols:
    if col not in S_df.columns:
        raise ValueError(f"Sfile 缺少必要列: {col}")

# 清洗
S_df['Gene_Symbol'] = S_df['Gene_Symbol'].astype(str).str.strip()
S_df = S_df[S_df['Gene_Symbol'] != '']
S_df['logFC'] = pd.to_numeric(S_df['logFC'], errors='coerce')
S_df['pv'] = pd.to_numeric(S_df['pv'], errors='coerce')
S_df = S_df.dropna(subset=['logFC', 'pv'])

# 筛选上下调
Sup_df = S_df[(S_df['logFC'] > 0) & (S_df['pv'] < 0.05)].copy()
Sdown_df = S_df[(S_df['logFC'] < 0) & (S_df['pv'] < 0.05)].copy()

# 去重，防止重复基因
Sup_genes = Sup_df['Gene_Symbol'].drop_duplicates().tolist()
Sdown_genes = Sdown_df['Gene_Symbol'].drop_duplicates().tolist()

# ========== 2. 构建合法 GMT ==========
gmt_lines = []

# GMT 标准格式: set_name \t description \t gene1 \t gene2 ...
# 只写非空集合，避免空 gene set 导致 prerank 崩溃
if len(Sdown_genes) > 0:
    gmt_lines.append('\t'.join(['down', 'na'] + Sdown_genes))

if len(Sup_genes) > 0:
    gmt_lines.append('\t'.join(['up', 'na'] + Sup_genes))

# 如果上下调都没有，直接输出 0
if len(gmt_lines) == 0:
    with open('out_lincs.txt', 'w') as file:
        file.write('0\n0')
    sys.exit(0)

gmt_file = 'gene_sets-defined.gmt'
with open(gmt_file, 'w', encoding='utf-8') as file:
    file.write('\n'.join(gmt_lines) + '\n')

# ========== 3. 读取药物 profile，整理成标准 prerank 输入 ==========
L_df = pd.read_csv(Lfile, header=None, sep='\t')

if L_df.shape[1] < 2:
    raise ValueError("Lfile 至少应有两列：gene 和 score")

# 只取前两列，明确为 gene + score
L_df = L_df.iloc[:, :2].copy()
L_df.columns = ['gene', 'score']

# 清洗
L_df['gene'] = L_df['gene'].astype(str).str.strip()
L_df = L_df[L_df['gene'] != '']
L_df['score'] = pd.to_numeric(L_df['score'], errors='coerce')
L_df = L_df.dropna(subset=['gene', 'score'])
L_df = L_df[~np.isinf(L_df['score'])]

# 去重（保留第一次出现）
L_df = L_df.drop_duplicates(subset='gene', keep='first')

# 按 score 排序
# 你原来是 ascending=True，这可能是你的设计（把更负的药物基因放前面）
# 如果你希望标准“高分在前”，改成 ascending=False
L_df_sorted = L_df.sort_values(by='score', ascending=True).reset_index(drop=True)

# 可选：检查与 disease gene 的交集；若完全没交集，直接给 0
all_sig_genes = set(Sup_genes) | set(Sdown_genes)
overlap = all_sig_genes & set(L_df_sorted['gene'])

if len(overlap) == 0:
    with open('out_lincs.txt', 'w') as file:
        file.write('0\n0')
    sys.exit(0)

# ========== 4. 跑 prerank ==========
# min_size 不建议用 0，至少设为 1 更稳妥
pre_res = gp.prerank(
    rnk=L_df_sorted,
    gene_sets=gmt_file,
    threads=4,
    min_size=1,
    max_size=7000,
    permutation_num=50,
    outdir=None,      # 调试时先别写磁盘
    seed=6,
    verbose=True
)

# ========== 5. 读取结果 ==========
ES_down = pre_res.results.get('down', {}).get('es', 0)
NES_down = pre_res.results.get('down', {}).get('nes', 0)

ES_up = pre_res.results.get('up', {}).get('es', 0)
NES_up = pre_res.results.get('up', {}).get('nes', 0)

# 综合得分
SSlincs = (ES_up - ES_down) / 2
SSlincs_Norm = (NES_up - NES_down) / 2

# ========== 6. 写出 ==========
with open('out_lincs.txt', 'w') as file:
    file.write(str(SSlincs))
    file.write('\n')
    file.write(str(SSlincs_Norm))