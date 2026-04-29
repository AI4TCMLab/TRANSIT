import os
import glob
import warnings
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, spearmanr
from scipy.spatial.distance import cosine


# =========================
# 1. 配置区
# =========================
# 药物表达谱目录（无表头：gene \t log2FC）
DRUG_DIR = r"G:\DOCTOR\1_心梗临界态课题记录\生信部分\4_药物预测算法代码及原始数据\图神经网络算法\调控网络\ITCM\DE_filtered"

# 疾病表达谱目录（有表头：Gene_Symbol \t logFC \t pv）
DISEASE_DIR = r"C:\Users\huihui1126\Desktop\药物预测方法测试\cmap\RRA测试\diseaselist-deg2function"

# 输出目录
OUT_DIR = r"G:\DOCTOR\1_心梗临界态课题记录\生信部分\4_药物预测算法代码及原始数据\图神经网络算法\调控网络\ITCM\Drug_Disease_Analysis_Output"

# 热图基因选择策略
STRICT_INTERSECTION = True
MIN_PRESENCE_RATIO = 0.3

# 相关性分析时，药物和疾病最少共同基因数
MIN_COMMON_GENES = 10

# 热图值截断
CLIP_VALUE = 3.0

# 是否画药物总热图
MAKE_DRUG_HEATMAP = True

# 组合图：自定义纳入绘制的药物
#CUSTOM_COMBO_DRUGS = ["S14S5", "S3S11", "S12S5", "S15S5", "S15S27", "S3S12","S3S6","S2S26","S12S29","S2S29",
#                      "S10S27", "S8S24", "S17S8", "S10S4", "S10S26", "S14S25","S10S28","S16S17","S17S5","S10S30"]
CUSTOM_COMBO_DRUGS = ["S14S5", "S15S26", "S2S26", "S6S29","S3S18","S6S26","S9S26", "S15S12", "S5S24","S10S26","S22", "S14S25"]


# 组合图中需要特别高亮的药物
HIGHLIGHT_DRUGS = ["S2S26", "S6S29"]

# 组合图排序方式
# "distance" 按距离从小到大
# "input"    按 CUSTOM_COMBO_DRUGS 里给的顺序
COMBO_ORDER_MODE = "distance"

# 散点图中标注距离最小的前多少个药物
SCATTER_LABEL_TOP_N = 15

# 共同基因热图显示多少个差异最大的基因
COMMON_HEATMAP_TOP_GENES = 50


# =========================
# 2. 通用读表函数
# =========================

def read_table_auto(path, sep="\t", header="infer", names=None):
    encodings = ["utf-8", "utf-8-sig", "gb18030", "gbk", "utf-16", "latin1"]
    last_err = None
    for enc in encodings:
        try:
            df = pd.read_csv(path, sep=sep, header=header, names=names, encoding=enc)
            return df
        except Exception as e:
            last_err = e
    raise last_err


# =========================
# 3. 读取药物表达谱
# =========================

def load_drug_signatures(drug_dir):
    """
    每个药物文件:
    gene \t log2FC
    无表头
    """
    drug_dict = {}
    files = glob.glob(os.path.join(drug_dir, "*.txt"))

    if not files:
        raise FileNotFoundError(f"药物目录下没有找到 txt 文件: {drug_dir}")

    for fp in files:
        drug_name = os.path.splitext(os.path.basename(fp))[0]

        df = read_table_auto(fp, sep="\t", header=None, names=["gene", "log2FC"])

        if df.shape[1] < 2:
            print(f"[跳过] {drug_name}: 列数不足 2")
            continue

        df = df.iloc[:, :2].copy()
        df.columns = ["gene", "log2FC"]

        df["gene"] = df["gene"].astype(str).str.strip()
        df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
        df = df.dropna(subset=["gene", "log2FC"])

        # 同一基因重复，保留绝对值最大的
        df["abs_fc"] = df["log2FC"].abs()
        df = df.sort_values("abs_fc", ascending=False).drop_duplicates("gene")
        df = df.drop(columns=["abs_fc"])

        if df.empty:
            print(f"[跳过] {drug_name}: 无有效数据")
            continue

        drug_dict[drug_name] = df.set_index("gene")["log2FC"]

    if not drug_dict:
        raise ValueError("没有成功读取任何药物表达谱")

    print(f"[OK] 成功读取药物表达谱: {len(drug_dict)} 个")
    return drug_dict


# =========================
# 4. 读取疾病表达谱
# =========================

def detect_disease_columns(df):
    cols = [c.strip() for c in df.columns]

    gene_candidates = ["Gene_Symbol", "gene", "Gene", "SYMBOL", "symbol"]
    fc_candidates = ["logFC", "log2FC", "logfc", "log2fc", "FC"]

    gene_col = None
    fc_col = None

    for c in gene_candidates:
        if c in cols:
            gene_col = c
            break

    for c in fc_candidates:
        if c in cols:
            fc_col = c
            break

    if gene_col is None:
        for c in cols:
            cl = c.lower()
            if "gene" in cl or "symbol" in cl:
                gene_col = c
                break

    if fc_col is None:
        for c in cols:
            cl = c.lower()
            if "logfc" in cl or "log2fc" in cl:
                fc_col = c
                break

    return gene_col, fc_col


def load_disease_signatures(disease_dir):
    """
    疾病文件有表头，典型列名:
    Gene_Symbol, logFC, pv
    """
    disease_dict = {}
    files = glob.glob(os.path.join(disease_dir, "*.txt"))

    if not files:
        raise FileNotFoundError(f"疾病目录下没有找到 txt 文件: {disease_dir}")

    for fp in files:
        disease_name = os.path.splitext(os.path.basename(fp))[0]

        df = read_table_auto(fp, sep="\t", header=0)
        df.columns = df.columns.astype(str).str.strip()

        gene_col, fc_col = detect_disease_columns(df)

        if gene_col is None or fc_col is None:
            print(f"[跳过] {disease_name}: 无法识别基因列/logFC列，列名为 {list(df.columns)}")
            continue

        tmp = df[[gene_col, fc_col]].copy()
        tmp.columns = ["gene", "logFC"]

        tmp["gene"] = tmp["gene"].astype(str).str.strip()
        tmp["logFC"] = pd.to_numeric(tmp["logFC"], errors="coerce")
        tmp = tmp.dropna(subset=["gene", "logFC"])

        # 同一基因重复，保留绝对值最大的
        tmp["abs_fc"] = tmp["logFC"].abs()
        tmp = tmp.sort_values("abs_fc", ascending=False).drop_duplicates("gene")
        tmp = tmp.drop(columns=["abs_fc"])

        if tmp.empty:
            print(f"[跳过] {disease_name}: 无有效数据")
            continue

        disease_dict[disease_name] = tmp.set_index("gene")["logFC"]

    if not disease_dict:
        raise ValueError("没有成功读取任何疾病表达谱")

    print(f"[OK] 成功读取疾病表达谱: {len(disease_dict)} 个")
    return disease_dict


# =========================
# 5. 药物热图
# =========================

def build_drug_matrix(drug_dict):
    mat = pd.concat(drug_dict, axis=1)
    mat.columns = mat.columns.get_level_values(0)
    return mat


def choose_heatmap_matrix(drug_matrix, strict_intersection=True, min_presence_ratio=0.3):
    n_drugs = drug_matrix.shape[1]
    strict_mat = drug_matrix.dropna(axis=0, how="any")

    if strict_intersection and strict_mat.shape[0] >= 2:
        print(f"[OK] 热图使用严格交集基因数: {strict_mat.shape[0]}")
        return strict_mat

    min_presence = max(2, int(np.ceil(n_drugs * min_presence_ratio)))
    keep = drug_matrix.notna().sum(axis=1) >= min_presence
    fallback_mat = drug_matrix.loc[keep].copy().fillna(0)

    print(f"[OK] 严格交集过少，热图回退为至少出现在 {min_presence}/{n_drugs} 个药物中的基因")
    print(f"[OK] 热图使用基因数: {fallback_mat.shape[0]}")
    return fallback_mat


def plot_drug_heatmap(mat, out_png):
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        print("[跳过] 热图矩阵太小，无法作图")
        return

    mat_plot = mat.clip(-CLIP_VALUE, CLIP_VALUE)

    sns.set(style="white")
    g = sns.clustermap(
        mat_plot,
        cmap="RdBu_r",
        center=0,
        figsize=(14, 16),
        row_cluster=True,
        col_cluster=True,
        xticklabels=True,
        yticklabels=False
    )
    g.fig.suptitle("Drug Signature Heatmap", y=1.02)
    g.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"[OK] 热图已保存: {out_png}")


# =========================
# 6. 核心指标计算
# =========================

def compute_one_disease_vs_all_drugs(disease_name, disease_ser, drug_dict, min_common_genes=10):
    """
    核心输出:
    - Pearson_r
    - Spearman_r
    - Cosine_similarity_common_genes_raw
    - Euclidean_distance_common_genes_raw

    Euclidean_distance_common_genes_raw:
    对 disease 和 drug 的共同基因取原始表达值，直接算欧氏距离
    """
    records = []

    for drug_name, drug_ser in drug_dict.items():
        common = disease_ser.index.intersection(drug_ser.index)

        if len(common) < min_common_genes:
            continue

        x = disease_ser.loc[common].astype(float).values
        y = drug_ser.loc[common].astype(float).values

        if np.std(x) == 0 or np.std(y) == 0:
            pear_r, pear_p = np.nan, np.nan
            spear_r, spear_p = np.nan, np.nan
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pear_r, pear_p = pearsonr(x, y)
                spear_r, spear_p = spearmanr(x, y)

        euc_raw = np.linalg.norm(x - y)

        try:
            cos_sim_raw = 1 - cosine(x, y)
        except Exception:
            cos_sim_raw = np.nan

        records.append({
            "Disease": disease_name,
            "Drug": drug_name,
            "n_common_genes": len(common),
            "Pearson_r": pear_r,
            "Pearson_p": pear_p,
            "Spearman_r": spear_r,
            "Spearman_p": spear_p,
            "Cosine_similarity_common_genes_raw": cos_sim_raw,
            "Euclidean_distance_common_genes_raw": euc_raw
        })

    if not records:
        return pd.DataFrame(columns=[
            "Disease", "Drug", "n_common_genes",
            "Pearson_r", "Pearson_p",
            "Spearman_r", "Spearman_p",
            "Cosine_similarity_common_genes_raw",
            "Euclidean_distance_common_genes_raw"
        ])

    res = pd.DataFrame(records)

    # 主排序仍按负相关优先
    res = res.sort_values(
        by=["Pearson_r", "Spearman_r", "Euclidean_distance_common_genes_raw"],
        ascending=[True, True, True],
        na_position="last"
    ).reset_index(drop=True)

    return res


# =========================
# 7. 组合图：自定义药物集合
# =========================

def plot_distance_combo_figure(
    score_df: pd.DataFrame,
    disease_name: str,
    out_prefix: str,
    selected_drugs: list = None,
    highlight_drugs: list = None,
    distance_col: str = "Euclidean_distance_common_genes_raw",
    order_mode: str = "distance"
):
    """
    左侧 = 疾病为中心的径向距离示意图
    右侧 = 横向条形图

    只绘制 selected_drugs 中指定的药物
    """
    if selected_drugs is None or len(selected_drugs) == 0:
        print(f"[跳过] {disease_name}: selected_drugs 为空")
        return

    if highlight_drugs is None:
        highlight_drugs = []

    if score_df is None or score_df.empty:
        print(f"[跳过] {disease_name}: score_df 为空")
        return

    required_cols = {"Drug", distance_col}
    if not required_cols.issubset(score_df.columns):
        raise ValueError(f"score_df 缺少必要列: {required_cols}")

    plot_df = score_df[score_df["Drug"].isin(selected_drugs)].copy()

    found_drugs = set(plot_df["Drug"].tolist())
    missing_drugs = [x for x in selected_drugs if x not in found_drugs]
    if missing_drugs:
        print(f"[提醒] {disease_name}: 以下药物未出现在结果表中，已跳过: {missing_drugs}")

    plot_df = plot_df[["Drug", distance_col]].dropna().copy()

    if plot_df.empty:
        print(f"[跳过] {disease_name}: 自定义药物在结果表中无可绘制数据")
        return

    if order_mode == "distance":
        plot_df = plot_df.sort_values(distance_col, ascending=True).reset_index(drop=True)
    elif order_mode == "input":
        plot_df["Drug"] = pd.Categorical(plot_df["Drug"], categories=selected_drugs, ordered=True)
        plot_df = plot_df.sort_values("Drug").reset_index(drop=True)
    else:
        raise ValueError("order_mode 只支持 'distance' 或 'input'")

    # 距离归一化到 [0.25, 1.0]
    dmin = plot_df[distance_col].min()
    dmax = plot_df[distance_col].max()

    if dmax == dmin:
        plot_df["radius"] = 0.6
    else:
        plot_df["radius"] = 0.25 + 0.75 * (plot_df[distance_col] - dmin) / (dmax - dmin)

    angles = np.linspace(20, 340, len(plot_df), endpoint=False)
    angles = np.deg2rad(angles)
    plot_df["angle"] = angles

    plot_df["x"] = plot_df["radius"] * np.cos(plot_df["angle"])
    plot_df["y"] = plot_df["radius"] * np.sin(plot_df["angle"])
    plot_df["is_highlight"] = plot_df["Drug"].isin(highlight_drugs)

    fig = plt.figure(figsize=(15, max(6, len(plot_df) * 0.45)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], wspace=0.25)

    # 左：径向距离图
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_aspect("equal")

    ax1.scatter(0, 0, s=320, marker="o", edgecolors="black", linewidths=1.5, zorder=5)
    ax1.text(0, 0, disease_name, ha="center", va="center", fontsize=10, fontweight="bold", zorder=6)

    for r in [0.25, 0.5, 0.75, 1.0]:
        circle = plt.Circle((0, 0), r, fill=False, linestyle="--", linewidth=0.6, alpha=0.3)
        ax1.add_patch(circle)

    for _, row in plot_df.iterrows():
        x, y = row["x"], row["y"]
        drug = row["Drug"]
        is_hi = row["is_highlight"]

        ax1.plot([0, x], [0, y], linestyle="--", linewidth=1.1)

        if is_hi:
            ax1.scatter(x, y, s=220, facecolors="none", edgecolors="red", linewidths=2.0, zorder=4)
        else:
            ax1.scatter(x, y, s=160, edgecolors="black", linewidths=0.8, zorder=4)

        ax1.text(
            x * 1.08, y * 1.08,
            drug if not is_hi else f"[H] {drug}",
            fontsize=9,
            ha="center", va="center"
        )

        midx, midy = x * 0.56, y * 0.56
        ax1.text(
            midx, midy,
            f"δ({disease_name},{drug})",
            fontsize=8,
            rotation=np.rad2deg(np.arctan2(y, x)),
            ha="center", va="center"
        )

    ax1.set_title(f"{disease_name}: radial distance map", fontsize=12)
    ax1.set_xlim(-1.28, 1.28)
    ax1.set_ylim(-1.28, 1.28)
    ax1.axis("off")

    # 右：横向条形图
    ax2 = fig.add_subplot(gs[0, 1])

    bar_df = plot_df.iloc[::-1].copy()

    bars = ax2.barh(
        y=bar_df["Drug"],
        width=bar_df[distance_col]
    )

    for i, drug in enumerate(bar_df["Drug"]):
        if drug in highlight_drugs:
            bars[i].set_edgecolor("red")
            bars[i].set_linewidth(2)

    ax2.set_xlabel(distance_col)
    ax2.set_ylabel("Drug")
    ax2.set_title(f"{disease_name}: selected drugs", fontsize=12)

    for i, (_, row) in enumerate(bar_df.iterrows()):
        ax2.text(
            row[distance_col],
            i,
            f" {row[distance_col]:.3f}",
            va="center",
            fontsize=8
        )

    plt.suptitle(
        f"Distance visualization based on {distance_col}",
        fontsize=14,
        y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    png_path = f"{out_prefix}.png"
    pdf_path = f"{out_prefix}.pdf"

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"[OK] 组合图已保存: {png_path}")
    print(f"[OK] 组合图PDF已保存: {pdf_path}")


# =========================
# 8. 距离 vs 相关性散点图
# =========================

def plot_distance_vs_correlation(
    score_df: pd.DataFrame,
    disease_name: str,
    out_prefix: str,
    highlight_drugs: list = None,
    label_top_n: int = 15
):
    if highlight_drugs is None:
        highlight_drugs = []

    if score_df is None or score_df.empty:
        print(f"[跳过] {disease_name}: score_df 为空")
        return

    plot_df = score_df.copy()
    plot_df["is_highlight"] = plot_df["Drug"].isin(highlight_drugs)

    plt.figure(figsize=(8, 7))
    plt.scatter(
        plot_df["Pearson_r"],
        plot_df["Euclidean_distance_common_genes_raw"],
        alpha=0.7
    )

    hi_df = plot_df[plot_df["is_highlight"]]
    if not hi_df.empty:
        plt.scatter(
            hi_df["Pearson_r"],
            hi_df["Euclidean_distance_common_genes_raw"],
            s=120,
            facecolors="none",
            edgecolors="red",
            linewidths=1.8
        )

    label_df = plot_df.nsmallest(label_top_n, "Euclidean_distance_common_genes_raw").copy()
    label_names = set(label_df["Drug"].tolist()) | set(highlight_drugs)

    for _, row in plot_df.iterrows():
        if row["Drug"] in label_names:
            txt = row["Drug"]
            if row["Drug"] in highlight_drugs:
                txt = f"[H] {txt}"
            plt.text(
                row["Pearson_r"],
                row["Euclidean_distance_common_genes_raw"],
                txt,
                fontsize=8
            )

    plt.xlabel("Pearson correlation")
    plt.ylabel("Euclidean distance on common genes (raw)")
    plt.title(f"{disease_name}: correlation vs raw common-gene distance")
    plt.tight_layout()

    png_path = f"{out_prefix}.png"
    pdf_path = f"{out_prefix}.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"[OK] 距离-相关性散点图已保存: {png_path}")
    print(f"[OK] 距离-相关性散点图PDF已保存: {pdf_path}")


# =========================
# 9. 单药物共同基因热图
# =========================

def plot_common_genes_heatmap_for_one_drug(
    disease_name: str,
    disease_ser: pd.Series,
    drug_name: str,
    drug_ser: pd.Series,
    out_prefix: str,
    top_n_genes: int = 50
):
    common = disease_ser.index.intersection(drug_ser.index)

    if len(common) < 2:
        print(f"[跳过] {disease_name} vs {drug_name}: 共同基因太少")
        return

    df = pd.DataFrame({
        disease_name: disease_ser.loc[common],
        drug_name: drug_ser.loc[common]
    })

    # 按两者差异绝对值排序，取最能贡献欧氏距离的前 top_n_genes
    df["abs_diff"] = (df[disease_name] - df[drug_name]).abs()
    df = df.sort_values("abs_diff", ascending=False).head(top_n_genes)
    df = df.drop(columns="abs_diff")

    plt.figure(figsize=(6, max(8, top_n_genes * 0.18)))
    sns.heatmap(
        df,
        cmap="RdBu_r",
        center=0,
        yticklabels=True
    )
    plt.title(f"{disease_name} vs {drug_name}\nCommon genes raw values")
    plt.tight_layout()

    png_path = f"{out_prefix}.png"
    pdf_path = f"{out_prefix}.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"[OK] 共同基因热图已保存: {png_path}")
    print(f"[OK] 共同基因热图PDF已保存: {pdf_path}")


# =========================
# 10. 主流程
# =========================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    score_dir = os.path.join(OUT_DIR, "Correlation_by_disease")
    combo_dir = os.path.join(OUT_DIR, "Distance_combo_plots")
    scatter_dir = os.path.join(OUT_DIR, "Distance_scatter_plots")
    heatmap_dir = os.path.join(OUT_DIR, "Common_gene_heatmaps")

    os.makedirs(score_dir, exist_ok=True)
    os.makedirs(combo_dir, exist_ok=True)
    os.makedirs(scatter_dir, exist_ok=True)
    os.makedirs(heatmap_dir, exist_ok=True)

    drug_dict = load_drug_signatures(DRUG_DIR)
    disease_dict = load_disease_signatures(DISEASE_DIR)

    if MAKE_DRUG_HEATMAP:
        drug_matrix = build_drug_matrix(drug_dict)
        heatmap_matrix = choose_heatmap_matrix(
            drug_matrix,
            strict_intersection=STRICT_INTERSECTION,
            min_presence_ratio=MIN_PRESENCE_RATIO
        )
        heatmap_path = os.path.join(OUT_DIR, "Drug_Cluster_Heatmap.png")
        plot_drug_heatmap(heatmap_matrix, heatmap_path)

    all_score_tables = []

    for disease_name, disease_ser in disease_dict.items():
        print(f"\n[分析疾病] {disease_name}")

        # 核心结果表
        score_df = compute_one_disease_vs_all_drugs(
            disease_name=disease_name,
            disease_ser=disease_ser,
            drug_dict=drug_dict,
            min_common_genes=MIN_COMMON_GENES
        )

        out_score = os.path.join(score_dir, f"{disease_name}_drug_correlation.tsv")
        score_df.to_csv(out_score, sep="\t", index=False)
        print(f"[OK] 相关性结果已保存: {out_score} (共 {score_df.shape[0]} 个药物)")

        all_score_tables.append(score_df)

        # 组合图：只绘制自定义药物集合
        combo_prefix = os.path.join(combo_dir, f"{disease_name}_distance_combo")
        plot_distance_combo_figure(
            score_df=score_df,
            disease_name=disease_name,
            out_prefix=combo_prefix,
            selected_drugs=CUSTOM_COMBO_DRUGS,
            highlight_drugs=HIGHLIGHT_DRUGS,
            distance_col="Euclidean_distance_common_genes_raw",
            order_mode=COMBO_ORDER_MODE
        )

        # 距离 vs 相关性散点图
        scatter_prefix = os.path.join(scatter_dir, f"{disease_name}_distance_vs_correlation")
        plot_distance_vs_correlation(
            score_df=score_df,
            disease_name=disease_name,
            out_prefix=scatter_prefix,
            highlight_drugs=HIGHLIGHT_DRUGS,
            label_top_n=SCATTER_LABEL_TOP_N
        )

        # 感兴趣药物共同基因热图
        for drug_name in HIGHLIGHT_DRUGS:
            if drug_name not in drug_dict:
                print(f"[提醒] {disease_name}: 高亮药物 {drug_name} 不在 drug_dict 中，跳过共同基因热图")
                continue

            one_heatmap_prefix = os.path.join(heatmap_dir, f"{disease_name}_vs_{drug_name}")
            plot_common_genes_heatmap_for_one_drug(
                disease_name=disease_name,
                disease_ser=disease_ser,
                drug_name=drug_name,
                drug_ser=drug_dict[drug_name],
                out_prefix=one_heatmap_prefix,
                top_n_genes=COMMON_HEATMAP_TOP_GENES
            )

    if all_score_tables:
        summary = pd.concat(all_score_tables, ignore_index=True)
        summary_path = os.path.join(OUT_DIR, "All_Disease_Drug_Correlation_Summary.tsv")
        summary.to_csv(summary_path, sep="\t", index=False)
        print(f"\n[OK] 总汇总表已保存: {summary_path}")

    print("\n全部完成。")


if __name__ == "__main__":
    main()