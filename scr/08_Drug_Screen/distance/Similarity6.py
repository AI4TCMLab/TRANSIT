import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# 自动避让标签
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("[提示] 未安装 adjustText，标签将不会自动避让。建议运行: python -m pip install adjustText")


# =========================
# 1. 配置区
# =========================

# 药物表达谱目录
DRUG_DIR = r"G:\DOCTOR\1_心梗临界态课题记录\生信部分\4_药物预测算法代码及原始数据\图神经网络算法\调控网络\ITCM\DE_filtered"

# DNB-TARGET 疾病谱目录
DISEASE_DIR = r"C:\Users\huihui1126\Desktop\药物预测方法测试\cmap\RRA测试\diseaselist-deg2function"

# 排名文件路径
RANK_FILE = r"G:\DOCTOR\1_心梗临界态课题记录\生信部分\4_药物预测算法代码及原始数据\图神经网络算法\调控网络\ITCM\排名文件.xlsx"  # 请修改为实际路径
RANK_SHEET = "early_macro_clip_norm_p1"

# 输出目录
OUT_DIR = r"G:\DOCTOR\1_心梗临界态课题记录\生信部分\4_药物预测算法代码及原始数据\图神经网络算法\调控网络\ITCM\Drug_DNB_RawDisease_TotalN_Weighted_Output"

# 药物筛选阈值
DRUG_ABS_LOGFC_CUTOFF = 0
DRUG_PVALUE_CUTOFF = 1


# 感兴趣药物 - 用特殊标记突出显示
INTEREST_DRUGS = ["S14S8", "S3S11", "S15S21", "S18" ,"S5S25", "S15S25",  "S11S25"]

# 不感兴趣药物 - 从图中删除
UNINTEREST_DRUGS = []  # 可以在这里添加不感兴趣的药物

# 是否标注所有点
LABEL_ALL_POINTS = False

# 如果不标所有点，每个象限最多标几个
TOP_LABEL_PER_QUADRANT = 5

# 离群值处理参数
REMOVE_OUTLIERS = False
OUTLIER_METHOD = "percentile"
OUTLIER_PERCENTILE_LOW = 0.5
OUTLIER_PERCENTILE_HIGH = 99.5
OUTLIER_IQR_MULTIPLIER = 1.5

# 标签参数
LABEL_FONT_SIZE = 8
INTEREST_DRUG_MARKER_SIZE = 100
INTEREST_DRUG_MARKER_COLOR = "gold"
INTEREST_DRUG_MARKER_EDGE_COLOR = "black"
INTEREST_DRUG_MARKER_SHAPE = "s"

# 排名相关参数
USE_RANK_INFO = True  # 是否使用排名信息
RANK_SIZE_MIN = 1  # 最小排名对应的点大小
RANK_SIZE_MAX = 200  # 最大排名对应的点大小
RANK_COLOR = "gray"  # 普通药物的颜色（统一为灰色）
RANK_ALPHA = 0.7  # 普通药物的透明度


# =========================
# 2. 基础函数
# =========================

def read_table_auto(path, sep="\t", header="infer", names=None):
    encodings = ["utf-8", "utf-8-sig", "gb18030", "gbk", "utf-16", "latin1"]
    last_err = None
    for enc in encodings:
        try:
            return pd.read_csv(path, sep=sep, header=header, names=names, encoding=enc)
        except Exception as e:
            last_err = e
    raise last_err


def zscore_series_custom(s, method="numpy"):
    """
    自定义z-score计算函数
    """
    s = pd.to_numeric(s, errors="coerce")
    mean_val = s.mean(skipna=True)
    
    if method == "numpy":
        std_val = np.std(s.dropna(), ddof=1)  # 分母n-1
    elif method == "pandas":
        std_val = s.std(skipna=True)  # pandas默认ddof=1，分母n-1
    elif method == "excel":
        std_val = np.std(s.dropna(), ddof=0)  # 分母n，类似于Excel
    else:
        raise ValueError(f"未知的z-score计算方法: {method}")
    
    if pd.isna(std_val) or std_val == 0:
        return pd.Series(np.nan, index=s.index)
    
    z_scores = (s - mean_val) / std_val
    
    return z_scores


def setup_plot_style():
    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "Arial",
        "font.size": 9
    })


# =========================
# 3. 读取药物表达谱
# =========================

def load_drug_signatures(drug_dir, abs_fc_cutoff=0.5, pvalue_cutoff=0.05, exclude_drugs=None):
    """
    返回:
        drug_dict: {drug_name: pd.Series(index=gene, values=log2FC)}
    """
    drug_dict = {}
    files = glob.glob(os.path.join(drug_dir, "*.txt"))

    if not files:
        raise FileNotFoundError(f"药物目录下没有找到 txt 文件: {drug_dir}")
    
    if exclude_drugs is None:
        exclude_drugs = []

    kept_gene_counts = []
    count_with_p = 0
    count_without_p = 0
    bad_files = []
    excluded_count = 0

    for fp in files:
        drug_name = os.path.splitext(os.path.basename(fp))[0]
        
        # 检查是否在排除列表中
        if drug_name in exclude_drugs:
            excluded_count += 1
            continue

        try:
            df = read_table_auto(fp, sep="\t", header=0)
            df.columns = df.columns.astype(str).str.strip()

            # 形式1：有 gene/log2FC/pvalue
            if {"gene", "log2FC", "pvalue"}.issubset(set(df.columns)):
                tmp = df[["gene", "log2FC", "pvalue"]].copy()
                tmp["gene"] = tmp["gene"].astype(str).str.strip()
                tmp["log2FC"] = pd.to_numeric(tmp["log2FC"], errors="coerce")
                tmp["pvalue"] = pd.to_numeric(tmp["pvalue"], errors="coerce")
                tmp = tmp.dropna(subset=["gene", "log2FC", "pvalue"])

                tmp["abs_fc"] = tmp["log2FC"].abs()
                tmp = tmp.sort_values("abs_fc", ascending=False).drop_duplicates("gene")
                tmp = tmp[
                    (tmp["abs_fc"] >= abs_fc_cutoff) &
                    (tmp["pvalue"] < pvalue_cutoff)
                ].copy()

                if tmp.empty:
                    print(f"[跳过] {drug_name}: 筛选后无有效基因")
                    continue

                ser = tmp.set_index("gene")["log2FC"]
                ser.index = ser.index.astype(str).str.strip()
                drug_dict[drug_name] = ser

                kept_gene_counts.append(len(tmp))
                count_with_p += 1
                continue

            # 形式2：无表头 gene/log2FC
            df2 = read_table_auto(fp, sep="\t", header=None, names=["gene", "log2FC"])
            df2["gene"] = df2["gene"].astype(str).str.strip()
            df2["log2FC"] = pd.to_numeric(df2["log2FC"], errors="coerce")
            df2 = df2.dropna(subset=["gene", "log2FC"])

            df2["abs_fc"] = df2["log2FC"].abs()
            df2 = df2.sort_values("abs_fc", ascending=False).drop_duplicates("gene")
            df2 = df2[df2["abs_fc"] >= abs_fc_cutoff].copy()

            if df2.empty:
                print(f"[跳过] {drug_name}: 筛选后无有效基因")
                continue

            ser = df2.set_index("gene")["log2FC"]
            ser.index = ser.index.astype(str).str.strip()
            drug_dict[drug_name] = ser

            kept_gene_counts.append(len(df2))
            count_without_p += 1

        except Exception as e:
            print(f"[跳过] {drug_name}: 读取失败，错误: {e}")
            bad_files.append(drug_name)

    if not drug_dict:
        raise ValueError("没有成功读取任何药物表达谱")

    print(f"[OK] 成功读取药物表达谱: {len(drug_dict)} 个")
    print(f"[OK] 排除不感兴趣药物: {excluded_count} 个")
    print(f"[OK] 其中带 pvalue 的药物文件: {count_with_p} 个；仅有 gene+log2FC 的药物文件: {count_without_p} 个")

    if kept_gene_counts:
        print(
            f"[OK] 药物筛选后保留基因数: 中位数={int(np.median(kept_gene_counts))}, "
            f"最小={np.min(kept_gene_counts)}, 最大={np.max(kept_gene_counts)}"
        )

    if bad_files:
        print(f"[提示] 读取失败的药物文件数: {len(bad_files)}")

    return drug_dict


# =========================
# 4. 读取疾病表
# =========================

def load_disease_tables(disease_dir):
    disease_table_dict = {}
    files = glob.glob(os.path.join(disease_dir, "*.txt"))

    if not files:
        raise FileNotFoundError(f"疾病目录下没有找到 txt 文件: {disease_dir}")

    for fp in files:
        disease_name = os.path.splitext(os.path.basename(fp))[0]
        try:
            df = read_table_auto(fp, sep="\t", header=0)
            df.columns = df.columns.astype(str).str.strip()
            disease_table_dict[disease_name] = df.copy()
        except Exception as e:
            print(f"[跳过] {disease_name}: 读取失败，错误: {e}")

    if not disease_table_dict:
        raise ValueError("没有成功读取任何疾病表")

    print(f"[OK] 成功读取疾病表: {len(disease_table_dict)} 个")
    return disease_table_dict


# =========================
# 5. 读取排名信息
# =========================

def load_rank_info(rank_file, sheet_name):
    """
    读取排名信息文件
    
    返回:
        rank_dict: {drug_name: rank_value}
    """
    try:
        # 尝试读取Excel文件
        rank_df = pd.read_excel(rank_file, sheet_name=sheet_name)
        print(f"[OK] 成功读取排名文件: {rank_file}")
        print(f"[OK] 工作表: {sheet_name}")
        
        # 检查列名
        print(f"[INFO] 排名文件列名: {rank_df.columns.tolist()}")
        
        # 假设列名为"NO"和"Rank"，但可能有不同的列名
        if "NO" in rank_df.columns and "Rank" in rank_df.columns:
            rank_df = rank_df[["NO", "Rank"]]
        elif len(rank_df.columns) >= 2:
            # 如果列名不同，使用前两列
            rank_df = rank_df.iloc[:, :2]
            rank_df.columns = ["NO", "Rank"]
        else:
            raise ValueError("排名文件格式不正确，需要至少两列")
        
        # 清理数据
        rank_df["NO"] = rank_df["NO"].astype(str).str.strip()
        rank_df["Rank"] = pd.to_numeric(rank_df["Rank"], errors="coerce")
        rank_df = rank_df.dropna(subset=["NO", "Rank"])
        
        # 转换为字典
        rank_dict = dict(zip(rank_df["NO"], rank_df["Rank"]))
        
        print(f"[INFO] 成功读取 {len(rank_dict)} 个药物的排名信息")
        print(f"[INFO] 排名范围: 最小={rank_df['Rank'].min()}, 最大={rank_df['Rank'].max()}")
        
        # 显示前10个药物的排名
        print("[INFO] 前10个药物排名:")
        for i, (drug, rank) in enumerate(list(rank_dict.items())[:10]):
            print(f"  {i+1}. {drug}: {rank}")
        
        return rank_dict
        
    except Exception as e:
        print(f"[警告] 读取排名文件失败: {e}")
        print(f"[提示] 将不使用排名信息")
        return {}


# =========================
# 6. 核心：基于原始疾病表，逐药物计算
# =========================

def compute_one_drug_scores_from_raw_disease(drug_name, drug_ser, disease_df):
    """
    基于原始疾病表计算单个药物的GOOD和BAD DNB基因的均值。
    """
    required_cols = {"DNB_Symbol", "DNB_log2fc"}
    if not required_cols.issubset(disease_df.columns):
        raise ValueError(f"疾病表缺少必要列: {required_cols}")

    tmp = disease_df[["DNB_Symbol", "DNB_log2fc"]].copy()
    tmp["DNB_Symbol"] = tmp["DNB_Symbol"].astype(str).str.strip()
    tmp["DNB_log2fc"] = pd.to_numeric(tmp["DNB_log2fc"], errors="coerce")
    tmp = tmp.dropna(subset=["DNB_Symbol", "DNB_log2fc"])

    # 药物映射到疾病 DNB
    drug_map = drug_ser.copy()
    drug_map.index = drug_map.index.astype(str).str.strip()

    tmp["drug_log2FC"] = tmp["DNB_Symbol"].map(drug_map)
    tmp = tmp.dropna(subset=["drug_log2FC"]).copy()

    # 分离 GOOD 和 BAD DNB
    good_dnb = tmp[tmp["DNB_log2fc"] > 0]
    bad_dnb = tmp[tmp["DNB_log2fc"] < 0]

    # 计算 GOOD 和 BAD DNB 上药物的 mean_log2FC
    good_mean_log2fc = good_dnb["drug_log2FC"].mean() if not good_dnb.empty else np.nan
    bad_mean_log2fc = bad_dnb["drug_log2FC"].mean() if not bad_dnb.empty else np.nan

    return {
        "drug": drug_name,
        "good_mean_log2fc": good_mean_log2fc,
        "bad_mean_log2fc": bad_mean_log2fc,
        "n_good_dnb": good_dnb["DNB_Symbol"].nunique(),
        "n_bad_dnb": bad_dnb["DNB_Symbol"].nunique()
    }


def compute_all_drug_scores(drug_dict, disease_df):
    """计算所有药物的分数"""
    rows = []
    for drug_name, drug_ser in drug_dict.items():
        row = compute_one_drug_scores_from_raw_disease(drug_name, drug_ser, disease_df)
        rows.append(row)

    res = pd.DataFrame(rows)
    
    # 计算总有效DNB数
    res["n_total_valid_dnb"] = res["n_good_dnb"].fillna(0) + res["n_bad_dnb"].fillna(0)
    
    # 对所有药物的mean_log2fc进行z-score标准化
    res["good_log2fc_z"] = zscore_series_custom(res["good_mean_log2fc"], method="numpy")
    res["bad_log2fc_z"] = zscore_series_custom(res["bad_mean_log2fc"], method="numpy")
    
    return res


# =========================
# 7. 四象限分类（使用log2fc_z）
# =========================

def classify_quadrants(score_df):
    """使用good_log2fc_z和bad_log2fc_z进行四象限分类"""
    out = score_df.copy()

    # 基于 good_log2fc_z 和 bad_log2fc_z 划分四象限
    out["quadrant"] = np.nan
    out["effect_type"] = "Unclassified"

    # 注意：x轴是bad_log2fc_z，y轴是good_log2fc_z
    # Q1: x>0, y>0
    # Q2: x>0, y<0  (reversing: 逆转疾病状态)
    # Q3: x<0, y<0
    # Q4: x<0, y>0  (aggravating: 加重疾病状态)
    
    out.loc[(out["good_log2fc_z"] > 0) & (out["bad_log2fc_z"] > 0), "quadrant"] = "Q1"
    out.loc[(out["good_log2fc_z"] < 0) & (out["bad_log2fc_z"] > 0), "quadrant"] = "Q2"
    out.loc[(out["good_log2fc_z"] < 0) & (out["bad_log2fc_z"] < 0), "quadrant"] = "Q3"
    out.loc[(out["good_log2fc_z"] > 0) & (out["bad_log2fc_z"] < 0), "quadrant"] = "Q4"

    out.loc[out["quadrant"] == "Q2", "effect_type"] = "Reversing"
    out.loc[out["quadrant"] == "Q4", "effect_type"] = "Aggravating"

    return out


# =========================
# 8. 生成详细四象限表格
# =========================

def create_quadrant_detail_table(quadrant_df, rank_dict=None):
    """
    生成详细四象限表格，包含以下列：
    drug, good_log2fc, bad_log2fc, n_total_valid_dnb, n_good_dnb, n_bad_dnb, 
    good_log2fc_z, bad_log2fc_z, quadrant, rank
    """
    detail_df = quadrant_df.copy()
    
    # 重命名列，使其符合要求的格式
    detail_df = detail_df.rename(columns={
        "good_mean_log2fc": "good_log2fc",
        "bad_mean_log2fc": "bad_log2fc"
    })
    
    # 添加排名信息
    if rank_dict:
        detail_df["rank"] = detail_df["drug"].map(rank_dict)
    
    # 按象限排序，然后按药物名称排序
    quadrant_order = {"Q1": 1, "Q2": 2, "Q3": 3, "Q4": 4}
    detail_df["quadrant_order"] = detail_df["quadrant"].map(quadrant_order)
    detail_df = detail_df.sort_values(["quadrant_order", "drug"]).reset_index(drop=True)
    detail_df = detail_df.drop(columns=["quadrant_order"])
    
    # 选择并重新排列列
    columns_order = [
        "drug", "good_log2fc", "bad_log2fc", 
        "n_total_valid_dnb", "n_good_dnb", "n_bad_dnb",
        "good_log2fc_z", "bad_log2fc_z", "quadrant", "effect_type"
    ]
    
    # 如果有序位列，添加到列顺序中
    if "rank" in detail_df.columns:
        columns_order.append("rank")
    
    # 确保所有需要的列都存在
    for col in columns_order:
        if col not in detail_df.columns:
            detail_df[col] = np.nan
    
    detail_df = detail_df[columns_order]
    
    return detail_df


# =========================
# 9. 离群值检测与处理
# =========================

def detect_and_remove_outliers(df, columns=["good_log2fc_z", "bad_log2fc_z"], 
                             method="percentile", lower_percentile=0.5, upper_percentile=99.5, 
                             iqr_multiplier=1.5, interest_drugs=None):
    """
    检测并标记离群值
    """
    df_clean = df.copy()
    outlier_mask = pd.Series(False, index=df.index)
    outlier_info = {
        "method": method,
        "total_outliers": 0,
        "outliers_by_column": {},
        "outlier_details": []
    }
    
    for col in columns:
        if col not in df.columns:
            continue
            
        # 跳过包含NaN的行
        col_data = df[col].dropna()
        if len(col_data) < 4:
            continue
            
        if method == "percentile":
            lower_bound = np.percentile(col_data, lower_percentile)
            upper_bound = np.percentile(col_data, upper_percentile)
            outlier_info["method_info"] = f"Percentile: {lower_percentile}% - {upper_percentile}%"
            
        elif method == "iqr":
            Q1 = col_data.quantile(0.25)
            Q3 = col_data.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - iqr_multiplier * IQR
            upper_bound = Q3 + iqr_multiplier * IQR
            outlier_info["method_info"] = f"IQR: multiplier={iqr_multiplier}"
        else:
            raise ValueError(f"未知的离群值检测方法: {method}")
        
        # 标记离群值
        col_outliers = (df[col] < lower_bound) | (df[col] > upper_bound)
        
        # 如果指定了感兴趣药物，则不将它们标记为离群值
        if interest_drugs and "drug" in df.columns:
            interest_mask = df["drug"].isin(interest_drugs)
            col_outliers = col_outliers & ~interest_mask
        
        outlier_mask = outlier_mask | col_outliers
        
        # 记录离群值信息
        n_outliers = col_outliers.sum()
        outlier_info["outliers_by_column"][col] = n_outliers
        
        # 记录每个离群值的详细信息
        for idx in df[col_outliers].index:
            outlier_info["outlier_details"].append({
                "drug": df.loc[idx, "drug"] if "drug" in df.columns else idx,
                "column": col,
                "value": df.loc[idx, col],
                "bounds": (lower_bound, upper_bound)
            })
    
    # 统计总离群值数量
    outlier_info["total_outliers"] = outlier_mask.sum()
    
    # 去除离群值
    df_outliers_removed = df_clean[~outlier_mask].copy()
    
    print(f"[INFO] 离群值检测结果 (方法: {method}):")
    if method == "percentile":
        print(f"  - 百分位数范围: {lower_percentile}% - {upper_percentile}%")
    elif method == "iqr":
        print(f"  - IQR倍数: {iqr_multiplier}")
    print(f"  - 检测列: {columns}")
    print(f"  - 总共离群值: {outlier_info['total_outliers']} 个")
    
    if interest_drugs:
        print(f"  - 感兴趣药物不参与离群值检测: {interest_drugs}")
    
    for col, count in outlier_info["outliers_by_column"].items():
        print(f"  - 列 '{col}' 离群值: {count} 个")
    
    if outlier_info["total_outliers"] > 0:
        print("  - 离群值药物示例 (前10个):")
        for detail in outlier_info["outlier_details"][:10]:
            print(f"    - 药物 '{detail['drug']}': {detail['column']}={detail['value']:.3f} (范围: [{detail['bounds'][0]:.3f}, {detail['bounds'][1]:.3f}])")
        if len(outlier_info["outlier_details"]) > 10:
            print(f"    ... 还有 {len(outlier_info['outlier_details']) - 10} 个离群值未显示")
    
    return df_outliers_removed, outlier_info


# =========================
# 10. 根据排名计算点的大小
# =========================

def calculate_point_size(rank, rank_dict, min_size=20, max_size=120):
    """
    根据排名计算点的大小
    
    参数:
    rank: 排名值
    rank_dict: 所有药物的排名字典
    min_size: 最小点大小
    max_size: 最大点大小
    
    返回:
    点的大小
    """
    if pd.isna(rank):
        return (min_size + max_size) / 2  # 默认大小
    
    # 获取所有排名值
    all_ranks = list(rank_dict.values())
    
    if not all_ranks:
        return (min_size + max_size) / 2
    
    # 转换为数值
    all_ranks = [float(r) for r in all_ranks if not pd.isna(r)]
    
    if not all_ranks:
        return (min_size + max_size) / 2
    
    # 排名越小，点越大
    min_rank = min(all_ranks)
    max_rank = max(all_ranks)
    
    # 归一化排名（排名越小，值越大）
    if max_rank == min_rank:
        normalized_rank = 0.5
    else:
        normalized_rank = 1.0 - ((float(rank) - min_rank) / (max_rank - min_rank))
    
    # 计算点大小
    size = min_size + normalized_rank * (max_size - min_size)
    
    return size


# =========================
# 11. 标注与绘图（使用排名信息）
# =========================

def plot_quadrant(score_df, out_pdf, title, outlier_info=None, 
                 interest_drugs=None, remove_uninterest_drugs=None, 
                 is_full_drugs_plot=False, rank_dict=None):
    """
    绘制四象限图，使用good_log2fc_z和bad_log2fc_z
    
    参数:
    score_df: 包含绘图数据的DataFrame
    out_pdf: 输出PDF文件路径
    title: 图标题
    outlier_info: 离群值信息字典
    interest_drugs: 感兴趣药物列表
    remove_uninterest_drugs: 要从图中删除的不感兴趣药物列表
    is_full_drugs_plot: 是否为全部药物图（不去离群值）
    rank_dict: 排名字典
    """
    # 使用log2fc_z进行绘图
    x_col = "bad_log2fc_z"  # x轴：bad DNB z-score
    y_col = "good_log2fc_z"  # y轴：good DNB z-score
    
    plot_df = score_df.dropna(subset=[x_col, y_col]).copy()
    
    # 如果不感兴趣药物列表不为空，从图中删除这些药物
    if remove_uninterest_drugs:
        before_count = len(plot_df)
        plot_df = plot_df[~plot_df["drug"].isin(remove_uninterest_drugs)].copy()
        after_count = len(plot_df)
        print(f"[INFO] 从不感兴趣药物列表中删除 {before_count - after_count} 个药物")
    
    if plot_df.empty:
        print("[跳过] 没有足够数据用于绘制四象限图")
        return

    setup_plot_style()
    fig, ax = plt.subplots(figsize=(9, 7))

    # 计算每个象限的药物数量
    quadrant_counts = {}
    for q in ["Q1", "Q2", "Q3", "Q4"]:
        sub = plot_df[plot_df["quadrant"] == q]
        quadrant_counts[q] = len(sub)
    
    # 分离感兴趣药物和普通药物
    interest_drugs_in_plot = []
    normal_drugs = []
    
    if interest_drugs:
        interest_drugs_in_plot = plot_df[plot_df["drug"].isin(interest_drugs)]
        normal_drugs = plot_df[~plot_df["drug"].isin(interest_drugs)]
    else:
        normal_drugs = plot_df
    
    # 绘制普通药物（统一为灰色，根据排名调整大小）
    if not normal_drugs.empty:
        sizes = []
        for _, row in normal_drugs.iterrows():
            drug_name = row["drug"]
            if rank_dict and drug_name in rank_dict:
                rank = rank_dict[drug_name]
                size = calculate_point_size(rank, rank_dict, RANK_SIZE_MIN, RANK_SIZE_MAX)
            else:
                size = (RANK_SIZE_MIN + RANK_SIZE_MAX) / 2
            sizes.append(size)
        
        # 使用统一的灰色绘制所有普通药物
        ax.scatter(
            normal_drugs[x_col], normal_drugs[y_col],
            s=sizes, alpha=RANK_ALPHA,
            color=RANK_COLOR,
            edgecolor="white", linewidth=0.5,
            label=f"All drugs (n={len(normal_drugs)})",
            zorder=2
        )
    
    # 绘制感兴趣药物（特殊标记，根据排名调整大小）
    if not interest_drugs_in_plot.empty:
        interest_sizes = []
        for _, row in interest_drugs_in_plot.iterrows():
            drug_name = row["drug"]
            if rank_dict and drug_name in rank_dict:
                rank = rank_dict[drug_name]
                size = calculate_point_size(rank, rank_dict, RANK_SIZE_MIN, RANK_SIZE_MAX)
            else:
                size = INTEREST_DRUG_MARKER_SIZE
            interest_sizes.append(size)
        
        ax.scatter(
            interest_drugs_in_plot[x_col], interest_drugs_in_plot[y_col],
            s=interest_sizes, alpha=0.9,
            color=INTEREST_DRUG_MARKER_COLOR,
            edgecolor=INTEREST_DRUG_MARKER_EDGE_COLOR, linewidth=1.5,
            marker=INTEREST_DRUG_MARKER_SHAPE,
            zorder=3
        )

    # 添加坐标轴
    ax.axhline(0, color="black", linewidth=0.9, linestyle="--", zorder=1)
    ax.axvline(0, color="black", linewidth=0.9, linestyle="--", zorder=1)

    # 设置标签
    ax.set_xlabel("Bad DNB z-score (log2FC)", fontsize=12)
    ax.set_ylabel("Good DNB z-score (log2FC)", fontsize=12)
    ax.set_title(title, fontsize=13, pad=12)
    
    # 添加离群值信息
    if outlier_info and outlier_info["total_outliers"] > 0 and not is_full_drugs_plot:
        if outlier_info.get("method") == "percentile":
            method_info = f"Percentile: {OUTLIER_PERCENTILE_LOW}% - {OUTLIER_PERCENTILE_HIGH}%"
        else:
            method_info = f"IQR: {OUTLIER_IQR_MULTIPLIER}x"
        ax.text(0.02, 0.02, f"Outliers removed: {outlier_info['total_outliers']}\n({method_info})",
                transform=ax.transAxes, fontsize=9, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.5))
    
    # 如果是全部药物图，添加标注
    if is_full_drugs_plot:
        ax.text(0.02, 0.02, "All drugs (no outlier removal)",
                transform=ax.transAxes, fontsize=9, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.5))
    
    # 添加排名信息说明
    if rank_dict and len(rank_dict) > 0:
        ax.text(0.02, 0.10, "Point size indicates rank\n(smaller rank = larger point)",
                transform=ax.transAxes, fontsize=8, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.5))

    # 添加标签
    texts = []
    
    if interest_drugs:
        for _, row in interest_drugs_in_plot.iterrows():
            texts.append(ax.text(row[x_col], row[y_col], 
                                 row["drug"], fontsize=LABEL_FONT_SIZE+2, 
                                 fontweight='bold', color='red', alpha=0.9))
    
    if LABEL_ALL_POINTS:
        for _, row in normal_drugs.iterrows():
            texts.append(ax.text(row[x_col], row[y_col], 
                                 row["drug"], fontsize=LABEL_FONT_SIZE, 
                                 alpha=0.7))
    else:
        for q in ["Q1", "Q2", "Q3", "Q4"]:
            sub = normal_drugs[normal_drugs["quadrant"] == q]
            if sub.empty:
                continue
                
            sub = sub.copy()
            sub["distance"] = np.sqrt(sub[x_col]**2 + sub[y_col]**2)
            sub_sorted = sub.sort_values("distance", ascending=False)
            
            count = 0
            for _, row in sub_sorted.iterrows():
                if count >= TOP_LABEL_PER_QUADRANT:
                    break
                texts.append(ax.text(row[x_col], row[y_col], 
                                     row["drug"], fontsize=LABEL_FONT_SIZE, 
                                     alpha=0.7))
                count += 1
    
    # 自动避让标签
    if HAS_ADJUST_TEXT and texts:
        try:
            adjust_text(texts, ax=ax, 
                       arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
                       expand_points=(1.5, 1.5), expand_text=(1.2, 1.2),
                       force_text=(0.5, 1.0), force_points=(0.5, 1.0))
        except Exception as e:
            print(f"[警告] adjust_text失败: {e}")
            pass
    
    # 添加图例
    from matplotlib.patches import Patch
    
    legend_elements = []
    
    if interest_drugs:
        legend_elements.append(
            Patch(facecolor=INTEREST_DRUG_MARKER_COLOR, edgecolor=INTEREST_DRUG_MARKER_EDGE_COLOR, 
                  linewidth=1.5, label='Interest Drugs')
        )
    
    if rank_dict and len(rank_dict) > 0:
        # 添加排名说明
        legend_elements.append(
            Patch(facecolor=RANK_COLOR, edgecolor='white', linewidth=0.5, 
                  alpha=RANK_ALPHA, label=f'All drugs (n={len(normal_drugs)})')
        )
    else:
        legend_elements.append(
            Patch(facecolor=RANK_COLOR, edgecolor='white', linewidth=0.5, 
                  alpha=RANK_ALPHA, label=f'All drugs (n={len(normal_drugs)})')
        )
    
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=9)

    # 美化图形
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    ax.tick_params(axis="both", labelsize=10, width=1.0)
    
    # 设置坐标轴范围
    x_min, x_max = plot_df[x_col].min(), plot_df[x_col].max()
    y_min, y_max = plot_df[y_col].min(), plot_df[y_col].max()
    
    x_margin = (x_max - x_min) * 0.1
    y_margin = (y_max - y_min) * 0.1
    ax.set_xlim(x_min - x_margin, x_max + x_margin)
    ax.set_ylim(y_min - y_margin, y_max + y_margin)

    plt.tight_layout()
    plt.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=300)
    plt.close()
    print(f"[OK] 四象限图已保存: {out_pdf}")


# =========================
# 12. 主程序
# =========================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # 创建输出目录
    quadrant_plot_dir = os.path.join(OUT_DIR, "Quadrant_Plots")
    data_table_dir = os.path.join(OUT_DIR, "Data_Tables")
    quadrant_table_dir = os.path.join(OUT_DIR, "Quadrant_Tables")
    outlier_info_dir = os.path.join(OUT_DIR, "Outlier_Info")
    comparison_dir = os.path.join(OUT_DIR, "Comparison_Plots")
    
    for dir_path in [quadrant_plot_dir, data_table_dir, quadrant_table_dir, 
                     outlier_info_dir, comparison_dir]:
        os.makedirs(dir_path, exist_ok=True)

    # 步骤1：读取排名信息
    rank_dict = {}
    if USE_RANK_INFO and os.path.exists(RANK_FILE):
        rank_dict = load_rank_info(RANK_FILE, RANK_SHEET)
    else:
        print("[INFO] 未使用排名信息或排名文件不存在")

    # 步骤2：读取药物和疾病文件
    # 注意：不感兴趣的药物在读取阶段就被排除
    drug_dict = load_drug_signatures(
        DRUG_DIR,
        abs_fc_cutoff=DRUG_ABS_LOGFC_CUTOFF,
        pvalue_cutoff=DRUG_PVALUE_CUTOFF,
        exclude_drugs=UNINTEREST_DRUGS
    )
    disease_table_dict = load_disease_tables(DISEASE_DIR)
    
    # 汇总所有疾病的结果
    all_disease_results = {}

    for disease_name, disease_df in disease_table_dict.items():
        print(f"\n" + "="*60)
        print(f"[分析] {disease_name}")
        print("="*60)

        try:
            # 步骤3：逐药物计算
            print(f"[INFO] 开始计算药物分数...")
            score_df = compute_all_drug_scores(drug_dict, disease_df)
            print(f"[INFO] 计算完成，共 {len(score_df)} 个药物")

            # 步骤4：计算四象限分类
            quadrant_df = classify_quadrants(score_df)
            
            # 步骤5：检查感兴趣药物是否存在
            if INTEREST_DRUGS:
                found_interest = [drug for drug in INTEREST_DRUGS if drug in quadrant_df["drug"].values]
                not_found = [drug for drug in INTEREST_DRUGS if drug not in quadrant_df["drug"].values]
                if found_interest:
                    print(f"[INFO] 找到 {len(found_interest)}/{len(INTEREST_DRUGS)} 个感兴趣药物: {found_interest}")
                if not_found:
                    print(f"[警告] 未找到 {len(not_found)} 个感兴趣药物: {not_found}")
            
            # 保存原始结果表格
            raw_data_path = os.path.join(data_table_dir, f"{disease_name}_raw_results.csv")
            quadrant_df.to_csv(raw_data_path, index=False, encoding="utf-8-sig")
            print(f"[OK] 原始数据表格已保存: {raw_data_path}")
            
            # 步骤6：生成详细四象限表格（包含排名）
            detail_table = create_quadrant_detail_table(quadrant_df, rank_dict)
            detail_table_path = os.path.join(quadrant_table_dir, f"{disease_name}_quadrant_detail.csv")
            detail_table.to_csv(detail_table_path, index=False, encoding="utf-8-sig")
            print(f"[OK] 详细四象限表格已保存: {detail_table_path}")
            
            # 记录原始数据信息
            n_total = len(quadrant_df)
            n_valid = quadrant_df.dropna(subset=["good_log2fc_z", "bad_log2fc_z"]).shape[0]
            print(f"[INFO] 药物总数: {n_total}, 有效药物数: {n_valid}")
            
            # 步骤7：绘制包含全部药物的图（不去离群值）
            print(f"\n[INFO] 生成包含全部药物的四象限图（不去离群值）...")
            full_drugs_plot_path = os.path.join(comparison_dir, f"{disease_name}_all_drugs_plot.pdf")
            plot_quadrant(
                quadrant_df,
                full_drugs_plot_path,
                f"{disease_name}: All Drugs (No Outlier Removal)",
                None,  # 不提供离群值信息
                INTEREST_DRUGS,
                None,  # 不删除不感兴趣药物（已经在读取时删除）
                is_full_drugs_plot=True,
                rank_dict=rank_dict
            )
            
            # 保存全部药物的数据表格
            all_drugs_table_path = os.path.join(comparison_dir, f"{disease_name}_all_drugs_data.csv")
            detail_table.to_csv(all_drugs_table_path, index=False, encoding="utf-8-sig")
            print(f"[OK] 全部药物数据表格已保存: {all_drugs_table_path}")
            
            # 步骤8：离群值检测与处理
            filtered_df = quadrant_df.copy()
            outlier_info = None
            
            if REMOVE_OUTLIERS and n_valid > 0:
                print(f"\n[INFO] 正在进行离群值检测 (方法: {OUTLIER_METHOD})...")
                filtered_df, outlier_info = detect_and_remove_outliers(
                    quadrant_df, 
                    columns=["good_log2fc_z", "bad_log2fc_z"],  # 使用log2fc_z进行离群值检测
                    method=OUTLIER_METHOD,
                    lower_percentile=OUTLIER_PERCENTILE_LOW,
                    upper_percentile=OUTLIER_PERCENTILE_HIGH,
                    iqr_multiplier=OUTLIER_IQR_MULTIPLIER,
                    interest_drugs=INTEREST_DRUGS
                )
                
                # 保存离群值信息
                if outlier_info and outlier_info["total_outliers"] > 0:
                    outlier_df = pd.DataFrame(outlier_info["outlier_details"])
                    outlier_path = os.path.join(outlier_info_dir, f"{disease_name}_outliers.csv")
                    outlier_df.to_csv(outlier_path, index=False, encoding="utf-8-sig")
                    print(f"[OK] 离群值信息已保存: {outlier_path}")
                    
                    # 保存去除离群值后的数据表格
                    filtered_path = os.path.join(data_table_dir, f"{disease_name}_filtered_results.csv")
                    filtered_df.to_csv(filtered_path, index=False, encoding="utf-8-sig")
                    print(f"[OK] 去离群值后数据表格已保存: {filtered_path}")
                    
                    # 生成去离群值后的详细四象限表格
                    filtered_detail_table = create_quadrant_detail_table(filtered_df, rank_dict)
                    filtered_detail_path = os.path.join(quadrant_table_dir, f"{disease_name}_quadrant_detail_filtered.csv")
                    filtered_detail_table.to_csv(filtered_detail_path, index=False, encoding="utf-8-sig")
                    print(f"[OK] 去离群值后详细四象限表格已保存: {filtered_detail_path}")
            else:
                print(f"[INFO] 跳过离群值检测")
            
            # 保存绘图数据
            plot_data_path = os.path.join(data_table_dir, f"{disease_name}_plot_data.csv")
            filtered_df.to_csv(plot_data_path, index=False, encoding="utf-8-sig")
            print(f"[OK] 绘图数据表格已保存: {plot_data_path}")
            
            # 记录最终数据统计
            n_final = filtered_df.dropna(subset=["good_log2fc_z", "bad_log2fc_z"]).shape[0]
            print(f"[INFO] 最终用于绘图的药物数: {n_final}")
            
            # 统计各象限药物数量
            quadrant_counts = filtered_df["quadrant"].value_counts().to_dict()
            for q in ["Q1", "Q2", "Q3", "Q4"]:
                count = quadrant_counts.get(q, 0)
                if q in ["Q2", "Q4"]:
                    effect = "(reversing)" if q == "Q2" else "(aggravating)" if q == "Q4" else ""
                    print(f"  - {q} {effect}: {count} 个药物")
                else:
                    print(f"  - {q}: {count} 个药物")
            
            # 步骤9：绘制去除离群值后的图
            quadrant_plot_fp = os.path.join(quadrant_plot_dir, f"{disease_name}_quadrant_plot.pdf")
            plot_quadrant(
                filtered_df,
                quadrant_plot_fp,
                f"{disease_name}: Filtered Drugs (Outliers Removed)",
                outlier_info,
                INTEREST_DRUGS,
                None,  # 不删除不感兴趣药物
                is_full_drugs_plot=False,
                rank_dict=rank_dict
            )
            
            # 保存到汇总结果
            all_disease_results[disease_name] = {
                "total_drugs": n_total,
                "valid_drugs": n_valid,
                "outliers_removed": outlier_info["total_outliers"] if outlier_info else 0,
                "final_drugs": n_final,
                "quadrant_counts": quadrant_counts
            }
            
        except Exception as e:
            print(f"[跳过] {disease_name}: 分析失败，错误: {e}")
            import traceback
            traceback.print_exc()
    
    # 保存所有疾病的汇总结果
    if all_disease_results:
        summary_df = pd.DataFrame.from_dict(all_disease_results, orient="index")
        summary_path = os.path.join(OUT_DIR, "disease_analysis_summary.csv")
        summary_df.to_csv(summary_path, encoding="utf-8-sig")
        print(f"\n[OK] 疾病分析汇总表已保存: {summary_path}")
        
        # 打印汇总信息
        print("\n" + "="*60)
        print("分析汇总:")
        print("="*60)
        for disease_name, stats in all_disease_results.items():
            print(f"\n{disease_name}:")
            print(f"  - 总药物数: {stats['total_drugs']}")
            print(f"  - 有效药物数: {stats['valid_drugs']}")
            if REMOVE_OUTLIERS:
                print(f"  - 去除离群值: {stats['outliers_removed']}")
            print(f"  - 最终药物数: {stats['final_drugs']}")
            for q in ["Q1", "Q2", "Q3", "Q4"]:
                count = stats['quadrant_counts'].get(q, 0)
                if count > 0:
                    effect = "(reversing)" if q == "Q2" else "(aggravating)" if q == "Q4" else ""
                    print(f"  - {q} {effect}: {count}")

    print("\n" + "="*60)
    print("全部完成。")
    print(f"结果保存在: {OUT_DIR}")
    print(f"  - 全部药物图（对比用）: {comparison_dir}")
    print(f"  - 四象限图（去除离群值）: {quadrant_plot_dir}")
    print(f"  - 详细四象限表格: {quadrant_table_dir}")
    print(f"  - 数据表格: {data_table_dir}")
    if REMOVE_OUTLIERS:
        print(f"  - 离群值信息: {outlier_info_dir}")


if __name__ == "__main__":
    main()