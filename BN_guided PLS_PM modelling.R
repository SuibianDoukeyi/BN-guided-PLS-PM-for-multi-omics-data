library(plspm)
library(tidyverse)

# ==== Step 1: 读取数据 ====
l2_df <- read.csv("D:/01文件/001论文/博士论文/05毒性转录代谢/PLSPM/Level 2/L2 plspm 82.csv", header = TRUE, check.names = FALSE)

# ==== Step 2: 提取映射表（潜变量→显变量） ====
l2_latent_df <- l2_df[, c("type2", "var_id")]  # type2 是潜变量，var_id 是显变量

# ==== Step 3: 提取表达量矩阵 ====
l2_expr_matrix <- l2_df[, !(colnames(l2_df) %in% c("type1", "type2", "var_id"))]
rownames(l2_expr_matrix) <- l2_df$var_id  # 行名设为显变量ID

# ==== Step 4: 标准化观测变量矩阵 ====
l2_measured_df_t <- as.data.frame(scale(t(l2_expr_matrix)))  # 行：样本，列：显变量

# ==== Step 5: 构建 blocks 列表（潜变量 → 显变量）====
l2_blocks <- split(l2_latent_df$var_id, l2_latent_df$type2)
l2_latent_names <- names(l2_blocks)

# ==== Step 6: 构建环路径矩阵 ====
generate_cycle_path_matrix <- function(latent_names) {
  n <- length(latent_names)
  mat <- matrix(0, n, n, dimnames = list(latent_names, latent_names))
  
  # 环形连接：每个 latent 指向下一个，最后一个指向第一个
  for (i in seq_len(n)) {
    next_i <- if (i == n) 1 else i + 1
    from <- latent_names[i]
    to <- latent_names[next_i]
    if (match(to, latent_names) > match(from, latent_names)) {
      mat[to, from] <- 1  # 满足下三角
    } else {
      mat[from, to] <- 1  # 若不满足，则转换方向为下三角
    }
  }
  
  # 删除对角线
  diag(mat) <- 0
  
  # 强制清除上三角（仅保留下三角）
  mat[upper.tri(mat)] <- 0
  
  return(mat)
}
l2_path_matrix <- generate_cycle_path_matrix(l2_latent_names)

# ==== Step 7: 执行 PLS-PM 初建模 ====
l2_pls_model <- plspm(
  Data = l2_measured_df_t,
  path_matrix = l2_path_matrix,
  blocks = l2_blocks,
  modes = rep("A", length(l2_blocks))
)


# ==== Step 8: 提取 outer model ====
l2_outer_df <- l2_pls_model$outer_model %>%
  tibble::rownames_to_column(var = "observed") %>%
  rename(
    latent = block,
    weight = weight,
    loading = loading,
    communality = communality,
    redundancy = redundancy
  )

# ==== Step 9: 筛选标准 ====
l2_selected_vars <- l2_outer_df %>%
  filter(abs(loading) >= 0.7, communality >= 0.5)

# ==== Step 10: 提取显变量表达矩阵 ====
l2_selected_ids <- as.character(l2_selected_vars$name)
l2_measured_selected_t <- l2_measured_df_t[, l2_selected_ids, drop = FALSE]

# ==== Step 12: 重构 blocks ====
l2_filtered_blocks <- l2_selected_vars %>%
  group_by(latent) %>%
  summarise(vars = list(as.character(name))) %>%
  deframe()

# ==== Step 13: 裁剪路径矩阵 ====
l2_filtered_latents <- names(l2_filtered_blocks)
l2_filtered_path_matrix <- generate_cycle_path_matrix(l2_filtered_latents)

l2_plspm_model_final <- plspm(
  Data = as.data.frame(l2_measured_selected_t),
  path_matrix = l2_filtered_path_matrix,
  blocks = l2_filtered_blocks,
  modes = rep("A", length(l2_filtered_blocks))
)

# ==== Step 15: 提取潜变量得分 ====
l2_lv_scores_final <- as.data.frame(t(l2_plspm_model_final$scores))

# ==== Step 16: 将潜变量得分映射到 L1 表 ====

# 读取 type1 → type2 对应关系表
l2_map_df <- read.csv("D:/01文件/001论文/博士论文/05毒性转录代谢/PLSPM/Level 2/L1 latent to L2.csv", stringsAsFactors = FALSE)

# 转换 lv_scores 行名为列，准备合并
l2_lv_scores_df <- l2_lv_scores_final %>%
  rownames_to_column(var = "type2")

# 合并映射
l2_mapped_scores <- l2_map_df %>%
  left_join(l2_lv_scores_df, by = "type2")


###1###
# 读取 L1 的原始补充数据
l2_l1_sub_df <- read.csv("D:/01文件/001论文/博士论文/05毒性转录代谢/PLSPM/Level 1/l1 sub l1012.csv", stringsAsFactors = FALSE)
# 提取前两列变量信息
var_info <- l2_l1_sub_df[, 1:2]

# 提取样本表达量数据（12列）
expr_mat <- l2_l1_sub_df[, 3:ncol(l2_l1_sub_df)]

# 对每一行进行标准化（按行 scale）
expr_scaled <- t(scale(t(as.matrix(expr_mat))))

# 合并变量信息和标准化结果
l2_l1_sub_scaled <- cbind(var_info, expr_scaled)
# 合并原始和 Level2 推导变量（按下方追加）
l2_final_l1_data <- bind_rows(l2_l1_sub_scaled, l2_mapped_scores)



###下面这几行代码在首次必须运行，避免不scale；然后后“注释”“隐藏”掉，避免二次scale

###2 先合并再scale
# l2_l1_sub_df <- read.csv("D:/01文件/001论文/博士论文/05毒性转录代谢/PLSPM/Level 1/l1 sub l1012.csv", stringsAsFactors = FALSE)
# l2_final_l1_data <- bind_rows(l2_l1_sub_df, l2_mapped_scores)
# var_info <- l2_final_l1_data[, 1:2]
# expr_mat <- l2_final_l1_data[, 3:ncol(l2_final_l1_data)]
# expr_scaled <- t(scale(t(as.matrix(expr_mat))))
# l2_final_l1_scaled <- cbind(var_info, expr_scaled)






####level 1####

# === 1. 加载必要包 ===
library(plspm)
library(bnlearn)
library(Rgraphviz)
library(plspm2)
library(igraph)
# === 2. 读取并预处理数据 ===
df <- l2_final_l1_data
latent_df <- df[, c("type1", "type2")]
data_matrix <- df[, !(colnames(df) %in% c("type1", "type2"))]
rownames(data_matrix) <- df$type2

# === 3. 标准化观测变量矩阵 ===
measured_df_t <- as.data.frame(t(data_matrix))

# === 4. 构建潜变量→显变量 blocks ===
blocks <- split(latent_df$type2, latent_df$type1)
latent_names <- names(blocks)




#################
#################
# === 5. 构建先验路径矩阵 ===
# 潜变量名称
latent_vars <- c(
  "Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae", "Biodegradation",
  "Drug_Metabolism", "Auxin_Metabolism", "Proliferation",
  "Carbohydrate_Metabolism", "Lipid_Metabolism", "Amino_Acid_Metabolism",
  "Oxidative_Response", "Photosystem",
  "Microbial_Function", "Microbial_Diversity", "Microbial_Composition"
)

# 初始化空矩阵
prior_path_matrix <- matrix(0, nrow = length(latent_vars), ncol = length(latent_vars))
colnames(prior_path_matrix) <- rownames(prior_path_matrix) <- latent_vars


## 初级生产与能量流动路径
prior_path_matrix["Photosystem", "Oxidative_Response"] <- 1  # 光合作用产生氧化压力
prior_path_matrix["Photosystem", "Growth_Endpoint"] <- 1      # 光合作用直接影响生长
prior_path_matrix["Photosystem", "Carbohydrate_Metabolism"] <- 1  # 光合产物驱动碳代谢
prior_path_matrix["Photosystem", "Microbial_Function"] <- 1   # 光合产物支持微生物功能


##核心代谢调控路径
prior_path_matrix["Carbohydrate_Metabolism", "Lipid_Metabolism"] <- 1  # 糖代谢提供脂质合成前体
prior_path_matrix["Carbohydrate_Metabolism", "Amino_Acid_Metabolism"] <- 1  # 糖代谢提供氨基酸合成碳骨架
prior_path_matrix["Amino_Acid_Metabolism", "Proliferation"] <- 1      # 氨基酸支持细胞增殖
prior_path_matrix["Oxidative_Response", "Lipid_Metabolism"] <- 1      # 氧化应激影响脂代谢
prior_path_matrix["Oxidative_Response", "Amino_Acid_Metabolism"] <- 1 # 氧化应激影响氨基酸代谢

#生长调控路径
prior_path_matrix["Auxin_Metabolism", "Proliferation"] <- 1    # 生长素代谢促进增殖
prior_path_matrix["Auxin_Metabolism", "Growth_Endpoint"] <- 1  # 生长素直接影响生长
prior_path_matrix["Proliferation", "Growth_Endpoint"] <- 1     # 细胞增殖促进生长
prior_path_matrix["Carbohydrate_Metabolism", "Growth_Endpoint"] <- 1  # 碳代谢支持生长
prior_path_matrix["Lipid_Metabolism", "Growth_Endpoint"] <- 1   # 脂代谢支持生长
prior_path_matrix["Amino_Acid_Metabolism", "Growth_Endpoint"] <- 1  # 氨基酸代谢支持生长
prior_path_matrix["Oxidative_Response", "Growth_Endpoint"] <- 1  # 氧化应激影响生长

#氧化应激与毒性响应路径
prior_path_matrix["Drug_Metabolism", "Oxidative_Response"] <- 1  # 药物代谢诱导氧化应激
prior_path_matrix["Oxidative_Response", "Toxicity_Endpoint"] <- 1  # 氧化应激导致毒性
prior_path_matrix["Drug_Metabolism", "Toxicity_Endpoint"] <- 1    # 药物代谢直接导致毒性

#微生物群落功能路径
prior_path_matrix["Microbial_Composition", "Microbial_Function"] <- 1  # 群落组成决定功能
prior_path_matrix["Microbial_Diversity", "Microbial_Function"] <- 1    # 多样性影响功能
prior_path_matrix["Microbial_Composition", "Biodegradation"] <- 1      # 群落组成影响降解

#抗生素残留与降解路径
prior_path_matrix["Microbial_Function", "Antibiotics_Residue_in_Algae"] <- 1  # 微生物功能影响残留
prior_path_matrix["Drug_Metabolism", "Antibiotics_Residue_in_Algae"] <- 1     # 药物代谢影响残留
prior_path_matrix["Microbial_Function", "Biodegradation"] <- 1      # 微生物功能驱动降解
prior_path_matrix["Drug_Metabolism", "Biodegradation"] <- 1         # 药物代谢影响降解

# === 终端变量约束 确保终端变量仅作为终点 ===
terminal_vars <- c("Toxicity_Endpoint", "Growth_Endpoint", 
                   "Antibiotics_Residue_in_Algae", "Biodegradation")
prior_path_matrix[terminal_vars, ] <- 0  # 清除所有出向路径


# === 自动生成拓扑排序顺序函数 ===
get_ordered_vars <- function(prior_matrix) {
  g <- igraph::graph_from_adjacency_matrix(prior_matrix, mode = "directed")
  top_order <- igraph::topo_sort(g, mode = "out")
  reversed_order <- rev(names(top_order))
  return(reversed_order)
}

# 应用拓扑排序
ordered_vars <- get_ordered_vars(prior_path_matrix)

# 重新排列路径矩阵为下三角结构
prior_path_matrix <- prior_path_matrix[ordered_vars, ordered_vars]

# 查看最终拓扑顺序和路径矩阵
ordered_vars
prior_path_matrix



########
# 建模
pls_model_sparse <- plspm(
  Data = measured_df_t,
  path_matrix = prior_path_matrix,
  blocks = blocks,
  modes = rep("A", length(blocks)),
  scheme = "path" 
)



summary(pls_model_sparse)
# === 7. 使用潜变量得分构建贝叶斯网络，强制结果变量 ===
latent_scores <- pls_model_sparse$scores

# 强制结果潜变量（不能作为其他变量的因变量）
result_vars <- c("Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae","Biodegradation")



# === 先验路径 ===
# 路径重要性评估函数
# 已有27路径的矩阵 prior_path_matrix
# 定义必须保留的核心路径（基于科学重要性）
core_paths <- list(
  c("Proliferation", "Growth_Endpoint"),
  c("Oxidative_Response", "Toxicity_Endpoint"),
  c("Microbial_Function", "Biodegradation"),
  c("Drug_Metabolism", "Antibiotics_Residue_in_Algae"),
  c("Photosystem", "Carbohydrate_Metabolism"),
  c("Carbohydrate_Metabolism", "Amino_Acid_Metabolism"),
  c("Amino_Acid_Metabolism", "Proliferation"),
  c("Drug_Metabolism", "Oxidative_Response"),
  c("Microbial_Composition", "Microbial_Function"),
  c("Microbial_Function", "Antibiotics_Residue_in_Algae"),
  c("Microbial_Diversity", "Growth_Endpoint"),  # 可选但建议保留
  c("Auxin_Metabolism", "Proliferation")       # 可选但建议保留
)

# 创建白名单数据框
white_list <- do.call(rbind, lapply(core_paths, function(p) {
  data.frame(from = p[1], to = p[2])
}))

# 定义终端变量
terminal_vars <- c("Toxicity_Endpoint", "Growth_Endpoint", 
                   "Antibiotics_Residue_in_Algae", "Biodegradation")

# 创建黑名单：禁止终端变量指向任何变量
black_list <- expand.grid(
  from = terminal_vars,
  to = latent_vars
)

# 添加自环禁止（通常BN不允许自环）
black_list <- rbind(black_list, data.frame(
  from = latent_vars,
  to = latent_vars
))

library(bnlearn)


# 执行带约束的贝叶斯网络结构学习
bn_model <- hc(
  x = as.data.frame(latent_scores),
  whitelist = white_list,  # 强制保留核心路径
  blacklist = black_list,  # 禁止终端变量作为原因
  score = "bic-g"           # 使用BIC分数（自动惩罚复杂度）
)

# 查看学习到的结构
graphviz.plot(bn_model, fontsize = 20)


#评估路径强度
# 执行自助法路径强度评估
arc_strength <- boot.strength(
  as.data.frame(latent_scores), 
  R = 200,                # 重复次数（根据数据量调整）
  algorithm = "hc",
  algorithm.args = list(
    whitelist = white_list,
    blacklist = black_list,
    score = "bic-g"
  )
)

# 筛选显著路径（强度>0.8且方向确定性>0.5）
significant_arcs <- arc_strength[
  arc_strength$strength > 0.8 & arc_strength$direction > 0.5, 
]

# 合并核心路径和显著发现的次要路径
final_arcs <- unique(rbind(
  white_list,
  significant_arcs[, c("from", "to")]
))

# 创建最终先验矩阵
final_prior_matrix <- matrix(0, nrow = length(latent_vars), ncol = length(latent_vars))
colnames(final_prior_matrix) <- rownames(final_prior_matrix) <- latent_vars

for(i in 1:nrow(final_arcs)) {
  from <- final_arcs$from[i]
  to <- final_arcs$to[i]
  final_prior_matrix[from, to] <- 1
}

# 应用拓扑排序
final_ordered_vars <- get_ordered_vars(final_prior_matrix)
final_prior_matrix <- final_prior_matrix[final_ordered_vars, final_ordered_vars]

#验证精简模型
# 比较模型复杂度
cat("原始路径数:", sum(prior_path_matrix), "\n")
cat("精简路径数:", sum(final_prior_matrix), "\n")

# 检查终端变量约束
terminal_out_degree <- rowSums(final_prior_matrix[terminal_vars, ])
cat("终端变量出度:", terminal_out_degree, "\n")  # 应全为0




# === 8. 从 BN 获取邻接矩阵并拓扑排序重排 ===
path_matrix_bn <- amat(bn_model)

diag(path_matrix_bn) <- 0
rownames(path_matrix_bn) <- colnames(path_matrix_bn) <- colnames(latent_scores)

# 拓扑排序（确保所有箭头都落在下三角）
bn_edges <- arcs(bn_model)
bn_graph <- igraph::graph_from_data_frame(bn_edges, directed = TRUE)

# 获取拓扑排序（解释变量在前，因变量在后）
topo_order <- rev(topo_sort(bn_graph, mode = "out")) |> names()

# 重新排序路径矩阵
path_matrix_ordered <- path_matrix_bn[topo_order, topo_order]

# 只保留下三角路径
path_matrix_final <- path_matrix_ordered
path_matrix_final[upper.tri(path_matrix_final)] <- 0

# === 对 blocks 重排序以匹配 path_matrix_final 的顺序 ===
latent_order <- colnames(path_matrix_final)

# 仅保留并重排 blocks 中与 path_matrix_final 匹配的潜变量
blocks_aligned <- blocks[latent_order]



# === 9. 使用 BN 路径矩阵构建最终 PLS-PM 模型 ===
pls_model_final <- plspm(
  Data = measured_df_t,
  path_matrix = path_matrix_final,
  blocks = blocks_aligned,
  modes = rep("A", length(blocks_aligned)),
  scheme = "path"
)

summary(pls_model_final)


# 按照 path_matrix 的顺序重新排列 blocks
blocks_aligned <- blocks[colnames(path_matrix_final)]
all(unlist(blocks) %in% colnames(measured_df_t))
rownames(data_matrix) <- df$type2 


names(blocks)
colnames(path_matrix_final)
head(colnames(measured_df_t))

# === 10. 使用 BN 结构坐标绘图（基于 igraph + layout_from_bn） ===
library(igraph)
library(qgraph)

# 将贝叶斯网络转换为 igraph 对象（方向：from → to）
bn_edges <- arcs(bn_model)
bn_graph <- graph_from_data_frame(bn_edges, directed = TRUE)

# 使用 igraph 自动布局（可选 layout 如：layout_with_sugiyama / layout_as_tree / layout_nicely）
coords <- layout_with_sugiyama(bn_graph)$layout
rownames(coords) <- V(bn_graph)$name

# 提取路径系数矩阵
path_coefs <- pls_model_final$path_coefs

# 构建路径边表（含系数）
edge_idx <- which(path_coefs != 0, arr.ind = TRUE)
edge_df <- data.frame(
  from = rownames(path_coefs)[edge_idx[, 1]],
  to = colnames(path_coefs)[edge_idx[, 2]],
  coef = path_coefs[edge_idx]
)

# 构建 igraph 对象用于绘图
g_pls <- graph_from_data_frame(edge_df, directed = TRUE, vertices = rownames(coords))

# === 确保坐标顺序一致 ===
coords_mat <- coords[match(V(g_pls)$name, rownames(coords)), ]

# === 11. ggraph 美化绘图 ===
library(igraph)
library(tidygraph)
library(ggraph)
library(dplyr)
library(grid)
library(scales)
library(ggforce)
# === 确保坐标是数值型，并绑定到 tidygraph 对象 ===
coords_df <- as.data.frame(coords_mat)
coords_df$name <- rownames(coords_mat)
colnames(coords_df)[1:2] <- c("x", "y")
coords_df$x <- as.numeric(coords_df$x)
coords_df$y <- as.numeric(coords_df$y)

# 合并坐标到节点数据
graph_tbl <- as_tbl_graph(g_pls) %>%
  activate(nodes) %>%
  left_join(coords_df, by = "name")


# === 自动生成 special_edges 模板 ===
special_edges_template <- edge_df %>%
  select(from, to) %>%
  mutate(
    curvature = 0.15,  # 默认值，可手动修改
    angle = 20         # 默认箭头锐角度数，可手动修改
  )

# 查看或导出以便编辑
print(special_edges_template)


# === 绘图 ===
p_plspm_final <- ggraph(graph_tbl, layout = "manual", x = x, y = y) +  # 使用 BN 提供的手动布局
  
  # === 边：弯曲边并显示路径系数 ===
  geom_edge_arc(
    aes(
      label = round(coef, 2.5),           # 路径系数文本
      color = coef,                     # 边的颜色代表系数正负
      width = abs(coef)                 # 边宽表示系数大小
    ),
    strength = 0.25,  # ← 控制弧度，0 是直线，1 是全弧（默认），建议 0.3–0.7 之间调节
    arrow = grid::arrow( length = unit(3, "mm"),
                         type = "closed",         # 保持实心
                         angle = 20,              # 尖锐角度（默认30）
                         ends = "last"),  # ← 添加箭头
    end_cap = circle(8, 'mm'),    # ← 提前终止边，留出箭头空间
    angle_calc = 'along',               # 标签沿边方向旋转
    label_dodge = unit(2.5, "mm"),        # 标签离边的距离
    label_size = 4,                     # 标签字体大小
    show.legend = TRUE ,                 # 显示图例
    lineend = "round"                   # 边线两端为圆形
  ) +
  
  # === 边颜色与粗细比例映射设置 ===
  scale_edge_color_gradient2(
    low = "blue", mid = "gray70", high = "red",  # 红→灰→蓝表示负→零→正系数
    midpoint = 0,
    name = "Path Coefficient"
  ) +
  scale_edge_width(range = c(0.5, 2)) +         # 控制边宽范围
  
  # === 节点：矩形标签，显示名称 ===
  geom_node_label(
    aes(label = name),
    family = "sans",
    size = 4,
    label.size = 0.25,
    label.padding = unit(2, "mm"),
    label.r = unit(2, "mm"),     # 圆角大幅增加，接近椭圆
    fill = "lightblue",             # 椭圆填色
    color = "black"
  )+
  
  # === 主题设置：去除背景、添加标题 ===
  theme_void() +
  ggtitle("PLS-PM BN Causal Network") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")  # 居中标题
  )
pls_model_final


summary( pls_model_final)

str(pls_model_final)

pls_model_final$outer_model




# ==== 计算 ΔBIC（放在 bn_model 与 latent_scores 已生成之后）====

library(bnlearn)

# 使用的评分类型：latent_scores 为连续型 → 高斯 BIC
bic_type <- "bic-g"   # 若是离散数据，用 "bic"

# 数据（用你的潜变量得分）
bn_data <- as.data.frame(latent_scores)

# --- 辅助函数：从邻接矩阵或白名单构建 BN 对象（去掉 igraph::is.dag 检查） ---
build_bn_from_amat <- function(amat_mat) {
  stopifnot(is.matrix(amat_mat))
  nodes <- colnames(amat_mat)
  net <- empty.graph(nodes)
  one_idx <- which(amat_mat == 1, arr.ind = TRUE)
  if (nrow(one_idx) > 0) {
    froms <- rownames(amat_mat)[one_idx[, 1]]
    tos   <- colnames(amat_mat)[one_idx[, 2]]
    for (i in seq_len(nrow(one_idx))) {
      net <- set.arc(net, from = froms[i], to = tos[i])
    }
  }
  # 不做 igraph::is.dag 检查；若结构有环，后面的 score() 会报错并被 tryCatch 捕获
  net
}

build_bn_from_whitelist <- function(nodes, wl_df) {
  net <- empty.graph(nodes)
  if (!is.null(wl_df) && nrow(wl_df) > 0) {
    for (i in seq_len(nrow(wl_df))) {
      net <- set.arc(net, from = wl_df$from[i], to = wl_df$to[i])
    }
  }
  # 同样不做 igraph::is.dag 检查
  net
}


#--- 构造候选模型清单（沿用你已有对象名） ---
# 1) 学得的结构（hc 结果）
bn_learned <- bn_model

# 2) 先验路径矩阵（prior_path_matrix）
bn_prior   <- build_bn_from_amat(prior_path_matrix)

# 3) 稳健边 + 白名单合并后的矩阵（final_prior_matrix）
bn_robust  <- build_bn_from_amat(final_prior_matrix)

# 4) 仅白名单核心路径
bn_white   <- build_bn_from_whitelist(latent_vars, white_list)

# 5) 独立网络（空图）
bn_indep   <- empty.graph(latent_vars)

candidates <- list(
  learned = bn_learned,
  prior   = bn_prior,
  robust  = bn_robust,
  whitelist_only = bn_white,
  independence  = bn_indep
)

#--- 计算各模型的 BIC ---
rows <- list()
for (nm in names(candidates)) {
  net <- candidates[[nm]]
  # 尝试评分，若结构非法会报错；用 tryCatch 保护一下
  sc <- tryCatch(score(net, data = bn_data, type = bic_type),
                 error = function(e) NA_real_)
  rows[[nm]] <- data.frame(model = nm, BIC = sc, stringsAsFactors = FALSE)
}
bic_tbl <- do.call(rbind, rows)

# 去除无法评分的模型（若有）
bic_tbl <- subset(bic_tbl, is.finite(BIC))

# ΔBIC（相对最优）
best_bic <- min(bic_tbl$BIC)
bic_tbl$DeltaBIC_to_best <- bic_tbl$BIC - best_bic

# ΔBIC（相对参考）：你可以选 "prior" 或 "independence" 或 "learned"
reference_name <- "prior"   # ← 如需改参照模型，改这里即可
ref_bic <- bic_tbl$BIC[bic_tbl$model == reference_name][1]
bic_tbl$DeltaBIC_vs_ref <- bic_tbl$BIC - ref_bic

# 排序并查看
bic_tbl <- bic_tbl[order(bic_tbl$BIC), ]
print(bic_tbl, row.names = FALSE)

# 如需导出：
# write.csv(bic_tbl, "BN_BIC_summary.csv", row.names = FALSE)









# === 找出所有节点的"源头“，并统计出现频率 ===
# === 将路径矩阵转为 igraph 对象 ===
plspm_graph <- graph_from_adjacency_matrix(
  path_matrix_final,
  mode = "directed",
  diag = FALSE
)

# === 统计每个节点可达的后继节点数量 ===
reach_counts <- sapply(V(plspm_graph)$name, function(node) {
  length(subcomponent(plspm_graph, node, mode = "out")) - 1
})
reach_sorted <- sort(reach_counts, decreasing = TRUE)
print(reach_sorted)


# === 定义终端节点 ===
terminal_nodes <- c("Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae","Biodegradation")

# === 识别所有终端节点的祖先节点 ===
all_ancestors <- unique(unlist(sapply(terminal_nodes, function(node) {
  ancestors(bn_model, node)
})))

# === 提取所有祖先到终端的路径 ===
all_paths <- list()
for (target in terminal_nodes) {
  for (source in all_ancestors) {
    paths <- all_simple_paths(bn_graph, from = source, to = target, mode = "out")
    if (length(paths) > 0) {
      all_paths <- c(all_paths, paths)
    }
  }
}

# === 合并路径中涉及的所有边 ===
ancestor_edges <- unique(do.call(rbind, lapply(all_paths, function(path) {
  path_nodes <- names(path)
  if (length(path_nodes) < 2) return(NULL)
  data.frame(
    from = head(path_nodes, -1),
    to = tail(path_nodes, -1),
    stringsAsFactors = FALSE
  )
})))

coef_matrix <- pls_model_final$path_coefs

# 提取路径系数
ancestor_edges$coef <- mapply(function(f, t) {
  coef_matrix[f, t]
}, ancestor_edges$from, ancestor_edges$to)


# === 构建 igraph 对象用于绘图（祖先网络子图） ===
g_ancestor <- graph_from_data_frame(ancestor_edges, directed = TRUE)

# === 自动布局或使用 sugiyama 布局 ===
layout_coords <- layout_with_sugiyama(g_ancestor)$layout
rownames(layout_coords) <- V(g_ancestor)$name
layout_df <- as.data.frame(layout_coords)
layout_df$name <- rownames(layout_coords)
colnames(layout_df)[1:2] <- c("x", "y")

# === tidygraph 转换并绘图 ===
library(tidygraph)
library(ggraph)
library(ggforce)
library(igraph)
g_tbl <- as_tbl_graph(g_ancestor) %>%
  activate(nodes) %>%
  left_join(layout_df, by = "name")

# === 计算边的中点位置，用于展示系数文本 ===
edge_label_df <- ancestor_edges %>%
  rowwise() %>%
  mutate(
    x = (layout_df$x[layout_df$name == from] + layout_df$x[layout_df$name == to]) / 2,
    y = (layout_df$y[layout_df$name == from] + layout_df$y[layout_df$name == to]) / 2
  )
# === 绘制最终图像 ===
p_ancestor <- ggraph(g_tbl, layout = "manual", x = x, y = y) +
  # 原始边
  geom_edge_link(
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(8, 'mm'),
    edge_width = 0.8,
    color = "darkblue"
  ) +
  
  # 显示边权值（路径系数）
  geom_text(
    data = edge_label_df,
    aes(x = x, y = y, label = round(coef, 2)),
    size = 3,
    color = "black",
    vjust = -0.5
  ) +
  
  # 节点标签
  geom_node_label(
    aes(label = name),
    size = 4,
    label.padding = unit(2, "mm"),
    label.r = unit(2, "mm"),
    fill = ifelse(V(g_ancestor)$name %in% terminal_nodes, "salmon", "lightblue"),
    color = "black"
  ) +
  
  theme_void() +
  ggtitle("祖先节点至终端变量路径图") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
# 显示图
print(p_ancestor)


# === 1. 获取所有祖先节点===
# terminal_nodes 已定义
all_ancestors <- unique(unlist(sapply(terminal_nodes, function(node) {
  ancestors(bn_model, node)
})))

# === 2. 获取每个祖先节点的入度（in-degree）===
# 从 bn_model 提取所有 arcs
arcs_df <- as.data.frame(arcs(bn_model), stringsAsFactors = FALSE)

# 计算每个节点的入度
in_degree <- table(arcs_df$to)

# 筛选出祖先节点中 in-degree == 0 的（无父节点）
root_ancestors <- setdiff(all_ancestors, names(in_degree))  # 没有出现在“被指向”的变量中

# 输出最源头祖先节点
print(root_ancestors)




# === 1. 模块定义 ===
module_def <- data.frame(
  name = c(
    "Microbial_Function", "Microbial_Diversity", "Microbial_Composition",
    "Carbohydrate_Metabolism", "Amino_Acid_Metabolism", "Lipid_Metabolism", "Proliferation",
    "Photosystem", "Oxidative_Response", "Auxin_Metabolism", "Drug_Metabolism"
  ),
  module = c(
    rep("Microbial", 3),
    rep("Primary metabolism", 4),
    rep("Secondary metabolism", 4)
  )
)

# === 2. 终端与祖先节点 ===
terminal_nodes <- c("Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae", "Biodegradation")


# 手动横向排列终端节点顺序（左→右）
# 改进后的横坐标：更紧凑更居中
terminal_x_map <- data.frame(
  name = terminal_nodes,
  x_terminal = c(0.5, 1.5, 2.5, 3.5)  # 统一左起展开
)

# === 3. Sugiyama 初始布局 ===
layout_coords <- layout_with_sugiyama(g_ancestor)$layout
rownames(layout_coords) <- V(g_ancestor)$name
layout_df <- as.data.frame(layout_coords)
colnames(layout_df)[1:2] <- c("x", "y")
layout_df$name <- rownames(layout_coords)

# Step 2: 坐标构建
layout_df_mod <- layout_df %>%
  mutate(name = as.character(name)) %>%
  left_join(module_def, by = "name") %>%
  mutate(
    role = case_when(
      name %in% terminal_nodes ~ "Terminal",
      name %in% root_ancestors ~ "Root",
      TRUE ~ "Intermediate"
    )
  ) %>%
  left_join(terminal_x_map, by = "name") %>%
  mutate(
    # 横向位置（终端单独排，模块 1-2-3）
    x = case_when(
      role == "Terminal" ~ x_terminal,
      module == "Primary metabolism" ~ 1,
      module == "Secondary metabolism" ~ 2,
      module == "Microbial" ~ 3,
      TRUE ~ 2
    ),
    # 纵向结构层级
    y = case_when(
      role == "Root" ~ 2.5,
      role == "Terminal" ~ 1.5,
      TRUE ~ 2
    )
  ) %>%
  group_by(x, y) %>%
  mutate(y_spread = y + (row_number() - median(row_number())) * 0.25) %>%
  ungroup()

# Step 3: 构建图
graph_tbl <- as_tbl_graph(g_pls) %>%
  activate(nodes) %>%
  left_join(layout_df_mod %>% select(name, x, y = y_spread), by = "name")

# Step 4: 绘图（已去除模块背景与图例）
p_plspm_final <- ggraph(graph_tbl, layout = "manual", x = x, y = y) +
  geom_edge_arc(
    aes(
      label = ifelse(abs(coef) < 0.1, "", round(coef, 2)),
      color = coef,
      alpha = ifelse(abs(coef) < 0.1, 0.2, 1),
      width = abs(coef)
    ),
    strength = 0.25,
    arrow = grid::arrow(length = unit(3, "mm"), type = "closed", angle = 20, ends = "last"),
    end_cap = circle(6, 'mm'),
    angle_calc = 'along',
    label_dodge = unit(3, "mm"),
    label_size = 3.8,
    show.legend = TRUE,
    lineend = "round"
  ) +
  scale_edge_color_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0,
    name = "Path Coef."
  ) +
  scale_edge_width(range = c(0.4, 2.2)) +
  scale_edge_alpha(guide = "none") +
  geom_node_label(
    aes(label = name),
    family = "Times New Roman",
    fontface = "bold",
    size = 5,
    label.size = 0,
    label.padding = unit(2, "mm"),
    label.r = unit(2.5, "mm"),
    fill = "lightblue",
    color = "black"
  ) +
  theme_void() +
  ggtitle("Causal network of latent variables") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"),
    legend.position = "right",
    legend.title = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 9, family = "Times New Roman")
  )

# Step 5: 显示图
p_plspm_final




### 源头节点重绘
# === 构建 igraph 子图 ===
g_ancestor <- graph_from_data_frame(ancestor_edges, directed = TRUE)

# === 模块定义 ===
module_def <- data.frame(
  name = c(
    "Microbial_Function", "Microbial_Diversity", "Microbial_Composition",
    "Carbohydrate_Metabolism", "Amino_Acid_Metabolism", "Lipid_Metabolism", "Proliferation",
    "Photosystem", "Oxidative_Response", "Auxin_Metabolism", "Drug_Metabolism"
  ),
  module = c(
    rep("Microbial", 3),
    rep("Primary metabolism", 4),
    rep("Secondary metabolism", 4)
  )
)

# === 定义终端节点、祖先节点、横坐标映射 ===
terminal_nodes <- c("Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae", "Biodegradation")
ancestors_p_ancestors <- root_ancestors

terminal_x_map <- data.frame(
  name = terminal_nodes,
  x_terminal = c(0.5, 1.5, 2.5, 3.5)
)

# === 构建节点坐标 layout_df ===
layout_coords <- layout_with_sugiyama(g_ancestor)$layout
rownames(layout_coords) <- V(g_ancestor)$name
layout_df <- as.data.frame(layout_coords)
colnames(layout_df)[1:2] <- c("x", "y")
layout_df$name <- rownames(layout_coords)

layout_df_mod <- layout_df %>%
  mutate(name = as.character(name)) %>%
  left_join(module_def, by = "name") %>%
  mutate(
    role = case_when(
      name %in% terminal_nodes ~ "Terminal",
      name %in% ancestors_p_ancestors ~ "Root",
      TRUE ~ "Intermediate"
    )
  ) %>%
  left_join(terminal_x_map, by = "name") %>%
  mutate(
    x = case_when(
      role == "Terminal" ~ x_terminal,
      module == "Primary metabolism" ~ 1,
      module == "Secondary metabolism" ~ 2,
      module == "Microbial" ~ 3,
      TRUE ~ 2
    ),
    y = case_when(
      role == "Root" ~ 2.5,
      role == "Terminal" ~ 1.5,
      TRUE ~ 2.0
    )
  ) %>%
  group_by(x, y) %>%
  mutate(y = y + (row_number() - median(row_number())) * 0.25) %>%
  ungroup()

layout_df_mod <- layout_df_mod %>%
  mutate(
    x = ifelse(name %in% ancestors_p_ancestors, x + 0.3, x)
  )


# === 颜色设置 ===
layout_df_mod$fill_color <- case_when(
  layout_df_mod$name %in% ancestors_p_ancestors ~ "#F4A582",  # root
  layout_df_mod$name %in% terminal_nodes ~ "#92C5DE",         # terminal
  TRUE ~ "#f0f8ff"                                             # other
)

# === 构建 tidygraph 对象 ===
graph_tbl <- as_tbl_graph(g_ancestor) %>%
  activate(nodes) %>%
  left_join(layout_df_mod %>% select(name, x, y, fill_color), by = "name")

# === 计算边标签位置 ===
edge_label_df <- ancestor_edges %>%
  rowwise() %>%
  mutate(
    x = mean(layout_df_mod$x[layout_df_mod$name %in% c(from, to)]),
    y = mean(layout_df_mod$y[layout_df_mod$name %in% c(from, to)])
  )

# === 绘图 ===
p_ancestor_final <- ggraph(graph_tbl, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(edge_width = abs(coef), edge_color = coef, alpha = ifelse(abs(coef) < 0.1, 0.2, 1)),
    arrow = arrow(length = unit(3, "mm"), type = "closed", ends = "last"),
    end_cap = circle(5, "mm"),
    lineend = "round"
  ) +
  scale_edge_color_gradient2(
    low = "#313695", mid = "#f7f7f7", high = "#a50026", midpoint = 0,
    name = "Path Coef."
  ) +
  scale_edge_width(range = c(0.4, 2)) +
  scale_edge_alpha(guide = "none") +
  
  geom_text(
    data = edge_label_df %>%
      mutate(
        nudge_x = runif(n(), -0.05, 0.05),
        nudge_y = runif(n(), -0.05, 0.05),
        color = scales::col_numeric(
          palette = c("#313695", "#f7f7f7", "#a50026"),
          domain = c(min(ancestor_edges$coef), max(ancestor_edges$coef))
        )(coef)
      ),
    aes(x = x + nudge_x, y = y + nudge_y, label = round(coef, 2), color = I(color)),
    size = 4,
    family = "Times New Roman",
    show.legend = FALSE
  ) +
  
  geom_node_label(
    aes(label = name, fill = I(fill_color)),
    size = 4.6,
    label.size = 0,
    label.padding = unit(2.2, "mm"),
    label.r = unit(2.2, "mm"),
    color = "black",
    fontface = "bold",
    family = "Times New Roman"
  ) +
  
  theme_void() +
  ggtitle("Causal Paths from Ancestors to Terminal Latent Variables") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, family = "Times New Roman", face = "bold"),
    legend.title = element_text(size = 10, family = "Times New Roman"),
    legend.text = element_text(size = 9, family = "Times New Roman")
  )

# === 输出图形 ===
print(p_ancestor_final)









###总图Flow Betweenness（中介流动中心性）


library(sna)
library(igraph)
library(dplyr)

# === 从 PLS 模型提取路径结构 ===
path_coefs <- pls_model_final$path_coefs
edge_idx <- which(path_coefs != 0, arr.ind = TRUE)

# 构建边集
pls_edges <- data.frame(
  from = rownames(path_coefs)[edge_idx[, 1]],
  to = colnames(path_coefs)[edge_idx[, 2]],
  coef = path_coefs[edge_idx]
)

# 构建 igraph 对象（整个 inner model 图）
g_full <- graph_from_data_frame(pls_edges, directed = TRUE)

# 提取邻接矩阵（flowbet 使用邻接矩阵形式）
adj_matrix <- as_adjacency_matrix(g_full, attr = "coef", sparse = FALSE)

# === 执行 flowbet 分析 ===
flow_scores <- flowbet(adj_matrix)

# === 输出中心变量评分 ===
flowbet_global <- data.frame(
  variable = rownames(adj_matrix),
  flow_betweenness = flow_scores
) %>%
  arrange(desc(abs(flow_betweenness)))  # 按绝对值降序排列

# 显示结果
print(flowbet_global)









# ---- 1) 抽取 B/A（尽量复用现有对象；局部变量，不影响外部） ----
.get_BA <- function() {
  if (exists("pls_model_final") && !is.null(pls_model_final$path_coefs)) {
    B <- as.matrix(pls_model_final$path_coefs)
  } else stop("未找到 pls_model_final$path_coefs；请在 PLS 拟合后再运行本追加块。")
  
  if (is.null(rownames(B)) || is.null(colnames(B))) stop("B 需要行列名（潜变量名）。")
  B[is.na(B)] <- 0
  diag(B) <- 0
  
  if (exists("path_matrix_final")) {
    A <- as.matrix(path_matrix_final)
    A <- A[rownames(B), colnames(B)]
  } else {
    A <- (B != 0) + 0
    dimnames(A) <- dimnames(B)
  }
  diag(A) <- 0
  list(B = B, A = A)
}

# ---- 2) 计算 (I-B)^(-1) 与 T=(I-B)^(-1)-I；不可逆时用 Katz 衰减兜底 ----
.compute_T_inv <- function(B, alpha_max = 0.9) {
  n <- nrow(B); I <- diag(n); rownames(I) <- colnames(I) <- rownames(B)
  INV <- tryCatch(solve(I - B), error = function(e) NULL)
  if (!is.null(INV)) return(list(INV = INV, T = INV - I, alpha = 1))
  rho <- tryCatch(max(Mod(eigen(abs(B), only.values = TRUE)$values)), error = function(e) 0)
  alpha <- if (is.finite(rho) && rho > 0) min(alpha_max / rho, 0.99) else 0.9
  INV <- solve(I - alpha * B)
  list(INV = INV, T = INV - I, alpha = alpha)
}

# ---- 3) 自动识别终端与源头（如未显式提供） ----
.infer_endpoints <- function(A) {                         # 出度=0 → 终端
  od <- rowSums(A); names(od)[od == 0]
}
.detect_sources <- function(A, endpoints) {               # 入度=0 且非终端 → 源头
  indeg <- colSums(A)
  S <- setdiff(names(indeg[indeg == 0]), endpoints)
  if (!length(S)) S <- setdiff(colnames(A), endpoints)    # 兜底：非终端全体
  S
}

# ---- 4) IFC 主函数 ----
ifc_scores <- function(B, A, endpoints = NULL, exclude_endpoints = TRUE, alpha_max = 0.9){
  tmp <- .compute_T_inv(B, alpha_max); T <- tmp$T; alpha_used <- tmp$alpha
  
  # 终端集合：优先用你代码里已定义的对象
  E <- if (!is.null(endpoints) && length(endpoints)) {
    intersect(endpoints, colnames(B))
  } else if (exists("terminal_nodes")) {
    intersect(get("terminal_nodes"), colnames(B))
  } else if (exists("terminal_vars")) {
    intersect(get("terminal_vars"), colnames(B))
  } else {
    .infer_endpoints(A)
  }
  
  # 源头集合
  S <- .detect_sources(A, E)
  
  # 易感度：max_{s∈S} |T_{s→j}|；终端影响：sum_{e∈E} |T_{j→e}|
  ifc_in  <- if (length(S)) apply(abs(T[S, , drop = FALSE]), 2, max) else rep(0, nrow(B))
  ifc_out <- if (length(E)) rowSums(abs(T[, E, drop = FALSE])) else rep(NA_real_, nrow(B))
  names(ifc_in) <- names(ifc_out) <- colnames(B)
  
  df <- data.frame(
    module       = colnames(B),
    IFC_in_Linf  = as.numeric(ifc_in[colnames(B)]),
    IFC_out_L1   = as.numeric(ifc_out[colnames(B)]),
    IFC_flow     = as.numeric(ifc_in[colnames(B)] * ifc_out[colnames(B)]),
    alpha_used   = alpha_used,
    row.names    = NULL
  )
  if (exclude_endpoints && length(E)) df <- df[!(df$module %in% E), , drop = FALSE]
  df[order(-df$IFC_flow), , drop = FALSE]
}

# ---- 5) 执行并产出 IFC_table（不覆盖你已有变量） ----
BA <- .get_BA()
endpoints_for_ifc <- if (exists("terminal_nodes")) {
  intersect(terminal_nodes, colnames(BA$B))
} else if (exists("terminal_vars")) {
  intersect(terminal_vars, colnames(BA$B))
} else NULL

IFC_table <- ifc_scores(BA$B, BA$A, endpoints = endpoints_for_ifc, exclude_endpoints = TRUE)
print(head(IFC_table, 15))










# === 每个潜变量的回归图 ===

##自动识别最优模型
library(ggplot2)
library(dplyr)
library(patchwork)

# === 1. 读取并标准化数据 ===
dose_df <- read.csv("D:/01文件/001论文/博士论文/05毒性转录代谢/PLSPM/Level 1/DoseM.csv", check.names = FALSE)
# 从 pls_plspm_final 提取潜变量得分
p_lv_scores <- as.data.frame(pls_model_final$scores)
latent_scores_df <- p_lv_scores

# 标准化 Dose 向量（用于拟合），保留原始 Dose 用于横坐标展示
dose_vector_raw <- dose_df$Dose
dose_vector_scaled <- dose_df$Dose  # 转为 numeric 向量

# 标准化所有潜变量
latent_scores_scaled <- as.data.frame(latent_scores_df)

# === 2. 拟合函数（多模型） ===
fit_models <- function(x, y) {
  models <- list()
  metrics <- data.frame(Model = character(), R2 = numeric(), Coef = numeric(), P = numeric(), ModelP = numeric(), stringsAsFactors = FALSE)
  
  # Linear
  fit1 <- lm(y ~ x)
  models$linear <- fit1
  f1 <- summary(fit1)$fstatistic
  model_p1 <- pf(f1[1], f1[2], f1[3], lower.tail = FALSE)
  metrics <- rbind(metrics, data.frame(Model = "linear", R2 = summary(fit1)$r.squared,
                                       Coef = coef(fit1)[2], P = coef(summary(fit1))[2,4], ModelP = model_p1))
  
  # Quadratic
  fit2 <- lm(y ~ poly(x, 2, raw = TRUE))
  models$quadratic <- fit2
  f2 <- summary(fit2)$fstatistic
  model_p2 <- pf(f2[1], f2[2], f2[3], lower.tail = FALSE)
  metrics <- rbind(metrics, data.frame(Model = "quadratic", R2 = summary(fit2)$r.squared,
                                       Coef = coef(fit2)[3], P = coef(summary(fit2))[3,4], ModelP = model_p2))
  
  # Logarithmic
  if (all(x > 0)) {
    fit3 <- lm(y ~ log(x))
    models$log <- fit3
    f3 <- summary(fit3)$fstatistic
    model_p3 <- pf(f3[1], f3[2], f3[3], lower.tail = FALSE)
    metrics <- rbind(metrics, data.frame(Model = "log", R2 = summary(fit3)$r.squared,
                                         Coef = coef(fit3)[2], P = coef(summary(fit3))[2,4], ModelP = model_p3))
  }
  
  # Exponential
  if (all(y > 0)) {
    fit4 <- lm(log(y) ~ x)
    models$exp <- fit4
    f4 <- summary(fit4)$fstatistic
    model_p4 <- pf(f4[1], f4[2], f4[3], lower.tail = FALSE)
    metrics <- rbind(metrics, data.frame(Model = "exp", R2 = summary(fit4)$r.squared,
                                         Coef = coef(fit4)[2], P = coef(summary(fit4))[2,4], ModelP = model_p4))
  }
  
  # Power
  if (all(x > 0) && all(y > 0)) {
    fit5 <- lm(log(y) ~ log(x))
    models$power <- fit5
    f5 <- summary(fit5)$fstatistic
    model_p5 <- pf(f5[1], f5[2], f5[3], lower.tail = FALSE)
    metrics <- rbind(metrics, data.frame(Model = "power", R2 = summary(fit5)$r.squared,
                                         Coef = coef(fit5)[2], P = coef(summary(fit5))[2,4], ModelP = model_p5))
  }
  
  best_idx <- which.max(metrics$R2)
  best_model <- models[[metrics$Model[best_idx]]]
  
  return(list(
    best_model = best_model,
    best_type = metrics$Model[best_idx],
    best_coef = metrics$Coef[best_idx],
    best_p = metrics$P[best_idx],
    best_r2 = metrics$R2[best_idx],
    best_model_p = metrics$ModelP[best_idx]
  ))
}

# === 3. 批量拟合所有潜变量 ===
results <- lapply(colnames(latent_scores_scaled), function(var) {
  y <- latent_scores_scaled[[var]]
  fit_result <- fit_models(dose_vector_scaled, y)
  data.frame(
    LatentVar = var,
    BestModel = fit_result$best_type,
    Coef = fit_result$best_coef,
    R2 = fit_result$best_r2,
    p = fit_result$best_p,
    model_p = fit_result$best_model_p
  )
})
regression_df <- do.call(rbind, results)

# === 4. 显著性标注 ===
regression_df$Signif <- cut(regression_df$p,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                            labels = c("***", "**", "*", ".", ""))

# === 5. 加入原始 Dose 并排序绘图 ===
latent_scores_df$Dose <- dose_vector_raw

desired_order <- c(
  "Growth_Endpoint", "Toxicity_Endpoint", "Antibiotics_Residue_in_Algae",
  "Biodegradation", "Drug_Metabolism",
  "Auxin_Metabolism", "Proliferation", "Carbohydrate_Metabolism",
  "Lipid_Metabolism", "Amino_Acid_Metabolism",
  "Oxidative_Response", "Photosystem",
  "Microbial_Function", "Microbial_Diversity", "Microbial_Composition"
)

regression_df <- regression_df %>%
  mutate(LatentVar = factor(LatentVar, levels = desired_order)) %>%
  arrange(LatentVar)

# === 6. 绘图 ===
plot_list <- lapply(as.character(regression_df$LatentVar), function(var) {
  df_plot <- data.frame(
    Dose = dose_vector_raw,
    Score = latent_scores_scaled[[var]]
  )
  
  info <- regression_df %>% filter(LatentVar == var)
  best_model <- as.character(info$BestModel)
  model_display <- ifelse(best_model == "quadratic", "U-shaped", best_model)
  p_text_model <- ifelse(info$model_p < 0.01, "Model P < 0.01",
                         ifelse(info$model_p < 0.05, "Model P < 0.05", "Model P > 0.05"))
  
  p1 <- ggplot(df_plot, aes(x = Dose, y = Score)) +
    geom_point(size = 2) +
    labs(
      title = var,
      subtitle = paste(
        paste0(model_display, ", Coef = ", round(info$Coef, 4)),
        paste0("R² = ", round(info$R2, 2), ", ", p_text_model),
        sep = "\n"
      ),
      y = "Latent Score"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 14, face = "bold"),  # ⬅️ 坐标轴标题 
      axis.text  = element_text(size = 12, face = "bold")   # ⬅️ 坐标轴刻度
    )
  
  if (best_model == "linear") {
    p2 <- p1 + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "blue")
  } else if (best_model == "quadratic") {
    model_quad <- lm(Score ~ poly(Dose, 2, raw = TRUE), data = df_plot)
    coef_quad <- coef(model_quad)
    a <- coef_quad[3]
    b <- coef_quad[2]
    x_peak <- -b / (2 * a)
    y_peak <- predict(model_quad, newdata = data.frame(Dose = x_peak))
    
    # 默认曲线
    p2 <- p1 +
      geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = TRUE, color = "#5066a1")
    
    # 非特殊变量才绘制 EDmax
    if (var != "Microbial_Function") {
      p2 <- p2 +
        geom_point(aes(x = x_peak, y = y_peak), color = "red", size = 2.5) +
        annotate("text", x = x_peak, y = y_peak,
                 label = paste0("EDmax = ", round(x_peak, 2)),
                 vjust = -1.2, size = 3, color = "red")
    }
  } else if (best_model == "log") {
    p2 <- p1 + geom_smooth(method = "lm", formula = y ~ log(x), se = TRUE, color = "orange")
  } else if (best_model == "exp") {
    p2 <- p1 + geom_smooth(method = "lm", formula = log(y) ~ x, se = FALSE, color = "purple")
  } else if (best_model == "power") {
    p2 <- p1 + geom_smooth(method = "lm", formula = log(y) ~ log(x), se = FALSE, color = "red")
  }
  
  return(p2)
})

# === 输出拼图，每行 5 张图 ===
final_plot <- wrap_plots(plot_list, ncol = 5)
print(final_plot)



###Dose 回归系数图

# === 数据准备（确保是字符）===
regression_df_plot <- regression_df %>%
  mutate(LatentVar = as.character(LatentVar)) %>%
  arrange(Coef) %>%  # 按系数值排序
  mutate(LatentVar = factor(LatentVar, levels = LatentVar))  # 保持排序顺序

# === 散点图绘制 ===
ggplot(regression_df_plot, aes(x = LatentVar, y = Coef)) +
  geom_point(aes(color = Coef, size = abs(Coef))) +
  geom_text(aes(label = round(Coef, 2)), vjust = -1.5, size = 3.5) +
  scale_color_gradient2(low = "red", mid = "gray80", high = "blue", midpoint = 0) +
  scale_size(range = c(1, 8)) +
  labs(
    title = "Dose → Latent Variables",
    subtitle = "Each point represents the best-fit regression coefficient",
    x = "Latent Variable",
    y = "Regression Coefficient"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

library(dplyr)
library(igraph)

# === 1. Dose → 潜变量的边（从 regression_df）===
dose_edges <- regression_df %>%
  transmute(from = "Dose", to = LatentVar, coef = Coef)

# === 2. 提取 PLS-PM 路径系数边 ===
path_coefs <- pls_model_final$path_coefs
edge_idx <- which(path_coefs != 0, arr.ind = TRUE)

pls_edges <- data.frame(
  from = rownames(path_coefs)[edge_idx[, 1]],
  to = colnames(path_coefs)[edge_idx[, 2]],
  coef = path_coefs[edge_idx]
)

# === 3. 合并网络边 ===
all_edges <- bind_rows(dose_edges, pls_edges)


# 确保 from/to 为字符型
all_edges$from <- as.character(all_edges$from)
all_edges$to <- as.character(all_edges$to)

# 构建所有节点
vertex_df <- data.frame(name = unique(c(all_edges$from, all_edges$to)), stringsAsFactors = FALSE)

# 构建 igraph 图
g_full <- graph_from_data_frame(all_edges, vertices = vertex_df, directed = TRUE)


calc_effects <- function(target) {
  paths <- all_simple_paths(g_full, from = "Dose", to = target, mode = "out")
  
  total_effect <- 0
  indirect_effect <- 0
  direct_effect <- 0
  
  for (path in paths) {
    node_names <- names(path)
    if (length(node_names) < 2) next
    
    edge_pairs <- paste(head(node_names, -1), tail(node_names, -1))
    coefs <- sapply(edge_pairs, function(pair) {
      edge <- all_edges %>% filter(paste(from, to) == pair)
      if (nrow(edge) > 0) return(edge$coef) else return(1)
    })
    
    path_effect <- prod(coefs)
    
    if (length(node_names) == 2) {
      direct_effect <- path_effect
    } else {
      indirect_effect <- indirect_effect + path_effect
    }
    
    total_effect <- total_effect + path_effect
  }
  
  return(data.frame(
    LatentVar = target,
    Direct = direct_effect,
    Indirect = indirect_effect,
    Total = total_effect
  ))
}

latent_vars <- colnames(pls_model_final$path_coefs)
effect_df <- do.call(rbind, lapply(latent_vars, calc_effects))


library(ggplot2)
effect_df_long <- tidyr::pivot_longer(effect_df, cols = c("Direct", "Indirect", "Total"),
                                      names_to = "EffectType", values_to = "Value")

ggplot(effect_df_long, aes(x = LatentVar, y = Value, fill = EffectType)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Direct" = "#E67E22", "Indirect" = "gray70", "Total" = "#5066a1")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Dose Effects on Latent Variables",
    x = "Latent Variable",
    y = "Effect Value",
    fill = "Effect Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )



###level 1→level 2→显变量；潜→显变量   观测模型


#########
library(tidyverse)
library(ggraph)
library(igraph)
library(dplyr)
library(tibble)

#提取观测模型
outer_L1 <- pls_model_final$outer_model %>%
  as.data.frame() %>%
  rownames_to_column(var = "observed") %>%     # 显变量名
  rename(latent = block, weight = weight) %>%  # 潜变量名和权重
  select(latent, observed, weight)

outer_L2 <- l2_pls_model$outer_model %>%
  as.data.frame() %>%
  rownames_to_column(var = "observed") %>%     # 显变量名
  rename(latent = block, weight = weight) %>%  # 潜变量名和权重
  select(latent, observed, weight)


# === 1. 读取 outer_model 数据 ===
df_l1 <- outer_L1
df_l2 <- outer_L2

# 标准化列名
colnames(df_l1) <- c("latent", "observed", "weight")
colnames(df_l2) <- c("latent", "observed", "weight")

# 构建边表
edges_l1_l2 <- df_l1 %>%
  rename(from = latent, to = observed, weight12 = weight)

# Level2 → 显变量（Top 5）
edges_l2_obs <- df_l2 %>%
  group_by(latent) %>%
  arrange(desc(abs(weight))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  rename(from = latent, to = observed, weight23 = weight)

# 所有的 level1 潜变量
level1_list <- unique(df_l1$latent)

# === 通用绘图函数 ===
plot_three_layer <- function(level1_latent) {
  # 获取 Level2 潜变量（即 Level1 的显变量）
  l2_nodes <- edges_l1_l2 %>% filter(from == level1_latent) %>% pull(to)
  
  # Level1 → Level2 边
  edges1 <- edges_l1_l2 %>% filter(from == level1_latent)
  
  # Level2 → 显变量边（Top 5）
  edges2 <- edges_l2_obs %>% filter(from %in% l2_nodes)
  
  # 合并所有边
  edges_all <- bind_rows(
    edges1 %>% select(from, to, weight = weight12),
    edges2 %>% select(from, to, weight = weight23)
  )
  
  # 构建节点并分层
  nodes <- data.frame(name = unique(c(edges_all$from, edges_all$to))) %>%
    mutate(layer = case_when(
      name == level1_latent ~ 1,
      name %in% edges1$to ~ 2,
      name %in% edges2$to ~ 3
    ))
  
  # 布局：压缩横向空间（可调）
  layer_spacing <- c(0.5, 1.5, 2.5)
  nodes <- nodes %>%
    group_by(layer) %>%
    mutate(
      x = layer_spacing[layer],
      y = seq(-1, -n(), length.out = n())
    ) %>%
    ungroup()
  
  # 连接边：默认中心点坐标
  edges_all <- edges_all %>%
    left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x_end = x, y_end = y)
  
  # ==== 加入横向边缘偏移，使线条连接变量框的边缘 ====
  x_offset <- 0.15  # 推荐范围 0.1–0.2，太大线会错开
  
  # 标红最大路径
  edges_all <- edges_all %>%
    group_by(from) %>%
    mutate(is_max = abs(weight) == max(abs(weight))) %>%
    ungroup()
  
  # Level2 潜变量的显变量数量注释
  latent_annot_df_l2 <- df_l2 %>%
    group_by(latent) %>%
    summarise(n_total = n()) %>%
    filter(latent %in% edges2$from) %>%
    left_join(nodes, by = c("latent" = "name")) %>%
    mutate(label = paste0("Total: ", n_total),
           x = x,
           y = y - 0.4)
  
  # === 绘图 ===
  ggplot() +
    # 边
    geom_segment(data = edges_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
                     color = is_max),
                 linewidth = 1.3) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "#999999")) +
    
    # 路径权重标签
    geom_text(data = edges_all,
              aes(x = (x_start + x_end)/2,
                  y = (y_start + y_end)/2,
                  label = sprintf("%.3f", weight)),
              size = 4.2, fontface = "italic", color = "black") +
    
    # Level2 注释
    geom_text(data = latent_annot_df_l2,
              aes(x = x, y = y, label = label),
              size = 4.2, fontface = "italic", color = "black") +
    
    # 节点
    geom_label(data = nodes,
               aes(x = x, y = y, label = name, fill = factor(layer)),
               color = "black", size = 4.2, label.size = 0.6) +
    
    # 配色与主题
    scale_fill_manual(values = c("lightblue", "#FF6633", "#a3d393")) +
    theme_void() +
    ggtitle(paste("Three-layer Structure for", level1_latent)) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Times New Roman"))
}



# 绘制所有 Level1 潜变量图
plots <- lapply(level1_list, plot_three_layer)
for (p in plots) print(p)

# 设置保存目录
#output_dir <- "D:/01文件/001论文/博士论文/05毒性转录代谢/Rplot/PLSPM 计算mol/Full Outer"

# 如果目录不存在，则创建
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 自定义保存函数（PDF）
save_three_layer_plot_pdf <- function(plot_obj, level1_latent, width_mm = 210, height_mm = 210, dpi = 900) {
  ggsave(
    filename = file.path(output_dir, paste0("three_layer_structure_", level1_latent, ".pdf")),
    plot = plot_obj,
    device = cairo_pdf,              # 使用高质量 PDF 引擎
    width = width_mm,
    height = height_mm,
    dpi = dpi,
    units = "mm",                   # 使用毫米作为单位
    limitsize = FALSE
  )
}

# 查看并保存所有潜变量图
plots <- lapply(level1_list, function(name) {
  p <- plot_three_layer(name)             # 绘图
  print(p)                                # 可视化
  save_three_layer_plot_pdf(p, name)      # 保存
  p
})



###重要因果 主链
library(igraph)
library(dplyr)




# 构建节点表（防止遗漏）
vertex_df <- data.frame(
  name = unique(c(all_edges$from, all_edges$to)),
  stringsAsFactors = FALSE
)

# 构建 igraph 图对象
g_all <- graph_from_data_frame(all_edges, vertices = vertex_df, directed = TRUE)

# 添加边属性
E(g_all)$weight <- abs(E(g_all)$coef)   # 正值表示路径强度

# 终端变量列表
terminal_vars <- c(
  "Growth_Endpoint",
  "Toxicity_Endpoint",
  "Antibiotics_Residue_in_Algae",
  "Biodegradation"
)

# 路径提取函数
find_dose_path <- function(target) {
  sp <- shortest_paths(
    g_all,
    from = "Dose",
    to = target,
    weights = -E(g_all)$weight,  # 最大累计影响力路径（绝对值最大）
    output = "both"
  )
  
  node_path <- names(sp$vpath[[1]])
  
  # 提取路径边：仅保留路径内相邻节点对
  path_edges <- all_edges %>%
    filter(paste(from, to) %in% paste(head(node_path, -1), tail(node_path, -1)))
  
  list(
    target = target,
    path_nodes = node_path,
    path_edges = path_edges
  )
}
# 执行路径提取
dose_paths <- lapply(terminal_vars, find_dose_path)
names(dose_paths) <- terminal_vars



library(dplyr)
library(igraph)

# == Step 1: 构建函数：提取 Dose → terminal 所有路径及乘积 ==
get_all_paths_with_effects <- function(g, all_edges, target) {
  paths <- all_simple_paths(g, from = "Dose", to = target, mode = "out")
  
  path_df <- data.frame()
  
  for (path in paths) {
    node_names <- names(path)
    if (length(node_names) < 2) next
    
    edge_pairs <- paste(head(node_names, -1), tail(node_names, -1))
    coefs <- sapply(edge_pairs, function(pair) {
      edge <- all_edges %>% filter(paste(from, to) == pair)
      if (nrow(edge) > 0) return(edge$coef) else return(1)
    })
    
    path_effect <- prod(coefs)
    
    path_df <- bind_rows(path_df, data.frame(
      target = target,
      path = paste(node_names, collapse = " → "),
      n_nodes = length(node_names),
      effect = path_effect
    ))
  }
  return(path_df)
}
# Step 2: 提取所有路径效果
all_path_results <- lapply(terminal_vars, function(tar) {
  get_all_paths_with_effects(g_all, all_edges, tar)
})
all_path_df <- bind_rows(all_path_results)

# Step 3: 找每个终端变量的最大路径效应
key_paths <- all_path_df %>%
  group_by(target) %>%
  filter(abs(effect) == max(abs(effect))) %>%
  ungroup()

print(key_paths)

# 计算并输出每个终端变量的总效应系数
main_path <- lapply(dose_paths, function(res) {
  total_effect <- prod(res$path_edges$coef)
  
  # 打印输出（可选）
  cat("\n=== Dose →", res$target, "的关键路径 ===\n")
  cat("路径节点：", paste(res$path_nodes, collapse = " → "), "\n")
  cat("路径边（含系数）：\n")
  print(res$path_edges)
  cat("总效应系数（Dose →", res$target, "）:", round(total_effect, 4), "\n")
  
  # 返回结构化结果
  list(
    target = res$target,
    path = paste(res$path_nodes, collapse = " → "),
    edges = res$path_edges,
    total_effect = total_effect
  )
})
main_path









##寻找最重要间接效应链
##提取所有 Dose → Tn 的“间接路径”
get_indirect_paths <- function(g, all_edges, target) {
  all_paths <- all_simple_paths(g, from = "Dose", to = target, mode = "out")
  
  indirect_df <- data.frame()
  
  for (p in all_paths) {
    nodes <- names(p)
    if (length(nodes) <= 2) next  # 跳过直接路径
    
    edge_pairs <- paste(head(nodes, -1), tail(nodes, -1))
    coefs <- sapply(edge_pairs, function(pair) {
      edge <- all_edges %>% filter(paste(from, to) == pair)
      if (nrow(edge) > 0) return(edge$coef) else return(1)
    })
    
    indirect_df <- bind_rows(indirect_df, data.frame(
      target = target,
      path = paste(nodes, collapse = " → "),
      n_nodes = length(nodes),
      effect = prod(coefs)
    ))
  }
  return(indirect_df)
}


##对每个终端变量提取最强间接路径
indirect_paths_all <- lapply(terminal_vars, function(tar) {
  get_indirect_paths(g_all, all_edges, tar)
})
indirect_df <- bind_rows(indirect_paths_all)

# 找最大效应路径
key_indirect_paths <- indirect_df %>%
  group_by(target) %>%
  filter(abs(effect) == max(abs(effect))) %>%
  ungroup()


##删除这条路径后，重新统计全部间接路径效应之和（总间接效应）
get_total_indirect_effect <- function(target, g, all_edges) {
  paths <- all_simple_paths(g, from = "Dose", to = target, mode = "out")
  effs <- sapply(paths, function(p) {
    nodes <- names(p)
    if (length(nodes) <= 2) return(0)
    edge_pairs <- paste(head(nodes, -1), tail(nodes, -1))
    coefs <- sapply(edge_pairs, function(pair) {
      edge <- all_edges %>% filter(paste(from, to) == pair)
      if (nrow(edge) > 0) return(edge$coef) else return(1)
    })
    prod(coefs)
  })
  sum(effs)
}


#模拟断链，检验间接效应是否明显下降
analyze_key_path_impact <- function(path_row) {
  nodes <- strsplit(path_row$path, " → ")[[1]]
  edge_list <- data.frame(
    from = head(nodes, -1),
    to = tail(nodes, -1)
  )
  
  # 删除路径边
  edges_removed <- anti_join(all_edges, edge_list, by = c("from", "to"))
  g_removed <- graph_from_data_frame(edges_removed, directed = TRUE)
  
  before <- get_total_indirect_effect(path_row$target, g_all, all_edges)
  after <- get_total_indirect_effect(path_row$target, g_removed, edges_removed)
  
  data.frame(
    target = path_row$target,
    key_path = path_row$path,
    key_path_effect = path_row$effect,
    total_indirect_before = before,
    total_indirect_after = after,
    percent_drop = round(1 - after / before, 3)
  )
}

impact_results <- key_indirect_paths %>%
  rowwise() %>%
  do(analyze_key_path_impact(.))

print(impact_results)






# === 构建 Dose → 中间潜变量 → 终端变量 所有有效路径结构 ===

# 设置终端变量
terminal_vars <- c(
  "Growth_Endpoint",
  "Toxicity_Endpoint",
  "Antibiotics_Residue_in_Algae",
  "Biodegradation"
)

# 初始化结果表
all_paths_df <- data.frame()

# 构建图对象
g_all <- igraph::graph_from_data_frame(all_edges, vertices = vertex_df, directed = TRUE)

# 递归函数：查找所有路径
for (target in terminal_vars) {
  all_paths <- igraph::all_simple_paths(g_all, from = "Dose", to = target)
  
  for (p in all_paths) {
    node_path <- names(p)
    
    # 提取路径边（边对）
    edge_table <- data.frame(
      From = head(node_path, -1),
      To = tail(node_path, -1),
      stringsAsFactors = FALSE
    )
    
    # 添加系数
    edge_table <- dplyr::left_join(edge_table, all_edges, by = c("From" = "from", "To" = "to"))
    
    # 增加路径整体信息
    edge_table$Target <- target
    edge_table$FullPath <- paste(node_path, collapse = " → ")
    edge_table$ProductContribution <- prod(edge_table$coef)
    edge_table$SumContribution <- sum(edge_table$coef)
    
    # 合并入总表
    all_paths_df <- dplyr::bind_rows(all_paths_df, edge_table)
  }
}

all_paths_df






###level 1→level 2→显变量；潜→显变量   观测模型


#########
library(tidyverse)
library(ggraph)
library(igraph)
library(dplyr)
library(tibble)

#提取观测模型
outer_L1 <- pls_model_final$outer_model %>%
  as.data.frame() %>%
  rownames_to_column(var = "observed") %>%     # 显变量名
  rename(latent = block, weight = weight) %>%  # 潜变量名和权重
  select(latent, observed, weight)

outer_L2 <- l2_pls_model$outer_model %>%
  as.data.frame() %>%
  rownames_to_column(var = "observed") %>%     # 显变量名
  rename(latent = block, weight = weight) %>%  # 潜变量名和权重
  select(latent, observed, weight)


# === 1. 读取 outer_model 数据 ===
df_l1 <- outer_L1
df_l2 <- outer_L2

# 标准化列名
colnames(df_l1) <- c("latent", "observed", "weight")
colnames(df_l2) <- c("latent", "observed", "weight")

# 构建边表
edges_l1_l2 <- df_l1 %>%
  rename(from = latent, to = observed, weight12 = weight)

# Level2 → 显变量（Top 5）
edges_l2_obs <- df_l2 %>%
  group_by(latent) %>%
  arrange(desc(abs(weight))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  rename(from = latent, to = observed, weight23 = weight)

# 所有的 level1 潜变量
level1_list <- unique(df_l1$latent)

# === 通用绘图函数 ===
plot_three_layer <- function(level1_latent) {
  # 获取所有 Level2 节点（即 Level1 的 observed）
  l2_nodes <- edges_l1_l2 %>% filter(from == level1_latent) %>% pull(to)
  
  # 构建完整的 Level1 → Level2
  edges1 <- edges_l1_l2 %>% filter(from == level1_latent)
  
  # 构建 Level2 → 显变量（限制 top5）
  edges2 <- edges_l2_obs %>% filter(from %in% l2_nodes)
  
  # 合并边
  edges_all <- bind_rows(
    edges1 %>% select(from, to, weight = weight12),
    edges2 %>% select(from, to, weight = weight23)
  )
  
  # 构建节点并分层
  nodes <- data.frame(name = unique(c(edges_all$from, edges_all$to))) %>%
    mutate(layer = case_when(
      name == level1_latent ~ 1,
      name %in% edges1$to ~ 2,
      name %in% edges2$to ~ 3
    ))
  
  # 保留所有节点
  nodes <- nodes %>%
    group_by(layer) %>%
    mutate(x = layer - 1, y = seq(-1, -n(), length.out = n()))
  
  # 位置连接
  edges_all <- edges_all %>%
    left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
    rename(x_start = x, y_start = y) %>%
    left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x_end = x, y_end = y)
  
  # 标红最大值
  edges_all <- edges_all %>%
    group_by(from) %>%
    mutate(is_max = abs(weight) == max(abs(weight))) %>%
    ungroup()
  
  latent_annot_df_l2 <- df_l2 %>%
    group_by(latent) %>%
    summarise(n_total = n()) %>%
    filter(latent %in% edges2$from) %>%  # 仅保留当前图中有展示的 level2 潜变量
    left_join(nodes, by = c("latent" = "name")) %>%
    mutate(label = paste0("Total: ", n_total),
           x = x,
           y = y - 0.4)  
  # 绘图
  ggplot() +
    # 边线（红色高亮 + 灰色）
    geom_segment(data = edges_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end,
                     color = is_max),
                 linewidth = 1.5) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "gray60")) +
    
    # 添加路径系数（保留2位小数）
    geom_text(data = edges_all,
              aes(x = (x_start + x_end)/2,
                  y = (y_start + y_end)/2,
                  label = sprintf("%.3f", weight)),
              size = 5, fontface = "italic", color = "black") +
    
    # 添加 Level2 的观测变量数量注释
    geom_text(data = latent_annot_df_l2,
              aes(x = x, y = y, label = label),
              size = 5, fontface = "italic", color = "black") +
    
    # 节点绘制
    geom_label(data = nodes,
               aes(x = x, y = y, label = name, fill = factor(layer)),
               color = "black", size = 5, label.size = 0.75) +
    
    # 配色与主题
    scale_fill_manual(values = c("lightblue", "#FF6633", "#a3d393")) +
    theme_void() +
    ggtitle(paste("Three-layer Structure for", level1_latent)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# 绘制所有 Level1 潜变量图
plots <- lapply(level1_list, plot_three_layer)
for (p in plots) print(p)




###子图Flow Betweenness（中介流动中心性）

# 终端变量列表
terminal_vars <-c(
  "Growth_Endpoint",
  "Toxicity_Endpoint",
  "Antibiotics_Residue_in_Algae",
  "Biodegradation"
)

collect_upstream <- function(node, g_rev, visited = character()) {
  if (node %in% visited) return(visited)
  parents <- names(neighbors(g_rev, node, mode = "out"))
  if (length(parents) == 0) return(c(visited, node))
  visited <- c(visited, node)
  for (p in parents) {
    visited <- union(visited, collect_upstream(p, g_rev, visited))
  }
  return(visited)
}

# 构建反向图一次即可
g_rev <- igraph::reverse_edges(g_pls)

# 提取子图 + 计算 flowbet() 的函数
# 修改函数：增加 g_rev 参数
get_flowbet_for_terminal <- function(terminal, graph, g_rev) {
  upstream_nodes <- collect_upstream(terminal, g_rev)  # 递归找上游
  g_sub <- induced_subgraph(graph, vids = upstream_nodes)
  adj <- as_adjacency_matrix(g_sub, attr = "coef", sparse = FALSE)
  flow_scores <- flowbet(adj)
  
  data.frame(
    terminal = terminal,
    variable = rownames(adj),
    flow_betweenness = flow_scores
  )
}

# 批量计算 flowbet
flowbet_all <- lapply(terminal_vars, function(tar) {
  get_flowbet_for_terminal(tar, g_pls, g_rev)
}) %>% bind_rows()

# 排序查看
flowbet_all <- flowbet_all %>%
  mutate(abs_flow = abs(flow_betweenness)) %>%
  arrange(terminal, desc(abs_flow)) %>%
  select(-abs_flow)  # 排序后移除临时列

print(flowbet_all)

###重要因果 结构
# === 只使用 PLS 路径（不包含 Dose） ===

# 加载包
library(igraph)
library(dplyr)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(purrr)

# 构建 PLS 路径边数据
path_coefs <- pls_model_final$path_coefs
edge_idx <- which(path_coefs != 0, arr.ind = TRUE)

pls_edges <- data.frame(
  from = rownames(path_coefs)[edge_idx[, 1]],
  to = colnames(path_coefs)[edge_idx[, 2]],
  coef = path_coefs[edge_idx]
)

# 构建网络图对象
vertex_df <- data.frame(name = unique(c(pls_edges$from, pls_edges$to)))
g_pls <- graph_from_data_frame(pls_edges, vertices = vertex_df, directed = TRUE)
E(g_pls)$weight <- abs(E(g_pls)$coef)

# 构建反向图，用于递归向上追踪因果路径
g_rev <- igraph::reverse_edges(g_pls)

# 递归函数：收集所有祖先节点
collect_upstream <- function(node, visited = character()) {
  if (node %in% visited) return(visited)
  parents <- names(neighbors(g_rev, node, mode = "out"))
  if (length(parents) == 0) return(c(visited, node))
  visited <- c(visited, node)
  for (p in parents) {
    visited <- union(visited, collect_upstream(p, visited))
  }
  return(visited)
}

# 终端变量列表
terminal_vars <- c(
  "Growth_Endpoint",
  "Toxicity_Endpoint",
  "Antibiotics_Residue_in_Algae",
  "Biodegradation"
)

# 绘图函数：递归可达的所有路径子图
plot_upstream_recursive <- function(target) {
  # === Step 1: 获取目标变量的所有上游节点 ===
  upstream_nodes <- collect_upstream(target)  # 递归收集能最终指向 target 的所有变量
  g_sub <- induced_subgraph(g_pls, vids = upstream_nodes)  # 提取子图
  g_tbl <- as_tbl_graph(g_sub)  # 转为 tidygraph 格式，便于 ggraph 可视化
  
  # === Step 2: 绘图开始 ===
  ggraph(g_tbl, layout = "fr") +  # 使用 fruchterman-reingold 力导布局，更适合区分交叉边
    
    # === 边绘制 ===
    geom_edge_link2(
      aes(
        label = round(coef, 2),      # 显示边上的路径系数，保留2位小数
        edge_width = abs(coef)       # 边的粗细与系数绝对值相关
      ),
      label_size = 3.8,              # 边上标签的字体大小
      label_dodge = unit(4.5, "mm"), # 标签与边之间的垂直偏移，避免重叠
      arrow = arrow(length = unit(3.5, "mm"), type = "closed"),  # 设定箭头长度
      end_cap = circle(8, 'pt'),     # 控制箭头与节点框之间的距离
      edge_colour = "grey40"         # 边的颜色
    ) +
    
    scale_edge_width(range = c(0.4, 1.2)) +  # 设置边粗细的范围（最细 ~ 最粗）
    
    # === 节点绘制 ===
    geom_node_label(
      aes(label = name),            # 节点标签使用变量名
      fill = "lightblue",           # 节点背景颜色
      color = "black",              # 节点文字颜色
      label.padding = unit(5, "pt"),# 标签内边距（控制文字与边框距离）
      label.size = 0.3,             # 节点边框粗细
      size = 5                      # 节点文字字体大小
    ) +
    
    # === 图整体样式 ===
    theme_void() +                  # 去除背景、坐标轴等
    labs(title = paste("Upstream Recursive →", target)) +  # 图标题
    coord_cartesian(clip = "off") + # 禁止坐标系裁剪，避免边缘被切掉
    theme(
      plot.margin = margin(50, 50, 50, 50),        # 图四周的边距，避免遮挡
      plot.title = element_text(hjust = 0.5, size = 16)  # 标题居中与字号
    )
}

plot_upstream_recursive <- function(target) {
  # 1. 从终变量递归收集所有上游变量（包含它们之间的路径）
  upstream_nodes <- collect_upstream(target)
  
  # 2. 从原图中提取对应子图（包含所有相关节点和边）
  g_sub <- induced_subgraph(g_pls, vids = upstream_nodes)
  g_tbl <- as_tbl_graph(g_sub)
  
  # 3. 绘图：使用 fr（force-directed）布局，更有利于节点自动分散
  ggraph(g_tbl, layout = "fr") +
    
    # 4. 边：路径方向+标签+粗细+颜色均反映路径系数大小
    geom_edge_link2(
      aes(
        label = round(coef, 2),         # 显示路径系数
        edge_width = abs(coef),         # 系数越大，边越粗
        edge_colour = abs(coef)         # 系数越大，颜色越深
      ),
      label_size = 3.8,                 # 标签字体大小
      label_dodge = unit(12, "mm"),    # 标签与边的距离
      arrow = arrow(length = unit(3.5, "mm"), type = "closed"),  # 实心箭头
      end_cap = circle(5, 'pt')         # 箭头和节点框之间留白（防止重叠）
    ) +
    
    # 5. 设置边宽范围
    scale_edge_width(range = c(0.4, 1.6)) +
    
    # 6. 设置边颜色渐变：颜色越深表示路径重要性越高
    scale_edge_colour_gradient(
      low = "lightblue", high = "red"
    ) +
    
    # 7. 节点标签样式
    geom_node_label(
      aes(label = name),
      fill = "lightblue",      # 标签背景色
      color = "black",         # 字体颜色
      label.padding = unit(5, "pt"),  # 标签边距
      label.size = 0.3,        # 标签边框线宽
      size = 5                 # 字体大小
    ) +
    
    # 8. 取消背景，保留最简风格
    theme_void() +
    
    # 9. 图标题
    labs(title = paste("Upstream Recursive →", target)) +
    
    # 10. 允许超出图区域的标签不被裁剪
    coord_cartesian(clip = "off") +
    
    # 11. 设置图边距和标题样式
    theme(
      plot.margin = margin(50, 50, 50, 50),         # 上右下左留白
      plot.title = element_text(hjust = 0.5, size = 16)  # 居中标题
    )
}

# 单独绘图：每个终端变量分别调用
plot_upstream_recursive("Growth_Endpoint")
plot_upstream_recursive("Toxicity_Endpoint")
plot_upstream_recursive("Antibiotics_Residue_in_Algae")


summary(pls_model_final)

