#MaxQuant 下游 R 分析

# 1. 加载包
library(tidyverse)
library(pheatmap)

# 2. 读取文件
mq <- read.delim(
  "C:\\Users\\叶子欣\\Desktop\\proteomics_assignment\\raw\\combined\\txt\\proteinGroups.txt",
  stringsAsFactors = FALSE,
  sep = "\t",
  check.names = FALSE
)

# 3. 数据清洗
mq_clean <- mq

# 自动识别蛋白ID列
protein_id_col <- ifelse("Protein IDs" %in% colnames(mq_clean), 
                         "Protein IDs", 
                         "Protein.IDs")

# 过滤反向库
if ("Reverse" %in% colnames(mq_clean)) {
  mq_clean <- mq_clean %>% filter(Reverse != "+")
} else if ("Reverse." %in% colnames(mq_clean)) {
  mq_clean <- mq_clean %>% filter(Reverse. != "+")
}

# 过滤污染蛋白
if ("Potential.contaminant" %in% colnames(mq_clean)) {
  mq_clean <- mq_clean %>% filter(Potential.contaminant != "+")
} else if ("Potential.contaminant." %in% colnames(mq_clean)) {
  mq_clean <- mq_clean %>% filter(Potential.contaminant. != "+")
}

# 过滤低可信度蛋白（肽段数≥2）
if ("Peptides" %in% colnames(mq_clean)) {
  mq_clean <- mq_clean %>% filter(Peptides >= 2)
}

# 把蛋白ID设为行名
mq_clean <- mq_clean %>% column_to_rownames(protein_id_col)

# 4. 提取并清洗 LFQ 定量信息
lfq_cols <- grep("LFQ.intensity", colnames(mq_clean), value = TRUE)
expr <- as.matrix(mq_clean[, lfq_cols])  # 强制转为矩阵，避免list类型报错
expr_log <- log2(expr + 1)

# 5. 缺失值处理
expr_log[is.na(expr_log)] <- 0
expr_log[is.infinite(expr_log)] <- 0

# 过滤掉无效蛋白
expr_log <- expr_log[rowSums(expr_log != 0) >= 2, ]

# 6. 画 PCA 图
pca <- prcomp(t(expr_log), scale = FALSE)
pc_df <- as.data.frame(pca$x[, 1:2])
pc_df$sample <- rownames(pc_df)

ggplot(pc_df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 4, color = "firebrick") +
  geom_text(vjust = -1, size = 3.5) +
  labs(title = "PCA Plot (Label-Free Quantification)",
       x = paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100, 1), "%)"),
       y = paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100, 1), "%)")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

# 7. 画热图
# 先不做行标准化，直接画原始对数数据的热图
pheatmap(
  top50,
  scale = "none",  # 关键：去掉行标准化，避免产生NaN
  show_rownames = FALSE,
  treeheight_row = 15,
  treeheight_col = 20,
  main = "Top 50 Most Variable Proteins (No Scaling)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

# 手动对top50做行标准化
top50_scaled <- t(scale(t(top50), center = TRUE, scale = TRUE))

# 检查一遍有没有NaN
any(is.na(top50_scaled))

# 画热图
pheatmap(
  top50_scaled,
  scale = "none",  # 已经手动标准化过了，这里不用再scale
  show_rownames = FALSE,
  treeheight_row = 15,
  treeheight_col = 20,
  main = "Top 50 Most Variable Proteins (Scaled)",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)

pheatmap(
  top50,
  scale = "none",         
  cluster_rows = TRUE,    
  cluster_cols = FALSE,  
  show_rownames = FALSE,
  treeheight_row = 15,
  main = "Top 50 Most Variable Proteins"
)


pheatmap(
  top50,
  scale = "none",
  cluster_rows = FALSE,   # 不聚类行
  cluster_cols = FALSE,   # 不聚类列
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = "Top 50 Most Variable Proteins"
)


library(tidyverse)

# 读取数据
mq <- read.delim(
  "C:\\Users\\叶子欣\\Desktop\\proteomics_assignment\\raw\\combined\\txt\\proteinGroups.txt",
  stringsAsFactors = FALSE,
  sep = "\t",
  check.names = FALSE
)

# 清洗：只保留肽段数≥2的蛋白（去掉污染/反向库的过滤，避免列名不匹配）
mq_clean <- mq %>% filter(Peptides >= 2)

# 提取3个样本的LFQ强度
expr <- mq_clean %>% select(starts_with("LFQ intensity"))
colnames(expr) <- c("Control", "Treat1", "Treat2")

# 对数化
expr_log <- log2(expr + 1)

# 计算处理组均值和log2FC
expr_log$Treat_mean <- rowMeans(expr_log[, 2:3], na.rm = TRUE)
expr_log$log2FC <- expr_log$Treat_mean - expr_log$Control

# 标记差异蛋白（|log2FC| ≥ 1）
expr_log$sig <- case_when(
  expr_log$log2FC >= 1 ~ "Up",
  expr_log$log2FC <= -1 ~ "Down",
  TRUE ~ "NS"
)

# 输出差异数量
cat("上调蛋白数：", sum(expr_log$sig == "Up"), "\n")
cat("下调蛋白数：", sum(expr_log$sig == "Down"), "\n")

# 画火山图
ggplot(expr_log, aes(x = log2FC, y = 0, color = sig)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  labs(title = "Volcano Plot (Control vs Treatments)", x = "log2(Fold Change)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# 导出差异蛋白表格
diff_result <- expr_log %>%
  filter(sig != "NS") %>%
  mutate(ProteinID = mq_clean$`Protein IDs`) %>%
  select(ProteinID, Control, Treat1, Treat2, Treat_mean, log2FC, sig)

write.csv(diff_result, "差异蛋白_最终结果.csv", row.names = FALSE)


# ===================== 差异分析+火山图 =====================
library(tidyverse)

# 1. 读取数据
mq <- read.delim(
  "C:\\Users\\叶子欣\\Desktop\\proteomics_assignment\\raw\\combined\\txt\\proteinGroups.txt",
  stringsAsFactors = FALSE,
  sep = "\t",
  check.names = FALSE
)

# 2. 清洗数据（仅保留肽段数≥2的蛋白，避免列名不匹配）
mq_clean <- mq %>% filter(Peptides >= 2)

# 3. 提取 3 个样本的 LFQ 强度
expr <- mq_clean %>% select(starts_with("LFQ intensity"))
colnames(expr) <- c("Control", "Treat1", "Treat2")

# 4. 对数化处理
expr_log <- log2(expr + 1)

# 5. 计算差异（处理组均值 vs 对照组）
expr_log$Treat_mean <- rowMeans(expr_log[, 2:3], na.rm = TRUE)
expr_log$log2FC <- expr_log$Treat_mean - expr_log$Control

# 6. 标记差异蛋白（|log2FC| ≥ 1）
expr_log$sig <- case_when(
  expr_log$log2FC >= 1 ~ "Up",
  expr_log$log2FC <= -1 ~ "Down",
  TRUE ~ "NS"
)

# 7. 合并蛋白ID到 expr_log
expr_log$ProteinID <- mq_clean$`Protein IDs`

# 8. 统计差异蛋白数量
cat("上调蛋白数：", sum(expr_log$sig == "Up"), "\n")
cat("下调蛋白数：", sum(expr_log$sig == "Down"), "\n")

# 9. 画火山图
ggplot(expr_log, aes(x = log2FC, y = 0, color = sig)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  labs(title = "Volcano Plot (Control vs Treatments)", x = "log2(Fold Change)") +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# 10. 导出差异蛋白表格
diff_result <- expr_log %>%
  filter(sig != "NS") %>%
  select(ProteinID, Control, Treat1, Treat2, Treat_mean, log2FC, sig)

write.csv(diff_result, "差异蛋白_最终结果.csv", row.names = FALSE)

getwd()

# 加载包
library(clusterProfiler)
library(org.Hs.eg.db)  # 人源蛋白注释库，如果是其他物种替换成对应的库
library(ggplot2)
library(dplyr)
library(stringr)

# 读取生成的差异蛋白表格
diff_proteins <- read.csv("差异蛋白_最终结果.csv", stringsAsFactors = FALSE)

# 处理蛋白ID：拆分Uniprot ID，取第一个有效ID
diff_proteins$Uniprot_ID <- str_split(diff_proteins$ProteinID, ";", simplify = TRUE)[, 1]

# 提取上调、下调蛋白的Uniprot ID
up_proteins <- diff_proteins %>% filter(sig == "Up") %>% pull(Uniprot_ID)
down_proteins <- diff_proteins %>% filter(sig == "Down") %>% pull(Uniprot_ID)
all_diff_proteins <- c(up_proteins, down_proteins)

# 打印蛋白数量，确认读取成功
cat("上调蛋白数：", length(up_proteins), "\n")
cat("下调蛋白数：", length(down_proteins), "\n")
cat("总差异蛋白数：", length(all_diff_proteins), "\n")



# ===================== 差异蛋白富集分析代码 =====================
# 1. 安装并加载所需包

library(clusterProfiler)
library(org.Hs.eg.db)  # 通用蛋白注释库，适配无脊椎动物蛋白
library(ggplot2)
library(dplyr)
library(stringr)
library(enrichplot)

# 2. 读取差异蛋白数据
diff_proteins <- read.csv("C:\\Users\\叶子欣\\Desktop\\差异蛋白_最终结果.csv", stringsAsFactors = FALSE)

# 3. 处理蛋白ID：拆分Uniprot ID，取第一个有效ID
diff_proteins$Uniprot_ID <- str_split(diff_proteins$ProteinID, ";", simplify = TRUE)[, 1]
up_proteins <- diff_proteins %>% filter(sig == "Up") %>% pull(Uniprot_ID)
down_proteins <- diff_proteins %>% filter(sig == "Down") %>% pull(Uniprot_ID)
all_diff_proteins <- c(up_proteins, down_proteins)

# 4. ID转换：Uniprot ID → ENTREZ ID
id_map <- mapIds(
  x = org.Hs.eg.db,
  keys = all_diff_proteins,
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)
id_map <- id_map[!is.na(id_map)]
cat("成功转换的蛋白数：", length(id_map), "\n")

# 5. GO功能富集分析（生物学过程BP、细胞组分CC、分子功能MF）
go_enrich <- enrichGO(
  gene = id_map,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# 导出GO富集结果
write.csv(go_enrich@result, "C:\\Users\\叶子欣\\Desktop\\GO功能富集分析结果.csv", row.names = FALSE)

# 6. KEGG通路富集分析
kegg_enrich <- enrichKEGG(
  gene = id_map,
  organism = "hsa",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# 导出KEGG富集结果
write.csv(kegg_enrich@result, "C:\\Users\\叶子欣\\Desktop\\KEGG通路富集分析结果.csv", row.names = FALSE)

# 7. 高分可视化：GO富集气泡图
go_top20 <- go_enrich@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

ggplot(go_top20, aes(x = Count, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue") +
  labs(
    title = "Top 20 GO Functional Enrichment",
    x = "Gene Number",
    y = "GO Term",
    color = "Adjusted p-value",
    size = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )

ggsave("C:\\Users\\叶子欣\\Desktop\\GO富集气泡图.png", width = 10, height = 8, dpi = 300)

# KEGG通路富集条形图
kegg_top15 <- kegg_enrich@result %>%
  arrange(p.adjust) %>%
  slice_head(n = 15)

ggplot(kegg_top15, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), fill = Count)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(
    title = "Top 15 KEGG Pathway Enrichment",
    x = "-log10(Adjusted p-value)",
    y = "KEGG Pathway",
    fill = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold")
  )

ggsave("C:\\Users\\叶子欣\\Desktop\\KEGG富集条形图.png", width = 10, height = 6, dpi = 300)

# 9. 输出核心结果
cat("\n===== 富集分析完成！=====\n")
cat("GO富集结果已保存到：GO功能富集分析结果.csv\n")
cat("KEGG富集结果已保存到：KEGG通路富集分析结果.csv\n")
cat("GO富集气泡图已保存到：GO富集气泡图.png\n")
cat("KEGG富集条形图已保存到：KEGG富集条形图.png\n")