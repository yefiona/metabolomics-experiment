# ======================
# 1. 清空环境 + 加载包
# ======================
rm(list = ls())
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

desktop <- file.path(Sys.getenv("HOME"), "Desktop")

# ======================
# 2. 读取数据
# ======================
counts <- read.csv("C:/Users/叶子欣/Desktop/exp2_example_CountMatrix.csv", row.names = 1)
sample <- read.csv("C:/Users/叶子欣/Desktop/exp2_example_sample.csv", row.names = 1)

# ======================
# 3. DESeq2 差异表达分析
# ======================
# 构建分组信息（C vs D）
group <- factor(c("C", "C", "C", "D", "D", "D"))
colData <- data.frame(group = group)

# 创建DESeq2对象
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ group
)

# 过滤低表达基因
dds <- dds[rowSums(counts(dds)) > 10, ]

# 差异分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "D", "C"))
res_df <- as.data.frame(res)

# ======================
# 4. 筛选差异基因
# ======================
padj_cut <- 0.05
logfc_cut <- 1

res_df$sig <- "ns"
res_df$sig[!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange > logfc_cut] <- "up"
res_df$sig[!is.na(res_df$padj) & res_df$padj < padj_cut & res_df$log2FoldChange < -logfc_cut] <- "down"

diff_genes <- rownames(subset(res_df, sig %in% c("up", "down")))

# 统计数量
up <- sum(res_df$sig == "up", na.rm = T)
down <- sum(res_df$sig == "down", na.rm = T)

cat("差异基因统计：\n上调：", up, "\n下调：", down, "\n总计：", length(diff_genes), "\n")

# ======================
# 5. 画火山图
# ======================
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.2) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_vline(xintercept = c(-1,1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "D vs C 差异基因火山图", x = "log2FC", y = "-log10(Adjusted P)")

ggsave(file.path(desktop, "火山图.pdf"), p_volcano, width = 8, height = 6)

# ======================
# 6. GO 富集分析（BP）
# ======================
go <- enrichGO(
  gene = diff_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "fdr",
  qvalueCutoff = 0.05
)

# ======================
# 7. GO 气泡图
# ======================
p_go_dot <- dotplot(go, showCategory = 10, title = "GO Biological Process") +
  theme_bw() + scale_color_gradient(low = "red", high = "blue")

ggsave(file.path(desktop, "GO气泡图.pdf"), p_go_dot, width = 10, height = 8)

# ======================
# 8. GO 圆形网络图
# ======================
p_go_circle <- cnetplot(go,
                        showCategory = 10,
                        layout = "circle",
                        colorEdge = TRUE)

ggsave(file.path(desktop, "GO圆形网络图.pdf"), p_go_circle, width = 10, height = 10)

# ======================
# 完成
# ======================
cat("\n✅ 全部分析完成！所有图片已保存到桌面！\n")