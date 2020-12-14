# Correlation analysis of METABRIC data
# Weihua Guo, Ph.D.
# 09/04/2020

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(viridis))
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(gridExtra))

data_dir <- "C:/Users/wguo/Documents/temp_work_athome/"
res_dir <- "C:/Users/wguo/Documents/temp_work_athome/metabric_corr/"
input_df <- read.csv(paste(data_dir, "correlation_matrix_complete.csv", sep = ""), header = T, row.names = 1)
subtype <- "ER"

cate_cols <- c("menopause_state","grade","stage","PR.expr", "pam50", "chemotherapy")
for (icc in cate_cols) {
  input_df[,icc] <- as.factor(input_df[,icc])
}
num_cols <- c("age","tumor_size", "Tex", "CD8A",  "PTPRC", "oncotype_dx", "proliferation","ost", "dmfs")
if (subtype == "TNBC") {
  cate_cols <- c("menopause_state","grade","stage","pam50", "chemotherapy")
  for (icc in cate_cols) {
    input_df[,icc] <- as.factor(input_df[,icc])
  }
  num_cols <- c("age","tumor_size", "Tex", "CD8A",  "PTPRC", "oncotype_dx", "proliferation","ost", "dmfs")
}

all_cols <- c(cate_cols, num_cols)
all_comb <- combn(all_cols, 2)

plot_m <- matrix(NA, length(all_cols), length(all_cols))
plot_m[lower.tri(plot_m, diag = F)] <- 1:ncol(all_comb)

all_cols <- c(cate_cols, num_cols)
all_comb <- combn(all_cols, 2)
plot_list <- vector(mode = "list", length = ncol(all_comb))
names(plot_list) <- 1:ncol(all_comb)

plot_df <- input_df[input_df$subtype == "ER+",]
plotPf <- paste(res_dir, "ER_", sep = "")

if (subtype == "TNBC") {
  cat("Forest plot...\n")
  os_fit <- coxph(Surv(ost, ose) ~ Tex + menopause_state + age + oncotype_dx + chemotherapy + CD8A + PTPRC + proliferation + tumor_size + grade, data = plot_df)
  os_ggf <- ggforest(os_fit)
  ggsave(paste(plotPf, "coxph_forest_os_all.png", sep = ""), os_ggf, dpi = 300, width = 6, height =6, limitsize = FALSE)
  
  dmfs_fit <- coxph(Surv(dmfs, dmfe) ~ Tex + menopause_state + age + oncotype_dx + chemotherapy + CD8A + PTPRC + proliferation + tumor_size + grade, data = plot_df)
  dmfs_ggf <- ggforest(dmfs_fit)
  ggsave(paste(plotPf, "coxph_forest_dmfs_all.png", sep = ""), dmfs_ggf, dpi = 300, width = 6, height =6, limitsize = FALSE)
} else {
  cat("Forest plot...\n")
  os_fit <- coxph(Surv(ost, ose) ~ Tex + menopause_state + age + PR.expr + oncotype_dx + chemotherapy + CD8A + PTPRC + proliferation + tumor_size + pam50 + grade, data = plot_df)
  os_ggf <- ggforest(os_fit)
  ggsave(paste(plotPf, "coxph_forest_os_all.png", sep = ""), os_ggf, dpi = 300, width = 6, height =6, limitsize = FALSE)
  
  dmfs_fit <- coxph(Surv(dmfs, dmfe) ~ Tex + menopause_state + age + PR.expr + oncotype_dx + chemotherapy + CD8A + PTPRC + proliferation + tumor_size + pam50 + grade, data = plot_df)
  dmfs_ggf <- ggforest(dmfs_fit)
  ggsave(paste(plotPf, "coxph_forest_dmfs_all.png", sep = ""), dmfs_ggf, dpi = 300, width = 6, height =6, limitsize = FALSE)
}

cor_res <- as.data.frame(t(all_comb))
colnames(cor_res) <- c("V1", "V2")
cor_res$r <- NA
cor_res$p <- NA
cor_res$method <- "NONE"


cat("One-to-one correlation...\n")
for (ic in 1:nrow(cor_res)) {
  x_name <- all_comb[1,ic]
  y_name <- all_comb[2, ic]
  cat("\tx:", x_name, "\ty:", y_name, "\n")
  if (x_name %in% num_cols & y_name %in% num_cols) {
    tmp_res <- cor.test(plot_df[,x_name], plot_df[, y_name])
    cor_res$r[ic] <- tmp_res$estimate
    cor_res$p[ic] <- tmp_res$p.value
    cor_res$method[ic] <- "Pearson"
    
  }
  if (x_name %in% cate_cols & y_name %in% cate_cols) {
    tmp_tbl <- as.data.frame.matrix(table(plot_df[,x_name], plot_df[,y_name]))
    tmp_res <- chisq.test(tmp_tbl, correct = FALSE)
    cor_res$r[ic] <- sqrt(tmp_res$statistic/sum(tmp_tbl))
    cor_res$p[ic] <- tmp_res$p.value
    cor_res$method[ic] <- "Crammer's V"
    
  }
  if (x_name %in% cate_cols & y_name %in% num_cols) {
    tmp_res <- aov(as.formula(paste(y_name, x_name, sep = "~")), data = plot_df)
    tmp_sum <- summary.lm(tmp_res)
    cor_res$r[ic] <- sqrt(tmp_sum$r.squared)
    cor_res$p[ic] <- tmp_sum$coefficients[2,ncol(tmp_sum$coefficients)]
    cor_res$method[ic] <- "ANOVA (Linear model)"
  }
  
  if (y_name %in% cate_cols & x_name %in% num_cols) {
    tmp_res <- aov(as.formula(paste(x_name, y_name, sep = "~")), data = plot_df)
    tmp_sum <- summary.lm(tmp_res)
    cor_res$r[ic] <- sqrt(tmp_sum$r.squared)
    cor_res$p[ic] <- tmp_sum$coefficients[2,ncol(tmp_sum$coefficients)]
    cor_res$method[ic] <- "ANOVA (Linear model)"
  }
}
cts <- as.data.frame(table(cor_res$V1))
cor_res$V1C <- 0
for (i in 1:nrow(cts)) {
  tmp_mask <- cor_res$V1 == cts[i,1]
  cor_res$V1C[tmp_mask] <- cts[i,2]
}

cts <- as.data.frame(table(cor_res$V2))
cor_res$V2C <- 0
for (i in 1:nrow(cts)) {
  tmp_mask <- cor_res$V2 == cts[i,1]
  cor_res$V2C[tmp_mask] <- cts[i,2]
}

cor_res$pflag <- "ns"
cor_res$pflag[cor_res$p <= 0.05] <- "*"
cor_res$pflag[cor_res$p <= 0.01] <- "**"
cor_res$pflag[cor_res$p <= 0.001] <- "***"
cor_res$pflag[cor_res$p <= 0.001] <- "****"
cor_gg <- ggplot(cor_res, aes(x = reorder(V1, -V1C), y = reorder(V2, -V2C))) +
  geom_point(aes(color = r, size = pflag)) +
  scale_color_viridis_c(option = "plasma", breaks = c(-0.2, 0.0, 0.3, 0.6, 0.9)) +
  scale_size_manual(breaks = c("ns", "*", "**", "***", "****"), values = c(1,4,5,7,9)) +
  labs(size = "Significance", color = "Correlation\ncoefficients") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        legend.position = "right")
ggsave(paste(plotPf, "correlation_dotplot_triangle.png", sep = ""), cor_gg, dpi = 300, width = 7.2, height = 6, limitsize = FALSE)


cat("One-to-one visualization...\n")
for (ic in 1:ncol(all_comb)) {
  x_name <- all_comb[1,ic]
  y_name <- all_comb[2, ic]
  cat("\tx:", x_name, "\ty:", y_name, "\n")
  if (x_name %in% cate_cols & y_name %in% num_cols) {
    cat("\t\tMixed situation...\n")
    tmp_gg <- ggplot(plot_df, aes_string(x = x_name, y = y_name)) +
      geom_boxplot(aes_string(color = x_name)) +
      geom_jitter(position = position_jitterdodge(0.2), alpha = 0.54, size = 2.4, aes_string(color = x_name)) +
      stat_compare_means(label = "p.signif", label.x = 1.35, method = "anova") +
      theme_classic() +
      theme(axis.text=element_text(size=5.67, color = "#000000"),
	    axis.title=element_text(size=5.67, color = "#000000", face='bold'))

  }
  if (x_name %in% num_cols & y_name %in% cate_cols) {
    tmp_gg <- ggplot(plot_df, aes_string(x = x_name, y = y_name)) +
      geom_boxplot(aes_string(color = x_name)) +
      geom_jitter(position = position_jitterdodge(0.2), alpha = 0.54, size = 2.4, aes_string(color = x_name)) +
      stat_compare_means(label = "p.signif", label.x = 1.35, method = "anova") +
      theme_classic() +
      theme(axis.text=element_text(size=5.67, color = "#000000"),
	    axis.title=element_text(size=5.67, color = "#000000", face='bold'))
  }
  if (x_name %in% cate_cols & y_name %in% cate_cols) {
    cat("Both categorical\n")
    tmp_table <- as.data.frame(table(plot_df[,x_name], plot_df[, y_name]))
    colnames(tmp_table) <- c(x_name, y_name, "num")
    tmp_test <- fisher.test(table(plot_df[,x_name], plot_df[, y_name]),simulate.p.value=TRUE,B=1e7)
    pflag <- "ns"
    if (tmp_test$p.value <= 0.05) {
      pflag <- "*"
    } else if (tmp_test$p.value <= 0.01) {
      pflag <- "**"
    } else if (tmp_test$p.value <= 0.001) {
      pflag <- "***"
    } else if (tmp_test$p.value <= 0.0001) {
      pflag <- "****"
    }
    tmp_gg <- ggplot(tmp_table, aes_string(x = x_name, y = "num")) +
      geom_bar(stat = "identity", position = position_stack(), aes_string(fill = y_name)) +
      scale_fill_brewer(palette = "Paired") + 
      labs(y = "Patient number") +
      annotate("text", x = 2, y = 2*max(tmp_table$num), label = pflag) + 
      theme_classic() +
      theme(axis.text=element_text(size=5.67, color = "#000000"),
	    axis.title=element_text(size=5.67, color = "#000000", face='bold'))

  }
  
  if (x_name %in% num_cols & y_name %in% num_cols) {
    cat("Both numericals\n")
    tmp_gg <- ggscatter(plot_df, x = x_name, y = y_name, color = "firebrick",
                        size = 0.6, add = "reg.line", conf.int = TRUE, alpha = 0.54,
                        add.params = list(color = "blue", fill = "lightgray")) +
      stat_cor(method = "pearson") +
      theme_classic() +
      theme(axis.text=element_text(size=5.67, color = "#000000"),
	    axis.title=element_text(size=5.67, color = "#000000", face='bold'))

  }
  plot_list[[ic]] <- tmp_gg
  ggsave(paste(plotPf, x_name, "_vs_",y_name, "_specific_correlation_triangle.png", sep = ""), tmp_gg, dpi = 300, width = 2, height = 1.5)
}

tri_gg <- grid.arrange(grobs = plot_list, layout_matrix = plot_m)
ggsave(paste(plotPf, "specific_correlation_triangle.pdf", sep = ""), tri_gg, dpi = 300, width = 27, height = 21, limitsize = FALSE)
ggsave(paste(plotPf, "specific_correlation_triangle.png", sep = ""), tri_gg, dpi = 300, width = 27, height = 21, limitsize = FALSE)
