# Correlation analysis of METABRIC data
# Weihua Guo, Ph.D.
# 09/04/2020

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(gridExtra))

data_dir <- "C:/Users/wguo/Documents/temp_work_athome/"
res_dir <- "C:/Users/wguo/Documents/temp_work_athome/metabric_corr/"
input_df <- read.csv(paste(data_dir, "correlation_matrix_complete.csv", sep = ""), header = T, row.names = 1)

cate_cols <- c("menopause_state","grade","stage","PR.expr", "pam50", "chemotherapy")
for (icc in cate_cols) {
  input_df[,icc] <- as.factor(input_df[,icc])
}
num_cols <- c("age","tumor_size", "Tex", "CD8A",  "PTPRC", "oncotype_dx", "proliferation","ost", "dmfs")

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
for (ic in 1:ncol(all_comb)) {
  x_name <- all_comb[1,ic]
  y_name <- all_comb[2, ic]
  cat("\tx\t", x_name, "\ty\t", y_name, "\n")
  if (x_name %in% cate_cols & y_name %in% num_cols) {
    cat("\t\tMixed situation...\n")
    tmp_gg <- ggplot(plot_df, aes_string(x = x_name, y = y_name)) +
      geom_boxplot(aes_string(color = x_name)) +
      geom_jitter(position = position_jitterdodge(0.2), alpha = 0.54, size = 2.4, aes_string(color = x_name)) +
      stat_compare_means(label = "p.signif", label.x = 1.35) +
      theme_classic()
  }
  if (x_name %in% num_cols & y_name %in% cate_cols) {
    tmp_gg <- ggplot(plot_df, aes_string(x = x_name, y = y_name)) +
      geom_boxplot(aes_string(color = x_name)) +
      geom_jitter(position = position_jitterdodge(0.2), alpha = 0.54, size = 2.4, aes_string(color = x_name)) +
      stat_compare_means(label = "p.signif", label.x = 1.35) +
      theme_classic()
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
      theme_classic()
  }
  
  if (x_name %in% num_cols & y_name %in% num_cols) {
    cat("Both numericals\n")
    tmp_gg <- ggscatter(plot_df, x = x_name, y = y_name, color = "firebrick",
                        size = 0.6, add = "reg.line", conf.int = TRUE, alpha = 0.54,
                        add.params = list(color = "blue", fill = "lightgray")) +
      stat_cor(method = "pearson") +
      theme_classic()
  }
  plot_list[[ic]] <- tmp_gg
}

tri_gg <- grid.arrange(grobs = plot_list, layout_matrix = plot_m)
ggsave(paste(plotPf, "specific_correlation_triangle.pdf", sep = ""), tri_gg, dpi = 300, width = 49, height = 42, limitsize = FALSE)
