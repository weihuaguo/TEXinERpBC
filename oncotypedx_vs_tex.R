# For oncotype DX vs Tex
# Weihua Guo, Ph.D
# 02/04/2021

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(stringr))
suppressMessages(library(cocor))
suppressMessages(library(scales))
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(survMisc))
suppressMessages(library(rstatix))

data_dir <- 'C:/Users/wguo/OneDrive - City of Hope National Medical Center/tmp_works/oncotype_dx_pre_vs_post/'
bc_type <- 'ER'
meno_status <- 'all'
onco_status <- 'int'
plot_title <- paste('Menopause', meno_status, '- Oncotype Dx', onco_status)
plot_pf <- paste(data_dir, 'final_05152021', sep = '')

df <- read.csv(paste(data_dir, 'metabric_', bc_type, '_tex_oncotype_sigscore.csv', sep = ''), row.names = 1)

df$oncotype_unscale <- df$oncotype
latest_onco <- read.csv(paste(data_dir, 'oncotype_function_sigscore_', bc_type, '_cbpt_allmatch.csv', sep = ''), row.names = 1)
df$oncotype_scale <- latest_onco[rownames(df),'sigscore']

clinical_df <- read.csv(paste(data_dir, 'brca_metabric_clinical_data.csv', sep = ''), header = T, row.names = 1,
                          check.names = T)
rownames(clinical_df) <- str_replace(rownames(clinical_df), '-', '.')

cat("Compare oncotype Dx and Tex...\n")
all_tex_oncodx <- ggscatter(df, x = 'oncotype_unscale', y = 'Tex', size = 1,
                            add = 'reg.line', add.params = list(color = 'blue', fill = 'lightgray'),
                            conf.int = TRUE) +
  stat_cor(method = 'pearson') + 
  labs(x = 'Oncotype Dx', y = 'Tex', title = 'All the samples')
ggsave(paste(plot_pf, 'oncotype_tex_scatter.png', sep = ''), all_tex_oncodx,
       dpi = 600, width = 4.5, height = 4.8)

pre_tex_oncodx <- ggscatter(df[df$menopausal_State=='pre',], x = 'oncotype_unscale', y = 'Tex', size = 1,
                            add = 'reg.line', add.params = list(color = 'blue', fill = 'lightgray'),
                            conf.int = TRUE) +
  stat_cor(method = 'pearson') + 
  labs(x = 'Oncotype Dx', y = 'Tex', title = "Pre-menopause")
ggsave(paste(plot_pf, 'oncotype_tex_scatter_pre.png', sep = ''), pre_tex_oncodx,
       dpi = 600, width = 4.5, height = 4.8)

post_tex_oncodx <- ggscatter(df[df$menopausal_State=='post',], x = 'oncotype_unscale', y = 'Tex', size = 1,
                            add = 'reg.line', add.params = list(color = 'blue', fill = 'lightgray'),
                            conf.int = TRUE) +
  stat_cor(method = 'pearson') + 
  labs(x = 'Oncotype Dx', y = 'Tex', title = "Post-menopause")
ggsave(paste(plot_pf, 'oncotype_tex_scatter_post.png', sep = ''), post_tex_oncodx,
       dpi = 600, width = 4.5, height = 4.8)


cat("Compare the actual and in-house oncotype Dx...\n")
onco_comp_sct <- ggscatter(df, x = 'oncotype_unscale', y = 'oncotype_scale',
                           add = 'reg.line', add.params = list(color = 'blue', fill = 'lightgray'),
                           conf.int = TRUE) +
  stat_cor(method = 'pearson') +
  geom_vline(xintercept = quantile(df$oncotype_unscale, c(0.15, 0.85))) +
  geom_hline(yintercept = c(11,25), color = 'red') +
  geom_hline(yintercept = quantile(df$oncotype_scale, c(0.15, 0.85)))

ggsave(paste(plot_pf, 'oncotype_compare_scatter.png', sep = ''), onco_comp_sct,
       dpi = 600, width = 6, height = 6)

onco_comp_hists <- ggplot(df, aes(x = oncotype_scale)) + 
  geom_histogram(fill = 'skyblue', color = 'blue', bins = 81) +
  geom_vline(xintercept = c(11,25), color = 'red') +
  geom_vline(xintercept = quantile(df$oncotype_scale, c(0.15, 0.85))) +
  labs(title = 'Potential actual oncotype DX') + theme_classic()

onco_comp_histu <- ggplot(df, aes(x = oncotype_unscale)) + 
  geom_histogram(fill = 'skyblue', color = 'blue', bins = 81) +
  geom_vline(xintercept = quantile(df$oncotype_unscale, c(0.15, 0.85))) +
  labs(title = 'Oncotype Dx from sig.score') + theme_classic()

onco_comp_hist <- ggarrange(onco_comp_hists, onco_comp_histu, ncol = 1, nrow = 2)

ggsave(paste(plot_pf, 'oncotype_compare_histogram.png', sep = ''), onco_comp_hist,
       dpi = 600, width = 6, height = 6)

plot_pf <- paste(plot_pf, bc_type, meno_status, 'onco', onco_status, sep = '_')

if (meno_status == 'all') {
  cat('INDEPENDENT!!!\n')
  fit <- coxph(Surv(ost, ose) ~ menopausal_State+grade+tsize+oncotype+Tex, 
               data = df)
  forest_gg <- ggforest(fit)
  ggsave(paste(plot_pf, 'multiv_all_tex_onco_meno_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)
  
  fit <- coxph(Surv(rfst, rfse) ~ menopausal_State+grade+tsize+oncotype+Tex, 
               data = df)
  forest_gg <- ggforest(fit)
  ggsave(paste(plot_pf, 'multiv_all_tex_onco_meno_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)
}
if (meno_status == 'pre') {
  use_df <- df[df$menopausal_State == 'pre',]
} else if (meno_status == 'post') {
  use_df <- df[df$menopausal_State == 'post',]  
} else {use_df <- df}

use_df$group_onco <- ifelse(use_df$oncotype > quantile(use_df$oncotype, 0.85), 'High', 
                            ifelse(use_df$oncotype < quantile(use_df$oncotype, 0.15), 'Low', 'Medium'))
use_df$group_tex <- ifelse(use_df$Tex > quantile(use_df$Tex, 0.75), 'High', 
                           ifelse(use_df$Tex < quantile(use_df$Tex, 0.25), 'Low', 'Medium'))

write.csv(use_df, paste(plot_pf, 'organized_df.csv', sep = ''))
stop("TEST")

cat("\tThree group comparison...\n")
m_use_df <- use_df
m_use_df$group_final <- ifelse(m_use_df$group_onco == 'Medium', 
                               str_c('Oncotype Medium - Tex ', m_use_df$group_tex), 
                               str_c('Oncotype ', m_use_df$group_onco))
m_use_df <- m_use_df[!str_detect(m_use_df$group_final, 'Tex Medium'),]

fit <- survfit(Surv(ost, ose) ~ group_final, data = m_use_df)
surv_gg <- ggsurvplot(fit, data = m_use_df, pval = TRUE,
                      title = paste('Menopause', meno_status), legend = 'right')
png(paste(plot_pf, 'os_onco_medium_three_group.png', sep = '_'), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()
ps_res <- pairwise_survdiff(Surv(ost, ose) ~ group_final, data = m_use_df)
write.csv(ps_res$p.value, paste(plot_pf, 'os_onco_medium_three_group_pairwise.csv', sep = '_'))

fit <- survfit(Surv(rfst, rfse) ~ group_final, data = m_use_df)
surv_gg <- ggsurvplot(fit, data = m_use_df, pval = TRUE,
                      title = paste('Menopause', meno_status), legend = 'right')
png(paste(plot_pf, 'rfs_onco_medium_three_group.png', sep = '_'), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()
ps_res <- pairwise_survdiff(Surv(rfst, rfse) ~ group_final, data = m_use_df)
write.csv(ps_res$p.value, paste(plot_pf, 'rfs_onco_medium_three_group_pairwise.csv', sep = '_'))

ml_use_df <- use_df
ml_use_df$group_final <- ifelse(ml_use_df$group_onco != 'High', 
                               str_c('Oncotype Low/Medium - Tex ', ml_use_df$group_tex), 
                               str_c('Oncotype ', ml_use_df$group_onco))
ml_use_df <- ml_use_df[!str_detect(ml_use_df$group_final, 'Tex Medium'),]

fit <- survfit(Surv(ost, ose) ~ group_final, data = ml_use_df)
surv_gg <- ggsurvplot(fit, data = ml_use_df, pval = TRUE,
                      title = paste('Menopause', meno_status), legend = 'right')
png(paste(plot_pf, 'os_onco_medium-low_three_group.png', sep = '_'), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()
ps_res <- pairwise_survdiff(Surv(ost, ose) ~ group_final, data = ml_use_df)
write.csv(ps_res$p.value, paste(plot_pf, 'os_onco_medium-low_three_group_pairwise.csv', sep = '_'))

fit <- survfit(Surv(rfst, rfse) ~ group_final, data = ml_use_df)
surv_gg <- ggsurvplot(fit, data = ml_use_df, pval = TRUE,
                      title = paste('Menopause', meno_status), legend = 'right')
png(paste(plot_pf, 'rfs_onco_medium-low_three_group.png', sep = '_'), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()
ps_res <- pairwise_survdiff(Surv(rfst, rfse) ~ group_final, data = ml_use_df)
write.csv(ps_res$p.value, paste(plot_pf, 'rfs_onco_medium-low_three_group_pairwise.csv', sep = '_'))

if (onco_status == 'high') {
  tmp_df <- use_df[use_df$group_onco == 'High',]
} else if (onco_status == 'low') {
  tmp_df <- use_df[use_df$group_onco == 'Low',]
} else if (onco_status == 'int') {
  tmp_df <- use_df[use_df$group_onco == 'Medium',]
} else {
  tmp_df <- use_df[use_df$group_onco != 'High',]
}

colnames(clinical_df)
clin_cols <- c("Age.at.Diagnosis", "Chemotherapy", "Lymph.nodes.examined.positive", 
               "Tumor.Size", "Tumor.Stage", "Neoplasm.Histologic.Grade")
use_clin <- clinical_df[,clin_cols]
tmp_df <- cbind(tmp_df, use_clin[rownames(tmp_df),])
tmp_df$Tumor.Stage <- as.factor(tmp_df$Tumor.Stage)
tmp_df$Neoplasm.Histologic.Grade <- as.factor(tmp_df$Neoplasm.Histologic.Grade)
tmp_df$Chemotherapy <- as.factor(tmp_df$Chemotherapy)

all_tmp_df <- tmp_df

onco_tex_sub <- ggscatter(all_tmp_df, x = 'Tex', y = 'oncotype', color = 'group_tex',
                           add = 'reg.line', add.params = list(color = 'blue', fill = 'lightgray'),
                           conf.int = TRUE) +
  stat_cor(method = 'pearson') +
  labs(title = plot_title, color = 'Tex groups', y = 'Oncotype DX')
ggsave(paste(plot_pf, '_oncotype_tex_corr.png', sep = ''), onco_tex_sub,
       dpi = 600, width = 6, height = 6)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(ost, ose) ~ oncotype+Tex, 
             data = all_tmp_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, 'multiv_sub_tex_onco_os.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_sub_tex_onco_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(rfst, rfse) ~ oncotype+Tex, 
             data = all_tmp_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, 'multiv_sub_tex_onco_rfs.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_sub_tex_onco_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)


cat('Hazard ratio -- Multivariant -- ALL\n')
fit <- coxph(Surv(ost, ose) ~ age+grade+tsize+oncotype+Tex,
             data = all_tmp_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, 'multiv_all_tex_onco_os.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_all_tex_onco_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant -- ALL\n')
fit <- coxph(Surv(rfst, rfse) ~ age+grade+tsize+oncotype+Tex, 
             data = all_tmp_df)
sum_fit <- summary(fit)
sum_fit_df <- cbind(sum_fit$coefficients, sum_fit$conf.int)
write.csv(sum_fit_df, paste(plot_pf, 'multiv_all_tex_onco_rfs.csv', sep = '_'))
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_all_tex_onco_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)


tmp_df <- tmp_df[tmp_df$group_tex != 'Medium',]
onco_tex_box <- ggplot(tmp_df, aes(x = group_tex, y = oncotype, color = group_tex)) +
  geom_boxplot() +
  geom_point() +
  labs(color = 'Tex groups', y = 'Oncotype DX', title = plot_title) +
  stat_compare_means(aes(group = group_tex), label.x = 1.25, method = 't.test') +
  theme_classic() +
  theme(axis.title.x = element_blank())
ggsave(paste(plot_pf, '_oncotype_tex_comp.png', sep = ''), onco_tex_box,
       dpi = 600, width = 4.5, height = 3)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(ost, ose) ~ Tex+oncotype, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_tex_noint_onco_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(rfst, rfse) ~ Tex+oncotype, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_tex_noint_onco_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(ost, ose) ~ group_tex+oncotype, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_tex_group_onco_os.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)

cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(rfst, rfse) ~ group_tex+oncotype, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_tex_group_onco_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 6, height = 4)


cat('Hazard ratio -- Multivariant...\n')
fit <- coxph(Surv(ost, ose) ~ Age.at.Diagnosis+Tumor.Size+Tumor.Stage+Neoplasm.Histologic.Grade+Tex, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_hr_os.png', sep = '_'), forest_gg, dpi = 600, width = 9, height = 6)

fit <- coxph(Surv(rfst, rfse) ~ Age.at.Diagnosis+Tumor.Size+Tumor.Stage+Neoplasm.Histologic.Grade+Tex, 
             data = tmp_df)
forest_gg <- ggforest(fit)
ggsave(paste(plot_pf, 'multiv_hr_rfs.png', sep = '_'), forest_gg, dpi = 600, width = 9, height = 6)

cat('\t\tForest plot...Univariant\n')
#https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/
covariates <- c("Tex", "Age.at.Diagnosis", "Tumor.Size", "Tumor.Stage", "Neoplasm.Histologic.Grade")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(ost, ose)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = tmp_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1,1], digits=2);#coeficient beta
                         HR <-signif(x$coef[1,2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"],2)
                         names <- rownames(x$coef)[1]
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value, names)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value", "sig")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$HR <- as.numeric(as.character(res$HR))
res$HRL <- as.numeric(as.character(res$HRL))
res$HRU <- as.numeric(as.character(res$HRU))
# res$sig <- c('Tex', 'OncotypeDX', clin_cols)
res$p_format <- str_c('p-val = ', res$p.value)

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = 'Hazard ratio (univariant tests, overall survival)', y = 'Variables', title = plot_title) +
  theme_classic()
ggsave(paste(plot_pf, 'univ_os_hr.png', sep = '_'), form_gg, dpi = 600, width = 9, height = 4.5)
write.csv(res, paste(plot_pf, 'univ_os_hr.csv', sep = '_'))

covariates <- c("Tex", "Age.at.Diagnosis", "Tumor.Size", "Tumor.Stage", "Neoplasm.Histologic.Grade")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(rfst, rfse)~', x)))
univ_models <- lapply(univ_formulas, function(x){coxph(x, data = tmp_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1,1], digits=2);#coeficient beta
                         HR <-signif(x$coef[1,2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[1,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[1,"upper .95"],2)
                         names <- rownames(x$coef)[1]
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value, names)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value", "sig")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$HR <- as.numeric(as.character(res$HR))
res$HRL <- as.numeric(as.character(res$HRL))
res$HRU <- as.numeric(as.character(res$HRU))
# res$sig <- c('Tex', 'OncotypeDX', clin_cols)
res$p_format <- str_c('p-val = ', res$p.value)

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  scale_x_continuous(trans = log2_trans(),
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) +
  labs(x = 'Hazard ratio (univariant tests, relapse-free survival)', y = 'Variables', title = plot_title) +
  theme_classic()
ggsave(paste(plot_pf, 'univ_rfs_hr.png', sep = '_'), form_gg, dpi = 600, width = 9, height = 4.5)
write.csv(res, paste(plot_pf, 'univ_rfs_hr.csv', sep = '_'))
stop("HERE")

cat('\t\tTex\n')
fit <- survfit(Surv(ost, ose) ~ group_tex, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(plot_pf, 'os_tex_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_tex, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(plot_pf, 'rfs_tex_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('\t\tTex high: Chemo\n')
tmp_hi_df <- tmp_df[tmp_df$group_tex == 'High',]
fit <- survfit(Surv(ost, ose) ~ Chemotherapy, data = tmp_hi_df)
surv_gg <- ggsurvplot(fit, data = tmp_hi_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Chemotherapy')
png(paste(plot_pf, 'os_chemo_texhi_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ Chemotherapy, data = tmp_hi_df)
surv_gg <- ggsurvplot(fit, data = tmp_hi_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Chemotherapy')
png(paste(plot_pf, 'rfs_chemo_texhi_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('\t\tTex int/low: Chemo\n')
tmp_lo_df <- all_tmp_df[all_tmp_df$group_tex != 'High',]
fit <- survfit(Surv(ost, ose) ~ Chemotherapy, data = tmp_lo_df)
surv_gg <- ggsurvplot(fit, data = tmp_lo_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Chemotherapy')
png(paste(plot_pf, 'os_chemo_tex_lo_int_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ Chemotherapy, data = tmp_lo_df)
surv_gg <- ggsurvplot(fit, data = tmp_lo_df, pval = TRUE,
                      title = plot_title,
                      legend.title = 'Chemotherapy')
png(paste(plot_pf, 'rfs_chemo_tex_lo_int_km.png', sep = '_'), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

stop("TEST")
#############################SEPARATION####





cat('Hazard ratio...\n')
fit <- coxph(Surv(ost, ose) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)

fit <- coxph(Surv(rfst, rfse) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)
stop("TEST")

pre_df <- df[df$menopausal_State == 'pre',]
post_df <- df[df$menopausal_State == 'post',]
cat('Hazard ratio...\n')
fit <- coxph(Surv(ost, ose) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)

fit <- coxph(Surv(rfst, rfse) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)

use_df <- df
meno_state <- 'all80'
cat('\t', meno_state, '\n')
use_df$group_onco <- ifelse(use_df$oncotype > quantile(use_df$oncotype, 0.80), 'High', 
                            ifelse(use_df$oncotype < quantile(use_df$oncotype, 0.15), 'Low', 'Medium'))
use_df$group_onco <- ifelse(use_df$oncotype > quantile(use_df$oncotype, 0.80), 'High', 'Medium')
use_df$group_tex <- ifelse(use_df$Tex > quantile(use_df$Tex, 0.75), 'High', 
                           ifelse(use_df$Tex < quantile(use_df$Tex, 0.25), 'Low', 'Medium'))

tmp_df <- use_df[use_df$group_onco == 'Medium',]
tmp_df <- tmp_df[tmp_df$group_tex != 'Medium',]

cat('\t\tTex\n')
fit <- survfit(Surv(ost, ose) ~ group_tex, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_os_tex_oncomedium.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_tex, data = tmp_df)
surv_gg <- ggsurvplot(fit, data = tmp_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_rfs_tex_oncomedium.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()
stop("TEST")
##########################################################

use_df$hhll <- str_c('Tex:', use_df$group_tex, '-Oncotype:', use_df$group_onco)

cat('\t\tFour group....\n')
hhll_df <- use_df[use_df$group_tex != 'Medium' & use_df$group_onco != 'Medium',]
fit <- survfit(Surv(ost, ose) ~ hhll, data = hhll_df)
surv_gg <- ggsurvplot(fit, data = hhll_df, pval = TRUE,
                      title = meno_state, legend = 'right', legend.title = 'Combined groups')
png(paste(data_dir, meno_state, '_os_hhll.png', sep = ''), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ hhll, data = hhll_df)
surv_gg <- ggsurvplot(fit, data = hhll_df, pval = TRUE,
                      title = meno_state, legend = 'right', legend.title = 'Combined groups')
png(paste(data_dir, meno_state, '_rfs_hhll.png', sep = ''), res = 600, width = 9, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('\t\tCutP\n')
tex_coxph <- coxph(Surv(ost, ose) ~ Tex, data = use_df)
tex_coxph <- cutp(tex_coxph)$Tex
merge_tex_cutp <- merge(tex_coxph, use_df, by = 'Tex')
merge_tex_cutp$sig <- 'Tex score'
merge_tex_cutp$stype <- 'OS'
onco_coxph <- coxph(Surv(ost, ose) ~ oncotype, data = use_df)
onco_coxph <- cutp(onco_coxph)$oncotype
merge_onco_cutp <- merge(onco_coxph, use_df, by = 'oncotype')
merge_onco_cutp$sig <- 'Oncotype DX'
merge_onco_cutp$stype <- 'OS'
merge_cutp <- rbind(merge_tex_cutp, merge_onco_cutp)

tex_mean <- mean(merge_cutp$U[merge_cutp$sig == 'Tex score'])

onco_time <- mean(use_df$Tex/use_df$oncotype)
onco_mean <- mean(merge_cutp$U[merge_cutp$sig == 'Oncotype DX'])

use_df$comb <- use_df$Tex*tex_mean + use_df$oncotype*onco_mean*onco_time

hybrid_coxph <- coxph(Surv(ost, ose) ~ comb, data = use_df)
hybrid_coxph <- cutp(hybrid_coxph)$comb
merge_hybrid_cutp <- merge(hybrid_coxph, use_df, by = 'comb')
merge_hybrid_cutp$sig <- 'Hybrid'
merge_hybrid_cutp$stype <- 'OS'
merge_hybrid_cutp$comb <- NULL
merge_cutp <- rbind(merge_cutp, merge_hybrid_cutp)

use_df$group_comb <- ifelse(use_df$comb > quantile(use_df$comb, 0.75), 'High', 
                           ifelse(use_df$comb < quantile(use_df$comb, 0.25), 'Low', 'Medium'))

spr_cutp <- spread(merge_cutp[, c('pid', 'U', 'sig', 'Tex', 'oncotype')], 'sig', 'U')
cor_gg <- ggscatter(spr_cutp, x='Tex score', y = 'Oncotype DX', 
                    size = 1, alpha = 0.3,
                    add = 'reg.line',
                    conf.int = TRUE, 
                    rug = TRUE) +
  stat_cor() + 
  labs(y = 'log-rank score\n(Oncotype DX)', x = 'log-rank score\n(Tex)') +
  theme_classic() +
  theme(legend.position = 'bottom')
ggsave(paste(data_dir, meno_state, '_os_cutp_correlation.png', sep = ''), cor_gg, dpi = 600, width = 6, height = 4)

comps <- list(c('Tex score', 'Oncotype DX'), c('Tex score', 'Hybrid'), c('Oncotype DX', 'Hybrid'))
box_gg <- ggplot(merge_cutp, aes(x = sig, y = U, color = Tex)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(0.3)) +
  stat_compare_means(label = 'p.signif', comparisons = comps) +
  scale_color_distiller(palette = 'RdYlBu') +
  labs(x = 'Signatures', y = 'log-rank score', fill = 'Tex', title = meno_state) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste(data_dir, meno_state, '_os_cutp_box.png', sep = ''), box_gg, dpi = 600, width = 3, height = 4)

cat('\t\tGroup counts...\n')
group_cts <- as.data.frame(table(use_df$group_onco, use_df$group_tex))
colnames(group_cts) <- c('Oncotype', 'Tex', 'Num')
cts_gg <- ggplot(group_cts, aes_string(x = 'Tex', y = 'Oncotype', color = 'Num', size = 'Num')) +
  geom_point() +
  scale_color_distiller(palette = 'Spectral') +
  labs(title = meno_state) +
  theme_bw() +
  theme(legend.position = 'top')
ggsave(paste(data_dir, meno_state, '_group_onco_vs_tex.png', sep = ''), cts_gg, dpi = 600, width = 6, height = 3)

cat('\t\tForest plot...MANUALLL\n')
#https://www.r-bloggers.com/2016/12/cox-proportional-hazards-model/
covariates <- c("Tex", "oncotype",  "comb")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(ost, ose)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = use_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$sig <- c('Tex', 'OncotypeDX', 'Hybrid')
res$p_format <- str_c('p-val = ', res$p.value)

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1.0, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  labs(x = 'Hazard ratio (univariant tests, overall survival)', y = 'Signatures') +
  theme_classic()
ggsave(paste(data_dir, meno_state, '_os_hr.png', sep = ''), form_gg, dpi = 600, width = 6, height = 3)

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(rfst, rfse)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = use_df)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$logtest["pvalue"], digits=2)
                         logrank.test<-signif(x$logtest["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         res<-c(beta, HR, HR.confint.lower,HR.confint.upper, logrank.test, p.value)
                         names(res)<-c("beta", "HR", "HRL", "HRU", "logrank", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$sig <- c('Tex', 'OncotypeDX', 'Hybrid')
res$p_format <- str_c('p-val = ', res$p.value)

form_gg <- ggplot(res, aes(x = HR, y = sig)) +
  geom_point() +
  geom_errorbar(aes(xmin = HRL, xmax = HRU), width = 0.2) +
  geom_vline(xintercept = 1.0, linetype = 'dashed', color = 'red') +
  geom_text(aes(label = p_format), hjust = -1, vjust = -1.5) +
  labs(x = 'Hazard ratio (univariant tests, relapse-free survival)', y = 'Signatures') +
  theme_classic()
ggsave(paste(data_dir, meno_state, '_rfs_hr.png', sep = ''), form_gg, dpi = 600, width = 6, height = 3)

if (FALSE) {
fit <- coxph(Surv(ost, ose) ~ oncotype, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)

fit <- coxph(Surv(rfst, rfse) ~ Tex + oncotype, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)
}

cat('\t\tTex\n')
tex_df <- use_df[use_df$group_tex != 'Medium',]
fit <- survfit(Surv(ost, ose) ~ group_tex, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_os_tex.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_tex, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Tex',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_rfs_tex.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('\t\tOncotype DX\n')
tex_df <- use_df[use_df$group_onco != 'Medium',]
fit <- survfit(Surv(ost, ose) ~ group_onco, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Oncotype DX',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_os_oncotype.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_onco, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Oncotype DX',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_rfs_oncotype.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat('\t\tComb\n')
tex_df <- use_df[use_df$group_comb != 'Medium',]
fit <- survfit(Surv(ost, ose) ~ group_comb, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Hybrid',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_os_hybrid.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

fit <- survfit(Surv(rfst, rfse) ~ group_comb, data = tex_df)
surv_gg <- ggsurvplot(fit, data = tex_df, pval = TRUE,
                      title = meno_state,
                      legend.title = 'Hybrid',
                      legend.labs = c('High', 'Low'))
png(paste(data_dir, meno_state, '_rfs_hybrid.png', sep = ''), res = 600, width = 6, height = 4, units = 'in')
print(surv_gg)
gar <- dev.off()

cat("Correlation...\n")
cor_gg <- ggscatter(df, x='Tex', y = 'oncotype', color = 'menopausal_State',
                    size = 1, alpha = 0.3,
                    add = 'reg.line',
                    conf.int = TRUE, 
                    palette = "jco",
                    rug = TRUE) +
  stat_cor(aes(color = menopausal_State), label.x = 6.4) + 
  labs(y = 'Oncotype DX', color = 'Menopausal states', fill = 'Menopausal states') +
  theme_classic() +
  theme(legend.position = 'bottom')
ggsave(paste(data_dir, 'colored_correlation.png', sep = ''), cor_gg, dpi = 600, width = 6, height = 4)

cor_gg <- ggscatter(df, x='Tex', y = 'oncotype', color = 'menopausal_State',
                    size = 1, alpha = 0.3,
                    facet.by = 'menopausal_State',
                    add = 'reg.line',
                    conf.int = TRUE, 
                    palette = "jco",
                    rug = TRUE) +
  stat_cor(aes(color = menopausal_State), label.x = 6.4) + 
  labs(y = 'Oncotype DX', color = 'Menopausal states', fill = 'Menopausal states') +
  theme_classic() +
  theme(legend.position = 'bottom')
ggsave(paste(data_dir, 'colored_facet_correlation.png', sep = ''), cor_gg, dpi = 600, width = 6, height = 4)
print(table(df$menopausal_State))

