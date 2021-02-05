# For oncotype DX vs Tex
# Weihua Guo, Ph.D
# 02/04/2021

rm(list = ls())
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(cocor))
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(survMisc))

data_dir <- 'C:/Users/wguo/OneDrive - City of Hope National Medical Center/tmp_works/oncotype_dx_pre_vs_post/'
bc_type <- 'ER'
df <- read.csv(paste(data_dir, 'metabric_', bc_type, '_tex_oncotype_sigscore.csv', sep = ''), row.names = 1)

pre_df <- df[df$menopausal_State == 'pre',]
post_df <- df[df$menopausal_State == 'post',]
cat('Hazard ratio...\n')
fit <- coxph(Surv(ost, ose) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)

fit <- coxph(Surv(rfst, rfse) ~ menopausal_State + Tex + oncotype, data = df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, 'all_together_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 4)


cat('\tPre\n')
use_df <- pre_df
meno_state <- 'pre'
use_df$group_onco <- ifelse(use_df$oncotype > quantile(use_df$oncotype, 0.75), 'High', 
                            ifelse(use_df$oncotype < quantile(use_df$oncotype, 0.25), 'Low', 'Medium'))
use_df$group_tex <- ifelse(use_df$Tex > quantile(use_df$Tex, 0.75), 'High', 
                           ifelse(use_df$Tex < quantile(use_df$Tex, 0.25), 'Low', 'Medium'))

tex_coxph <- coxph(Surv(ost, ose) ~ Tex, data = use_df)
tex_coxph <- cutp(tex_coxph)$Tex
colnames(tex_coxph)[2:ncol(tex_coxph)] <- paste('tex_', colnames(tex_coxph)[2:ncol(tex_coxph)], sep = '')
merge_tex_cutp <- merge(tex_coxph, use_df, by = 'Tex')
onco_coxph <- coxph(Surv(ost, ose) ~ oncotype, data = use_df)
onco_coxph <- cutp(onco_coxph)$oncotype
colnames(onco_coxph)[2:ncol(onco_coxph)] <- paste('onco_', colnames(onco_coxph)[2:ncol(onco_coxph)], sep = '')
merge_onco_cutp <- merge(onco_coxph, use_df, by = 'oncotype')
merge_cutp <- merge(merge_tex_cutp, merge_onco_cutp, by = 'pid')

cor_gg <- ggscatter(merge_cutp, x='tex_U', y = 'onco_U', 
                    size = 1, alpha = 0.3,
                    add = 'reg.line',
                    conf.int = TRUE, 
                    rug = TRUE) +
  stat_cor() + 
  labs(y = 'log-rank score\n(Oncotype DX)', x = 'log-rank score\n(Tex)') +
  theme_classic() +
  theme(legend.position = 'bottom')
ggsave(paste(data_dir, meno_state, '_cutp_correlation.png', sep = ''), cor_gg, dpi = 600, width = 6, height = 4)
stop('TEST')

group_cts <- as.data.frame(table(use_df$group_onco, use_df$group_tex))
colnames(group_cts) <- c('Oncotype', 'Tex', 'Num')
cts_gg <- ggplot(group_cts, aes_string(x = 'Tex', y = 'Oncotype', color = 'Num', size = 'Num')) +
  geom_point() +
  scale_color_distiller(palette = 'Spectral') +
  labs(title = meno_state) +
  theme_bw() +
  theme(legend.position = 'top')
ggsave(paste(data_dir, meno_state, '_group_onco_vs_tex.png', sep = ''), cts_gg, dpi = 600, width = 6, height = 3)

fit <- coxph(Surv(ost, ose) ~ group_tex + group_onco, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)

fit <- coxph(Surv(rfst, rfse) ~ group_tex + group_onco, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)

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

cat('\tPost\n')
use_df <- post_df
meno_state <- 'post'
use_df$group_onco <- ifelse(use_df$oncotype > quantile(use_df$oncotype, 0.75), 'High', 
                            ifelse(use_df$oncotype < quantile(use_df$oncotype, 0.25), 'Low', 'Medium'))
use_df$group_tex <- ifelse(use_df$Tex > quantile(use_df$Tex, 0.75), 'High', 
                           ifelse(use_df$Tex < quantile(use_df$Tex, 0.25), 'Low', 'Medium'))

group_cts <- as.data.frame(table(use_df$group_onco, use_df$group_tex))
colnames(group_cts) <- c('Oncotype', 'Tex', 'Num')
cts_gg <- ggplot(group_cts, aes_string(x = 'Tex', y = 'Oncotype', color = 'Num', size = 'Num')) +
  geom_point() +
  scale_color_distiller(palette = 'Spectral') +
  labs(title = meno_state) +
  theme_bw() +
  theme(legend.position = 'top')
ggsave(paste(data_dir, meno_state, '_group_onco_vs_tex.png', sep = ''), cts_gg, dpi = 600, width = 6, height = 3)

fit <- coxph(Surv(ost, ose) ~ group_tex + group_onco, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_os_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)

fit <- coxph(Surv(rfst, rfse) ~ group_tex + group_onco, data = use_df)
forest_gg <- ggforest(fit)
ggsave(paste(data_dir, meno_state, '_rfs_hr.png', sep = ''), forest_gg, dpi = 600, width = 6, height = 3)

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

