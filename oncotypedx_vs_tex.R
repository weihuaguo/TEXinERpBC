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
suppressMessages(library(survminer))
suppressMessages(library(survival))
suppressMessages(library(survMisc))
suppressMessages(library(rstatix))

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

