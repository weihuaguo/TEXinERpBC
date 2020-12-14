library(fgsea)
library(ggplot2)
library(dplyr)
library(tidyr)

gmt_dir <- '/home/weihua/mnts/group_plee/Weihua/PD1_CD39_Public_Data/final_used_data'
gmt_file <- system.file(gmt_dir, 'signature_v4.gmt', package='fgsea')
tex_pub_sigs <- gmtPathways(paste(gmt_dir, "/signature_v5.gmt", sep = ""))
print(str(tex_pub_sigs))
tex_sig <- read.table(paste(gmt_dir, "/ranked_tex_sig_long_adjp010_pos_neg.rnk", sep = ""),
		      header=F, colClasses = c('character', 'numeric'))
print(head(tex_sig))
# tex_sig <- setNames(tex_sig$V1, tex_sig$V2)
colnames(tex_sig) <- c('ID', 't')
tex_sig <- setNames(tex_sig$t, tex_sig$ID)

fgseaRes <- fgseaMultilevel(pathways=tex_pub_sigs, stats=tex_sig, minSize=15, maxSize=500)
print(head(fgseaRes))
fgsea_tidy <- as.data.frame(fgseaRes %>% arrange(padj))
# write.csv(as.data.frame(fgsea_tidy), paste(gmt_dir, "/fgsea_public_signatures_results.csv", sep = ""))

print(fgsea_tidy)
for (gs in names(tex_pub_sigs)) {
	en_gg <- plotEnrichment(tex_pub_sigs[[gs]], tex_sig) + labs(title=gs)
	ggsave(paste(gmt_dir, "/", gs, "_enrichment_plot.png", sep = ""), en_gg,
	       dpi=300, width=6, height=4)
}
