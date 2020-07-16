#!/bin/env Rscript

#args: normal_stats out_prefix
args <- commandArgs(trailingOnly=T)

if (length(commandArgs) != 2){
	message("Usage: Rscript process_propca_stat.R python_output outprefix")
	q(save='no',status=1)
}

library(data.table)

stats <- fread(args[1],header=F,stringsAsFactors=F,data.table=F)
rownames(stats) <- stats[,1]
stats[,1] <- NULL

colnames(stats) <- paste0("PC",1:ncol(stats))

sdevs <- as.numeric(apply(stats,2,sd))

stats <- t(apply(stats,1,function(x){x / sdevs})) ** 2

write.table(stats,paste0(args[2],".chisq"),quote=F,row.names=T,col.names=F,sep="\t")

lambda_gc <- apply(stats,2,function(x){median(x)/qchisq(0.5,df=1)})

write.table(lambda_gc,paste0(args[2],".lambda_gc"),quote=F,row.names=T,col.names=F,sep="\t")

combined <- data.frame(rowSums(stats))
write.table(combined,paste0(args[2],".combined.chisq"),quote=F,row.names=T,col.names=F,sep="\t")

combined_lgc <- median(combined[,1])/qchisq(0.5,df=ncol(stats))

write.table(combined_lgc,paste0(args[2],".combined.lambda_gc"),quote=F,row.names=F,col.names=F)

probs <- ppoints(nrow(stats))
chisq_quants <- qchisq(probs,df=1)

qq <- sort(chisq_quants)
for (i in 1:ncol(stats)){
        qq <- cbind(qq,sort(stats[,i]))
}
qq <- data.frame(qq)
colnames(qq) <- c("chi-square",paste0("PC",1:ncol(stats)))

combined_quants <- qchisq(probs,df=5)
combined_df <- data.frame(cbind(sort(combined_quants),sort(combined[,1])))

library(tibble)

combined_df <- add_column(combined_df, PC = "Combined",.after=1)
colnames(combined_df) <- c("chi-square","PC","Magnitude")

library(ggplot2)
library(dplyr)
library(tidyr)

png(paste0(args[2],"_qqplot.png"),width=1920,height=1080,res=300)
first <- qq %>% gather(PC, Magnitude, starts_with("PC")) 
first <- rbind(first,data.frame(combined_df,check.names=F))
first %>% ggplot(aes(x=`chi-square`, y=Magnitude, colour=PC)) + geom_point() + theme_bw() + geom_abline(slope=1,intercept=0) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title=element_blank()) + xlab("Theoretical Quantiles") + ylab("Empirical Quantiles")
dev.off()

penalty <- nrow(stats) * (ncol(stats)+1)

stats_pval <- pchisq(stats,df=1,lower.tail=F)
combined_pval <- pchisq(combined[,1],df=5,lower.tail=F)
names(combined_pval) <- rownames(combined)

write.table(stats_pval,paste0(args[2],"_pval.no_penalty"),quote=F,row.names=T,col.names=F,sep="\t")
write.table(combined_pval,paste0(args[2],"_combined_pval.no_penalty"),quote=F,row.names=T,col.names=F,sep="\t")

stats_pval_adj <- stats_pval * penalty
combined_pval_adj <- combined_pval * penalty

stats_pval_adj[stats_pval_adj > 1] <- 1
combined_pval_adj[combined_pval_adj > 1] <- 1

stats_sig_rows <- apply(stats_pval_adj,1,function(x){any(x<0.05)})
#combined_pval_rows <- apply(combined_pval_adj,1,function(x){any(x<0.05)})

stats_sig <- stats_pval_adj[stats_sig_rows,]
combined_sig <- combined_pval_adj[combined_pval_adj < 0.05]

write.table(stats_sig,paste0(args[2],"_pval.sig"),quote=F,row.names=T,col.names=F,sep="\t")
write.table(combined_sig,paste0(args[2],"_combined_pval.sig"),quote=F,row.names=T,col.names=F,sep="\t")

