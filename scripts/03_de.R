## 03_de.R -- clean DE for both cohorts, cross-cohort consistency, cache objects
source("scripts/lib.R")
r36 <- run_de("GSE36791")
r73 <- run_de("GSE73378")
for(r in list(r36, r73)){
  cat(sprintf("%s: %d genes tested | DEGs(FDR<.05,|lfc|>1)=%d (up=%d down=%d)\n",
      r$tag, nrow(r$tt), nrow(r$deg), sum(r$deg$logFC>0), sum(r$deg$logFC<0)))
  write.csv(r$tt, file.path(RES, sprintf("DEG_%s.csv", r$tag)), row.names=FALSE)
  write.csv(r$deg, file.path(RES, sprintf("DEGsig_%s.csv", r$tag)), row.names=FALSE)
}
cat("sample GSE36791 top DEGs:\n"); print(head(r36$tt[,c("gene","logFC","P.Value","adj.P.Val")],8))

# ---- cross-cohort consistency ----
g36 <- r36$tt; g73 <- r73$tt
idx <- match(g36$gene, g73$gene); ok <- !is.na(idx)
lfc36 <- g36$logFC[ok]; lfc73 <- g73$logFC[idx[ok]]; p73 <- g73$P.Value[idx[ok]]
cat(sprintf("\ncommon genes=%d | logFC Pearson r=%.3f Spearman rho=%.3f\n",
    sum(ok), cor(lfc36,lfc73), cor(lfc36,lfc73,method="spearman")))
d36 <- r36$deg
i2 <- match(d36$gene, g73$gene); k2 <- !is.na(i2)
conc <- sum(p73[i2[k2]]<0.05 & sign(g73$logFC[i2[k2]])==sign(d36$logFC[k2]))
cat(sprintf("GSE36791 DEGs testable in GSE73378=%d, concordant nominal-significant=%d (%.1f%%)\n",
    sum(k2), conc, 100*conc/sum(k2)))
write.csv(data.frame(metric=c("n_common","pearson_r","spearman_rho","n_deg36_testable","n_concordant","pct_concordant"),
    value=c(sum(ok),cor(lfc36,lfc73),cor(lfc36,lfc73,method="spearman"),sum(k2),conc,100*conc/sum(k2))),
    file.path(RES,"cross_cohort_consistency.csv"), row.names=FALSE)

saveRDS(list(GSE36791=list(tt=r36$tt,deg=r36$deg,lg=r36$lg,pheno=r36$pheno,design=r36$design,fit=r36$fit)),
        file.path(SCRATCH,"r36.rds"))
saveRDS(list(GSE73378=list(tt=r73$tt,deg=r73$deg,lg=r73$lg,pheno=r73$pheno,design=r73$design,fit=r73$fit)),
        file.path(SCRATCH,"r73.rds"))
cat("\nDONE 03_de\n")
