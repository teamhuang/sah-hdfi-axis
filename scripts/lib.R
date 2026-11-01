## lib.R -- shared helpers. Run all scripts with cwd = project dir (CJK paths break literal strings).
suppressWarnings(suppressMessages({library(limma)}))

BASE    <- getwd()
RES     <- file.path(BASE,"results")
FIG     <- file.path(BASE,"figures")
DAT     <- file.path(BASE,"data")
SCRATCH <- "C:/Users/Admin/.qoderworkcn/workspace/muycr5n2rnyoywup/sah_rds"  # ASCII: saveRDS/gzfile OK
for(d in c(RES,FIG,SCRATCH)) dir.create(d, showWarnings=FALSE, recursive=TRUE)

# Okabe-Ito colorblind-safe palette + Arial base theme
OI <- c(blue="#0072B2", orange="#E69F00", sky="#56B4E9", green="#009E73",
        yellow="#F0E442", black="#000000", vermil="#D55E00", purple="#CC79A7")
theme_pub <- function(base=10){ theme_classic(base_size=base, base_family="Arial") +
  theme(axis.text=element_text(colour="black"), axis.title=element_text(colour="black"),
        legend.key=element_blank(), plot.title=element_text(face="bold")) }

save_png_pdf <- function(p, name, w, h){
  png(file.path(FIG, paste0(name,".png")), width=w, height=h, units="cm", res=300, bg="white")
  print(p); dev.off()
  tryCatch({ pdf(file.path(FIG, paste0(name,".pdf")), width=w/2.54, height=h/2.54);
             print(p); dev.off() }, error=function(e) cat("[pdf skip]",name,conditionMessage(e),"\n"))
}

load_cohort <- function(tag){
  e <- as.matrix(read.delim(file.path(DAT,sprintf("expr_linear_%s.tsv",tag)),
                            row.names=1, check.names=FALSE))
  p <- read.delim(file.path(DAT,sprintf("pheno_%s.tsv",tag)), check.names=FALSE)
  p <- p[p$gsm %in% colnames(e),]; rownames(p) <- p$gsm; p <- p[colnames(e),]
  list(expr=e, pheno=p)
}

run_de <- function(tag){
  co <- load_cohort(tag); e <- co$expr; p <- co$pheno
  lg <- normalizeBetweenArrays(log2(e+1), method="quantile")
  keep <- rowSums(lg>1) >= max(2, ceiling(0.1*ncol(lg)))
  lg <- lg[keep,,drop=FALSE]
  gname <- rownames(lg)                            # captured by position
  grp <- factor(p$group, levels=c("Control","SAH"))
  sex <- factor(ifelse(is.na(p$sex),"NA",p$sex))
  design <- model.matrix(~ sex + grp)
  fit <- eBayes(lmFit(lg, design))
  # manual statistics (lmFit/eBayes preserve feature ORDER -> safe against name mangling)
  tt <- data.frame(
    gene = gname,
    logFC = as.numeric(fit$coefficients[, "grpSAH"]),
    AveExpr = as.numeric(fit$Amean),
    t = as.numeric(fit$t[, "grpSAH"]),
    P.Value = as.numeric(fit$p.value[, "grpSAH"]),
    stringsAsFactors = FALSE)
  tt$adj.P.Val <- p.adjust(tt$P.Value, method="BH")
  tt <- tt[order(tt$P.Value),]
  stopifnot(nrow(tt)==length(gname), !any(duplicated(tt$gene)))
  list(gene=gname, lg=lg, design=design, fit=fit, tt=tt, pheno=p,
       deg=tt[tt$adj.P.Val<0.05 & abs(tt$logFC)>1,], tag=tag)
}
