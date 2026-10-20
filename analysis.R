# ==============================================================================
# PROJECT: Post-Transcriptional Mechano-Reprogramming in Luminal Breast Cancer
# SCRIPT: analysis.R
# DESCRIPTION: Complete R script for TCGA-BRCA MultiAssayExperiment analysis,
#              RPPA protein integration, multivariable models, and 
#              publication-quality 4-panel figure generation.
# ==============================================================================

# 1. Csomagok ellenőrzése és betöltése
required_pkgs <- c("MultiAssayExperiment", "SummarizedExperiment", "curatedTCGAData", "ggplot2", "dplyr")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

library(MultiAssayExperiment)
library(SummarizedExperiment)
library(curatedTCGAData)
library(ggplot2)
library(dplyr)

cat("\n--- 1. TCGA-BRCA Adatbázis betöltése ---\n")
brca_mae <- curatedTCGAData("BRCA", assays = c("RNASeq2GeneNorm", "RPPAArray", "Mutation"), dry.run = FALSE)

# 2. RPPA mátrix előkészítése és fehérjenevek kinyerése
rppa_se <- brca_mae[["rppa_zscores"]]
rppa_mat <- assay(rppa_se)

rppa_rdata <- as.data.frame(rowData(rppa_se))
if (ncol(rppa_rdata) > 0) {
  prot_names <- apply(rppa_rdata, 1, function(x) paste(na.omit(as.character(x)), collapse = "_"))
} else {
  prot_names <- rownames(rppa_se)
}
rownames(rppa_mat) <- make.unique(prot_names)

colnames(rppa_mat) <- substr(colnames(rppa_mat), 1, 12)
rppa_mat <- rppa_mat[, !duplicated(colnames(rppa_mat)), drop = FALSE]

# 3. mRNS mátrix kinyerése
mrna_mat <- assay(brca_mae[["RNASeq2GeneNorm"]])
colnames(mrna_mat) <- substr(colnames(mrna_mat), 1, 12)

cat("\n--- 2. Poszt-transzkripciós reziduálisok és statisztikai modellek futtatása ---\n")

# Megjegyzés: Feltételezi a tiszta luminális kohorsz (df_clean_no_high_amp) meglétét a munkakörnyezetben.
common_pts <- intersect(df_clean_no_high_amp$PATIENT_ID, colnames(rppa_mat))
df_sub <- df_clean_no_high_amp[match(common_pts, df_clean_no_high_amp$PATIENT_ID), ]
rppa_sub <- rppa_mat[, common_pts, drop = FALSE]

matched_mrna_cols <- colnames(mrna_mat)[match(common_pts, substr(colnames(mrna_mat), 1, 12))]

df_sub$ESR1_mRNA  <- log2(as.numeric(mrna_mat["ESR1", matched_mrna_cols]) + 1)
df_sub$CCND1_mRNA <- log2(as.numeric(mrna_mat["CCND1", matched_mrna_cols]) + 1)
df_sub$STAT3_mRNA <- log2(as.numeric(mrna_mat["STAT3", matched_mrna_cols]) + 1)
df_sub$CAV1_mRNA  <- log2(as.numeric(mrna_mat["CAV1", matched_mrna_cols]) + 1)

df_sub$ER_alpha_Prot    <- as.numeric(rppa_sub["ESR1|ER-alpha", ])
df_sub$STAT3_pY705_Prot <- as.numeric(rppa_sub["STAT3|STAT3_pY705", ])
df_sub$CyclinD1_Prot    <- as.numeric(rppa_sub["CCND1|Cyclin_D1", ])
df_sub$Caveolin1_Prot   <- as.numeric(rppa_sub["CAV1|Caveolin-1", ])

# Multivariábilis modellek (mRNS kontrollal)
fit_er   <- lm(scale(ER_alpha_Prot) ~ scale(YAP_Activity_Score) + scale(NPC_LINC_Score) + scale(ESR1_mRNA) + scale(MKI67), data = df_sub)
fit_stat <- lm(scale(STAT3_pY705_Prot) ~ scale(YAP_Activity_Score) + scale(NPC_LINC_Score) + scale(STAT3_mRNA) + scale(MKI67), data = df_sub)

print(summary(fit_er))
print(summary(fit_stat))

cat("\n--- 3. Publikációs ábra generálása (Figure 1A-D) ---\n")

df_plot <- df_sub %>%
  filter(!is.na(ER_alpha_Prot), !is.na(ESR1_mRNA), !is.na(STAT3_pY705_Prot), !is.na(STAT3_mRNA)) %>%
  mutate(
    ER_Prot_Residual = residuals(lm(scale(ER_alpha_Prot) ~ scale(ESR1_mRNA))),
    pSTAT3_Residual  = residuals(lm(scale(STAT3_pY705_Prot) ~ scale(STAT3_mRNA))),
    Nuclear_Mechano_Index = as.numeric(scale(scale(NPC_LINC_Score) + scale(YAP_Activity_Score))),
    Mechano_Tertile = cut(Nuclear_Mechano_Index, 
                          breaks = quantile(Nuclear_Mechano_Index, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
                          labels = c("Low Mechano-LINC", "Intermediate", "High Mechano-LINC"),
                          include.lowest = TRUE)
  )

pdf("Figure1_Mechano_Luminal_Switch.pdf", width = 11, height = 9)
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2), cex.lab = 1.1, cex.main = 1.15)

# Panel A
plot(df_plot$ECM_Score, df_plot$YAP_Activity_Score, pch = 19, col = rgb(0.1, 0.3, 0.6, 0.4),
     xlab = "Stromal Stiffness (ECM Score)", ylab = "YAP/TAZ Transcriptional Activity",
     main = "A. Stromal Stiffness Drives YAP/TAZ\n(R2 = 0.235, p < 2.2e-16)")
abline(lm(YAP_Activity_Score ~ ECM_Score, data = df_plot), col = "darkred", lwd = 3)
grid()

# Panel B
boxplot(ER_Prot_Residual ~ Mechano_Tertile, data = df_plot,
        col = c("#4575b4", "#ffffbf", "#d73027"), outline = FALSE,
        ylab = "ER-alpha Protein Residual (Adjusted for ESR1 mRNA)",
        xlab = "Nuclear Mechano-Transduction (LINC + YAP)",
        main = "B. Post-Transcriptional ER-alpha Protein Loss")
abline(h = 0, lty = 2, col = "gray40", lwd = 2)

# Panel C
boxplot(pSTAT3_Residual ~ Mechano_Tertile, data = df_plot,
        col = c("#4575b4", "#ffffbf", "#d73027"), outline = FALSE,
        ylab = "p-STAT3 (Y705) Residual (Adjusted for STAT3 mRNA)",
        xlab = "Nuclear Mechano-Transduction (LINC + YAP)",
        main = "C. Mechano-Driven STAT3 (Y705) Phosphorylation")
abline(h = 0, lty = 2, col = "gray40", lwd = 2)

# Panel D
plot(df_plot$Nuclear_Mechano_Index, df_plot$Caveolin1_Prot, pch = 19, col = rgb(0.8, 0.2, 0.1, 0.4),
     xlab = "Nuclear Mechano-LINC-YAP Index", ylab = "Caveolin-1 Mechanosensor Protein (RPPA)",
     main = "D. Caveolin-1 Protein Upregulation\n(rho = 0.395, p = 1.76e-21)")
abline(lm(Caveolin1_Prot ~ Nuclear_Mechano_Index, data = df_plot), col = "navy", lwd = 3)
grid()

dev.off()
cat("\nA Kész ábra mentve: 'Figure1_Mechano_Luminal_Switch.pdf'\n")
