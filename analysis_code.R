############################################
library(data.table)
library(dplyr)
library(TEMR)

############################################
## 1. EUR tool SNP (P<5e‑8 + clump)
############################################
eur_iv <- fread("193_MDD_EUR.csv", sep = ",") |> as.data.frame()
IV_list <- eur_iv$SNP
cat("EUR IV =", length(IV_list), "\n")      # Recommendations ≥ 10

# EUR β/SE
eur_exp <- eur_iv |> transmute(SNP,
                               beta_ex = beta,
                               se_ex   = se)

############################################
## 2. Read the remaining three exposure files (EAS / AFR / HIS/SAS)
############################################
eas_exp_full <- fread("MDD_EAS_Neff.csv", sep = ",",
                      select = c("SNP","beta","se"))
afr_exp_full <- fread("MDD_AFR_Neff.csv", sep = ",",
                      select = c("SNP","beta","se"))
his_exp_full <- fread("MDD_HIS_Neff.csv", sep = ",",
                      select = c("SNP","beta","se"))
sas_exp_full <- fread("MDD_SAS_Neff.csv", sep = ",",
                      select = c("SNP","beta","se"))

eas_exp <- eas_exp_full[SNP %in% IV_list] |>
  transmute(SNP, beta_ex = beta, se_ex = se)
afr_exp <- afr_exp_full[SNP %in% IV_list] |>
  transmute(SNP, beta_ex = beta, se_ex = se)
his_exp <- his_exp_full[SNP %in% IV_list] |>
  transmute(SNP, beta_ex = beta, se_ex = se)
sas_exp <- sas_exp_full[SNP %in% IV_list] |>
  transmute(SNP, beta_ex = beta, se_ex = se)


############################################
## 3. Read 5 ending files (EUR / EAS / AFR / HIS); Here we shall take asthma as an example.
############################################
eur_out_full <- fread("GCST90399686.tsv", sep = ",",
                      select = c("SNP","beta","se"))
eas_out_full <- fread("GCST90399685.txt", sep = ",",
                      select = c("SNP","beta","se"))
afr_out_full <- fread("GCST90399683.tsv", sep = ",",
                      select = c("SNP","beta","se"))
his_out_full <- fread("GCST90399684.tsv", sep = ",",
                      select = c("SNP","beta","se"))
sas_out_full <- fread("GCST90399688.tsv ", sep = ",",
                      select = c("SNP","beta","se"))


eur_out <- eur_out_full[SNP %in% IV_list] |>
  transmute(SNP, beta_out = beta, se_out = se)
eas_out <- eas_out_full[SNP %in% IV_list] |>
  transmute(SNP, beta_out = beta, se_out = se)
afr_out <- afr_out_full[SNP %in% IV_list] |>
  transmute(SNP, beta_out = beta, se_out = se)
his_out <- his_out_full[SNP %in% IV_list] |>
  transmute(SNP, beta_out = beta, se_out = se)
sas_out <- sas_out_full[SNP %in% IV_list] |>
  transmute(SNP, beta_out = beta, se_out = se)


########################################################################

## 4. Retain only rows where SE > 0   (all data frames)
########################################################################
keep_se_pos <- function(df, beta_col, se_col){
  df[df[[se_col]] > 0 , ]
}

eur_exp <- keep_se_pos(eur_exp,"beta_ex","se_ex")
eas_exp <- keep_se_pos(eas_exp,"beta_ex","se_ex")
afr_exp <- keep_se_pos(afr_exp,"beta_ex","se_ex")
his_exp <- keep_se_pos(his_exp,"beta_ex","se_ex")
sas_exp <- keep_se_pos(sas_exp,"beta_ex","se_ex")

eur_out <- keep_se_pos(eur_out,"beta_out","se_out")
eas_out <- keep_se_pos(eas_out,"beta_out","se_out")
afr_out <- keep_se_pos(afr_out,"beta_out","se_out")
his_out <- keep_se_pos(his_out,"beta_out","se_out")
sas_out <- keep_se_pos(sas_out,"beta_out","se_out")

########################################################################
## 5. Intersection & Row Alignment
########################################################################
common <- Reduce(intersect,
                 list(eur_exp$SNP, eas_exp$SNP, afr_exp$SNP,
                      his_exp$SNP, sas_exp$SNP,
                      eur_out$SNP, eas_out$SNP, afr_out$SNP,
                      his_out$SNP, sas_out$SNP))
stopifnot(length(common) >= 3)

to_match <- function(df) df[match(common, df$SNP), ]
eur_exp <- to_match(eur_exp); eas_exp <- to_match(eas_exp)
afr_exp <- to_match(afr_exp); his_exp <- to_match(his_exp); sas_exp <- to_match(sas_exp)
eur_out <- to_match(eur_out); eas_out <- to_match(eas_out)
afr_out <- to_match(afr_out); his_out <- to_match(his_out); sas_out <- to_match(sas_out)

########################################################################
## 6. Constructing the β/SE matrix
########################################################################
betaXG <- cbind(EAS=eas_exp$beta_ex,
                AFR=afr_exp$beta_ex,
                HIS=his_exp$beta_ex,
                SAS=sas_exp$beta_ex,
                auxiliary=eur_exp$beta_ex)

betaYG <- cbind(EAS=eas_out$beta_out,
                AFR=afr_out$beta_out,
                HIS=his_out$beta_out,
                SAS=sas_out$beta_out,
                auxiliary=eur_out$beta_out)

sebetaXG <- cbind(EAS=eas_exp$se_ex,
                  AFR=afr_exp$se_ex,
                  HIS=his_exp$se_ex,
                  SAS=sas_exp$se_ex,
                  auxiliary=eur_exp$se_ex)

sebetaYG <- cbind(EAS=eas_out$se_out,
                  AFR=afr_out$se_out,
                  HIS=his_out$se_out,
                  SAS=sas_out$se_out,
                  auxiliary=eur_out$se_out)

rownames(betaXG) <- rownames(betaYG) <- common

########################################################################
## 7. Minimum QC: SE lower bound, delete βₓG = 0 row
########################################################################
floor_se <- 1e-4
sebetaXG[sebetaXG < floor_se] <- floor_se
sebetaYG[sebetaYG < floor_se] <- floor_se

flag_zero <- apply(betaXG == 0, 1, any)    # 
if(any(flag_zero)){
  cat("delete β_XG = 0 SNP =", sum(flag_zero), "\n")
  betaXG   <- betaXG [!flag_zero, ]
  betaYG   <- betaYG [!flag_zero, ]
  sebetaXG <- sebetaXG[!flag_zero, ]
  sebetaYG <- sebetaYG[!flag_zero, ]
}
cat("Number of SNPs included in TEMR =", nrow(betaXG), "\n")
stopifnot(nrow(betaXG) >= 3)

## ================================================================
## 8 Export SNP β/SE (exposure end) to a wide table
##  ---------------------------------------------------------------
snp_table <- data.frame(
  SNP        = rownames(betaXG),
  beta_EAS   = betaXG[,"EAS"],   se_EAS   = sebetaXG[,"EAS"],
  beta_AFR   = betaXG[,"AFR"],   se_AFR   = sebetaXG[,"AFR"],
  beta_HIS   = betaXG[,"HIS"],   se_HIS   = sebetaXG[,"HIS"],
  beta_SAS   = betaXG[,"SAS"],   se_SAS   = sebetaXG[,"SAS"],
  beta_EUR   = betaXG[,"auxiliary"],       # EUR＝auxiliary
  se_EUR     = sebetaXG[,"auxiliary"]
)

# Save to working directory; filename or path may be modified
write.csv(snp_table,
          file = "TEMR_input_SNPs_beta_SE.csv",
          row.names = FALSE)
## ================================================================


########################################################################
## 9. ρ (rg) Here, the rg parameter in SAS is set to NULL. Please run the analysis with rho=NULL to obtain the final SAS values.
########################################################################
rg_list <- c(EAS=0.526, AFR=0.619, HIS=0.660, SAS=0.650)

########################################################################
## 10. TEMR Multi-Optimizer Iteration
########################################################################
opt_vec <- c("BFGS","Nelder-Mead","CG","L-BFGS-B")

for (mm in opt_vec) {
  message(">>> Optimizer: ", mm)
  dir.create(mm, showWarnings = FALSE)
  
  core_all   <- list()
  method_all <- list()
  
  for (tg in names(rg_list)) {
    
    bxg <- betaXG   [, c(tg,"auxiliary")]
    byg <- betaYG   [, c(tg,"auxiliary")]
    sxg <- sebetaXG [, c(tg,"auxiliary")]
    syg <- sebetaYG [, c(tg,"auxiliary")]
    
    rho2 <- matrix(c(1, rg_list[tg],
                     rg_list[tg], 1),
                   2, 2,
                   dimnames=list(c(tg,"auxiliary"),
                                 c(tg,"auxiliary")))
    
    res <- tryCatch(
      TEMR(bxg, byg, sxg, syg, rho=rho2, meth=mm),
      error=function(e){
        message("  ", tg, ": ρ failed to be set; switching to rho=NULL ")
        TEMR(bxg, byg, sxg, syg, rho=NULL, meth=mm)
      })
    
    core_all[[tg]] <- data.frame(population=tg,
                                 beta = res$resultall[[tg]]$beta,
                                 stat = res$resultall[[tg]]$stat,
                                 pval = res$resultall[[tg]]$pvalue)
    
    method_all[[tg]] <- rbind(
      data.frame(population=tg, method="IVW",
                 t(res$betaXY_IVW   [tg, ])),
      data.frame(population=tg, method="Egger",
                 t(res$betaXY_Egger[tg, ])),
      data.frame(population=tg, method="Final",
                 t(res$betaXY_Final[tg, ]))
    )
  }
  
  write.csv(do.call(rbind, core_all),
            file.path(mm, "TEMR_pairwise.csv"),    row.names=FALSE)
  write.csv(do.call(rbind, method_all),
            file.path(mm, "MRmethods_pairwise.csv"),row.names=FALSE)
  
  print(do.call(rbind, core_all))
}

message("All optimisation processes have been completed, and the results have been written to the corresponding folders.")



