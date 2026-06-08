# Supplementary: Sensitivity analysis for reviewer response

# Requires: 01_data_cleaning_processing.R to have been run first
# Dependencies: utils_ecology_myxo.R, vegan, phyloseq

## 0. Preparation
### 0-1. Load library 
source("utils_ecology_myxo.R")
load_ecology_packages()
library(vegan)

# ======================================================================
# PART 1: Japan TNS exclusion sensitivity
# ======================================================================

cat("=== PART 1: Japan TNS sensitivity ===\n\n")

## 1-1. Load data
DF_seq          <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds")
community_RT_orig <- readRDS("../_data/_processed_data/community_RT.rds")
community_RG_orig <- readRDS("../_data/_processed_data/community_RG.rds")

## 1-2. Identify TNS specimens within Japan
Jp_all   <- DF_seq %>% filter(Code == "Jp")
Jp_noTNS <- Jp_all %>% filter(!grepl("^TNS", final_colno))

cat("Jp specimens (total)      :", nrow(Jp_all), "\n")
cat("Jp specimens (TNS origin) :", nrow(Jp_all) - nrow(Jp_noTNS), "\n")
cat("Jp specimens (excl. TNS)  :", nrow(Jp_noTNS), "\n\n")

## 1-3. Rebuild community matrices with TNS-excluded Jp
DF_noTNS <- DF_seq %>%
  filter(Code != "Jp") %>%
  bind_rows(Jp_noTNS)

community_RT_noTNS <- create_community_matrix(DF_noTNS, "Code", "RT")
community_RG_noTNS <- create_community_matrix(DF_noTNS, "Code", "RG")

## 1-4. Calculate endemicity (endemic RT/RG count and rate) per region
# An RT/RG is "endemic" if it occurs in exactly one region.
calc_endemicity <- function(community, level = "RT") {
  cm <- as.data.frame(community)
  # Identify endemic RT/RGs: those present in exactly one region
  region_count <- colSums(cm > 0)
  is_endemic   <- (region_count == 1)
  
  result <- data.frame(
    Code              = rownames(cm),
    Total_specimens   = rowSums(cm),                                  # total barcoded specimens per region
    Endemic_specimens = rowSums(cm[, is_endemic, drop = FALSE]),      # specimens with endemic RT/RGs
    stringsAsFactors  = FALSE
  )
  result$Endemic_pct_spec <- round(result$Endemic_specimens / result$Total_specimens * 100, 1)
  # Rename columns to include level (RT or RG)
  colnames(result)[2:4] <- paste0(c("Total_spec_", "Endemic_spec_", "Endemic_pct_spec_"), level)
  return(result)
}

end_RT_orig  <- calc_endemicity(community_RT_orig,  "RT")
end_RT_noTNS <- calc_endemicity(community_RT_noTNS, "RT")
end_RG_orig  <- calc_endemicity(community_RG_orig,  "RG")
end_RG_noTNS <- calc_endemicity(community_RG_noTNS, "RG")

## 1-5. Print comparison for Jp
cat("--- RT endemicity: Jp ---\n")
cat("Original:\n")
print(end_RT_orig[end_RT_orig$Code == "Jp", ])
cat("TNS excluded:\n")
print(end_RT_noTNS[end_RT_noTNS$Code == "Jp", ])

cat("\n--- RG endemicity: Jp ---\n")
cat("Original:\n")
print(end_RG_orig[end_RG_orig$Code == "Jp", ])
cat("TNS excluded:\n")
print(end_RG_noTNS[end_RG_noTNS$Code == "Jp", ])

## 1-6. Save full comparison table (all regions, for completeness)
out_jp <- bind_rows(
  cbind(Analysis = "Original", end_RT_orig) %>%
    left_join(cbind(Analysis = "Original", end_RG_orig), by = c("Analysis", "Code")),
  cbind(Analysis = "No_TNS", end_RT_noTNS) %>%
    left_join(cbind(Analysis = "No_TNS", end_RG_noTNS), by = c("Analysis", "Code"))
)

write.csv(out_jp, "../_results/SA_Japan_TNS_endemicity.csv", row.names = FALSE)
cat("\nSaved: ../_results/SA_Japan_TNS_endemicity.csv\n")

# ======================================================================
# PART 2: GAP subsampling sensitivity (999 iterations, n = 300)
# ======================================================================
cat("\n=== PART 2: GAP subsampling sensitivity ===\n\n")

## 2-1. Parameters
N_ITER        <- 999
SUBSAMPLE_SIZE <- 300
set.seed(42)   # reproducibility

## 2-2. Load data (reload to ensure clean state)
DF_seq          <- readRDS("../_data/_processed_data/DF_biogeo_seq.rds")
community_RT_orig <- readRDS("../_data/_processed_data/community_RT.rds")
community_RG_orig <- readRDS("../_data/_processed_data/community_RG.rds")

GAP_all  <- DF_seq %>% filter(Code == "GAP")
DF_nonGAP <- DF_seq %>% filter(Code != "GAP")

cat("GAP specimens (original)   :", nrow(GAP_all), "\n")
cat("Subsampling to             :", SUBSAMPLE_SIZE, "specimens\n")
cat("Iterations                 :", N_ITER, "\n\n")

## 2-3. Rarefaction target: minimum specimen count across all regions
#  Computed from original matrices so the reference depth is consistent.
min_n_RT <- min(rowSums(community_RT_orig))
min_n_RG <- min(rowSums(community_RG_orig))
cat("Rarefaction depth (RT)     :", min_n_RT, "\n")
cat("Rarefaction depth (RG)     :", min_n_RG, "\n\n")

## 2-4. Compute reference values from original (full) dataset
d4_orig     <- as.matrix(vegdist(community_RT_orig, method = "jaccard", binary = FALSE))
d5_orig     <- as.matrix(vegdist(community_RT_orig, method = "horn"))
rar_RT_orig <- rarefy(community_RT_orig, sample = min_n_RT)
rar_RG_orig <- rarefy(community_RG_orig, sample = min_n_RG)

original_values <- data.frame(
  Metric   = c("Rarefied_RT", "Rarefied_RG", "D4_jaccard_GAP_mean", "D5_horn_GAP_mean"),
  Original = c(
    round(rar_RT_orig["GAP"], 2),
    round(rar_RG_orig["GAP"], 2),
    round(mean(d4_orig["GAP", rownames(d4_orig) != "GAP"]), 4),
    round(mean(d5_orig["GAP", rownames(d5_orig) != "GAP"]), 4)
  )
)
cat("Original GAP values:\n")
print(original_values)
cat("\n")

## 2-5. Subsampling loop
results <- data.frame(
  iter              = seq_len(N_ITER),
  Rarefied_RT       = NA_real_,
  Rarefied_RG       = NA_real_,
  D4_jaccard_GAP_mean = NA_real_,
  D5_horn_GAP_mean  = NA_real_
)

cat("Running iterations (this may take a few minutes)...\n")
pb <- txtProgressBar(min = 0, max = N_ITER, style = 3)

for (i in seq_len(N_ITER)) {
  
  # Subsample GAP specimens
  GAP_sub <- GAP_all %>% slice_sample(n = SUBSAMPLE_SIZE)
  
  # Rebuild full specimen data with subsampled GAP
  DF_sub <- bind_rows(DF_nonGAP, GAP_sub)
  
  # Rebuild community matrices
  cm_RT <- create_community_matrix(DF_sub, "Code", "RT")
  cm_RG <- create_community_matrix(DF_sub, "Code", "RG")
  
  # Rarefied richness for GAP at fixed depth
  rar_RT_i <- rarefy(cm_RT, sample = min_n_RT)
  rar_RG_i <- rarefy(cm_RG, sample = min_n_RG)
  results$Rarefied_RT[i] <- rar_RT_i["GAP"]
  results$Rarefied_RG[i] <- rar_RG_i["GAP"]
  
  # D4: weighted Jaccard (GAP vs all other 10 regions)
  d4_i <- as.matrix(vegdist(cm_RT, method = "jaccard", binary = FALSE))
  results$D4_jaccard_GAP_mean[i] <- mean(d4_i["GAP", rownames(d4_i) != "GAP"])
  
  # D5: Horn index (GAP vs all other 10 regions)
  d5_i <- as.matrix(vegdist(cm_RT, method = "horn"))
  results$D5_horn_GAP_mean[i] <- mean(d5_i["GAP", rownames(d5_i) != "GAP"])
  
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\nIterations complete.\n\n")

## 2-6. Summarize with 95% CI
summarize_ci <- function(x, label) {
  data.frame(
    Metric  = label,
    Mean    = round(mean(x), 4),
    SD      = round(sd(x),   4),
    CI_2.5  = round(quantile(x, 0.025), 4),
    CI_97.5 = round(quantile(x, 0.975), 4)
  )
}

subsampled_summary <- bind_rows(
  summarize_ci(results$Rarefied_RT,          "Rarefied_RT"),
  summarize_ci(results$Rarefied_RG,          "Rarefied_RG"),
  summarize_ci(results$D4_jaccard_GAP_mean,  "D4_jaccard_GAP_mean"),
  summarize_ci(results$D5_horn_GAP_mean,     "D5_horn_GAP_mean")
)

## 2-7. Merge original values with subsampled CI
final_table <- left_join(original_values, subsampled_summary, by = "Metric")

cat("=== Summary: original vs subsampled GAP (n=300, 999 iter) ===\n")
print(final_table, row.names = FALSE)

## 2-8. Save results
saveRDS(results,     "../_results/SA_GAP_subsampling_raw.rds")
write.csv(final_table, "../_results/SA_GAP_subsampling_summary.csv", row.names = FALSE)

cat("\nSaved:\n")
cat("  ../_results/SA_GAP_subsampling_summary.csv\n")
cat("  ../_results/SA_GAP_subsampling_raw.rds\n")
cat("\n=== SA_sensitivity_analysis.R complete ===\n")
