library(methylclock)
Train_data <- readRDS("~/Documents/Forensic Project/preExistingClocks.Train.RDS")
head(Train_data)

# Run Horvath clock specifically
# If your data has actual chronological ages
# (e.g., in a row called "age" or as a separate vector)

# ===== LOAD DATA =====

cat("Data dimensions:", dim(Train_data), "\n")
cat("Samples (columns):", ncol(Train_data), "\n")


# ===== RUN HORVATH CLOCK (WITHOUT AGE) =====
cat("\nRunning Horvath clock...\n")

horvath_results <- DNAmAge(Train_data, clocks = "Horvath")

cat("✓ Done!\n")
cat("Results dimensions:", dim(horvath_results), "\n")

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(horvath_results)) {
  actual_ages <- actual_ages[1:nrow(horvath_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
horvath_results$actual_age <- actual_ages

# Check
head(horvath_results)

# ===== CALCULATE METRICS =====
horvath_results$error <- horvath_results$Horvath - horvath_results$actual_age
horvath_results$abs_error <- abs(horvath_results$error)

mae <- mean(horvath_results$abs_error, na.rm = TRUE)
rmse <- sqrt(mean(horvath_results$error^2, na.rm = TRUE))
correlation <- cor(horvath_results$Horvath, horvath_results$actual_age, 
                   use = "complete.obs")

cat("\n=== HORVATH CLOCK PERFORMANCE ===\n")
cat("MAE:", round(mae, 2), "years\n")
cat("RMSE:", round(rmse, 2), "years\n")
cat("Correlation:", round(correlation, 3), "\n")

# ===== CREATE SCATTER PLOT =====
plot(horvath_results$actual_age, 
     horvath_results$Horvath,
     pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.6),
     cex = 1.2,
     xlab = "Chronological Age (years)",
     ylab = "Horvath Predicted Age (years)",
     main = "Horvath Clock Performance")

abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
abline(lm(horvath_results$Horvath ~ horvath_results$actual_age), 
       col = "darkblue", lwd = 2)
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("RMSE = ", round(rmse, 2), " years")
       ),
       bty = "n",
       cex = 1)

legend("bottomright",
       legend = c("Perfect prediction", "Actual fit"),
       col = c("red", "darkblue"),
       lty = c(2, 1),
       lwd = 2,
       bty = "n",
       cex = 0.9)




# ===== RUN ALL 4 CLOCKS (WITHOUT AGE) =====
cat("\nRunning all 4 methylation clocks...\n")

clock_results <- DNAmAge(Train_data, 
                         clocks = c("Horvath", "skinHorvath", "Hannum", "EN"))

cat("✓ Done!\n")
cat("Results dimensions:", dim(clock_results), "\n")
head(clock_results)

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(clock_results)) {
  actual_ages <- actual_ages[1:nrow(clock_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
clock_results$actual_age <- actual_ages

# Check
head(clock_results)

# ===== CHECK FOR MISSING VALUES =====
cat("\n=== Checking for missing values ===\n")
cat("NAs in Horvath:", sum(is.na(clock_results$Horvath)), "\n")
cat("NAs in SkinHorvath:", sum(is.na(clock_results$skinHorvath)), "\n")
cat("NAs in Hannum:", sum(is.na(clock_results$Hannum)), "\n")
cat("NAs in EN:", sum(is.na(clock_results$EN)), "\n")
cat("NAs in actual_age:", sum(is.na(clock_results$actual_age)), "\n")

# ===== IDENTIFY WORKING CLOCKS =====
working_clocks <- c()
if(sum(!is.na(clock_results$Horvath)) > 0) working_clocks <- c(working_clocks, "Horvath")
if(sum(!is.na(clock_results$skinHorvath)) > 0) working_clocks <- c(working_clocks, "skinHorvath")
if(sum(!is.na(clock_results$Hannum)) > 0) working_clocks <- c(working_clocks, "Hannum")
if(sum(!is.na(clock_results$EN)) > 0) working_clocks <- c(working_clocks, "EN")

cat("\n=== Working clocks ===\n")
print(working_clocks)

# ===== CALCULATE METRICS (ONLY FOR WORKING CLOCKS) =====
calculate_metrics <- function(predicted, actual, clock_name) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    cat("WARNING:", clock_name, "has no complete predictions\n")
    return(list(MAE = NA, RMSE = NA, Correlation = NA, N = 0))
  }
  
  predicted_clean <- predicted[complete_idx]
  actual_clean <- actual[complete_idx]
  error <- predicted_clean - actual_clean
  mae <- mean(abs(error))
  rmse <- sqrt(mean(error^2))
  correlation <- cor(predicted_clean, actual_clean)
  
  return(list(
    MAE = round(mae, 2),
    RMSE = round(rmse, 2),
    Correlation = round(correlation, 3),
    N = length(predicted_clean)
  ))
}

# Calculate for each clock
metrics_horvath <- calculate_metrics(clock_results$Horvath, clock_results$actual_age, "Horvath")
metrics_skin <- calculate_metrics(clock_results$skinHorvath, clock_results$actual_age, "skinHorvath")
metrics_hannum <- calculate_metrics(clock_results$Hannum, clock_results$actual_age, "Hannum")
metrics_en <- calculate_metrics(clock_results$EN, clock_results$actual_age, "EN")

# Create metrics table
metrics_table <- data.frame(
  Clock = c("Horvath", "skinHorvath", "Hannum", "EN"),
  N_samples = c(metrics_horvath$N, metrics_skin$N, 
                metrics_hannum$N, metrics_en$N),
  MAE = c(metrics_horvath$MAE, metrics_skin$MAE, 
          metrics_hannum$MAE, metrics_en$MAE),
  RMSE = c(metrics_horvath$RMSE, metrics_skin$RMSE, 
           metrics_hannum$RMSE, metrics_en$RMSE),
  Correlation = c(metrics_horvath$Correlation, metrics_skin$Correlation,
                  metrics_hannum$Correlation, metrics_en$Correlation)
)

# Print table
cat("\n=== PERFORMANCE METRICS ===\n")
print(metrics_table)

# Identify failed clocks
failed_clocks <- metrics_table$Clock[metrics_table$N_samples == 0]
if(length(failed_clocks) > 0) {
  cat("\n⚠️  Failed clocks (all NA):", paste(failed_clocks, collapse = ", "), "\n")
  cat("This is expected for semen data - blood-specific clocks may not work\n")
}

# ===== CREATE 4-PANEL FIGURE =====
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Function to create scatter plot
plot_clock <- function(actual, predicted, clock_name, mae_val, rmse_val, cor_val, n_val) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    # If no data, show message
    plot.new()
    text(0.5, 0.5, paste(clock_name, "\nNo predictions available\n(All NA)"), 
         cex = 1.2, col = "red")
    return()
  }
  actual_clean <- actual[complete_idx]
  predicted_clean <- predicted[complete_idx]
  
  plot(actual_clean, predicted_clean,
       pch = 19,
       col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.2,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = clock_name,
       xlim = range(c(actual_clean, predicted_clean)),
       ylim = range(c(actual_clean, predicted_clean)))
  
  # Perfect prediction line
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  
  # Regression line
  abline(lm(predicted_clean ~ actual_clean), col = "darkblue", lwd = 2)
  
  # Metrics
  if(!is.na(cor_val)) {
    legend("topleft",
           legend = c(
             paste0("r = ", cor_val),
             paste0("MAE = ", mae_val, " yrs"),
             paste0("RMSE = ", rmse_val, " yrs"),
             paste0("n = ", n_val)
           ),
           bty = "n",
           cex = 0.75)
  }
}

# Create all 4 panels
plot_clock(clock_results$actual_age, clock_results$Horvath,
           "Horvath Clock",
           metrics_table$MAE[1], metrics_table$RMSE[1], 
           metrics_table$Correlation[1], metrics_table$N_samples[1])

plot_clock(clock_results$actual_age, clock_results$skinHorvath,
           "Skin+Blood Horvath Clock",
           metrics_table$MAE[2], metrics_table$RMSE[2], 
           metrics_table$Correlation[2], metrics_table$N_samples[2])

plot_clock(clock_results$actual_age, clock_results$Hannum,
           "Hannum Clock",
           metrics_table$MAE[3], metrics_table$RMSE[3], 
           metrics_table$Correlation[3], metrics_table$N_samples[3])

plot_clock(clock_results$actual_age, clock_results$EN, 
           "Elastic Net Clock",
           metrics_table$MAE[4], metrics_table$RMSE[4], 
           metrics_table$Correlation[4], metrics_table$N_samples[4])
par(mfrow = c(1, 1))


#Hannum did not give a plot as there was not enough data available

Train50k <-- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_top_50K.RDS")
Train <--readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")

#cg16867657 vs clock results

dim(Train)

which(rownames(Train)== "cg16867657")

plot(clock_results$actual_age, Train[99804,])

#cg17147820

which(rownames(Train) == "cg17147820")
plot(clock_results$actual_age, Train[394604,])

#cg10528482
which(rownames(Train) == "cg10528482")
plot(clock_results$actual_age, Train[39955,])


# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Example: Your data frame with chromosomes and positions
my_data <- data.frame(
  chr = c("chr1", "chr2", "chr3"),
  pos = c(10497, 10525, 10542)
)

# Get annotation from 450k array
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Create a key for matching
anno$key <- paste0(anno$chr, ":", anno$pos)
my_data$key <- paste0(my_data$chr, ":", my_data$pos)

# Match and extract CpG IDs
my_data$cpg_id <- anno$Name[match(my_data$key, anno$key)]

# View results
print(my_data)



#VISAGE LOCI
# ===== INSTALL AND LOAD ANNOTATION PACKAGE =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load library
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# ===== CREATE VISAGE DATA FRAME FROM YOUR IMAGE =====
# This is YOUR data (from the table in your image)
visage_data <- data.frame(
  Gene = c("SYT7", "TUBB3", "SH2B2", "ARHGEF17", "EXOC3", "GALR2", "PPP2R2C"),
  CpG_ID = c("cg17147820", "cg18701351", "cg00018181", "cg09855959", 
             "cg10528482", "cg07909178", "cg02766173"),
  chr = c("chr11", "chr16", "chr7", "chr11", "chr5", "chr17", "chr4"),
  pos = c(61554783, 89921897, 102288444, 73311506, 525656, 76077795, 6473455)
)

cat("=== VISAGE Enhanced Tool Markers ===\n")
print(visage_data)


# ===== GET ANNOTATION FROM 450k ARRAY =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation database loaded\n")
cat("Total CpGs in 450k array:", nrow(anno), "\n\n")

# ===== CREATE MATCHING KEYS =====
# Combine chromosome + position into unique key
anno$key <- paste0(anno$chr, ":", anno$pos)
visage_data$key <- paste0(visage_data$chr, ":", visage_data$pos)


# ===== VERIFY CpG IDs MATCH POSITIONS =====
# Match positions to get CpG IDs from annotation
visage_data$cpg_verified <- anno$Name[match(visage_data$key, anno$key)]

# Check if CpG IDs match
visage_data$match <- visage_data$CpG_ID == visage_data$cpg_verified
cat("=== Verification Results ===\n")
print(visage_data[, c("Gene", "CpG_ID", "cpg_verified", "match")])

if(all(visage_data$match, na.rm = TRUE)) {
  cat("\n✓ All CpG IDs verified correctly!\n\n")
} else {
  cat("\n⚠️  Some mismatches detected\n")
  cat("Check the rows where match = FALSE\n\n")
}


# ===== GET ADDITIONAL ANNOTATION INFO =====
# Add gene names and other info from annotation
visage_data$anno_gene <- anno$UCSC_RefGene_Name[match(visage_data$key, anno$key)]
visage_data$cpg_island <- anno$Relation_to_Island[match(visage_data$key, anno$key)]


cat("=== Additional Annotation Info ===\n")
print(visage_data[, c("Gene", "CpG_ID", "anno_gene", "cpg_island")])

# ===== NOW LOAD YOUR METHYLATION DATA =====
library(readr)

Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age_data <- read.csv("Training data set.csv")
train_age <- train_age_data$Age

# ===== FIX: MATCH LENGTHS =====
cat("Data samples:", ncol(Train), "\n")
cat("Ages:", length(train_age), "\n")

# Trim ages to match data
if(length(train_age) > ncol(Train)) {
  train_age <- train_age[1:ncol(Train)]
  cat("Trimmed ages to:", length(train_age), "\n")
}

cat("\n=== Checking VISAGE Markers in Your Data ===\n")

# Check which VISAGE markers are in your data
visage_data$in_your_data <- visage_data$CpG_ID %in% rownames(Train)

print(visage_data[, c("Gene", "CpG_ID", "in_your_data")])

cat("\nAvailable in your data:", sum(visage_data$in_your_data), "/ 7\n\n")


  # ===== BUILD LINEAR MODEL WITH AVAILABLE MARKERS =====
available_markers <- visage_data$CpG_ID[visage_data$in_your_data]
available_genes <- visage_data$Gene[visage_data$in_your_data]

if(length(available_markers) >= 3) {
  
  cat("✓ Building model with", length(available_markers), "markers:\n")
  cat("  ", paste(available_genes, collapse = ", "), "\n\n")
  
  # Extract methylation data for these markers
  train_visage <- t(Train[available_markers, ])
  
  # Create modeling data frame
  train_df <- data.frame(
    age = train_age,
    train_visage
  )
  # Use gene names as column names
  colnames(train_df)[-1] <- available_genes
  
  
  # ===== TRAIN LINEAR MODEL =====
  visage_model <- lm(age ~ ., data = train_df)
  
  cat("=== VISAGE Linear Model Summary ===\n")
  print(summary(visage_model))
  # ===== MODEL EQUATION =====
  cat("\n=== Prediction Equation ===\n")
  coefs <- coef(visage_model)
  cat("Age = ", round(coefs[1], 2), "\n", sep = "")
  for(i in 2:length(coefs)) {
    cat("    ", ifelse(coefs[i] >= 0, "+", ""), " ", 
        round(coefs[i], 2), " × ", names(coefs)[i], "\n", sep = "")
  }
  
  # ===== PREDICTIONS =====
  predictions <- predict(visage_model)
  errors <- predictions - train_age
  
  
  # ===== PERFORMANCE METRICS =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, train_age)
  r2 <- summary(visage_model)$r.squared
  adj_r2 <- summary(visage_model)$adj.r.squared
  
  cat("\n========================================\n")
  cat("VISAGE MODEL PERFORMANCE\n")
  cat("========================================\n")
  cat("Markers:", length(available_markers), "/ 7\n")
  cat("Samples:", length(train_age), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("Adjusted R²:", round(adj_r2, 3), "\n")
  
  
  # ===== SCATTER PLOT =====
  plot(train_age, predictions,
       pch = 19, col = rgb(0.2, 0.6, 0.4, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = paste0("VISAGE Enhanced Tool (", length(available_markers), " CpGs)"))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ train_age), col = "darkgreen", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3))
         ),
         bty = "n", cex = 1.1)
  
  # ===== SAVE EVERYTHING =====
  results_df <- data.frame(
    sample_id = rownames(train_visage),
    actual_age = train_age,
    predicted_age = predictions,
    error = errors,
    abs_error = abs(errors)
  )
  
  write.csv(results_df, "VISAGE_results.csv", row.names = FALSE)
  write.csv(visage_data, "VISAGE_marker_info.csv", row.names = FALSE)
  saveRDS(visage_model, "VISAGE_model.RDS")
  
  cat("\n✓ All files saved!\n")
} else {
  cat("❌ Not enough markers available (need at least 3)\n")
  cat("Available:", sum(visage_data$in_your_data), "/ 7\n")
}

#Colouring by study#
# After you've made predictions with VISAGE model

# Get study labels
train_age_data <- read_csv("Training data set.csv")
study_labels <- train_age_data$GSE[1:length(predictions)]

# Define colors
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

colors_by_study <- study_colors[study_labels]
# Plot
plot(train_age, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "VISAGE Predicted Age (years)",
     main = "VISAGE Model - By Study")

abline(0, 1, col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       title = "Study")
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1.1)
# ===== FINAL SUMMARY TABLE =====
cat("\n========================================\n")
cat("SUMMARY TABLE\n")
cat("========================================\n")

summary_table <- data.frame(
  Gene = visage_data$Gene,
  CpG_ID = visage_data$CpG_ID,
  Chr = visage_data$chr,
  Position = visage_data$pos,
  In_Data = visage_data$in_your_data,
  Verified = ifelse(is.na(visage_data$match), "N/A", visage_data$match)
)

print(summary_table)



#JENKINS

library(readr)

# ===== INSTALL AND LOAD ANNOTATION =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ===== CREATE JENKINS DATA  =====
# Combining data 
jenkins_data <- data.frame(
  Gene = c(
    # Image 1
    "ADAMTS8", "ARC", "ARGHGEF10", "BCL11A", "C1ORF122", "C7ORF50",
    "CCDC144NL", "CLIC1", "DMPK", "FAM86C1", "FAM86JP", "FOXK1",
    "FSCN", "GAPDH", "GET4", "GNB2", "GPANK1", "GPR45", "KCNQ1",
    # Image 2
    "LDLRAD4", "LMO3", "LOC100133461", "MIR22HG", "MTMR8", "N10",
    "N12", "N22", "N23", "N24", "N27", "N30", "N8", "N9",
    "NCOR2", "NONE", "NSG1", "PAX2", "PITX1", "PRSS22",
    "PTPRN2.3", "PTPRN2.4", "PURA", "PYY2", "SECTM1", "SEMA6B",
    "SEZ6", "SLC22A18AS", "SOHLH1", "THBS3", "TNXB"
  ),
  chr = c(
    # Image 1
    "chr11", "chr8", "chr8", "chr2", "chr1", "chr7",
    "chr17", "chr6", "chr19", "chr11", "chr3", "chr7",
    "chr7", "chr12", "chr7", "chr7", "chr6", "chr2", "chr11",
    # Image 2
    "chr18", "chr12", "chr4", "chr17", "chrX", "chr1",
    "chr5", "chr19", "chr14", "chr6", "chr6", "chr15", "chr11", "chr7",
    "chr12", "chr10", "chr4", "chr10", "chr5", "chr16",
    "chr7", "chr7", "chr5", "chr17", "chr17", "chr19",
    "chr17", "chr11", "chr9", "chr1", "chr6"
  ),
  start = c(
    # Image 1
    130299298, 143694010, 1877888, 60680616, 38272200, 1083209,
    20798895, 31698492, 46282571, 71498202, 125634060, 4722778,
    5635134, 6641602, 914964, 100274361, 31630819, 105857809, 2554562,
    # Image 2
    13611370, 16760040, 3680721, 1617363, 63614857, 28423399,
    3593413, 4579481, 106004434, 170449417, 30432200, 27959473, 69260136, 35300077,
    124990897, 17347047, 4386726, 102509693, 134365728, 2908157,
    157523356, 158109339, 139492535, 26553567, 80278592, 4555999,
    27330794, 2909690, 138590204, 155176868, 32064146
  ),
  stop = c(
    # Image 1
    130299948, 143694548, 1878324, 60680762, 38273057, 1084163,
    20799770, 31699299, 46283081, 71499118, 125634453, 4723928,
    5635954, 6642355, 915832, 100275305, 31632542, 105859084, 2555577,
    # Image 2
    13611825, 16761003, 3681760, 1618296, 63615496, 28424202,
    3594276, 4580471, 106004608, 170450804, 30433944, 27960032, 69261045, 35301070,
    124991140, 17347392, 4387698, 102510569, 134366535, 2908935,
    157524159, 158110153, 139493491, 26554908, 80280331, 4556983,
    27332647, 2909716, 138590996, 155177784, 32065891
  )
)

cat("=== Jenkins et al. 50 Genomic Regions ===\n")
cat("Total regions:", nrow(jenkins_data), "\n\n")
print(head(jenkins_data, 50))

# ===== LOAD ANNOTATION DATABASE =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation loaded:", nrow(anno), "CpGs\n\n")

# ===== FIND CpGs WITHIN EACH JENKINS REGION =====
# Jenkins regions are genomic windows, not single CpGs
# Need to find ALL CpGs within each region

jenkins_cpgs <- list()

cat("=== Finding CpGs in Jenkins regions ===\n")

for(i in 1:nrow(jenkins_data)) {
  region_chr <- jenkins_data$chr[i]
  region_start <- jenkins_data$start[i]
  region_stop <- jenkins_data$stop[i]
  region_name <- jenkins_data$Gene[i]
  
  # Find CpGs in this region
  cpgs_in_region <- anno$Name[
    anno$chr == region_chr & 
      anno$pos >= region_start & 
      anno$pos <= region_stop
  ]
  
  jenkins_cpgs[[i]] <- cpgs_in_region
  
  cat(i, ". ", region_name, " (", region_chr, ":", region_start, "-", region_stop, "): ",
      length(cpgs_in_region), " CpGs\n", sep = "")
}
# Total CpGs across all regions
all_jenkins_cpgs <- unique(unlist(jenkins_cpgs))
cat("\nTotal unique CpGs across 50 regions:", length(all_jenkins_cpgs), "\n\n")

# ===== CHECK AVAILABILITY IN YOUR DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age <- read.csv("Training data set.csv")$Age
train_age <- train_age[1:ncol(Train)]  # Fix length

available_jenkins <- all_jenkins_cpgs %in% rownames(Train)
cat("=== Jenkins CpGs in Your Data ===\n")
cat("Total Jenkins CpGs:", length(all_jenkins_cpgs), "\n")
cat("Available in your data:", sum(available_jenkins), "\n")
cat("Percentage:", round(sum(available_jenkins)/length(all_jenkins_cpgs)*100, 1), "%\n\n")

# ===== BUILD MODEL WITH AVAILABLE JENKINS CpGs =====
if(sum(available_jenkins) >= 50) {
  
  jenkins_use <- all_jenkins_cpgs[available_jenkins]
  
  cat("✓ Building Jenkins model with", length(jenkins_use), "CpGs\n\n")
  
  # Extract data
  train_jenkins <- t(Train[jenkins_use, ])
  
  # Create data frame
  train_df <- data.frame(
    age = train_age,
    train_jenkins
  )
  
  # ===== TRAIN MODEL =====
  # Use elastic net due to large number of features
  library(glmnet)
  
  X_train <- as.matrix(train_df[, -1])
  y_train <- train_df$age
  
  complete_idx <- complete.cases(X_train, y_train)
  X_train <- X_train[complete_idx, ]
  y_train <- y_train[complete_idx]
  
  
  # Cross-validated elastic net
  jenkins_model <- cv.glmnet(X_train, y_train, 
                             alpha = 0.5,  # Elastic net
                             nfolds = 10)
  
  cat("=== Jenkins Model Trained ===\n")
  cat("Lambda min:", jenkins_model$lambda.min, "\n")
  cat("Lambda 1se:", jenkins_model$lambda.1se, "\n\n")
  

  # ===== PREDICTIONS =====
  predictions <- predict(jenkins_model, newx = X_train, s = "lambda.min")[,1]
  errors <- predictions - y_train
  
  
  # ===== PERFORMANCE =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, y_train)
  r2 <- summary(lm(predictions ~ y_train))$r.squared
  
  cat("=== Jenkins Model Performance ===\n")
  cat("CpGs used:", length(jenkins_use), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("\nPublished Jenkins performance: MAE = 2.04-2.37 years\n")
  
  
  # ===== PLOT =====
  plot(y_train, predictions,
       pch = 19, col = rgb(0.4, 0.2, 0.8, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Jenkins Predicted Age (years)",
       main = paste0("Jenkins Germ Line Calculator (", length(jenkins_use), " CpGs)"))
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ y_train), col = "purple", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3)),
           paste0("n = ", length(jenkins_use), " CpGs")
         ),
         bty = "n", cex = 1)
  
  # ===== SAVE =====
  results <- data.frame(
    actual_age = y_train,
    predicted_age = predictions,
    error = errors
  )
  
  write.csv(results, "Jenkins_predictions.csv", row.names = FALSE)
  saveRDS(jenkins_model, "Jenkins_model.RDS")
  # Save CpG list
  cpg_list <- data.frame(CpG_ID = jenkins_use)
  write.csv(cpg_list, "Jenkins_CpGs_used.csv", row.names = FALSE)
  
  cat("\n✓ Results saved!\n")
  
} else {
  cat("❌ Not enough Jenkins CpGs available\n")
  cat("Need ~200+ CpGs, have:", sum(available_jenkins), "\n")
  cat("Try loading full unfiltered data\n")
}


# ===== SUMMARY OF CpGs PER REGION =====
cat("\n=== CpGs per Jenkins Region ===\n")

region_summary <- data.frame(
  Region = jenkins_data$Gene,
  Chr = jenkins_data$chr,
  N_CpGs_total = sapply(jenkins_cpgs, length),
  N_CpGs_available = sapply(jenkins_cpgs, function(cpgs) sum(cpgs %in% rownames(Train)))
)

print(head(region_summary, 20))

cat("\nRegions with most CpGs:\n")
top_regions <- region_summary[order(region_summary$N_CpGs_total, decreasing = TRUE), ][1:10, ]
print(top_regions)

write.csv(region_summary, "Jenkins_region_summary.csv", row.names = FALSE)

#Colour by datasets
#Tidy up cg plot
#Calc correlation 50k variably methylated lociwith age 
#plot top 6 most correlating/ anticorrelating- absolute vals

# JENKINS COLOURED BY STUDY #

# ===== CREATE STUDY LABELS =====
# Load study information from your CSV
train_age_data <- read_csv("Training data set.csv")
train_study <- train_age_data$GSE

# Match to the samples you're actually using (after removing NAs)
study_labels <- train_study[complete_idx]

# Check
cat("Study distribution:\n")
print(table(study_labels))

# ===== DEFINE COLORS FOR 3 STUDIES =====
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

# Create color vector
colors_by_study <- study_colors[study_labels]

# ===== PLOT COLORED BY STUDY =====
plot(y_train, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "Jenkins Predicted Age (years)",
     main = "Jenkins Model - Colored by Study")

# Perfect prediction line
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Regression line
abline(lm(predictions ~ y_train), col = "black", lwd = 2)

# Legend
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " yrs"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study")



#glmnet package 
install.packages("glmnet", repos = "https://cran.us.r-project.org")


library(methylclock)
Train_data <- readRDS("~/Documents/Forensic Project/preExistingClocks.Train.RDS")
head(Train_data)

# Run Horvath clock specifically
# If your data has actual chronological ages
# (e.g., in a row called "age" or as a separate vector)

# ===== LOAD DATA =====

cat("Data dimensions:", dim(Train_data), "\n")
cat("Samples (columns):", ncol(Train_data), "\n")


# ===== RUN HORVATH CLOCK (WITHOUT AGE) =====
cat("\nRunning Horvath clock...\n")

horvath_results <- DNAmAge(Train_data, clocks = "Horvath")

cat("✓ Done!\n")
cat("Results dimensions:", dim(horvath_results), "\n")

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(horvath_results)) {
  actual_ages <- actual_ages[1:nrow(horvath_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
horvath_results$actual_age <- actual_ages

# Check
head(horvath_results)

# ===== CALCULATE METRICS =====
horvath_results$error <- horvath_results$Horvath - horvath_results$actual_age
horvath_results$abs_error <- abs(horvath_results$error)

mae <- mean(horvath_results$abs_error, na.rm = TRUE)
rmse <- sqrt(mean(horvath_results$error^2, na.rm = TRUE))
correlation <- cor(horvath_results$Horvath, horvath_results$actual_age, 
                   use = "complete.obs")

cat("\n=== HORVATH CLOCK PERFORMANCE ===\n")
cat("MAE:", round(mae, 2), "years\n")
cat("RMSE:", round(rmse, 2), "years\n")
cat("Correlation:", round(correlation, 3), "\n")

# ===== CREATE SCATTER PLOT =====
plot(horvath_results$actual_age, 
     horvath_results$Horvath,
     pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.6),
     cex = 1.2,
     xlab = "Chronological Age (years)",
     ylab = "Horvath Predicted Age (years)",
     main = "Horvath Clock Performance")

abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
abline(lm(horvath_results$Horvath ~ horvath_results$actual_age), 
       col = "darkblue", lwd = 2)
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("RMSE = ", round(rmse, 2), " years")
       ),
       bty = "n",
       cex = 1)

legend("bottomright",
       legend = c("Perfect prediction", "Actual fit"),
       col = c("red", "darkblue"),
       lty = c(2, 1),
       lwd = 2,
       bty = "n",
       cex = 0.9)




# ===== RUN ALL 4 CLOCKS (WITHOUT AGE) =====
cat("\nRunning all 4 methylation clocks...\n")

clock_results <- DNAmAge(Train_data, 
                         clocks = c("Horvath", "skinHorvath", "Hannum", "EN"))

cat("✓ Done!\n")
cat("Results dimensions:", dim(clock_results), "\n")
head(clock_results)

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(clock_results)) {
  actual_ages <- actual_ages[1:nrow(clock_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
clock_results$actual_age <- actual_ages

# Check
head(clock_results)

# ===== CHECK FOR MISSING VALUES =====
cat("\n=== Checking for missing values ===\n")
cat("NAs in Horvath:", sum(is.na(clock_results$Horvath)), "\n")
cat("NAs in SkinHorvath:", sum(is.na(clock_results$skinHorvath)), "\n")
cat("NAs in Hannum:", sum(is.na(clock_results$Hannum)), "\n")
cat("NAs in EN:", sum(is.na(clock_results$EN)), "\n")
cat("NAs in actual_age:", sum(is.na(clock_results$actual_age)), "\n")

# ===== IDENTIFY WORKING CLOCKS =====
working_clocks <- c()
if(sum(!is.na(clock_results$Horvath)) > 0) working_clocks <- c(working_clocks, "Horvath")
if(sum(!is.na(clock_results$skinHorvath)) > 0) working_clocks <- c(working_clocks, "skinHorvath")
if(sum(!is.na(clock_results$Hannum)) > 0) working_clocks <- c(working_clocks, "Hannum")
if(sum(!is.na(clock_results$EN)) > 0) working_clocks <- c(working_clocks, "EN")

cat("\n=== Working clocks ===\n")
print(working_clocks)

# ===== CALCULATE METRICS (ONLY FOR WORKING CLOCKS) =====
calculate_metrics <- function(predicted, actual, clock_name) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    cat("WARNING:", clock_name, "has no complete predictions\n")
    return(list(MAE = NA, RMSE = NA, Correlation = NA, N = 0))
  }
  
  predicted_clean <- predicted[complete_idx]
  actual_clean <- actual[complete_idx]
  error <- predicted_clean - actual_clean
  mae <- mean(abs(error))
  rmse <- sqrt(mean(error^2))
  correlation <- cor(predicted_clean, actual_clean)
  
  return(list(
    MAE = round(mae, 2),
    RMSE = round(rmse, 2),
    Correlation = round(correlation, 3),
    N = length(predicted_clean)
  ))
}

# Calculate for each clock
metrics_horvath <- calculate_metrics(clock_results$Horvath, clock_results$actual_age, "Horvath")
metrics_skin <- calculate_metrics(clock_results$skinHorvath, clock_results$actual_age, "skinHorvath")
metrics_hannum <- calculate_metrics(clock_results$Hannum, clock_results$actual_age, "Hannum")
metrics_en <- calculate_metrics(clock_results$EN, clock_results$actual_age, "EN")

# Create metrics table
metrics_table <- data.frame(
  Clock = c("Horvath", "skinHorvath", "Hannum", "EN"),
  N_samples = c(metrics_horvath$N, metrics_skin$N, 
                metrics_hannum$N, metrics_en$N),
  MAE = c(metrics_horvath$MAE, metrics_skin$MAE, 
          metrics_hannum$MAE, metrics_en$MAE),
  RMSE = c(metrics_horvath$RMSE, metrics_skin$RMSE, 
           metrics_hannum$RMSE, metrics_en$RMSE),
  Correlation = c(metrics_horvath$Correlation, metrics_skin$Correlation,
                  metrics_hannum$Correlation, metrics_en$Correlation)
)

# Print table
cat("\n=== PERFORMANCE METRICS ===\n")
print(metrics_table)

# Identify failed clocks
failed_clocks <- metrics_table$Clock[metrics_table$N_samples == 0]
if(length(failed_clocks) > 0) {
  cat("\n⚠️  Failed clocks (all NA):", paste(failed_clocks, collapse = ", "), "\n")
  cat("This is expected for semen data - blood-specific clocks may not work\n")
}

# ===== CREATE 4-PANEL FIGURE =====
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Function to create scatter plot
plot_clock <- function(actual, predicted, clock_name, mae_val, rmse_val, cor_val, n_val) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    # If no data, show message
    plot.new()
    text(0.5, 0.5, paste(clock_name, "\nNo predictions available\n(All NA)"), 
         cex = 1.2, col = "red")
    return()
  }
  actual_clean <- actual[complete_idx]
  predicted_clean <- predicted[complete_idx]
  
  plot(actual_clean, predicted_clean,
       pch = 19,
       col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.2,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = clock_name,
       xlim = range(c(actual_clean, predicted_clean)),
       ylim = range(c(actual_clean, predicted_clean)))
  
  # Perfect prediction line
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  
  # Regression line
  abline(lm(predicted_clean ~ actual_clean), col = "darkblue", lwd = 2)
  
  # Metrics
  if(!is.na(cor_val)) {
    legend("topleft",
           legend = c(
             paste0("r = ", cor_val),
             paste0("MAE = ", mae_val, " yrs"),
             paste0("RMSE = ", rmse_val, " yrs"),
             paste0("n = ", n_val)
           ),
           bty = "n",
           cex = 0.75)
  }
}

# Create all 4 panels
plot_clock(clock_results$actual_age, clock_results$Horvath,
           "Horvath Clock",
           metrics_table$MAE[1], metrics_table$RMSE[1], 
           metrics_table$Correlation[1], metrics_table$N_samples[1])

plot_clock(clock_results$actual_age, clock_results$skinHorvath,
           "Skin+Blood Horvath Clock",
           metrics_table$MAE[2], metrics_table$RMSE[2], 
           metrics_table$Correlation[2], metrics_table$N_samples[2])

plot_clock(clock_results$actual_age, clock_results$Hannum,
           "Hannum Clock",
           metrics_table$MAE[3], metrics_table$RMSE[3], 
           metrics_table$Correlation[3], metrics_table$N_samples[3])

plot_clock(clock_results$actual_age, clock_results$EN, 
           "Elastic Net Clock",
           metrics_table$MAE[4], metrics_table$RMSE[4], 
           metrics_table$Correlation[4], metrics_table$N_samples[4])
par(mfrow = c(1, 1))


#Hannum did not give a plot as there was not enough data available

Train50k <-- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_top_50K.RDS")
Train <--readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")

#cg16867657 vs clock results

dim(Train)

which(rownames(Train)== "cg16867657")

plot(clock_results$actual_age, Train[99804,])

#cg17147820

which(rownames(Train) == "cg17147820")
plot(clock_results$actual_age, Train[394604,])

#cg10528482
which(rownames(Train) == "cg10528482")
plot(clock_results$actual_age, Train[39955,])


# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Example: Your data frame with chromosomes and positions
my_data <- data.frame(
  chr = c("chr1", "chr2", "chr3"),
  pos = c(10497, 10525, 10542)
)

# Get annotation from 450k array
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Create a key for matching
anno$key <- paste0(anno$chr, ":", anno$pos)
my_data$key <- paste0(my_data$chr, ":", my_data$pos)

# Match and extract CpG IDs
my_data$cpg_id <- anno$Name[match(my_data$key, anno$key)]

# View results
print(my_data)



#VISAGE LOCI
# ===== INSTALL AND LOAD ANNOTATION PACKAGE =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load library
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# ===== CREATE VISAGE DATA FRAME FROM YOUR IMAGE =====
# This is YOUR data (from the table in your image)
visage_data <- data.frame(
  Gene = c("SYT7", "TUBB3", "SH2B2", "ARHGEF17", "EXOC3", "GALR2", "PPP2R2C"),
  CpG_ID = c("cg17147820", "cg18701351", "cg00018181", "cg09855959", 
             "cg10528482", "cg07909178", "cg02766173"),
  chr = c("chr11", "chr16", "chr7", "chr11", "chr5", "chr17", "chr4"),
  pos = c(61554783, 89921897, 102288444, 73311506, 525656, 76077795, 6473455)
)

cat("=== VISAGE Enhanced Tool Markers ===\n")
print(visage_data)


# ===== GET ANNOTATION FROM 450k ARRAY =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation database loaded\n")
cat("Total CpGs in 450k array:", nrow(anno), "\n\n")

# ===== CREATE MATCHING KEYS =====
# Combine chromosome + position into unique key
anno$key <- paste0(anno$chr, ":", anno$pos)
visage_data$key <- paste0(visage_data$chr, ":", visage_data$pos)


# ===== VERIFY CpG IDs MATCH POSITIONS =====
# Match positions to get CpG IDs from annotation
visage_data$cpg_verified <- anno$Name[match(visage_data$key, anno$key)]

# Check if CpG IDs match
visage_data$match <- visage_data$CpG_ID == visage_data$cpg_verified
cat("=== Verification Results ===\n")
print(visage_data[, c("Gene", "CpG_ID", "cpg_verified", "match")])

if(all(visage_data$match, na.rm = TRUE)) {
  cat("\n✓ All CpG IDs verified correctly!\n\n")
} else {
  cat("\n⚠️  Some mismatches detected\n")
  cat("Check the rows where match = FALSE\n\n")
}


# ===== GET ADDITIONAL ANNOTATION INFO =====
# Add gene names and other info from annotation
visage_data$anno_gene <- anno$UCSC_RefGene_Name[match(visage_data$key, anno$key)]
visage_data$cpg_island <- anno$Relation_to_Island[match(visage_data$key, anno$key)]


cat("=== Additional Annotation Info ===\n")
print(visage_data[, c("Gene", "CpG_ID", "anno_gene", "cpg_island")])

# ===== NOW LOAD YOUR METHYLATION DATA =====
library(readr)

Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age_data <- read.csv("Training data set.csv")
train_age <- train_age_data$Age

# ===== FIX: MATCH LENGTHS =====
cat("Data samples:", ncol(Train), "\n")
cat("Ages:", length(train_age), "\n")

# Trim ages to match data
if(length(train_age) > ncol(Train)) {
  train_age <- train_age[1:ncol(Train)]
  cat("Trimmed ages to:", length(train_age), "\n")
}

cat("\n=== Checking VISAGE Markers in Your Data ===\n")

# Check which VISAGE markers are in your data
visage_data$in_your_data <- visage_data$CpG_ID %in% rownames(Train)

print(visage_data[, c("Gene", "CpG_ID", "in_your_data")])

cat("\nAvailable in your data:", sum(visage_data$in_your_data), "/ 7\n\n")


  # ===== BUILD LINEAR MODEL WITH AVAILABLE MARKERS =====
available_markers <- visage_data$CpG_ID[visage_data$in_your_data]
available_genes <- visage_data$Gene[visage_data$in_your_data]

if(length(available_markers) >= 3) {
  
  cat("✓ Building model with", length(available_markers), "markers:\n")
  cat("  ", paste(available_genes, collapse = ", "), "\n\n")
  
  # Extract methylation data for these markers
  train_visage <- t(Train[available_markers, ])
  
  # Create modeling data frame
  train_df <- data.frame(
    age = train_age,
    train_visage
  )
  # Use gene names as column names
  colnames(train_df)[-1] <- available_genes
  
  
  # ===== TRAIN LINEAR MODEL =====
  visage_model <- lm(age ~ ., data = train_df)
  
  cat("=== VISAGE Linear Model Summary ===\n")
  print(summary(visage_model))
  # ===== MODEL EQUATION =====
  cat("\n=== Prediction Equation ===\n")
  coefs <- coef(visage_model)
  cat("Age = ", round(coefs[1], 2), "\n", sep = "")
  for(i in 2:length(coefs)) {
    cat("    ", ifelse(coefs[i] >= 0, "+", ""), " ", 
        round(coefs[i], 2), " × ", names(coefs)[i], "\n", sep = "")
  }
  
  # ===== PREDICTIONS =====
  predictions <- predict(visage_model)
  errors <- predictions - train_age
  
  
  # ===== PERFORMANCE METRICS =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, train_age)
  r2 <- summary(visage_model)$r.squared
  adj_r2 <- summary(visage_model)$adj.r.squared
  
  cat("\n========================================\n")
  cat("VISAGE MODEL PERFORMANCE\n")
  cat("========================================\n")
  cat("Markers:", length(available_markers), "/ 7\n")
  cat("Samples:", length(train_age), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("Adjusted R²:", round(adj_r2, 3), "\n")
  
  
  # ===== SCATTER PLOT =====
  plot(train_age, predictions,
       pch = 19, col = rgb(0.2, 0.6, 0.4, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = paste0("VISAGE Enhanced Tool (", length(available_markers), " CpGs)"))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ train_age), col = "darkgreen", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3))
         ),
         bty = "n", cex = 1.1)
  
  # ===== SAVE EVERYTHING =====
  results_df <- data.frame(
    sample_id = rownames(train_visage),
    actual_age = train_age,
    predicted_age = predictions,
    error = errors,
    abs_error = abs(errors)
  )
  
  write.csv(results_df, "VISAGE_results.csv", row.names = FALSE)
  write.csv(visage_data, "VISAGE_marker_info.csv", row.names = FALSE)
  saveRDS(visage_model, "VISAGE_model.RDS")
  
  cat("\n✓ All files saved!\n")
} else {
  cat("❌ Not enough markers available (need at least 3)\n")
  cat("Available:", sum(visage_data$in_your_data), "/ 7\n")
}

#Colouring by study#
# After you've made predictions with VISAGE model

# Get study labels
train_age_data <- read_csv("Training data set.csv")
study_labels <- train_age_data$GSE[1:length(predictions)]

# Define colors
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

colors_by_study <- study_colors[study_labels]
# Plot
plot(train_age, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "VISAGE Predicted Age (years)",
     main = "VISAGE Model - By Study")

abline(0, 1, col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       title = "Study")
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1.1)
# ===== FINAL SUMMARY TABLE =====
cat("\n========================================\n")
cat("SUMMARY TABLE\n")
cat("========================================\n")

summary_table <- data.frame(
  Gene = visage_data$Gene,
  CpG_ID = visage_data$CpG_ID,
  Chr = visage_data$chr,
  Position = visage_data$pos,
  In_Data = visage_data$in_your_data,
  Verified = ifelse(is.na(visage_data$match), "N/A", visage_data$match)
)

print(summary_table)



#JENKINS

library(readr)

# ===== INSTALL AND LOAD ANNOTATION =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ===== CREATE JENKINS DATA  =====
# Combining data 
jenkins_data <- data.frame(
  Gene = c(
    # Image 1
    "ADAMTS8", "ARC", "ARGHGEF10", "BCL11A", "C1ORF122", "C7ORF50",
    "CCDC144NL", "CLIC1", "DMPK", "FAM86C1", "FAM86JP", "FOXK1",
    "FSCN", "GAPDH", "GET4", "GNB2", "GPANK1", "GPR45", "KCNQ1",
    # Image 2
    "LDLRAD4", "LMO3", "LOC100133461", "MIR22HG", "MTMR8", "N10",
    "N12", "N22", "N23", "N24", "N27", "N30", "N8", "N9",
    "NCOR2", "NONE", "NSG1", "PAX2", "PITX1", "PRSS22",
    "PTPRN2.3", "PTPRN2.4", "PURA", "PYY2", "SECTM1", "SEMA6B",
    "SEZ6", "SLC22A18AS", "SOHLH1", "THBS3", "TNXB"
  ),
  chr = c(
    # Image 1
    "chr11", "chr8", "chr8", "chr2", "chr1", "chr7",
    "chr17", "chr6", "chr19", "chr11", "chr3", "chr7",
    "chr7", "chr12", "chr7", "chr7", "chr6", "chr2", "chr11",
    # Image 2
    "chr18", "chr12", "chr4", "chr17", "chrX", "chr1",
    "chr5", "chr19", "chr14", "chr6", "chr6", "chr15", "chr11", "chr7",
    "chr12", "chr10", "chr4", "chr10", "chr5", "chr16",
    "chr7", "chr7", "chr5", "chr17", "chr17", "chr19",
    "chr17", "chr11", "chr9", "chr1", "chr6"
  ),
  start = c(
    # Image 1
    130299298, 143694010, 1877888, 60680616, 38272200, 1083209,
    20798895, 31698492, 46282571, 71498202, 125634060, 4722778,
    5635134, 6641602, 914964, 100274361, 31630819, 105857809, 2554562,
    # Image 2
    13611370, 16760040, 3680721, 1617363, 63614857, 28423399,
    3593413, 4579481, 106004434, 170449417, 30432200, 27959473, 69260136, 35300077,
    124990897, 17347047, 4386726, 102509693, 134365728, 2908157,
    157523356, 158109339, 139492535, 26553567, 80278592, 4555999,
    27330794, 2909690, 138590204, 155176868, 32064146
  ),
  stop = c(
    # Image 1
    130299948, 143694548, 1878324, 60680762, 38273057, 1084163,
    20799770, 31699299, 46283081, 71499118, 125634453, 4723928,
    5635954, 6642355, 915832, 100275305, 31632542, 105859084, 2555577,
    # Image 2
    13611825, 16761003, 3681760, 1618296, 63615496, 28424202,
    3594276, 4580471, 106004608, 170450804, 30433944, 27960032, 69261045, 35301070,
    124991140, 17347392, 4387698, 102510569, 134366535, 2908935,
    157524159, 158110153, 139493491, 26554908, 80280331, 4556983,
    27332647, 2909716, 138590996, 155177784, 32065891
  )
)

cat("=== Jenkins et al. 50 Genomic Regions ===\n")
cat("Total regions:", nrow(jenkins_data), "\n\n")
print(head(jenkins_data, 50))

# ===== LOAD ANNOTATION DATABASE =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation loaded:", nrow(anno), "CpGs\n\n")

# ===== FIND CpGs WITHIN EACH JENKINS REGION =====
# Jenkins regions are genomic windows, not single CpGs
# Need to find ALL CpGs within each region

jenkins_cpgs <- list()

cat("=== Finding CpGs in Jenkins regions ===\n")

for(i in 1:nrow(jenkins_data)) {
  region_chr <- jenkins_data$chr[i]
  region_start <- jenkins_data$start[i]
  region_stop <- jenkins_data$stop[i]
  region_name <- jenkins_data$Gene[i]
  
  # Find CpGs in this region
  cpgs_in_region <- anno$Name[
    anno$chr == region_chr & 
      anno$pos >= region_start & 
      anno$pos <= region_stop
  ]
  
  jenkins_cpgs[[i]] <- cpgs_in_region
  
  cat(i, ". ", region_name, " (", region_chr, ":", region_start, "-", region_stop, "): ",
      length(cpgs_in_region), " CpGs\n", sep = "")
}
# Total CpGs across all regions
all_jenkins_cpgs <- unique(unlist(jenkins_cpgs))
cat("\nTotal unique CpGs across 50 regions:", length(all_jenkins_cpgs), "\n\n")

# ===== CHECK AVAILABILITY IN YOUR DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age <- read.csv("Training data set.csv")$Age
train_age <- train_age[1:ncol(Train)]  # Fix length

available_jenkins <- all_jenkins_cpgs %in% rownames(Train)
cat("=== Jenkins CpGs in Your Data ===\n")
cat("Total Jenkins CpGs:", length(all_jenkins_cpgs), "\n")
cat("Available in your data:", sum(available_jenkins), "\n")
cat("Percentage:", round(sum(available_jenkins)/length(all_jenkins_cpgs)*100, 1), "%\n\n")

# ===== BUILD MODEL WITH AVAILABLE JENKINS CpGs =====
if(sum(available_jenkins) >= 50) {
  
  jenkins_use <- all_jenkins_cpgs[available_jenkins]
  
  cat("✓ Building Jenkins model with", length(jenkins_use), "CpGs\n\n")
  
  # Extract data
  train_jenkins <- t(Train[jenkins_use, ])
  
  # Create data frame
  train_df <- data.frame(
    age = train_age,
    train_jenkins
  )
  
  # ===== TRAIN MODEL =====
  # Use elastic net due to large number of features
  library(glmnet)
  
  X_train <- as.matrix(train_df[, -1])
  y_train <- train_df$age
  
  complete_idx <- complete.cases(X_train, y_train)
  X_train <- X_train[complete_idx, ]
  y_train <- y_train[complete_idx]
  
  
  # Cross-validated elastic net
  jenkins_model <- cv.glmnet(X_train, y_train, 
                             alpha = 0.5,  # Elastic net
                             nfolds = 10)
  
  cat("=== Jenkins Model Trained ===\n")
  cat("Lambda min:", jenkins_model$lambda.min, "\n")
  cat("Lambda 1se:", jenkins_model$lambda.1se, "\n\n")
  

  # ===== PREDICTIONS =====
  predictions <- predict(jenkins_model, newx = X_train, s = "lambda.min")[,1]
  errors <- predictions - y_train
  
  
  # ===== PERFORMANCE =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, y_train)
  r2 <- summary(lm(predictions ~ y_train))$r.squared
  
  cat("=== Jenkins Model Performance ===\n")
  cat("CpGs used:", length(jenkins_use), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("\nPublished Jenkins performance: MAE = 2.04-2.37 years\n")
  
  
  # ===== PLOT =====
  plot(y_train, predictions,
       pch = 19, col = rgb(0.4, 0.2, 0.8, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Jenkins Predicted Age (years)",
       main = paste0("Jenkins Germ Line Calculator (", length(jenkins_use), " CpGs)"))
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ y_train), col = "purple", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3)),
           paste0("n = ", length(jenkins_use), " CpGs")
         ),
         bty = "n", cex = 1)
  
  # ===== SAVE =====
  results <- data.frame(
    actual_age = y_train,
    predicted_age = predictions,
    error = errors
  )
  
  write.csv(results, "Jenkins_predictions.csv", row.names = FALSE)
  saveRDS(jenkins_model, "Jenkins_model.RDS")
  # Save CpG list
  cpg_list <- data.frame(CpG_ID = jenkins_use)
  write.csv(cpg_list, "Jenkins_CpGs_used.csv", row.names = FALSE)
  
  cat("\n✓ Results saved!\n")
  
} else {
  cat("❌ Not enough Jenkins CpGs available\n")
  cat("Need ~200+ CpGs, have:", sum(available_jenkins), "\n")
  cat("Try loading full unfiltered data\n")
}


# ===== SUMMARY OF CpGs PER REGION =====
cat("\n=== CpGs per Jenkins Region ===\n")

region_summary <- data.frame(
  Region = jenkins_data$Gene,
  Chr = jenkins_data$chr,
  N_CpGs_total = sapply(jenkins_cpgs, length),
  N_CpGs_available = sapply(jenkins_cpgs, function(cpgs) sum(cpgs %in% rownames(Train)))
)

print(head(region_summary, 20))

cat("\nRegions with most CpGs:\n")
top_regions <- region_summary[order(region_summary$N_CpGs_total, decreasing = TRUE), ][1:10, ]
print(top_regions)

write.csv(region_summary, "Jenkins_region_summary.csv", row.names = FALSE)

#Colour by datasets
#Tidy up cg plot
#Calc correlation 50k variably methylated lociwith age 
#plot top 6 most correlating/ anticorrelating- absolute vals

# JENKINS COLOURED BY STUDY #

# ===== CREATE STUDY LABELS =====
# Load study information from your CSV
train_age_data <- read_csv("Training data set.csv")
train_study <- train_age_data$GSE

# Match to the samples you're actually using (after removing NAs)
study_labels <- train_study[complete_idx]

# Check
cat("Study distribution:\n")
print(table(study_labels))

# ===== DEFINE COLORS FOR 3 STUDIES =====
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

# Create color vector
colors_by_study <- study_colors[study_labels]

# ===== PLOT COLORED BY STUDY =====
plot(y_train, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "Jenkins Predicted Age (years)",
     main = "Jenkins Model - Colored by Study")

# Perfect prediction line
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Regression line
abline(lm(predictions ~ y_train), col = "black", lwd = 2)

# Legend
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " yrs"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study")



#glmnet package 
install.packages("glmnet", repos = "https://cran.us.r-project.org")


library(methylclock)
Train_data <- readRDS("~/Documents/Forensic Project/preExistingClocks.Train.RDS")
head(Train_data)

# Run Horvath clock specifically
# If your data has actual chronological ages
# (e.g., in a row called "age" or as a separate vector)

# ===== LOAD DATA =====

cat("Data dimensions:", dim(Train_data), "\n")
cat("Samples (columns):", ncol(Train_data), "\n")


# ===== RUN HORVATH CLOCK (WITHOUT AGE) =====
cat("\nRunning Horvath clock...\n")

horvath_results <- DNAmAge(Train_data, clocks = "Horvath")

cat("✓ Done!\n")
cat("Results dimensions:", dim(horvath_results), "\n")

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(horvath_results)) {
  actual_ages <- actual_ages[1:nrow(horvath_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
horvath_results$actual_age <- actual_ages

# Check
head(horvath_results)

# ===== CALCULATE METRICS =====
horvath_results$error <- horvath_results$Horvath - horvath_results$actual_age
horvath_results$abs_error <- abs(horvath_results$error)

mae <- mean(horvath_results$abs_error, na.rm = TRUE)
rmse <- sqrt(mean(horvath_results$error^2, na.rm = TRUE))
correlation <- cor(horvath_results$Horvath, horvath_results$actual_age, 
                   use = "complete.obs")

cat("\n=== HORVATH CLOCK PERFORMANCE ===\n")
cat("MAE:", round(mae, 2), "years\n")
cat("RMSE:", round(rmse, 2), "years\n")
cat("Correlation:", round(correlation, 3), "\n")

# ===== CREATE SCATTER PLOT =====
plot(horvath_results$actual_age, 
     horvath_results$Horvath,
     pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.6),
     cex = 1.2,
     xlab = "Chronological Age (years)",
     ylab = "Horvath Predicted Age (years)",
     main = "Horvath Clock Performance")

abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
abline(lm(horvath_results$Horvath ~ horvath_results$actual_age), 
       col = "darkblue", lwd = 2)
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("RMSE = ", round(rmse, 2), " years")
       ),
       bty = "n",
       cex = 1)

legend("bottomright",
       legend = c("Perfect prediction", "Actual fit"),
       col = c("red", "darkblue"),
       lty = c(2, 1),
       lwd = 2,
       bty = "n",
       cex = 0.9)




# ===== RUN ALL 4 CLOCKS (WITHOUT AGE) =====
cat("\nRunning all 4 methylation clocks...\n")

clock_results <- DNAmAge(Train_data, 
                         clocks = c("Horvath", "skinHorvath", "Hannum", "EN"))

cat("✓ Done!\n")
cat("Results dimensions:", dim(clock_results), "\n")
head(clock_results)

# ===== LOAD ACTUAL AGES SEPARATELY =====
train_age_data <- read.csv("Training data set.csv")
actual_ages <- train_age_data$Age

cat("\nActual ages loaded:", length(actual_ages), "\n")

# Match to results length
if(length(actual_ages) > nrow(clock_results)) {
  actual_ages <- actual_ages[1:nrow(clock_results)]
  cat("Trimmed ages to:", length(actual_ages), "\n")
}

# ===== ADD ACTUAL AGES TO RESULTS =====
clock_results$actual_age <- actual_ages

# Check
head(clock_results)

# ===== CHECK FOR MISSING VALUES =====
cat("\n=== Checking for missing values ===\n")
cat("NAs in Horvath:", sum(is.na(clock_results$Horvath)), "\n")
cat("NAs in SkinHorvath:", sum(is.na(clock_results$skinHorvath)), "\n")
cat("NAs in Hannum:", sum(is.na(clock_results$Hannum)), "\n")
cat("NAs in EN:", sum(is.na(clock_results$EN)), "\n")
cat("NAs in actual_age:", sum(is.na(clock_results$actual_age)), "\n")

# ===== IDENTIFY WORKING CLOCKS =====
working_clocks <- c()
if(sum(!is.na(clock_results$Horvath)) > 0) working_clocks <- c(working_clocks, "Horvath")
if(sum(!is.na(clock_results$skinHorvath)) > 0) working_clocks <- c(working_clocks, "skinHorvath")
if(sum(!is.na(clock_results$Hannum)) > 0) working_clocks <- c(working_clocks, "Hannum")
if(sum(!is.na(clock_results$EN)) > 0) working_clocks <- c(working_clocks, "EN")

cat("\n=== Working clocks ===\n")
print(working_clocks)

# ===== CALCULATE METRICS (ONLY FOR WORKING CLOCKS) =====
calculate_metrics <- function(predicted, actual, clock_name) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    cat("WARNING:", clock_name, "has no complete predictions\n")
    return(list(MAE = NA, RMSE = NA, Correlation = NA, N = 0))
  }
  
  predicted_clean <- predicted[complete_idx]
  actual_clean <- actual[complete_idx]
  error <- predicted_clean - actual_clean
  mae <- mean(abs(error))
  rmse <- sqrt(mean(error^2))
  correlation <- cor(predicted_clean, actual_clean)
  
  return(list(
    MAE = round(mae, 2),
    RMSE = round(rmse, 2),
    Correlation = round(correlation, 3),
    N = length(predicted_clean)
  ))
}

# Calculate for each clock
metrics_horvath <- calculate_metrics(clock_results$Horvath, clock_results$actual_age, "Horvath")
metrics_skin <- calculate_metrics(clock_results$skinHorvath, clock_results$actual_age, "skinHorvath")
metrics_hannum <- calculate_metrics(clock_results$Hannum, clock_results$actual_age, "Hannum")
metrics_en <- calculate_metrics(clock_results$EN, clock_results$actual_age, "EN")

# Create metrics table
metrics_table <- data.frame(
  Clock = c("Horvath", "skinHorvath", "Hannum", "EN"),
  N_samples = c(metrics_horvath$N, metrics_skin$N, 
                metrics_hannum$N, metrics_en$N),
  MAE = c(metrics_horvath$MAE, metrics_skin$MAE, 
          metrics_hannum$MAE, metrics_en$MAE),
  RMSE = c(metrics_horvath$RMSE, metrics_skin$RMSE, 
           metrics_hannum$RMSE, metrics_en$RMSE),
  Correlation = c(metrics_horvath$Correlation, metrics_skin$Correlation,
                  metrics_hannum$Correlation, metrics_en$Correlation)
)

# Print table
cat("\n=== PERFORMANCE METRICS ===\n")
print(metrics_table)

# Identify failed clocks
failed_clocks <- metrics_table$Clock[metrics_table$N_samples == 0]
if(length(failed_clocks) > 0) {
  cat("\n⚠️  Failed clocks (all NA):", paste(failed_clocks, collapse = ", "), "\n")
  cat("This is expected for semen data - blood-specific clocks may not work\n")
}

# ===== CREATE 4-PANEL FIGURE =====
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Function to create scatter plot
plot_clock <- function(actual, predicted, clock_name, mae_val, rmse_val, cor_val, n_val) {
  # Remove NAs
  complete_idx <- complete.cases(predicted, actual)
  
  if(sum(complete_idx) == 0) {
    # If no data, show message
    plot.new()
    text(0.5, 0.5, paste(clock_name, "\nNo predictions available\n(All NA)"), 
         cex = 1.2, col = "red")
    return()
  }
  actual_clean <- actual[complete_idx]
  predicted_clean <- predicted[complete_idx]
  
  plot(actual_clean, predicted_clean,
       pch = 19,
       col = rgb(0.2, 0.4, 0.8, 0.6),
       cex = 1.2,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = clock_name,
       xlim = range(c(actual_clean, predicted_clean)),
       ylim = range(c(actual_clean, predicted_clean)))
  
  # Perfect prediction line
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  
  # Regression line
  abline(lm(predicted_clean ~ actual_clean), col = "darkblue", lwd = 2)
  
  # Metrics
  if(!is.na(cor_val)) {
    legend("topleft",
           legend = c(
             paste0("r = ", cor_val),
             paste0("MAE = ", mae_val, " yrs"),
             paste0("RMSE = ", rmse_val, " yrs"),
             paste0("n = ", n_val)
           ),
           bty = "n",
           cex = 0.75)
  }
}

# Create all 4 panels
plot_clock(clock_results$actual_age, clock_results$Horvath,
           "Horvath Clock",
           metrics_table$MAE[1], metrics_table$RMSE[1], 
           metrics_table$Correlation[1], metrics_table$N_samples[1])

plot_clock(clock_results$actual_age, clock_results$skinHorvath,
           "Skin+Blood Horvath Clock",
           metrics_table$MAE[2], metrics_table$RMSE[2], 
           metrics_table$Correlation[2], metrics_table$N_samples[2])

plot_clock(clock_results$actual_age, clock_results$Hannum,
           "Hannum Clock",
           metrics_table$MAE[3], metrics_table$RMSE[3], 
           metrics_table$Correlation[3], metrics_table$N_samples[3])

plot_clock(clock_results$actual_age, clock_results$EN, 
           "Elastic Net Clock",
           metrics_table$MAE[4], metrics_table$RMSE[4], 
           metrics_table$Correlation[4], metrics_table$N_samples[4])
par(mfrow = c(1, 1))


#Hannum did not give a plot as there was not enough data available

Train50k <-- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_top_50K.RDS")
Train <--readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")

#cg16867657 vs clock results

dim(Train)

which(rownames(Train)== "cg16867657")

plot(clock_results$actual_age, Train[99804,])

#cg17147820

which(rownames(Train) == "cg17147820")
plot(clock_results$actual_age, Train[394604,])

#cg10528482
which(rownames(Train) == "cg10528482")
plot(clock_results$actual_age, Train[39955,])


# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load libraries
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Example: Your data frame with chromosomes and positions
my_data <- data.frame(
  chr = c("chr1", "chr2", "chr3"),
  pos = c(10497, 10525, 10542)
)

# Get annotation from 450k array
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Create a key for matching
anno$key <- paste0(anno$chr, ":", anno$pos)
my_data$key <- paste0(my_data$chr, ":", my_data$pos)

# Match and extract CpG IDs
my_data$cpg_id <- anno$Name[match(my_data$key, anno$key)]

# View results
print(my_data)



#VISAGE LOCI
# ===== INSTALL AND LOAD ANNOTATION PACKAGE =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# For 450k array
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

# Load library
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


# ===== CREATE VISAGE DATA FRAME FROM YOUR IMAGE =====
# This is YOUR data (from the table in your image)
visage_data <- data.frame(
  Gene = c("SYT7", "TUBB3", "SH2B2", "ARHGEF17", "EXOC3", "GALR2", "PPP2R2C"),
  CpG_ID = c("cg17147820", "cg18701351", "cg00018181", "cg09855959", 
             "cg10528482", "cg07909178", "cg02766173"),
  chr = c("chr11", "chr16", "chr7", "chr11", "chr5", "chr17", "chr4"),
  pos = c(61554783, 89921897, 102288444, 73311506, 525656, 76077795, 6473455)
)

cat("=== VISAGE Enhanced Tool Markers ===\n")
print(visage_data)


# ===== GET ANNOTATION FROM 450k ARRAY =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation database loaded\n")
cat("Total CpGs in 450k array:", nrow(anno), "\n\n")

# ===== CREATE MATCHING KEYS =====
# Combine chromosome + position into unique key
anno$key <- paste0(anno$chr, ":", anno$pos)
visage_data$key <- paste0(visage_data$chr, ":", visage_data$pos)


# ===== VERIFY CpG IDs MATCH POSITIONS =====
# Match positions to get CpG IDs from annotation
visage_data$cpg_verified <- anno$Name[match(visage_data$key, anno$key)]

# Check if CpG IDs match
visage_data$match <- visage_data$CpG_ID == visage_data$cpg_verified
cat("=== Verification Results ===\n")
print(visage_data[, c("Gene", "CpG_ID", "cpg_verified", "match")])

if(all(visage_data$match, na.rm = TRUE)) {
  cat("\n✓ All CpG IDs verified correctly!\n\n")
} else {
  cat("\n⚠️  Some mismatches detected\n")
  cat("Check the rows where match = FALSE\n\n")
}


# ===== GET ADDITIONAL ANNOTATION INFO =====
# Add gene names and other info from annotation
visage_data$anno_gene <- anno$UCSC_RefGene_Name[match(visage_data$key, anno$key)]
visage_data$cpg_island <- anno$Relation_to_Island[match(visage_data$key, anno$key)]


cat("=== Additional Annotation Info ===\n")
print(visage_data[, c("Gene", "CpG_ID", "anno_gene", "cpg_island")])

# ===== NOW LOAD YOUR METHYLATION DATA =====
library(readr)

Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age_data <- read.csv("Training data set.csv")
train_age <- train_age_data$Age

# ===== FIX: MATCH LENGTHS =====
cat("Data samples:", ncol(Train), "\n")
cat("Ages:", length(train_age), "\n")

# Trim ages to match data
if(length(train_age) > ncol(Train)) {
  train_age <- train_age[1:ncol(Train)]
  cat("Trimmed ages to:", length(train_age), "\n")
}

cat("\n=== Checking VISAGE Markers in Your Data ===\n")

# Check which VISAGE markers are in your data
visage_data$in_your_data <- visage_data$CpG_ID %in% rownames(Train)

print(visage_data[, c("Gene", "CpG_ID", "in_your_data")])

cat("\nAvailable in your data:", sum(visage_data$in_your_data), "/ 7\n\n")


  # ===== BUILD LINEAR MODEL WITH AVAILABLE MARKERS =====
available_markers <- visage_data$CpG_ID[visage_data$in_your_data]
available_genes <- visage_data$Gene[visage_data$in_your_data]

if(length(available_markers) >= 3) {
  
  cat("✓ Building model with", length(available_markers), "markers:\n")
  cat("  ", paste(available_genes, collapse = ", "), "\n\n")
  
  # Extract methylation data for these markers
  train_visage <- t(Train[available_markers, ])
  
  # Create modeling data frame
  train_df <- data.frame(
    age = train_age,
    train_visage
  )
  # Use gene names as column names
  colnames(train_df)[-1] <- available_genes
  
  
  # ===== TRAIN LINEAR MODEL =====
  visage_model <- lm(age ~ ., data = train_df)
  
  cat("=== VISAGE Linear Model Summary ===\n")
  print(summary(visage_model))
  # ===== MODEL EQUATION =====
  cat("\n=== Prediction Equation ===\n")
  coefs <- coef(visage_model)
  cat("Age = ", round(coefs[1], 2), "\n", sep = "")
  for(i in 2:length(coefs)) {
    cat("    ", ifelse(coefs[i] >= 0, "+", ""), " ", 
        round(coefs[i], 2), " × ", names(coefs)[i], "\n", sep = "")
  }
  
  # ===== PREDICTIONS =====
  predictions <- predict(visage_model)
  errors <- predictions - train_age
  
  
  # ===== PERFORMANCE METRICS =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, train_age)
  r2 <- summary(visage_model)$r.squared
  adj_r2 <- summary(visage_model)$adj.r.squared
  
  cat("\n========================================\n")
  cat("VISAGE MODEL PERFORMANCE\n")
  cat("========================================\n")
  cat("Markers:", length(available_markers), "/ 7\n")
  cat("Samples:", length(train_age), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("Adjusted R²:", round(adj_r2, 3), "\n")
  
  
  # ===== SCATTER PLOT =====
  plot(train_age, predictions,
       pch = 19, col = rgb(0.2, 0.6, 0.4, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = paste0("VISAGE Enhanced Tool (", length(available_markers), " CpGs)"))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ train_age), col = "darkgreen", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3))
         ),
         bty = "n", cex = 1.1)
  
  # ===== SAVE EVERYTHING =====
  results_df <- data.frame(
    sample_id = rownames(train_visage),
    actual_age = train_age,
    predicted_age = predictions,
    error = errors,
    abs_error = abs(errors)
  )
  
  write.csv(results_df, "VISAGE_results.csv", row.names = FALSE)
  write.csv(visage_data, "VISAGE_marker_info.csv", row.names = FALSE)
  saveRDS(visage_model, "VISAGE_model.RDS")
  
  cat("\n✓ All files saved!\n")
} else {
  cat("❌ Not enough markers available (need at least 3)\n")
  cat("Available:", sum(visage_data$in_your_data), "/ 7\n")
}

#Colouring by study#
# After you've made predictions with VISAGE model

# Get study labels
train_age_data <- read_csv("Training data set.csv")
study_labels <- train_age_data$GSE[1:length(predictions)]

# Define colors
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

colors_by_study <- study_colors[study_labels]
# Plot
plot(train_age, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "VISAGE Predicted Age (years)",
     main = "VISAGE Model - By Study")

abline(0, 1, col = "red", lwd = 2, lty = 2)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       title = "Study")
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " years"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1.1)
# ===== FINAL SUMMARY TABLE =====
cat("\n========================================\n")
cat("SUMMARY TABLE\n")
cat("========================================\n")

summary_table <- data.frame(
  Gene = visage_data$Gene,
  CpG_ID = visage_data$CpG_ID,
  Chr = visage_data$chr,
  Position = visage_data$pos,
  In_Data = visage_data$in_your_data,
  Verified = ifelse(is.na(visage_data$match), "N/A", visage_data$match)
)

print(summary_table)



#JENKINS

library(readr)

# ===== INSTALL AND LOAD ANNOTATION =====
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ===== CREATE JENKINS DATA  =====
# Combining data 
jenkins_data <- data.frame(
  Gene = c(
    # Image 1
    "ADAMTS8", "ARC", "ARGHGEF10", "BCL11A", "C1ORF122", "C7ORF50",
    "CCDC144NL", "CLIC1", "DMPK", "FAM86C1", "FAM86JP", "FOXK1",
    "FSCN", "GAPDH", "GET4", "GNB2", "GPANK1", "GPR45", "KCNQ1",
    # Image 2
    "LDLRAD4", "LMO3", "LOC100133461", "MIR22HG", "MTMR8", "N10",
    "N12", "N22", "N23", "N24", "N27", "N30", "N8", "N9",
    "NCOR2", "NONE", "NSG1", "PAX2", "PITX1", "PRSS22",
    "PTPRN2.3", "PTPRN2.4", "PURA", "PYY2", "SECTM1", "SEMA6B",
    "SEZ6", "SLC22A18AS", "SOHLH1", "THBS3", "TNXB"
  ),
  chr = c(
    # Image 1
    "chr11", "chr8", "chr8", "chr2", "chr1", "chr7",
    "chr17", "chr6", "chr19", "chr11", "chr3", "chr7",
    "chr7", "chr12", "chr7", "chr7", "chr6", "chr2", "chr11",
    # Image 2
    "chr18", "chr12", "chr4", "chr17", "chrX", "chr1",
    "chr5", "chr19", "chr14", "chr6", "chr6", "chr15", "chr11", "chr7",
    "chr12", "chr10", "chr4", "chr10", "chr5", "chr16",
    "chr7", "chr7", "chr5", "chr17", "chr17", "chr19",
    "chr17", "chr11", "chr9", "chr1", "chr6"
  ),
  start = c(
    # Image 1
    130299298, 143694010, 1877888, 60680616, 38272200, 1083209,
    20798895, 31698492, 46282571, 71498202, 125634060, 4722778,
    5635134, 6641602, 914964, 100274361, 31630819, 105857809, 2554562,
    # Image 2
    13611370, 16760040, 3680721, 1617363, 63614857, 28423399,
    3593413, 4579481, 106004434, 170449417, 30432200, 27959473, 69260136, 35300077,
    124990897, 17347047, 4386726, 102509693, 134365728, 2908157,
    157523356, 158109339, 139492535, 26553567, 80278592, 4555999,
    27330794, 2909690, 138590204, 155176868, 32064146
  ),
  stop = c(
    # Image 1
    130299948, 143694548, 1878324, 60680762, 38273057, 1084163,
    20799770, 31699299, 46283081, 71499118, 125634453, 4723928,
    5635954, 6642355, 915832, 100275305, 31632542, 105859084, 2555577,
    # Image 2
    13611825, 16761003, 3681760, 1618296, 63615496, 28424202,
    3594276, 4580471, 106004608, 170450804, 30433944, 27960032, 69261045, 35301070,
    124991140, 17347392, 4387698, 102510569, 134366535, 2908935,
    157524159, 158110153, 139493491, 26554908, 80280331, 4556983,
    27332647, 2909716, 138590996, 155177784, 32065891
  )
)

cat("=== Jenkins et al. 50 Genomic Regions ===\n")
cat("Total regions:", nrow(jenkins_data), "\n\n")
print(head(jenkins_data, 50))

# ===== LOAD ANNOTATION DATABASE =====
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("\nAnnotation loaded:", nrow(anno), "CpGs\n\n")

# ===== FIND CpGs WITHIN EACH JENKINS REGION =====
# Jenkins regions are genomic windows, not single CpGs
# Need to find ALL CpGs within each region

jenkins_cpgs <- list()

cat("=== Finding CpGs in Jenkins regions ===\n")

for(i in 1:nrow(jenkins_data)) {
  region_chr <- jenkins_data$chr[i]
  region_start <- jenkins_data$start[i]
  region_stop <- jenkins_data$stop[i]
  region_name <- jenkins_data$Gene[i]
  
  # Find CpGs in this region
  cpgs_in_region <- anno$Name[
    anno$chr == region_chr & 
      anno$pos >= region_start & 
      anno$pos <= region_stop
  ]
  
  jenkins_cpgs[[i]] <- cpgs_in_region
  
  cat(i, ". ", region_name, " (", region_chr, ":", region_start, "-", region_stop, "): ",
      length(cpgs_in_region), " CpGs\n", sep = "")
}
# Total CpGs across all regions
all_jenkins_cpgs <- unique(unlist(jenkins_cpgs))
cat("\nTotal unique CpGs across 50 regions:", length(all_jenkins_cpgs), "\n\n")

# ===== CHECK AVAILABILITY IN YOUR DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Training_Matrix.RDS")
train_age <- read.csv("Training data set.csv")$Age
train_age <- train_age[1:ncol(Train)]  # Fix length

available_jenkins <- all_jenkins_cpgs %in% rownames(Train)
cat("=== Jenkins CpGs in Your Data ===\n")
cat("Total Jenkins CpGs:", length(all_jenkins_cpgs), "\n")
cat("Available in your data:", sum(available_jenkins), "\n")
cat("Percentage:", round(sum(available_jenkins)/length(all_jenkins_cpgs)*100, 1), "%\n\n")

# ===== BUILD MODEL WITH AVAILABLE JENKINS CpGs =====
if(sum(available_jenkins) >= 50) {
  
  jenkins_use <- all_jenkins_cpgs[available_jenkins]
  
  cat("✓ Building Jenkins model with", length(jenkins_use), "CpGs\n\n")
  
  # Extract data
  train_jenkins <- t(Train[jenkins_use, ])
  
  # Create data frame
  train_df <- data.frame(
    age = train_age,
    train_jenkins
  )
  
  # ===== TRAIN MODEL =====
  # Use elastic net due to large number of features
  library(glmnet)
  
  X_train <- as.matrix(train_df[, -1])
  y_train <- train_df$age
  
  complete_idx <- complete.cases(X_train, y_train)
  X_train <- X_train[complete_idx, ]
  y_train <- y_train[complete_idx]
  
  
  # Cross-validated elastic net
  jenkins_model <- cv.glmnet(X_train, y_train, 
                             alpha = 0.5,  # Elastic net
                             nfolds = 10)
  
  cat("=== Jenkins Model Trained ===\n")
  cat("Lambda min:", jenkins_model$lambda.min, "\n")
  cat("Lambda 1se:", jenkins_model$lambda.1se, "\n\n")
  

  # ===== PREDICTIONS =====
  predictions <- predict(jenkins_model, newx = X_train, s = "lambda.min")[,1]
  errors <- predictions - y_train
  
  
  # ===== PERFORMANCE =====
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  correlation <- cor(predictions, y_train)
  r2 <- summary(lm(predictions ~ y_train))$r.squared
  
  cat("=== Jenkins Model Performance ===\n")
  cat("CpGs used:", length(jenkins_use), "\n")
  cat("MAE:", round(mae, 2), "years\n")
  cat("RMSE:", round(rmse, 2), "years\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R²:", round(r2, 3), "\n")
  cat("\nPublished Jenkins performance: MAE = 2.04-2.37 years\n")
  
  
  # ===== PLOT =====
  plot(y_train, predictions,
       pch = 19, col = rgb(0.4, 0.2, 0.8, 0.6), cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Jenkins Predicted Age (years)",
       main = paste0("Jenkins Germ Line Calculator (", length(jenkins_use), " CpGs)"))
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(predictions ~ y_train), col = "purple", lwd = 2)
  
  legend("topleft",
         legend = c(
           paste0("r = ", round(correlation, 3)),
           paste0("MAE = ", round(mae, 2), " years"),
           paste0("R² = ", round(r2, 3)),
           paste0("n = ", length(jenkins_use), " CpGs")
         ),
         bty = "n", cex = 1)
  
  # ===== SAVE =====
  results <- data.frame(
    actual_age = y_train,
    predicted_age = predictions,
    error = errors
  )
  
  write.csv(results, "Jenkins_predictions.csv", row.names = FALSE)
  saveRDS(jenkins_model, "Jenkins_model.RDS")
  # Save CpG list
  cpg_list <- data.frame(CpG_ID = jenkins_use)
  write.csv(cpg_list, "Jenkins_CpGs_used.csv", row.names = FALSE)
  
  cat("\n✓ Results saved!\n")
  
} else {
  cat("❌ Not enough Jenkins CpGs available\n")
  cat("Need ~200+ CpGs, have:", sum(available_jenkins), "\n")
  cat("Try loading full unfiltered data\n")
}


# ===== SUMMARY OF CpGs PER REGION =====
cat("\n=== CpGs per Jenkins Region ===\n")

region_summary <- data.frame(
  Region = jenkins_data$Gene,
  Chr = jenkins_data$chr,
  N_CpGs_total = sapply(jenkins_cpgs, length),
  N_CpGs_available = sapply(jenkins_cpgs, function(cpgs) sum(cpgs %in% rownames(Train)))
)

print(head(region_summary, 20))

cat("\nRegions with most CpGs:\n")
top_regions <- region_summary[order(region_summary$N_CpGs_total, decreasing = TRUE), ][1:10, ]
print(top_regions)

write.csv(region_summary, "Jenkins_region_summary.csv", row.names = FALSE)

#Colour by datasets
#Tidy up cg plot
#Calc correlation 50k variably methylated lociwith age 
#plot top 6 most correlating/ anticorrelating- absolute vals

# JENKINS COLOURED BY STUDY #

# ===== CREATE STUDY LABELS =====
# Load study information from your CSV
train_age_data <- read_csv("Training data set.csv")
train_study <- train_age_data$GSE

# Match to the samples you're actually using (after removing NAs)
study_labels <- train_study[complete_idx]

# Check
cat("Study distribution:\n")
print(table(study_labels))

# ===== DEFINE COLORS FOR 3 STUDIES =====
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

# Create color vector
colors_by_study <- study_colors[study_labels]

# ===== PLOT COLORED BY STUDY =====
plot(y_train, predictions,
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "Jenkins Predicted Age (years)",
     main = "Jenkins Model - Colored by Study")

# Perfect prediction line
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Regression line
abline(lm(predictions ~ y_train), col = "black", lwd = 2)

# Legend
legend("topleft",
       legend = c(
         paste0("r = ", round(correlation, 3)),
         paste0("MAE = ", round(mae, 2), " yrs"),
         paste0("R² = ", round(r2, 3))
       ),
       bty = "n", cex = 1)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study")



#glmnet package 
install.packages("glmnet", repos = "https://cran.us.r-project.org")



# ===== SIDE BY SIDE VISAGE AND JENKINS =====
png("visage_jenkins_comparison.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# ===== LEFT PANEL: VISAGE =====
plot(train_age, visage_predictions,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = "VISAGE Enhanced Tool (7 CpGs)")

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(visage_predictions ~ train_age), col = "black", lwd = 2)
for(study_name in plot_order) {
  idx <- study_labels_visage == study_name
  if(sum(idx) > 0) {
    points(train_age[idx], visage_predictions[idx],
           pch = 19, col = study_colors[study_name], cex = 1.5)
  }
}

legend("topleft",
       legend = c(
         paste0("r = ", round(visage_cor, 3)),
         paste0("MAE = ", round(visage_mae, 2), " years"),
         paste0("R² = ", round(visage_r2, 3))
       ),
       bty = "n", cex = 0.9)
# ===== RIGHT PANEL: JENKINS =====
plot(y_train, jenkins_predictions,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = "Jenkins Germline Calculator (~200 CpGs)")

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(jenkins_predictions ~ y_train), col = "black", lwd = 2)

for(study_name in plot_order) {
  idx <- study_labels_jenkins == study_name
  if(sum(idx) > 0) {
    points(y_train[idx], jenkins_predictions[idx],
           pch = 19, col = study_colors[study_name], cex = 1.5)
  }
}
legend("topleft",
       legend = c(
         paste0("r = ", round(jenkins_cor, 3)),
         paste0("MAE = ", round(jenkins_mae, 2), " years"),
         paste0("R² = ", round(jenkins_r2, 3))
       ),
       bty = "n", cex = 0.9)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study",
       cex = 0.9)

par(mfrow = c(1, 1))
dev.off()
