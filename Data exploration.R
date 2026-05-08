#Importing training and test sets 
train <- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_Top_50K.RDS")
test <- readRDS("~/Documents/Forensic Project/Jamie/Final_Test_Matrix_Top_50K.RDS")

#join together 
all(rownames(train) == rownames(test))
join <- cbind(train, test)

#calculate standard deviation
join.sd <- apply(join, 1, sd)

head(join.sd)
join.10k <- join[1:10000,]

#perform PCA
join.PCA <- prcomp(t(join.10k), scale=TRUE)

#plot PCA
plot(join.PCA$x, pch=19)

#colour by test/training 
colour <- c(rep("blue", ncol(train)), rep("red", ncol(test)))
plot(join.PCA$x, pch=19, col= colour)

# LABEL BY AGE 

# Your existing code (keep this)
train <- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_Top_50K.RDS")
test <- readRDS("~/Documents/Forensic Project/Jamie/Final_Test_Matrix_Top_50K.RDS")

# Check the structure of your data
str(train)
str(test)

# Extract age information
# Option C: If age is stored separately, you'll need to load it
train_age_data <- read.csv("Training data set.csv") 
test_age_data <- read.csv("Test data set.csv")

# Extract ONLY the age column (adjust "Age" to match your actual column name)
train_age <- train_age_data$Age  
test_age <- test_age_data$Age    

# Combine into one vector
all_ages <- c(train_age, test_age)

# Verify ages match data
cat("Train samples:", ncol(train), "| Train ages:", length(train_age), "\n")
cat("Test samples:", ncol(test), "| Test ages:", length(test_age), "\n")
cat("Age range:", range(all_ages), "\n")

# ===== PLOT 1: BY TRAIN/TEST (your existing plot) =====
colour <- c(rep("steelblue", ncol(train)), rep("red", ncol(test)))

# Calculate variance explained
pca_var <- join.PCA$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 1)

plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colour,
     cex = 1.5,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "PCA - Training vs Test Sets")
legend("topright",
       legend = c("Training", "Test"),
       col = c("steelblue", "red"),
       pch = 19,
       pt.cex = 1.5)

# ===== PLOT 2: BY AGE (NEW!) =====
library(RColorBrewer)

# Create color gradient from blue (young) to red (old)
age_colors <- colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))(100)

# Assign each sample to a color bin based on age
age_bins <- cut(all_ages, breaks = 100, labels = FALSE)
colors_by_age <- age_colors[age_bins]

# Create the age-colored plot
plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colors_by_age,
     cex = 1.5,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "PCA - Colored by Age")
# Add age legend
legend_ages <- seq(min(all_ages), max(all_ages), length.out = 5)
legend_colors <- age_colors[seq(1, 100, length.out = 5)]

legend("topright",
       legend = paste(round(legend_ages, 0), "years"),
       col = legend_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Age",
       bg = "white")

# ===== PLOT 3: SIDE BY SIDE COMPARISON =====
par(mfrow = c(1, 2))

# Left plot: Train/Test
plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colour,
     cex = 1.2,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "By Dataset")
legend("topright", 
       legend = c("Train", "Test"),
       col = c("steelblue", "red"), 
       pch = 19,
       cex = 0.8)

# Right plot: Age
plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colors_by_age,
     cex = 1.2,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "By Age")
legend("topright",
       legend = paste(round(legend_ages, 0), "yrs"),
       col = legend_colors,
       pch = 19,
       cex = 0.8,
       title = "Age")

par(mfrow = c(1, 1))  # Reset to single plot


# LABELLING BY STUDY

# Check if your CSV has a study column
names(train_age_data)
names(test_age_data)

# If there's a 'study' or 'dataset' column:
train_study <- train_age_data$GSE  # adjust column name
test_study <- test_age_data$GSE

# Combine
study_labels <- c(train_study, test_study)

# Check
table(study_labels)

# Or manually define colors for specific studies:
study_colors <- c("GSE149318" = "green",
 "GSE185445" = "pink",
  "GSE185920" = "lightblue")

# Create color vector for each sample
colors_by_study <- study_colors[study_labels]

# Plot
plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "PCA - Labeled by Study")
legend("topright",
       legend = names(study_colors),          # Study names from your color vector
       col = study_colors,                     # Colors you defined
       pch = 19,
       pt.cex = 1.5,
       title = "Study",
       bg = "white")

# Your existing code (keep this)
train_study <- train_age_data$GSE
test_study <- test_age_data$GSE
study_labels <- c(train_study, test_study)

study_colors <- c(
  "GSE149318" = "green",
  "GSE185445" = "pink",
  "GSE185920" = "lightblue"
)

colors_by_study <- study_colors[study_labels]

# Plot everything first
plot(join.PCA$x[, 1], join.PCA$x[, 2],
     pch = 19,
     col = colors_by_study,
     cex = 1.5,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "PCA - Labeled by Study")

# ===== ADD THESE LINES =====
# Re-plot GSE185445 on top (bigger)
idx_185445 <- which(study_labels == "GSE185445")
points(join.PCA$x[idx_185445, 1], 
       join.PCA$x[idx_185445, 2],
       pch = 19,
       col = "pink",
       cex = 1.5)

# Re-plot GSE149318 on top (biggest)
idx_149318 <- which(study_labels == "GSE149318")
points(join.PCA$x[idx_149318, 1], 
       join.PCA$x[idx_149318, 2],
       pch = 19,
       col = "green",
       cex = 1.5)
# ===== END =====
legend("topright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study",
       bg = "white")  


#PCA OF GSE185920 ONLY - COLORED BY AGE GROUPS

# ===== JOIN TRAIN AND TEST =====
join <- cbind(Train, Test)
all_ages <- c(train_ages_all, test_ages_all)
all_studies <- c(train_studies_all, test_studies_all)

cat("Total samples:", ncol(join), "\n")
cat("Studies:", table(all_studies), "\n\n")


# ===== FILTER TO GSE185920 ONLY =====
gse185920_idx <- all_studies == "GSE185920"

join_gse185920 <- join[, gse185920_idx]
ages_gse185920 <- all_ages[gse185920_idx]

cat("GSE185920 samples:", ncol(join_gse185920), "\n")
cat("Age range:", range(ages_gse185920), "\n\n")

# ===== SELECT TOP 10K VARIABLE CpGs (FROM GSE185920 ONLY) =====
join.sd <- apply(join_gse185920, 1, sd)
top_indices <- order(join.sd, decreasing = TRUE)[1:10000]
join.10k <- join_gse185920[top_indices, ]

cat("Selected top 10,000 variable CpGs\n\n")


# ===== PERFORM PCA =====
cat("Running PCA...\n")
pca_gse <- prcomp(t(join.10k), scale = TRUE)

# Calculate variance
pca_var <- pca_gse$sdev^2
pca_var_percent <- round(pca_var / sum(pca_var) * 100, 1)

cat("✓ PCA complete\n")
cat("PC1:", pca_var_percent[1], "%\n")
cat("PC2:", pca_var_percent[2], "%\n\n")

# ===== CREATE AGE GROUPS  =====
# Based on the "steps" in age prediction plot
age_groups <- cut(ages_gse185920,
                  breaks = c(0, 35, 43, 100),
                  labels = c("<35 years", "35-43 years", "44+ years"),
                  include.lowest = TRUE)

cat("Age group distribution:\n")
print(table(age_groups))
cat("\n")


# ===== DEFINE COLORS FOR AGE GROUPS =====
age_group_colors <- c(
  "<35 years" = "#3498db",     # Blue (young)
  "35-43 years" = "#f39c12",   # Orange (middle)
  "44+ years" = "#e74c3c"      # Red (older)
)

colors_by_age_group <- age_group_colors[age_groups]

# ===== CREATE PCA PLOT BY AGE GROUPS =====
png("PCA_GSE185920_age_groups.png", width = 1000, height = 800, res = 120)

plot(pca_gse$x[, 1], pca_gse$x[, 2],
     pch = 19,
     col = colors_by_age_group,
     cex = 1.8,
     xlab = paste0("PC1 (", pca_var_percent[1], "% variance)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "% variance)"),
     main = "PCA - GSE185920 Only (Colored by Age Group)",
     cex.lab = 1.2,
     cex.main = 1.3)

# Add legend
legend("bottomleft",
       legend = names(age_group_colors),
       col = age_group_colors,
       pch = 19,
       pt.cex = 2,
       cex = 1.2,
       title = "Age Group",
       bg = "white",
       box.lwd = 2)


dev.off()

cat("✓ PCA plot saved as 'PCA_GSE185920_age_groups.png'\n")

# ===== CHECK PC1 CORRELATION WITH AGE =====
cor_pc1_age <- cor(pca_gse$x[, 1], ages_gse185920)
cor_pc2_age <- cor(pca_gse$x[, 2], ages_gse185920)

cat("\n=== Correlation Analysis ===\n")
cat("PC1 vs Age: r =", round(cor_pc1_age, 3), "\n")
cat("PC2 vs Age: r =", round(cor_pc2_age, 3), "\n\n")

# ===== ALTERNATIVE: LARGER POINTS FOR CLARITY =====
png("PCA_GSE185920_age_groups_large.png", width = 1200, height = 900, res = 120)

par(mar = c(5, 5, 4, 2))

plot(pca_gse$x[, 1], pca_gse$x[, 2],
     type = "n",
     xlab = paste0("PC1 (", pca_var_percent[1], "% variance)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "% variance)"),
     main = "PCA of GSE185920 Dataset by Age Category",
     cex.lab = 1.3,
     cex.main = 1.4)

# Plot each age group separately for better control
for(group_name in names(age_group_colors)) {
  idx <- age_groups == group_name
  
  points(pca_gse$x[idx, 1], 
         pca_gse$x[idx, 2],
         pch = 19,
         col = age_group_colors[group_name],
         cex = 2)
}

# Add legend with counts
group_counts <- table(age_groups)
legend_labels <- paste0(names(age_group_colors), " (n=", group_counts, ")")

legend("bottomleft",
       legend = legend_labels,
       col = age_group_colors,
       pch = 19,
       pt.cex = 2.5,
       cex = 1.2,
       title = "Age Group",
       bg = "white",
       box.lwd = 2)

dev.off()

cat("✓ Large point PCA saved as 'PCA_GSE185920_age_groups_large.png'\n")

# ===== PC1 vs AGE SCATTER PLOT =====
png("PC1_vs_age_GSE185920.png", width = 1000, height = 700, res = 120)

plot(ages_gse185920, pca_gse$x[, 1],
     pch = 19,
     col = colors_by_age_group,
     cex = 1.8,
     xlab = "Chronological Age (years)",
     ylab = "PC1 Score",
     main = paste0("PC1 vs Age - GSE185920 Only\n(r = ", round(cor_pc1_age, 3), ")"),
     cex.lab = 1.2,
     cex.main = 1.3)

# Regression line
abline(lm(pca_gse$x[, 1] ~ ages_gse185920), 
       col = "black", 
       lwd = 3)
legend("topright",
       legend = legend_labels,
       col = age_group_colors,
       pch = 19,
       pt.cex = 2,
       cex = 1.1,
       title = "Age Group",
       bg = "white")

dev.off()

cat("✓ PC1 vs Age plot saved as 'PC1_vs_age_GSE185920.png'\n")

# ===== 3-PANEL FIGURE: PCA + PC1 vs AGE + AGE DISTRIBUTION =====
png("GSE185920_comprehensive_PCA.png", width = 1800, height = 600, res = 120)

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

# Panel 1: PCA
plot(pca_gse$x[, 1], pca_gse$x[, 2],
     pch = 19,
     col = colors_by_age_group,
     cex = 1.6,
     xlab = paste0("PC1 (", pca_var_percent[1], "%)"),
     ylab = paste0("PC2 (", pca_var_percent[2], "%)"),
     main = "A. PCA by Age Group",
     cex.lab = 1.2,
     cex.main = 1.3)
legend("bottomleft",
       legend = names(age_group_colors),
       col = age_group_colors,
       pch = 19,
       pt.cex = 1.8,
       cex = 1,
       title = "Age")

# Panel 2: PC1 vs Age
plot(ages_gse185920, pca_gse$x[, 1],
     pch = 19,
     col = colors_by_age_group,
     cex = 1.6,
     xlab = "Age (years)",
     ylab = "PC1 Score",
     main = paste0("B. PC1 vs Age (r=", round(cor_pc1_age, 3), ")"),
     cex.lab = 1.2,
     cex.main = 1.3)

abline(lm(pca_gse$x[, 1] ~ ages_gse185920), col = "black", lwd = 3)
# Panel 3: Age distribution
hist(ages_gse185920,
     breaks = 20,
     col = "steelblue",
     border = "white",
     xlab = "Age (years)",
     ylab = "Frequency",
     main = paste0("C. Age Distribution (n=", length(ages_gse185920), ")"),
     cex.lab = 1.2,
     cex.main = 1.3)

# Add vertical lines for age group boundaries
abline(v = 35, col = "red", lwd = 3, lty = 2)
abline(v = 43, col = "red", lwd = 3, lty = 2)
# Label groups
text(30, max(hist(ages_gse185920, plot=FALSE)$counts) * 0.9, 
     "<35", col = age_group_colors[1], cex = 1.2, font = 2)
text(39, max(hist(ages_gse185920, plot=FALSE)$counts) * 0.9, 
     "35-43", col = age_group_colors[2], cex = 1.2, font = 2)
text(48, max(hist(ages_gse185920, plot=FALSE)$counts) * 0.9, 
     "44+", col = age_group_colors[3], cex = 1.2, font = 2)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Comprehensive 3-panel saved as 'GSE185920_comprehensive_PCA.png'\n")


