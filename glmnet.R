install.packages("glmnet", repos = "https://cran.us.r-project.org")
library(glmnet)
library(readr)

# ===== LOAD YOUR DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_Top_50K.RDS")
train_age_data <- read_csv("Training data set.csv")
train_age <- train_age_data$Age[1:ncol(Train)]
train_age <- train_age_data$Age[1:ncol(Train)]
cat("Data loaded:\n")
cat("CpGs:", nrow(Train), "\n")
cat("Samples:", ncol(Train), "\n")
cat("Ages:", length(train_age), "\n\n")

# ===== PREPARE DATA FOR GLMNET =====
# glmnet needs:
# - x: matrix of predictors (samples as rows, features as columns)
# - y: vector of outcomes (ages)
# Transpose so samples are rows
x <- t(Train)  # Now: samples × CpGs
y <- train_age

cat("x dimensions:", dim(x), "\n")
cat("y length:", length(y), "\n")

# ===== CHECK AND REMOVE NAs =====
cat("\nChecking for NAs...\n")
cat("NAs in x:", sum(is.na(x)), "\n")
cat("NAs in y:", sum(is.na(y)), "\n")

# Remove incomplete cases
complete_idx <- complete.cases(x, y)
cat("Complete samples:", sum(complete_idx), "\n")

x <- x[complete_idx, ]
y <- y[complete_idx]

cat("After cleaning:\n")
cat("x:", dim(x), "\n")
cat("y:", length(y), "\n")
cat("NAs remaining:", sum(is.na(x)), "\n\n")

# ===== FIT LASSO MODEL (alpha = 1) =====
cat("Fitting Lasso model with cross-validation...\n")

# alpha = 1 for Lasso (alpha = 0 for Ridge, alpha = 0.5 for Elastic Net)
lasso_model <- cv.glmnet(x, y, alpha = 1, nfolds = 10)

cat("✓ Model fitted!\n\n")

# ===== MODEL SUMMARY =====
cat("=== Lasso Model Summary ===\n")
cat("Lambda min:", lasso_model$lambda.min, "\n")
cat("Lambda 1se:", lasso_model$lambda.1se, "\n")

# Number of features selected
coef_min <- coef(lasso_model, s = "lambda.min")
n_features_min <- sum(coef_min != 0) - 1  # Exclude intercept

coef_1se <- coef(lasso_model, s = "lambda.1se")
n_features_1se <- sum(coef_1se != 0) - 1

cat("Features selected (lambda.min):", n_features_min, "\n")
cat("Features selected (lambda.1se):", n_features_1se, "\n\n")

# ===== PLOT 1: CROSS-VALIDATION CURVE =====
plot(lasso_model)
title("Lasso Model - Cross-Validation", line = 2.5)

cat("✓ CV plot created\n")
# ===== PLOT 2: COEFFICIENT PATH =====
# Plot how coefficients change with lambda
plot(lasso_model$glmnet.fit, xvar = "lambda", label = TRUE)
title("Lasso Coefficient Paths", line = 2.5)

cat("✓ Coefficient path plot created\n")

# ===== MAKE PREDICTIONS =====
# Using lambda.min (best CV performance)
predictions_min <- predict(lasso_model, newx = x, s = "lambda.min")[,1]

# Using lambda.1se (more regularized, simpler model)
predictions_1se <- predict(lasso_model, newx = x, s = "lambda.1se")[,1]

# ===== CALCULATE PERFORMANCE =====
# Lambda.min performance
errors_min <- predictions_min - y
mae_min <- mean(abs(errors_min))
rmse_min <- sqrt(mean(errors_min^2))
cor_min <- cor(predictions_min, y)
r2_min <- summary(lm(predictions_min ~ y))$r.squared

# Lambda.1se performance
errors_1se <- predictions_1se - y
mae_1se <- mean(abs(errors_1se))
rmse_1se <- sqrt(mean(errors_1se^2))
cor_1se <- cor(predictions_1se, y)
r2_1se <- summary(lm(predictions_1se ~ y))$r.squared

cat("\n========================================\n")
cat("LASSO MODEL PERFORMANCE\n")
cat("========================================\n")
cat("\nLambda.min (best CV):\n")
cat("  Features:", n_features_min, "\n")
cat("  MAE:", round(mae_min, 2), "years\n")
cat("  RMSE:", round(rmse_min, 2), "years\n")
cat("  Correlation:", round(cor_min, 3), "\n")
cat("  R²:", round(r2_min, 3), "\n")

cat("\nLambda.1se (simpler model):\n")
cat("  Features:", n_features_1se, "\n")
cat("  MAE:", round(mae_1se, 2), "years\n")
cat("  RMSE:", round(rmse_1se, 2), "years\n")
cat("  Correlation:", round(cor_1se, 3), "\n")
cat("  R²:", round(r2_1se, 3), "\n")

# ===== GET STUDY LABELS =====
study_labels <- train_age_data$GSE[complete_idx]

study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

colors_by_study <- study_colors[study_labels]

# ===== PLOT 3: PREDICTIONS COLORED BY STUDY (lambda.min) =====
plot_order <- c("GSE185920", "GSE185445", "GSE149318")

plot(y, predictions_min,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Lasso Predicted Age (years)",
     main = "Lasso Model (lambda.min) - By Study")

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(predictions_min ~ y), col = "black", lwd = 2)

# Plot each study
for(study_name in plot_order) {
  idx <- study_labels == study_name
  
  if(study_name == "GSE185920") {
    points(y[idx], predictions_min[idx],
           pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1.2)
  } else {
    points(y[idx], predictions_min[idx],
           pch = 19, col = study_colors[study_name], cex = 1.5)
  }
}
legend("topleft",
       legend = c(
         paste0("r = ", round(cor_min, 3)),
         paste0("MAE = ", round(mae_min, 2), " yrs"),
         paste0(n_features_min, " CpGs")
       ),
       bty = "n", cex = 1)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study")
cat("\n✓ Prediction plot created\n")

# ===== PLOT 4: COMPARISON OF LAMBDA.MIN vs LAMBDA.1SE =====
par(mfrow = c(1, 2))

# Lambda.min
plot(y, predictions_min,
     pch = 19, col = colors_by_study, cex = 1.2,
     xlab = "Actual Age", ylab = "Predicted Age",
     main = paste0("Lambda.min (", n_features_min, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = paste0("MAE = ", round(mae_min, 2)), bty = "n")

# Lambda.1se
plot(y, predictions_1se,
     pch = 19, col = colors_by_study, cex = 1.2,
     xlab = "Actual Age", ylab = "Predicted Age",
     main = paste0("Lambda.1se (", n_features_1se, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = paste0("MAE = ", round(mae_1se, 2)), bty = "n")

par(mfrow = c(1, 1))


# ===== GET SELECTED FEATURES =====
cat("\n=== Selected CpG Sites (lambda.min) ===\n")

# Get non-zero coefficients
selected_coefs <- coef(lasso_model, s = "lambda.min")
selected_cpgs <- rownames(selected_coefs)[selected_coefs[,1] != 0]
selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]

cat("Top 20 selected CpGs:\n")
print(head(selected_cpgs, 20))

# Save selected features
selected_features <- data.frame(
  CpG = selected_cpgs,
  Coefficient = selected_coefs[selected_cpgs, 1]
)
selected_features <- selected_features[order(abs(selected_features$Coefficient), decreasing = TRUE), ]

cat("\nTop 10 by absolute coefficient:\n")
print(head(selected_features, 10))
write.csv(selected_features, "lasso_selected_features.csv", row.names = FALSE)

# ===== SAVE ALL RESULTS =====
results <- data.frame(
  actual_age = y,
  predicted_age_min = predictions_min,
  predicted_age_1se = predictions_1se,
  error_min = errors_min,
  error_1se = errors_1se,
  study = study_labels
)

write.csv(results, "lasso_predictions.csv", row.names = FALSE)
saveRDS(lasso_model, "lasso_model.RDS")


# ===== SUMMARY TABLE =====
performance_table <- data.frame(
  Model = c("Lambda.min", "Lambda.1se"),
  N_Features = c(n_features_min, n_features_1se),
  MAE = c(round(mae_min, 2), round(mae_1se, 2)),
  RMSE = c(round(rmse_min, 2), round(rmse_1se, 2)),
  Correlation = c(round(cor_min, 3), round(cor_1se, 3)),
  R2 = c(round(r2_min, 3), round(r2_1se, 3))
)

cat("\n========================================\n")
cat("LASSO MODEL COMPARISON\n")
cat("========================================\n")
print(performance_table)

write.csv(performance_table, "lasso_performance.csv", row.names = FALSE)
cat("\n========================================\n")
cat("FILES CREATED\n")
cat("========================================\n")
cat("1. lasso_predictions.csv - All predictions\n")
cat("2. lasso_selected_features.csv - Selected CpGs\n")
cat("3. lasso_performance.csv - Performance metrics\n")
cat("4. lasso_model.RDS - Trained model\n")
cat("5. VISAGE_by_study.png - Prediction plot\n")

# Create 4-panel figure with all glmnet diagnostic plots
png("lasso_diagnostic_plots.png", width = 1600, height = 1200, res = 120)

par(mfrow = c(2, 2))

# Plot 1: CV curve
plot(lasso_model, main = "Cross-Validation Curve")

# Plot 2: Coefficient paths
plot(lasso_model$glmnet.fit, xvar = "lambda", label = TRUE,
     main = "Coefficient Paths vs Lambda")
abline(v = log(lasso_model$lambda.min), col = "red", lty = 2)
abline(v = log(lasso_model$lambda.1se), col = "blue", lty = 2)

# Plot 3: Predictions (lambda.min)
plot(y, predictions_min,
     pch = 19, col = colors_by_study, cex = 1.2,
     xlab = "Actual Age", ylab = "Predicted Age",
     main = paste0("Predictions - lambda.min (", n_features_min, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = paste0("MAE = ", round(mae_min, 2)), bty = "n")

# Plot 4: Predictions (lambda.1se)
plot(y, predictions_1se,
     pch = 19, col = colors_by_study, cex = 1.2,
     xlab = "Actual Age", ylab = "Predicted Age",
     main = paste0("Predictions - lambda.1se (", n_features_1se, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = paste0("MAE = ", round(mae_1se, 2)), bty = "n")
par(mfrow = c(1, 1))
dev.off()

cat("✓ Diagnostic plots saved as 'lasso_diagnostic_plots.png'\n")

# ===== ALPHA PARAMETER =====
# alpha = 1  → Lasso (L1 penalty, feature selection)
# alpha = 0  → Ridge (L2 penalty, shrinkage)
# alpha = 0.5 → Elastic Net (mix of both)

# Lasso example
lasso <- cv.glmnet(x, y, alpha = 1)

# Ridge example
ridge <- cv.glmnet(x, y, alpha = 0)

# Elastic Net example
elastic <- cv.glmnet(x, y, alpha = 0.5)

# ===== FIT ALL THREE MODELS =====
cat("Fitting models...\n")

# Lasso (alpha = 1)
lasso_model <- cv.glmnet(x, y, alpha = 1, nfolds = 10)

# Ridge (alpha = 0)
ridge_model <- cv.glmnet(x, y, alpha = 0, nfolds = 10)

# Elastic Net (alpha = 0.5)
elastic_model <- cv.glmnet(x, y, alpha = 0.5, nfolds = 10)

cat("✓ All models fitted!\n\n")

# ===== PREDICTIONS =====
pred_lasso <- predict(lasso_model, newx = x, s = "lambda.min")[,1]
pred_ridge <- predict(ridge_model, newx = x, s = "lambda.min")[,1]
pred_elastic <- predict(elastic_model, newx = x, s = "lambda.min")[,1]


# ===== CALCULATE METRICS =====
calc_metrics <- function(pred, actual) {
  err <- pred - actual
  list(
    MAE = round(mean(abs(err)), 2),
    RMSE = round(sqrt(mean(err^2)), 2),
    Cor = round(cor(pred, actual), 3),
    R2 = round(summary(lm(pred ~ actual))$r.squared, 3)
  )
}
metrics_lasso <- calc_metrics(pred_lasso, y)
metrics_ridge <- calc_metrics(pred_ridge, y)
metrics_elastic <- calc_metrics(pred_elastic, y)

# Features selected
n_lasso <- sum(coef(lasso_model, s = "lambda.min") != 0) - 1
n_ridge <- sum(coef(ridge_model, s = "lambda.min") != 0) - 1
n_elastic <- sum(coef(elastic_model, s = "lambda.min") != 0) - 1

# ===== COMPARISON TABLE =====
comparison <- data.frame(
  Model = c("Lasso", "Ridge", "Elastic Net"),
  N_Features = c(n_lasso, n_ridge, n_elastic),
  MAE = c(metrics_lasso$MAE, metrics_ridge$MAE, metrics_elastic$MAE),
  RMSE = c(metrics_lasso$RMSE, metrics_ridge$RMSE, metrics_elastic$RMSE),
  Correlation = c(metrics_lasso$Cor, metrics_ridge$Cor, metrics_elastic$Cor),
  R2 = c(metrics_lasso$R2, metrics_ridge$R2, metrics_elastic$R2)
)

cat("========================================\n")
cat("MODEL COMPARISON\n")
cat("========================================\n")
print(comparison)

write.csv(comparison, "model_comparison.csv", row.names = FALSE)

# ===== GET STUDY LABELS =====
study_labels <- read_csv("Training data set.csv")$GSE[complete_idx]

study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

colors_by_study <- study_colors[study_labels]

# ===== CREATE 3-PANEL COMPARISON PLOT =====
png("lasso_ridge_elastic_comparison.png", width = 1800, height = 600, res = 120)

par(mfrow = c(1, 3))

plot_order <- c("GSE185920", "GSE185445", "GSE149318")

# Panel 1: Lasso
plot(y, pred_lasso, type = "n",
     xlab = "Actual Age", ylab = "Predicted Age",
     main = paste0("Lasso (", n_lasso, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
for(s in plot_order) {
  idx <- study_labels == s
  size <- ifelse(s == "GSE185920", 1.2, 2)
  points(y[idx], pred_lasso[idx], pch = 19, 
         col = study_colors[s], cex = size)
  
  # Panel 2: Ridge
  plot(y, pred_ridge, type = "n",
       xlab = "Actual Age", ylab = "Predicted Age",
       main = paste0("Ridge (", n_ridge, " CpGs)"))
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  for(s in plot_order) {
    idx <- study_labels == s
    size <- ifelse(s == "GSE185920", 1.2, 2)
    points(y[idx], pred_ridge[idx], pch = 19, 
           col = study_colors[s], cex = size)
  }
  legend("topleft", 
         legend = paste0("MAE = ", metrics_ridge$MAE, " yrs"), 
         bty = "n")
  # Panel 3: Elastic Net
  plot(y, pred_elastic, type = "n",
       xlab = "Actual Age", ylab = "Predicted Age",
       main = paste0("Elastic Net (", n_elastic, " CpGs)"))
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  for(s in plot_order) {
    idx <- study_labels == s
    size <- ifelse(s == "GSE185920", 1.2, 2)
    points(y[idx], pred_elastic[idx], pch = 19, 
           col = study_colors[s], cex = size)
  }
  legend("topleft", 
         legend = paste0("MAE = ", metrics_elastic$MAE, " yrs"), 
         bty = "n")
  
  # Add study legend to last panel
  legend("bottomright",
         legend = names(study_colors),
         col = study_colors,
         pch = 19,
         title = "Study",
         cex = 0.9)
  
  par(mfrow = c(1, 1))
  dev.off()
  
  cat("✓ Comparison plot saved!\n")
  
}
legend("topleft", 
       legend = paste0("MAE = ", metrics_lasso$MAE, " yrs"), 
       bty = "n")

#comparing alpha values

library(glmnet)
library(readr)

# ===== LOAD AND PREPARE DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_Top_50K.RDS")
train_age_data <- read_csv("Training data set.csv")

ages_all <- train_age_data$Age[1:ncol(Train)]
studies_all <- train_age_data$GSE[1:ncol(Train)]

x_all <- t(Train)
y_all <- ages_all
# Remove NAs
complete_idx <- complete.cases(x_all, y_all)
x <- x_all[complete_idx, ]
y <- y_all[complete_idx]
study_labels <- studies_all[complete_idx]

cat("Clean data:\n")
cat("  Samples:", nrow(x), "\n")
cat("  Features:", ncol(x), "\n\n")

# ===== FIT THREE MODELS =====
cat("Fitting models...\n")

# Ridge (alpha = 0)
ridge_model <- cv.glmnet(x, y, alpha = 0, nfolds = 10)
cat("  ✓ Ridge\n")

# Elastic Net (alpha = 0.5)
elastic_model <- cv.glmnet(x, y, alpha = 0.5, nfolds = 10)
cat("  ✓ Elastic Net\n")

# Lasso (alpha = 1)
lasso_model <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
cat("  ✓ Lasso\n\n")

# ===== PREDICTIONS =====
pred_ridge <- predict(ridge_model, newx = x, s = "lambda.min")[,1]
pred_elastic <- predict(elastic_model, newx = x, s = "lambda.min")[,1]
pred_lasso <- predict(lasso_model, newx = x, s = "lambda.min")[,1]

# ===== CALCULATE METRICS =====
calc_metrics <- function(pred, actual) {
  err <- pred - actual
  list(
    MAE = round(mean(abs(err)), 2),
    RMSE = round(sqrt(mean(err^2)), 2),
    Cor = round(cor(pred, actual), 3),
    R2 = round(summary(lm(pred ~ actual))$r.squared, 3),
    N_features = NA  # Will add separately
  )
}

metrics_ridge <- calc_metrics(pred_ridge, y)
metrics_elastic <- calc_metrics(pred_elastic, y)
metrics_lasso <- calc_metrics(pred_lasso, y)

# Add feature counts
metrics_ridge$N_features <- sum(coef(ridge_model, s = "lambda.min") != 0) - 1
metrics_elastic$N_features <- sum(coef(elastic_model, s = "lambda.min") != 0) - 1
metrics_lasso$N_features <- sum(coef(lasso_model, s = "lambda.min") != 0) - 1

# ===== SUMMARY TABLE =====
summary_table <- data.frame(
  Model = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"),
  Alpha = c(0, 0.5, 1),
  N_Features = c(metrics_ridge$N_features, metrics_elastic$N_features, metrics_lasso$N_features),
  MAE = c(metrics_ridge$MAE, metrics_elastic$MAE, metrics_lasso$MAE),
  RMSE = c(metrics_ridge$RMSE, metrics_elastic$RMSE, metrics_lasso$RMSE),
  Correlation = c(metrics_ridge$Cor, metrics_elastic$Cor, metrics_lasso$Cor),
  R2 = c(metrics_ridge$R2, metrics_elastic$R2, metrics_lasso$R2)
)

cat("========================================\n")
cat("MODEL COMPARISON\n")
cat("========================================\n")
print(summary_table)

write.csv(summary_table, "alpha_comparison_table.csv", row.names = FALSE)

# ===== DEFINE STUDY COLORS =====
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

plot_order <- c("GSE185920", "GSE185445", "GSE149318")

# ===== CREATE 3-PANEL FIGURE =====
png("alpha_comparison_3panel.png", width = 1800, height = 600, res = 120)

par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

# ===== PANEL 1: RIDGE (α = 0) =====
plot(y, pred_ridge,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Ridge (α=0)\n", metrics_ridge$N_features, " features"))

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(pred_ridge ~ y), col = "darkblue", lwd = 2)

for(study_name in plot_order) {
  idx <- study_labels == study_name
  if(study_name == "GSE185920") {
    points(y[idx], pred_ridge[idx],
           pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
  } else {
    points(y[idx], pred_ridge[idx],
           pch = 19, col = study_colors[study_name], cex = 1.8)
  }
}

legend("topleft",
       legend = c(
         paste0("r = ", metrics_ridge$Cor),
         paste0("MAE = ", metrics_ridge$MAE, " yrs"),
         paste0("R² = ", metrics_ridge$R2)
       ),
       bty = "n", cex = 0.9)

# ===== PANEL 2: ELASTIC NET (α = 0.5) =====
plot(y, pred_elastic,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Elastic Net (α=0.5)\n", metrics_elastic$N_features, " features"))

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(pred_elastic ~ y), col = "purple", lwd = 2)

for(study_name in plot_order) {
  idx <- study_labels == study_name
  if(study_name == "GSE185920") {
    points(y[idx], pred_elastic[idx],
           pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
  } else {
    points(y[idx], pred_elastic[idx],
           pch = 19, col = study_colors[study_name], cex = 1.8)
  }
}
legend("topleft",
       legend = c(
         paste0("r = ", metrics_elastic$Cor),
         paste0("MAE = ", metrics_elastic$MAE, " yrs"),
         paste0("R² = ", metrics_elastic$R2)
       ),
       bty = "n", cex = 0.9)

# ===== PANEL 3: LASSO (α = 1) =====
plot(y, pred_lasso,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Lasso (α=1)\n", metrics_lasso$N_features, " features"))

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(pred_lasso ~ y), col = "darkgreen", lwd = 2)

for(study_name in plot_order) {
  idx <- study_labels == study_name
  if(study_name == "GSE185920") {
    points(y[idx], pred_lasso[idx],
           pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
  } else {
    points(y[idx], pred_lasso[idx],
           pch = 19, col = study_colors[study_name], cex = 1.8)
  }
}
legend("topleft",
       legend = c(
         paste0("r = ", metrics_lasso$Cor),
         paste0("MAE = ", metrics_lasso$MAE, " yrs"),
         paste0("R² = ", metrics_lasso$R2)
       ),
       bty = "n", cex = 0.9)

# Add study legend to last panel
legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.3,
       title = "Study",
       cex = 0.8)


par(mfrow = c(1, 1))
dev.off()

cat("\n✓ 3-panel figure saved as 'alpha_comparison_3panel.png'\n")

# ===== SAVE MODELS =====
saveRDS(ridge_model, "ridge_model.RDS")
saveRDS(elastic_model, "elastic_model.RDS")
saveRDS(lasso_model, "lasso_model.RDS")

cat("✓ All models saved!\n")


# ===== BEST MODEL =====
cat("\n========================================\n")
cat("BEST MODEL\n")
cat("========================================\n")

best_idx <- which.min(summary_table$MAE)
cat("Best MAE:", summary_table$Model[best_idx], "\n")
cat("  MAE:", summary_table$MAE[best_idx], "years\n")
cat("  R²:", summary_table$R2[best_idx], "\n")
cat("  Features:", summary_table$N_Features[best_idx], "\n")

#TRAIN VS TEST 1se

library(glmnet)
library(readr)

# ===== LOAD DATA =====
Train <- readRDS("~/Documents/Forensic Project/Jamie/Final_Training_Matrix_Top_50K.RDS")
Test <- readRDS("~/Documents/Forensic Project/Jamie/Final_Test_Matrix_Top_50K.RDS")

train_age_data <- read_csv("Training data set.csv")
test_age_data <- read_csv("Test data set.csv")

# Get ages and studies
train_ages_all <- train_age_data$Age[1:ncol(Train)]
test_ages_all <- test_age_data$Age[1:ncol(Test)]

train_studies_all <- train_age_data$GSE[1:ncol(Train)]
test_studies_all <- test_age_data$GSE[1:ncol(Test)]
# ===== PREPARE TRAINING DATA =====
x_train_all <- t(Train)
y_train_all <- train_ages_all

# Remove NAs from training
complete_train <- complete.cases(x_train_all, y_train_all)

x_train <- x_train_all[complete_train, ]
y_train <- y_train_all[complete_train]
study_train <- train_studies_all[complete_train]

cat("Training data:\n")
cat("  Samples:", nrow(x_train), "\n")
cat("  Complete:", sum(complete_train), "\n\n")

# ===== PREPARE TEST DATA =====
x_test_all <- t(Test)
y_test_all <- test_ages_all

# Remove NAs from test
complete_test <- complete.cases(x_test_all, y_test_all)

x_test <- x_test_all[complete_test, ]
y_test <- y_test_all[complete_test]
study_test <- test_studies_all[complete_test]

cat("Test data:\n")
cat("  Samples:", nrow(x_test), "\n")
cat("  Complete:", sum(complete_test), "\n\n")

# ===== FIT MODELS ON TRAINING DATA =====
cat("Fitting models on training data...\n")

ridge_model <- cv.glmnet(x_train, y_train, alpha = 0, nfolds = 10)
elastic_model <- cv.glmnet(x_train, y_train, alpha = 0.5, nfolds = 10)
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)

cat("✓ All models fitted!\n\n")

# ===== PREDICTIONS (LAMBDA.1SE) =====
# Training set predictions
pred_train_ridge <- predict(ridge_model, newx = x_train, s = "lambda.1se")[,1]
pred_train_elastic <- predict(elastic_model, newx = x_train, s = "lambda.1se")[,1]
pred_train_lasso <- predict(lasso_model, newx = x_train, s = "lambda.1se")[,1]

# Test set predictions
pred_test_ridge <- predict(ridge_model, newx = x_test, s = "lambda.1se")[,1]
pred_test_elastic <- predict(elastic_model, newx = x_test, s = "lambda.1se")[,1]
pred_test_lasso <- predict(lasso_model, newx = x_test, s = "lambda.1se")[,1]

# ===== CALCULATE METRICS =====
calc_metrics <- function(pred, actual, name) {
  err <- pred - actual
  mae <- round(mean(abs(err)), 2)
  rmse <- round(sqrt(mean(err^2)), 2)
  cor_val <- round(cor(pred, actual), 3)
  r2 <- round(summary(lm(pred ~ actual))$r.squared, 3)
  
  cat(name, ":\n")
  cat("  MAE:", mae, "years\n")
  cat("  RMSE:", rmse, "years\n")
  cat("  r:", cor_val, "\n")
  cat("  R²:", r2, "\n\n")
  
  list(MAE = mae, RMSE = rmse, Cor = cor_val, R2 = r2)
}

cat("=== RIDGE (α=0) ===\n")
train_ridge_metrics <- calc_metrics(pred_train_ridge, y_train, "Training")
test_ridge_metrics <- calc_metrics(pred_test_ridge, y_test, "Test")

cat("=== ELASTIC NET (α=0.5) ===\n")
train_elastic_metrics <- calc_metrics(pred_train_elastic, y_train, "Training")
test_elastic_metrics <- calc_metrics(pred_test_elastic, y_test, "Test")

cat("=== LASSO (α=1) ===\n")
train_lasso_metrics <- calc_metrics(pred_train_lasso, y_train, "Training")
test_lasso_metrics <- calc_metrics(pred_test_lasso, y_test, "Test")

# ===== GET FEATURE COUNTS =====
n_ridge <- sum(coef(ridge_model, s = "lambda.1se") != 0) - 1
n_elastic <- sum(coef(elastic_model, s = "lambda.1se") != 0) - 1
n_lasso <- sum(coef(lasso_model, s = "lambda.1se") != 0) - 1


# ===== STUDY COLORS =====
study_colors <- c(
  "GSE185920" = "skyblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

plot_order <- c("GSE185920", "GSE185445", "GSE149318")

# ===== PLOTTING FUNCTION =====
plot_predictions <- function(y_data, pred, study_labels, title_text, metrics, study_colors) {
  
  plot(y_data, pred,
       type = "n",
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = title_text,
       xlim = range(c(y_data, pred)),
       ylim = range(c(y_data, pred)))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(pred ~ y_data), col = "black", lwd = 2)
  for(study_name in plot_order) {
    idx <- study_labels == study_name
    if(sum(idx) > 0) {
      if(study_name == "GSE185920") {
        points(y_data[idx], pred[idx],
               pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
      } else {
        points(y_data[idx], pred[idx],
               pch = 19, col = study_colors[study_name], cex = 1.8)
      }
    }
  }
  legend("topleft",
         legend = c(
           paste0("r = ", metrics$Cor),
           paste0("MAE = ", metrics$MAE, " yrs"),
           paste0("R² = ", metrics$R2)
         ),
         bty = "n", cex = 0.9)
  
  legend("bottomright",
         legend = names(study_colors),
         col = study_colors,
         pch = 19,
         pt.cex = 1.2,
         title = "Study",
         cex = 0.8)
}

# ===== FIGURE 1: RIDGE (α = 0) =====
png("Ridge_alpha0_train_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

plot_predictions(y_train, pred_train_ridge, study_train,
                 paste0("Ridge (α=0) - Training Set\n", n_ridge, " features"),
                 train_ridge_metrics, study_colors)

plot_predictions(y_test, pred_test_ridge, study_test,
                 paste0("Ridge (α=0) - Test Set\n", n_ridge, " features"),
                 test_ridge_metrics, study_colors)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Figure 1 saved: Ridge_alpha0_train_test.png\n")

# ===== FIGURE 2: ELASTIC NET (α = 0.5) =====
png("ElasticNet_alpha0.5_train_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

plot_predictions(y_train, pred_train_elastic, study_train,
                 paste0("Elastic Net (α=0.5) - Training Set\n", n_elastic, " features"),
                 train_elastic_metrics, study_colors)

plot_predictions(y_test, pred_test_elastic, study_test,
                 paste0("Elastic Net (α=0.5) - Test Set\n", n_elastic, " features"),
                 test_elastic_metrics, study_colors)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Figure 2 saved: ElasticNet_alpha0.5_train_test.png\n")

# ===== FIGURE 3: LASSO (α = 1) =====
png("Lasso_alpha1_train_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

plot_predictions(y_train, pred_train_lasso, study_train,
                 paste0("Lasso (α=1) - Training Set\n", n_lasso, " features"),
                 train_lasso_metrics, study_colors)

plot_predictions(y_test, pred_test_lasso, study_test,
                 paste0("Lasso (α=1) - Test Set\n", n_lasso, " features"),
                 test_lasso_metrics, study_colors)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Figure 3 saved: Lasso_alpha1_train_test.png\n")

# ===== SUMMARY TABLE (TRAIN VS TEST) =====
summary_table <- data.frame(
  Model = rep(c("Ridge", "Elastic Net", "Lasso"), each = 2),
  Alpha = rep(c(0, 0.5, 1), each = 2),
  Dataset = rep(c("Training", "Test"), 3),
  N_Features = rep(c(n_ridge, n_elastic, n_lasso), each = 2),
  MAE = c(train_ridge_metrics$MAE, test_ridge_metrics$MAE,
          train_elastic_metrics$MAE, test_elastic_metrics$MAE,
          train_lasso_metrics$MAE, test_lasso_metrics$MAE),
  RMSE = c(train_ridge_metrics$RMSE, test_ridge_metrics$RMSE,
           train_elastic_metrics$RMSE, test_elastic_metrics$RMSE,
           train_lasso_metrics$RMSE, test_lasso_metrics$RMSE),
  Correlation = c(train_ridge_metrics$Cor, test_ridge_metrics$Cor,
                  train_elastic_metrics$Cor, test_elastic_metrics$Cor,
                  train_lasso_metrics$Cor, test_lasso_metrics$Cor),
  R2 = c(train_ridge_metrics$R2, test_ridge_metrics$R2,
         train_elastic_metrics$R2, test_elastic_metrics$R2,
         train_lasso_metrics$R2, test_lasso_metrics$R2)
)
cat("\n========================================\n")
cat("SUMMARY: TRAINING vs TEST PERFORMANCE\n")
cat("========================================\n")
print(summary_table)

write.csv(summary_table, "train_test_comparison.csv", row.names = FALSE)

# ===== BEST MODEL =====
cat("\n=== Best Test Set Performance ===\n")
test_only <- summary_table[summary_table$Dataset == "Test", ]
best_idx <- which.min(test_only$MAE)

cat("Best model:", test_only$Model[best_idx], "(α =", test_only$Alpha[best_idx], ")\n")
cat("  MAE:", test_only$MAE[best_idx], "years\n")
cat("  R²:", test_only$R2[best_idx], "\n")
cat("  Features:", test_only$N_Features[best_idx], "\n")


cat("\n========================================\n")
cat("FILES CREATED\n")
cat("========================================\n")
cat("1. Ridge_alpha0_train_test.png\n")
cat("2. ElasticNet_alpha0.5_train_test.png\n")
cat("3. Lasso_alpha1_train_test.png\n")
cat("4. train_test_comparison.csv\n")
cat("5. ridge_model.RDS\n")
cat("6. elastic_model.RDS\n")
cat("7. lasso_model.RDS\n")

# ===== COMBINE INTO MASTER TABLE =====

# ===== FUNCTION TO CALCULATE METRICS =====
calculate_all_metrics <- function(model, x_train, y_train, x_test, y_test, alpha_val) {
  
  # Lambda.min results
  pred_train_min <- predict(model, newx = x_train, s = "lambda.min")[,1]
  pred_test_min <- predict(model, newx = x_test, s = "lambda.min")[,1]
  
  train_mae_min <- round(mean(abs(pred_train_min - y_train)), 2)
  train_rmse_min <- round(sqrt(mean((pred_train_min - y_train)^2)), 2)
  train_r2_min <- round(summary(lm(pred_train_min ~ y_train))$r.squared, 3)
  
  test_mae_min <- round(mean(abs(pred_test_min - y_test)), 2)
  test_rmse_min <- round(sqrt(mean((pred_test_min - y_test)^2)), 2)
  test_r2_min <- round(summary(lm(pred_test_min ~ y_test))$r.squared, 3)
  
  n_features_min <- sum(coef(model, s = "lambda.min") != 0) - 1
  
  # Lambda.1se results
  pred_train_1se <- predict(model, newx = x_train, s = "lambda.1se")[,1]
  pred_test_1se <- predict(model, newx = x_test, s = "lambda.1se")[,1]
  train_mae_1se <- round(mean(abs(pred_train_1se - y_train)), 2)
  train_rmse_1se <- round(sqrt(mean((pred_train_1se - y_train)^2)), 2)
  train_r2_1se <- round(summary(lm(pred_train_1se ~ y_train))$r.squared, 3)
  
  test_mae_1se <- round(mean(abs(pred_test_1se - y_test)), 2)
  test_rmse_1se <- round(sqrt(mean((pred_test_1se - y_test)^2)), 2)
  test_r2_1se <- round(summary(lm(pred_test_1se ~ y_test))$r.squared, 3)
  
  n_features_1se <- sum(coef(model, s = "lambda.1se") != 0) - 1
  
  # Return as data frame (2 rows: lambda.min and lambda.1se)
  data.frame(
    Alpha = rep(alpha_val, 2),
    Lambda = c("lambda.min", "lambda.1se"),
    N_Features = c(n_features_min, n_features_1se),
    Train_MAE = c(train_mae_min, train_mae_1se),
    Train_RMSE = c(train_rmse_min, train_rmse_1se),
    Train_R2 = c(train_r2_min, train_r2_1se),
    Test_MAE = c(test_mae_min, test_mae_1se),
    Test_RMSE = c(test_rmse_min, test_rmse_1se),
    Test_R2 = c(test_r2_min, test_r2_1se)
  )
}

# ===== CALCULATE FOR ALL THREE ALPHAS =====
cat("Calculating metrics...\n")

results_ridge <- calculate_all_metrics(ridge_model, x_train, y_train, x_test, y_test, 0)
results_elastic <- calculate_all_metrics(elastic_model, x_train, y_train, x_test, y_test, 0.5)
results_lasso <- calculate_all_metrics(lasso_model, x_train, y_train, x_test, y_test, 1)

master_table <- rbind(results_ridge, results_elastic, results_lasso)

# Reorder columns for clarity
master_table <- master_table[, c("Alpha", "Lambda", "N_Features", 
                                 "Train_MAE", "Train_RMSE", "Train_R2",
                                 "Test_MAE", "Test_RMSE", "Test_R2")]

cat("\n")
cat("================================================================================\n")
cat("                    COMPREHENSIVE MODEL COMPARISON                             \n")
cat("================================================================================\n\n")

print(master_table, row.names = FALSE)

# ===== SAVE TABLE =====
write.csv(master_table, "comprehensive_model_comparison.csv", row.names = FALSE)
cat("\n✓ Table saved as 'comprehensive_model_comparison.csv'\n")


# ===== CALCULATE OVERFITTING =====
master_table$Overfit_MAE <- master_table$Test_MAE - master_table$Train_MAE
master_table$Overfit_R2 <- master_table$Train_R2 - master_table$Test_R2

cat("\n=== Overfitting Analysis ===\n")
print(master_table[, c("Alpha", "Lambda", "Train_MAE", "Test_MAE", "Overfit_MAE")])

# ===== IDENTIFY BEST MODELS =====
cat("\n========================================\n")
cat("BEST MODELS\n")
cat("========================================\n")

# Best test MAE overall
best_overall <- which.min(master_table$Test_MAE)
cat("\n1. BEST TEST MAE (Overall):\n")
cat("   Alpha:", master_table$Alpha[best_overall], "\n")
cat("   Lambda:", master_table$Lambda[best_overall], "\n")
cat("   Features:", master_table$N_Features[best_overall], "\n")
cat("   Test MAE:", master_table$Test_MAE[best_overall], "years\n")
# Best lambda.1se (recommended for generalization)
best_1se <- master_table[master_table$Lambda == "lambda.1se", ]
best_1se_idx <- which.min(best_1se$Test_MAE)
cat("\n2. BEST LAMBDA.1SE (Best Generalization):\n")
cat("   Alpha:", best_1se$Alpha[best_1se_idx], "\n")
cat("   Features:", best_1se$N_Features[best_1se_idx], "\n")
cat("   Test MAE:", best_1se$Test_MAE[best_1se_idx], "years\n")

# Best balance (accuracy + simplicity)
# Score: normalized MAE + normalized features
mae_norm <- (master_table$Test_MAE - min(master_table$Test_MAE)) / 
  (max(master_table$Test_MAE) - min(master_table$Test_MAE))
feat_norm <- (master_table$N_Features - min(master_table$N_Features)) / 
  (max(master_table$N_Features) - min(master_table$N_Features))
balance_score <- 0.6 * mae_norm + 0.4 * feat_norm
# Best lambda.1se (recommended for generalization)
best_1se <- master_table[master_table$Lambda == "lambda.1se", ]
best_1se_idx <- which.min(best_1se$Test_MAE)
cat("\n2. BEST LAMBDA.1SE (Best Generalization):\n")
cat("   Alpha:", best_1se$Alpha[best_1se_idx], "\n")
cat("   Features:", best_1se$N_Features[best_1se_idx], "\n")
cat("   Test MAE:", best_1se$Test_MAE[best_1se_idx], "years\n")

# Best balance (accuracy + simplicity)
# Score: normalized MAE + normalized features
mae_norm <- (master_table$Test_MAE - min(master_table$Test_MAE)) / 
  (max(master_table$Test_MAE) - min(master_table$Test_MAE))
feat_norm <- (master_table$N_Features - min(master_table$N_Features)) / 
  (max(master_table$N_Features) - min(master_table$N_Features))
balance_score <- 0.6 * mae_norm + 0.4 * feat_norm


# ===== CREATE VISUALIZATION OF TABLE =====
png("model_comparison_table_viz.png", width = 1400, height = 1000, res = 120)

par(mfrow = c(2, 2))

# Plot 1: Test MAE comparison
mae_data <- master_table[, c("Alpha", "Lambda", "Test_MAE")]
mae_matrix <- matrix(mae_data$Test_MAE, nrow = 3, ncol = 2, byrow = FALSE)
colnames(mae_matrix) <- c("lambda.min", "lambda.1se")
rownames(mae_matrix) <- c("α=0 (Ridge)", "α=0.5 (Elastic)", "α=1 (Lasso)")

barplot(mae_matrix, beside = TRUE, 
        col = c("steelblue", "purple", "darkgreen"),
        main = "Test Set MAE",
        ylab = "MAE (years)",
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8))
# Plot 2: Number of features
feat_matrix <- matrix(master_table$N_Features, nrow = 3, ncol = 2, byrow = FALSE)
colnames(feat_matrix) <- c("lambda.min", "lambda.1se")
rownames(feat_matrix) <- c("α=0", "α=0.5", "α=1")

barplot(feat_matrix, beside = TRUE,
        col = c("steelblue", "purple", "darkgreen"),
        main = "Number of Features",
        ylab = "N Features",
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8))
# Plot 3: Test R²
r2_matrix <- matrix(master_table$Test_R2, nrow = 3, ncol = 2, byrow = FALSE)
colnames(r2_matrix) <- c("lambda.min", "lambda.1se")
rownames(r2_matrix) <- c("α=0", "α=0.5", "α=1")

barplot(r2_matrix, beside = TRUE,
        col = c("steelblue", "purple", "darkgreen"),
        main = "Test Set R²",
        ylab = "R²",
        ylim = c(0, 1),
        legend.text = TRUE,
        args.legend = list(x = "bottomright", cex = 0.8))
# Plot 4: Overfitting (Train - Test MAE)
overfit_matrix <- matrix(master_table$Overfit_MAE, nrow = 3, ncol = 2, byrow = FALSE)
colnames(overfit_matrix) <- c("lambda.min", "lambda.1se")
rownames(overfit_matrix) <- c("α=0", "α=0.5", "α=1")

barplot(overfit_matrix, beside = TRUE,
        col = c("steelblue", "purple", "darkgreen"),
        main = "Overfitting (Test MAE - Train MAE)",
        ylab = "MAE Difference (years)",
        legend.text = TRUE,
        args.legend = list(x = "topright", cex = 0.8))
abline(h = 0, col = "red", lty = 2)

par(mfrow = c(1, 1))
dev.off()

cat("\n✓ Visualization saved as 'model_comparison_table_viz.png'\n")

# ===== PRETTY PRINT TABLE =====
cat("\n\n")
cat("================================================================================\n")
cat("                         MODEL PERFORMANCE SUMMARY                              \n")
cat("================================================================================\n\n")

cat(sprintf("%-6s %-12s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n",
            "Alpha", "Lambda", "N_Feat", "Tr_MAE", "Tr_R2", "Te_MAE", "Te_R2", "Overfit", ""))
cat(strrep("-", 100), "\n")

for(i in 1:nrow(master_table)) {
  cat(sprintf("%-6.1f %-12s %-10d %-10.2f %-10.3f %-10.2f %-10.3f %-10.2f\n",
              master_table$Alpha[i],
              master_table$Lambda[i],
              master_table$N_Features[i],
              master_table$Train_MAE[i],
              master_table$Train_R2[i],
              master_table$Test_MAE[i],
              master_table$Test_R2[i],
              master_table$Overfit_MAE[i]))
}
cat("\n")

# ===== RECOMMENDATIONS BASED ON RESULTS =====
cat("========================================\n")
cat("RECOMMENDATIONS\n")
cat("========================================\n\n")

# Best accuracy
best_test_mae <- which.min(master_table$Test_MAE)
cat("🎯 BEST ACCURACY:\n")
cat("   Model: Alpha =", master_table$Alpha[best_test_mae], 
    ", Lambda =", master_table$Lambda[best_test_mae], "\n")
cat("   Test MAE:", master_table$Test_MAE[best_test_mae], "years\n")
cat("   Features:", master_table$N_Features[best_test_mae], "\n\n")
# Best generalization (lowest overfitting)
best_generalize <- which.min(abs(master_table$Overfit_MAE))
cat("🔄 BEST GENERALIZATION (Lowest Overfitting):\n")
cat("   Model: Alpha =", master_table$Alpha[best_generalize], 
    ", Lambda =", master_table$Lambda[best_generalize], "\n")
cat("   Overfit:", master_table$Overfit_MAE[best_generalize], "years\n")
cat("   Test MAE:", master_table$Test_MAE[best_generalize], "years\n")
cat("   Features:", master_table$N_Features[best_generalize], "\n\n")
# Most practical (fewest features, acceptable accuracy)
lambda_1se_only <- master_table[master_table$Lambda == "lambda.1se", ]
best_practical <- which.min(lambda_1se_only$N_Features)
cat("🔧 MOST PRACTICAL (Fewest Features, lambda.1se):\n")
cat("   Model: Alpha =", lambda_1se_only$Alpha[best_practical], "\n")
cat("   Features:", lambda_1se_only$N_Features[best_practical], "\n")
cat("   Test MAE:", lambda_1se_only$Test_MAE[best_practical], "years\n")
cat("   Test R²:", lambda_1se_only$Test_R2[best_practical], "\n\n")
# Recommended
cat("⭐ RECOMMENDED FOR FORENSIC USE:\n")
lasso_1se <- master_table[master_table$Alpha == 1 & master_table$Lambda == "lambda.1se", ]
cat("   Model: Lasso (α=1), lambda.1se\n")
cat("   Features:", lasso_1se$N_Features, "\n")
cat("   Test MAE:", lasso_1se$Test_MAE, "years\n")
cat("   Test R²:", lasso_1se$Test_R2, "\n")
cat("   Reason: Best balance of accuracy, simplicity, and generalization\n")

# ===== EXPORT TO FORMATTED TABLE =====
# Add descriptive model names
master_table$Model_Name <- paste0(
  ifelse(master_table$Alpha == 0, "Ridge",
         ifelse(master_table$Alpha == 0.5, "Elastic Net", "Lasso")),
  " (", master_table$Lambda, ")"
)

# Reorder columns
final_table <- master_table[, c("Model_Name", "Alpha", "Lambda", "N_Features",
                                "Train_MAE", "Train_RMSE", "Train_R2",
                                "Test_MAE", "Test_RMSE", "Test_R2",
                                "Overfit_MAE")]

write.csv(final_table, "complete_model_comparison.csv", row.names = FALSE)

# ===== CREATE PUBLICATION-READY TABLE =====
# Simplified version for presentation
pub_table <- data.frame(
  Model = c("Ridge", "Ridge", "Elastic Net", "Elastic Net", "Lasso", "Lasso"),
  Alpha = rep(c(0, 0.5, 1), each = 2),
  Lambda = rep(c("min", "1se"), 3),
  Features = master_table$N_Features,
  `Training MAE` = master_table$Train_MAE,
  `Test MAE` = master_table$Test_MAE,
  `Test R²` = master_table$Test_R2,
  check.names = FALSE
)

cat("\n\n")
cat("================================================================================\n")
cat("                    PUBLICATION-READY TABLE                                    \n")
cat("================================================================================\n\n")
print(pub_table, row.names = FALSE)

write.csv(pub_table, "publication_table.csv", row.names = FALSE)


# ===== SUMMARY STATISTICS =====
cat("\n\n")
cat("========================================\n")
cat("KEY FINDINGS\n")
cat("========================================\n\n")

cat("1. FEATURE SELECTION:\n")
cat("   Ridge uses ~", master_table$N_Features[master_table$Alpha == 0 & 
                                                 master_table$Lambda == "lambda.1se"], 
    " features\n", sep = "")
cat("   Elastic Net uses ~", master_table$N_Features[master_table$Alpha == 0.5 & 
                                                       master_table$Lambda == "lambda.1se"], 
    " features\n", sep = "")
cat("   Lasso uses ~", master_table$N_Features[master_table$Alpha == 1 & 
                                                 master_table$Lambda == "lambda.1se"], 
    " features\n\n", sep = "")

cat("2. ACCURACY:\n")
cat("   Best training MAE:", min(master_table$Train_MAE), "years\n")
cat("   Best test MAE:", min(master_table$Test_MAE), "years\n")
cat("   Range:", max(master_table$Test_MAE) - min(master_table$Test_MAE), "years difference\n\n")

cat("3. GENERALIZATION:\n")
avg_overfit_min <- mean(master_table$Overfit_MAE[master_table$Lambda == "lambda.min"])
avg_overfit_1se <- mean(master_table$Overfit_MAE[master_table$Lambda == "lambda.1se"])
cat("   Avg overfitting (lambda.min):", round(avg_overfit_min, 2), "years\n")
cat("   Avg overfitting (lambda.1se):", round(avg_overfit_1se, 2), "years\n")
cat("   → Lambda.1se generalizes better!\n\n")


# ==FIND OPTIMAL LAMBDA FOR <20 FEATURES==
  
  # ===== FIT LASSO MODEL =====
cat("Fitting Lasso model...\n")
lasso_model <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
cat("✓ Done!\n\n")


# ===== EXPLORE LAMBDA PATH =====
cat("=== Lambda Exploration ===\n")

# Get all lambdas tested
lambdas <- lasso_model$lambda

# Count features at each lambda
n_features <- sapply(lambdas, function(lambda) {
  sum(coef(lasso_model, s = lambda) != 0) - 1
})
# Create data frame
lambda_table <- data.frame(
  Lambda = round(lambdas, 4),
  N_Features = n_features
)

cat("Lambda range:", min(lambdas), "to", max(lambdas), "\n")
cat("Features range:", min(n_features), "to", max(n_features), "\n\n")

# ===== FIND LAMBDAS THAT GIVE <20 FEATURES =====
under_20 <- lambda_table[lambda_table$N_Features < 20, ]
cat("Lambdas with <20 features:\n")
print(head(under_20, 10))

# ===== TEST SPECIFIC LAMBDAS FOR <20 FEATURES =====
# Try a range of lambdas
test_lambdas <- c(0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0)

cat("\n=== Testing Specific Lambdas ===\n")

lambda_results <- data.frame(
  Lambda = numeric(),
  N_Features = numeric(),
  Train_MAE = numeric(),
  Train_R2 = numeric()
)
for(lambda_val in test_lambdas) {
  # Predict with this lambda
  pred <- predict(lasso_model, newx = x, s = lambda_val)[,1]
  
  # Count features
  n_feat <- sum(coef(lasso_model, s = lambda_val) != 0) - 1
  
  # Calculate metrics
  mae <- mean(abs(pred - y))
  r2 <- summary(lm(pred ~ y))$r.squared
  
  # Only keep if <20 features
  if(n_feat < 20) {
    lambda_results <- rbind(lambda_results,
                            data.frame(
                              Lambda = lambda_val,
                              N_Features = n_feat,
                              Train_MAE = round(mae, 2),
                              Train_R2 = round(r2, 3)
                            ))
  }
}
cat("\nResults for lambdas with <20 features:\n")
print(lambda_results)

# ===== FIND OPTIMAL LAMBDA =====
# Optimal = best MAE with <20 features
if(nrow(lambda_results) > 0) {
  
  optimal_idx <- which.min(lambda_results$Train_MAE)
  optimal_lambda <- lambda_results$Lambda[optimal_idx]
  optimal_features <- lambda_results$N_Features[optimal_idx]
  optimal_mae <- lambda_results$Train_MAE[optimal_idx]
  optimal_r2 <- lambda_results$Train_R2[optimal_idx]
  
  cat("\n========================================\n")
  cat("⭐ OPTIMAL MODEL (<20 features)\n")
  cat("========================================\n")
  cat("Lambda:", optimal_lambda, "\n")
  cat("Features:", optimal_features, "\n")
  cat("Training MAE:", optimal_mae, "years\n")
  cat("Training R²:", optimal_r2, "\n")
  
  # ===== GET SELECTED FEATURES =====
  optimal_coefs <- coef(lasso_model, s = optimal_lambda)
  selected_cpgs <- rownames(optimal_coefs)[optimal_coefs[,1] != 0]
  selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]
  
  cat("\n=== Selected CpGs (", length(selected_cpgs), ") ===\n", sep = "")
  print(selected_cpgs)
  
  # Get coefficients
  selected_features <- data.frame(
    CpG = selected_cpgs,
    Coefficient = optimal_coefs[selected_cpgs, 1]
  )
  selected_features <- selected_features[order(abs(selected_features$Coefficient), 
                                               decreasing = TRUE), ]
  cat("\nTop 10 features by importance:\n")
  print(head(selected_features, 10))
  
  
  # ===== PLOT LAMBDA vs FEATURES =====
  plot(log(lambdas), n_features,
       type = "l",
       lwd = 2,
       col = "steelblue",
       xlab = "log(Lambda)",
       ylab = "Number of Features",
       main = "Feature Selection Path")
  
  # Add horizontal line at 20
  abline(h = 20, col = "red", lwd = 2, lty = 2)
  
  # Mark optimal lambda
  abline(v = log(optimal_lambda), col = "darkgreen", lwd = 2, lty = 2)
  # Add lambda.min and lambda.1se
  abline(v = log(lasso_model$lambda.min), col = "blue", lwd = 1, lty = 3)
  abline(v = log(lasso_model$lambda.1se), col = "purple", lwd = 1, lty = 3)
  
  legend("topright",
         legend = c("< 20 features threshold", 
                    paste0("Optimal (λ=", optimal_lambda, ", n=", optimal_features, ")"),
                    "lambda.min",
                    "lambda.1se"),
         col = c("red", "darkgreen", "blue", "purple"),
         lwd = 2,
         lty = c(2, 2, 3, 3),
         cex = 0.8)
  
  
  # ===== PLOT OPTIMAL MODEL PREDICTIONS =====
  pred_optimal <- predict(lasso_model, newx = x, s = optimal_lambda)[,1]
  # Get study labels
  study_labels <- train_age_data$GSE[complete_idx]
  study_colors <- c("GSE185920" = "lightblue", 
                    "GSE185445" = "pink", 
                    "GSE149318" = "green")
  colors <- study_colors[study_labels]
  
  plot(y, pred_optimal,
       pch = 19,
       col = colors,
       cex = 1.5,
       xlab = "Chronological Age (years)",
       ylab = "Predicted Age (years)",
       main = paste0("Optimal Model\n(λ=", optimal_lambda, 
                     ", ", optimal_features, " features)"))
  
  abline(0, 1, col = "red", lwd = 2, lty = 2)
  abline(lm(pred_optimal ~ y), col = "black", lwd = 2)
  legend("topleft",
         legend = c(
           paste0("MAE = ", optimal_mae, " yrs"),
           paste0("R² = ", optimal_r2),
           paste0("n = ", optimal_features, " CpGs")
         ),
         bty = "n", cex = 1)
  
  legend("bottomright",
         legend = names(study_colors),
         col = study_colors,
         pch = 19,
         title = "Study",
         cex = 0.9)
  
  # ===== SAVE OPTIMAL MODEL RESULTS =====
  optimal_results <- data.frame(
    Lambda = optimal_lambda,
    N_Features = optimal_features,
    Train_MAE = optimal_mae,
    Train_R2 = optimal_r2,
    Selected_CpGs = paste(selected_cpgs, collapse = ", ")
  )
  
  write.csv(selected_features, "optimal_model_features.csv", row.names = FALSE)
  write.csv(optimal_results, "optimal_model_summary.csv", row.names = FALSE)
  
  cat("\n✓ Optimal model results saved!\n")
  
} else {
  cat("\n⚠️  No lambda values give <20 features\n")
  cat("Try increasing lambda manually\n")
}
  

# OPTIMAL CROSS VALIDATION 

# ===== FIT LASSO MODEL =====
cat("Fitting Lasso model...\n")
lasso_model <- cv.glmnet(x, y, alpha = 1, nfolds = 10)
cat("✓ Done!\n\n")


# ===== DEFINE YOUR OPTIMAL LAMBDA =====
# Replace these with YOUR actual values
optimal_lambda <- 0.7  # YOUR lambda
optimal_n_features <- 12  # YOUR feature count

# Verify
n_feat_verify <- sum(coef(lasso_model, s = optimal_lambda) != 0) - 1
cat("Your optimal lambda:", optimal_lambda, "\n")
cat("Features:", n_feat_verify, "\n\n")

# ===== CREATE CV PLOT =====
png("optimal_CV_plot.png", width = 1000, height = 800, res = 120)

# Standard CV plot from glmnet
plot(lasso_model,
     main = paste0("Lasso Cross-Validation\n",
                   "Optimal Model: λ = ", optimal_lambda, 
                   " (", optimal_n_features, " features)"),
     cex.main = 1.3)

# Add vertical line for YOUR optimal lambda (THICK RED)
abline(v = log(optimal_lambda), col = "red", lwd = 5, lty = 1)

# Add other reference lines
abline(v = log(lasso_model$lambda.min), col = "blue", lwd = 2, lty = 2)
abline(v = log(lasso_model$lambda.1se), col = "purple", lwd = 2, lty = 2)

# Get feature counts for comparison
n_min <- sum(coef(lasso_model, s = "lambda.min") != 0) - 1
n_1se <- sum(coef(lasso_model, s = "lambda.1se") != 0) - 1

# Legend
legend("topleft",
       legend = c(
         paste0("lambda.min (", n_min, " features)"),
         paste0("lambda.1se (", n_1se, " features)"),
         paste0("OPTIMAL: λ=", optimal_lambda, " (", optimal_n_features, " features)")
       ),
       col = c("blue", "purple", "red"),
       lwd = c(2, 2, 5),
       lty = c(2, 2, 1),
       cex = 1,
       bg = "white",
       box.lwd = 2)
# Add annotation pointing to optimal
lambda_idx <- which.min(abs(lasso_model$lambda - optimal_lambda))
optimal_mse <- lasso_model$cvm[lambda_idx]

text(log(optimal_lambda), optimal_mse,
     labels = paste0("  ← Optimal\n  ", optimal_n_features, " CpGs"),
     pos = 4,
     col = "red",
     cex = 1.1,
     font = 2)

dev.off()

cat("✓ CV plot saved as 'optimal_CV_plot.png'\n")

# ===== ALSO CREATE: NUMBER OF FEATURES PATH =====
png("lambda_features_path.png", width = 1000, height = 700, res = 120)

# Calculate features at each lambda
lambdas <- lasso_model$lambda
n_features <- sapply(lambdas, function(l) {
  sum(coef(lasso_model, s = l) != 0) - 1
})

# Plot
plot(log(lambdas), n_features,
     type = "l",
     lwd = 3,
     col = "steelblue",
     xlab = "log(Lambda)",
     ylab = "Number of Features",
     main = "Feature Selection Path")
# Add horizontal line at 20
abline(h = 20, col = "orange", lwd = 2, lty = 2)
text(max(log(lambdas)) * 0.7, 22, 
     "20 feature threshold", col = "orange", cex = 1)

# Mark YOUR optimal
abline(v = log(optimal_lambda), col = "red", lwd = 5)
abline(h = optimal_n_features, col = "red", lwd = 2, lty = 3)

points(log(optimal_lambda), optimal_n_features,
       pch = 19, col = "red", cex = 3)

text(log(optimal_lambda), optimal_n_features,
     labels = paste0("  Optimal\n  λ=", optimal_lambda, "\n  ", 
                     optimal_n_features, " features"),
     pos = 4,
     col = "red",
     font = 2,
     cex = 1.1)
# Other references
abline(v = log(lasso_model$lambda.min), col = "blue", lwd = 2, lty = 2)
abline(v = log(lasso_model$lambda.1se), col = "purple", lwd = 2, lty = 2)

legend("topright",
       legend = c("lambda.min", "lambda.1se", "Optimal", "20-feature limit"),
       col = c("blue", "purple", "red", "orange"),
       lwd = c(2, 2, 5, 2),
       lty = c(2, 2, 1, 2),
       cex = 0.9)

dev.off()

cat("✓ Feature path saved as 'lambda_features_path.png'\n")


#OPTIMAL VS TEST SET 

# ===== PREPARE TRAINING DATA =====
train_ages <- train_age_data$Age[1:ncol(Train)]
x_train_all <- t(Train)

complete_train <- complete.cases(x_train_all, train_ages)
x_train <- x_train_all[complete_train, ]
y_train <- train_ages[complete_train]


# ===== PREPARE TEST DATA =====
test_ages <- test_age_data$Age[1:ncol(Test)]
test_studies <- test_age_data$GSE[1:ncol(Test)]
x_test_all <- t(Test)

complete_test <- complete.cases(x_test_all, test_ages)
x_test <- x_test_all[complete_test, ]
y_test <- test_ages[complete_test]
study_test <- test_studies[complete_test]

# ===== PREPARE TRAINING DATA =====
train_ages <- train_age_data$Age[1:ncol(Train)]
x_train_all <- t(Train)

complete_train <- complete.cases(x_train_all, train_ages)
x_train <- x_train_all[complete_train, ]
y_train <- train_ages[complete_train]


# ===== PREPARE TEST DATA =====
test_ages <- test_age_data$Age[1:ncol(Test)]
test_studies <- test_age_data$GSE[1:ncol(Test)]
x_test_all <- t(Test)

complete_test <- complete.cases(x_test_all, test_ages)
x_test <- x_test_all[complete_test, ]
y_test <- test_ages[complete_test]
study_test <- test_studies[complete_test]

cat("Training samples:", nrow(x_train), "\n")
cat("Test samples:", nrow(x_test), "\n\n")


# ===== FIT LASSO MODEL ON TRAINING DATA =====
cat("Fitting Lasso model...\n")
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)
cat("✓ Done!\n\n")


# ===== YOUR OPTIMAL LAMBDA =====
optimal_lambda <- 0.7
optimal_n_features <- 12

cat("Optimal lambda:", optimal_lambda, "\n")
cat("Features:", optimal_n_features, "\n\n")

# ===== PREDICT ON TEST SET =====
pred_test <- predict(lasso_model, newx = x_test, s = optimal_lambda)[,1]


# ===== CALCULATE TEST METRICS =====
test_errors <- pred_test - y_test
test_mae <- mean(abs(test_errors))
test_rmse <- sqrt(mean(test_errors^2))
test_cor <- cor(pred_test, y_test)
test_r2 <- summary(lm(pred_test ~ y_test))$r.squared

cat("=== Test Set Performance ===\n")
cat("MAE:", round(test_mae, 2), "years\n")
cat("RMSE:", round(test_rmse, 2), "years\n")
cat("Correlation:", round(test_cor, 3), "\n")
cat("R²:", round(test_r2, 3), "\n\n")

# ===== STUDY COLORS =====
study_colors <- c(
  "GSE185920" = "lightblue",
  "GSE185445" = "pink",
  "GSE149318" = "green"
)

plot_order <- c("GSE185920", "GSE185445", "GSE149318")


# ===== CREATE TEST SET PLOT =====
png("optimal_model_test_set.png", width = 1000, height = 800, res = 120)

# Empty plot
plot(y_test, pred_test,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Optimal Lasso Model - Test Set Performance\n",
                   "λ = ", optimal_lambda, " (", optimal_n_features, " CpGs)"),
     cex.main = 1.2)

# Add reference lines
abline(0, 1, col = "red", lwd = 3, lty = 2)
abline(lm(pred_test ~ y_test), col = "black", lwd = 2.5)

# Plot by study (layered)
for(study_name in plot_order) {
  idx <- study_test == study_name
  if(sum(idx) > 0) {
    if(study_name == "GSE185920") {
      # Largest study: smaller, semi-transparent
      points(y_test[idx], pred_test[idx],
             pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1.2)
    } else {
      # Smaller studies: larger, opaque
      points(y_test[idx], pred_test[idx],
             pch = 19, col = study_colors[study_name], cex = 2)
    }
  }
}

# Metrics legend
legend("topleft",
       legend = c(
         paste0("r = ", round(test_cor, 3)),
         paste0("MAE = ", round(test_mae, 2), " years"),
         paste0("RMSE = ", round(test_rmse, 2), " years"),
         paste0("R² = ", round(test_r2, 3)),
         paste0("n = ", optimal_n_features, " CpGs")
       ),
       bty = "n",
       cex = 1.1)

# Metrics legend
legend("topleft",
       legend = c(
         paste0("r = ", round(test_cor, 3)),
         paste0("MAE = ", round(test_mae, 2), " years"),
         paste0("RMSE = ", round(test_rmse, 2), " years"),
         paste0("R² = ", round(test_r2, 3)),
         paste0("n = ", optimal_n_features, " CpGs")
       ),
       bty = "n",
       cex = 1.1)

# Study legend
legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.5,
       title = "Study",
       bg = "white")

# Lines legend
legend("bottom",
       legend = c("Perfect prediction", "Fitted line"),
       col = c("red", "black"),
       lwd = c(3, 2.5),
       lty = c(2, 1),
       bty = "n",
       cex = 0.9)

dev.off()
cat("\n✓ Test set plot saved as 'optimal_model_test_set.png'\n")

# ===== ALSO CREATE: TRAIN VS TEST COMPARISON =====
png("optimal_model_train_vs_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

# Get training predictions
pred_train <- predict(lasso_model, newx = x_train, s = optimal_lambda)[,1]
train_mae <- mean(abs(pred_train - y_train))
train_r2 <- summary(lm(pred_train ~ y_train))$r.squared

train_studies <- train_age_data$GSE[complete_train]
colors_train <- study_colors[train_studies]

# LEFT: Training Set
plot(y_train, pred_train,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Training Set\n(λ=", optimal_lambda, ", ", 
                   optimal_n_features, " CpGs)"))
# ===== ALSO CREATE: TRAIN VS TEST COMPARISON =====
png("optimal_model_train_vs_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

# Get training predictions
pred_train <- predict(lasso_model, newx = x_train, s = optimal_lambda)[,1]
train_mae <- mean(abs(pred_train - y_train))
train_r2 <- summary(lm(pred_train ~ y_train))$r.squared

train_studies <- train_age_data$GSE[complete_train]
colors_train <- study_colors[train_studies]

# LEFT: Training Set
plot(y_train, pred_train,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Training Set\n(λ=", optimal_lambda, ", ", 
                  optimal_n_features, " CpGs)"))



# ===== ALSO CREATE: TRAIN VS TEST COMPARISON =====
png("optimal_model_train_vs_test.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

# Get training predictions
pred_train <- predict(lasso_model, newx = x_train, s = optimal_lambda)[,1]
train_mae <- mean(abs(pred_train - y_train))
train_r2 <- summary(lm(pred_train ~ y_train))$r.squared

train_studies <- train_age_data$GSE[complete_train]
colors_train <- study_colors[train_studies]

# LEFT: Training Set
plot(y_train, pred_train,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Training Set\n(λ=", optimal_lambda, ", ", 
                  optimal_n_features, " CpGs)"))
abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(pred_train ~ y_train), col = "black", lwd = 2)

for(study_name in plot_order) {
  idx <- train_studies == study_name
  if(sum(idx) > 0) {
    if(study_name == "GSE185920") {
      points(y_train[idx], pred_train[idx],
             pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
    } else {
      points(y_train[idx], pred_train[idx],
             pch = 19, col = study_colors[study_name], cex = 1.8)
    }
  }
}

legend("topleft",
       legend = c(
         paste0("MAE = ", round(train_mae, 2), " yrs"),
         paste0("R² = ", round(train_r2, 3))
       ),
       bty = "n", cex = 1.1)


# RIGHT: Test Set
plot(y_test, pred_test,
     type = "n",
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Test Set\n(λ=", optimal_lambda, ", ", 
                   optimal_n_features, " CpGs)"))

abline(0, 1, col = "red", lwd = 2, lty = 2)
abline(lm(pred_test ~ y_test), col = "black", lwd = 2)

for(study_name in plot_order) {
  idx <- study_test == study_name
  if(sum(idx) > 0) {
    if(study_name == "GSE185920") {
      points(y_test[idx], pred_test[idx],
             pch = 19, col = rgb(0.68, 0.85, 0.90, 0.5), cex = 1)
    } else {
      points(y_test[idx], pred_test[idx],
             pch = 19, col = study_colors[study_name], cex = 1.8)
    }
  }
}

legend("topleft",
       legend = c(
         paste0("MAE = ", round(test_mae, 2), " yrs"),
         paste0("R² = ", round(test_r2, 3))
       ),
       bty = "n", cex = 1.1)

legend("bottomright",
       legend = names(study_colors),
       col = study_colors,
       pch = 19,
       pt.cex = 1.3,
       title = "Study",
       cex = 0.9)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Train vs Test comparison saved as 'optimal_model_train_vs_test.png'\n")

# ===== RESIDUAL ANALYSIS FOR TEST SET =====
png("optimal_model_test_residuals.png", width = 1200, height = 400, res = 120)

par(mfrow = c(1, 3))

# Residuals vs Predicted
plot(pred_test, test_errors,
     pch = 19, col = study_colors[study_test], cex = 1.3,
     xlab = "Predicted Age", ylab = "Residual (Predicted - Actual)",
     main = "Residuals vs Predicted")
abline(h = 0, col = "red", lwd = 2, lty = 2)

# Histogram of errors
hist(test_errors,
     breaks = 20,
     col = "lightblue",
     border = "white",
     main = "Distribution of Errors",
     xlab = "Error (years)")
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = mean(test_errors), col = "blue", lwd = 2, lty = 1)

# Absolute error by age
plot(y_test, abs(test_errors),
     pch = 19, col = study_colors[study_test], cex = 1.3,
     xlab = "Actual Age", ylab = "Absolute Error",
     main = "Error by Age")
abline(h = test_mae, col = "red", lwd = 2, lty = 2)
text(min(y_test) + 5, test_mae + 0.5,
     labels = paste0("Mean = ", round(test_mae, 2)),
     col = "red", pos = 3)

par(mfrow = c(1, 1))
dev.off()

cat("✓ Residual analysis saved as 'optimal_model_test_residuals.png'\n")


# ===== SAVE TEST SET RESULTS =====
test_results <- data.frame(
  sample_id = rownames(x_test),
  actual_age = y_test,
  predicted_age = pred_test,
  error = test_errors,
  abs_error = abs(test_errors),
  study = study_test
)

write.csv(test_results, "optimal_model_test_results.csv", row.names = FALSE)

cat("✓ Test results saved as 'optimal_model_test_results.csv'\n")



# model vs GSE185920

# ===== PREPARE TRAINING DATA =====
train_ages_all <- train_age_data$Age[1:ncol(Train)]
train_studies_all <- train_age_data$GSE[1:ncol(Train)]

x_train_all <- t(Train)

# Filter to complete cases
complete_train <- complete.cases(x_train_all, train_ages_all)
x_train_complete <- x_train_all[complete_train, ]
y_train_complete <- train_ages_all[complete_train]
study_train_complete <- train_studies_all[complete_train]

# ===== FILTER TO GSE185920 ONLY (TRAINING) =====
gse185920_train_idx <- study_train_complete == "GSE185920"

x_train_gse <- x_train_complete[gse185920_train_idx, ]
y_train_gse <- y_train_complete[gse185920_train_idx]

cat("=== GSE185920 Training Data ===\n")
cat("Total training samples:", nrow(x_train_complete), "\n")
cat("GSE185920 samples:", nrow(x_train_gse), "\n\n")


# ===== PREPARE TEST DATA =====
test_ages_all <- test_age_data$Age[1:ncol(Test)]
test_studies_all <- test_age_data$GSE[1:ncol(Test)]

x_test_all <- t(Test)

# Filter to complete cases
complete_test <- complete.cases(x_test_all, test_ages_all)
x_test_complete <- x_test_all[complete_test, ]
y_test_complete <- test_ages_all[complete_test]
study_test_complete <- test_studies_all[complete_test]

# ===== FILTER TO GSE185920 ONLY (TEST) =====
gse185920_test_idx <- study_test_complete == "GSE185920"

x_test_gse <- x_test_complete[gse185920_test_idx, ]
y_test_gse <- y_test_complete[gse185920_test_idx]

cat("=== GSE185920 Test Data ===\n")
cat("Total test samples:", nrow(x_test_complete), "\n")
cat("GSE185920 samples:", nrow(x_test_gse), "\n\n")

# ===== FIT MODEL ON GSE185920 TRAINING DATA ONLY =====
cat("Fitting Lasso on GSE185920 training data...\n")
lasso_gse <- cv.glmnet(x_train_gse, y_train_gse, alpha = 1, nfolds = 10)
cat("✓ Done!\n\n")


# ===== USE YOUR OPTIMAL LAMBDA =====
optimal_lambda <- 0.7
optimal_n_features <- 12

cat("Using optimal lambda:", optimal_lambda, "\n")
cat("Features:", optimal_n_features, "\n\n")

# ===== PREDICTIONS =====
# Training (GSE185920 only)
pred_train_gse <- predict(lasso_gse, newx = x_train_gse, s = optimal_lambda)[,1]

# Test (GSE185920 only)
pred_test_gse <- predict(lasso_gse, newx = x_test_gse, s = optimal_lambda)[,1]

# ===== CALCULATE METRICS =====
# Training
train_errors <- pred_train_gse - y_train_gse
train_mae <- mean(abs(train_errors))
train_rmse <- sqrt(mean(train_errors^2))
train_cor <- cor(pred_train_gse, y_train_gse)
train_r2 <- summary(lm(pred_train_gse ~ y_train_gse))$r.squared

# Test
test_errors <- pred_test_gse - y_test_gse
test_mae <- mean(abs(test_errors))
test_rmse <- sqrt(mean(test_errors^2))
test_cor <- cor(pred_test_gse, y_test_gse)
test_r2 <- summary(lm(pred_test_gse ~ y_test_gse))$r.squared

cat("=== GSE185920 ONLY - Performance ===\n")
cat("\nTraining:\n")
cat("  Samples:", length(y_train_gse), "\n")
cat("  MAE:", round(train_mae, 2), "years\n")
cat("  R²:", round(train_r2, 3), "\n")

cat("\nTest:\n")
cat("  Samples:", length(y_test_gse), "\n")
cat("  MAE:", round(test_mae, 2), "years\n")
cat("  R²:", round(test_r2, 3), "\n\n")


# ===== CREATE 2-PANEL PLOT (TRAIN & TEST) =====
png("optimal_model_GSE185920_only.png", width = 1400, height = 700, res = 120)

par(mfrow = c(1, 2))

# LEFT: Training Set (GSE185920)
plot(y_train_gse, pred_train_gse,
     pch = 19,
     col = rgb(0.3, 0.5, 0.8, 0.7),
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Training Set - GSE185920 Only\n",
                   "λ = ", optimal_lambda, " (", optimal_n_features, " CpGs)"))

abline(0, 1, col = "red", lwd = 3, lty = 2)
abline(lm(pred_train_gse ~ y_train_gse), col = "darkblue", lwd = 2.5)

legend("topleft",
       legend = c(
         paste0("n = ", length(y_train_gse), " samples"),
         paste0("r = ", round(train_cor, 3)),
         paste0("MAE = ", round(train_mae, 2), " years"),
         paste0("R² = ", round(train_r2, 3))
       ),
       bty = "n",
       cex = 1.1)

legend("bottomright",
       legend = c("Perfect prediction", "Fitted line"),
       col = c("red", "darkblue"),
       lwd = c(3, 2.5),
       lty = c(2, 1),
       bty = "n",
       cex = 0.9)

# RIGHT: Test Set (GSE185920)
plot(y_test_gse, pred_test_gse,
     pch = 19,
     col = rgb(0.8, 0.3, 0.5, 0.7),
     cex = 1.5,
     xlab = "Chronological Age (years)",
     ylab = "Predicted Age (years)",
     main = paste0("Test Set - GSE185920 Only\n",
                   "λ = ", optimal_lambda, " (", optimal_n_features, " CpGs)"))

abline(0, 1, col = "red", lwd = 3, lty = 2)
abline(lm(pred_test_gse ~ y_test_gse), col = "darkgreen", lwd = 2.5)

legend("topleft",
       legend = c(
         paste0("n = ", length(y_test_gse), " samples"),
         paste0("r = ", round(test_cor, 3)),
         paste0("MAE = ", round(test_mae, 2), " years"),
         paste0("R² = ", round(test_r2, 3))
       ),
       bty = "n",
       cex = 1.1)

legend("bottomright",
       legend = c("Perfect prediction", "Fitted line"),
       col = c("red", "darkgreen"),
       lwd = c(3, 2.5),
       lty = c(2, 1),
       bty = "n",
       cex = 0.9)

par(mfrow = c(1, 1))
dev.off()

cat("\n✓ Plot saved as 'optimal_model_GSE185920_only.png'\n")

# ===== SAVE RESULTS =====
# Training results
train_results_gse <- data.frame(
  dataset = "GSE185920",
  set = "Training",
  actual_age = y_train_gse,
  predicted_age = pred_train_gse,
  error = train_errors,
  abs_error = abs(train_errors)
)

# Test results
test_results_gse <- data.frame(
  dataset = "GSE185920",
  set = "Test",
  actual_age = y_test_gse,
  predicted_age = pred_test_gse,
  error = test_errors,
  abs_error = abs(test_errors)
)

# Combine
all_results_gse <- rbind(train_results_gse, test_results_gse)

write.csv(all_results_gse, "GSE185920_predictions.csv", row.names = FALSE)

cat("✓ Results saved as 'GSE185920_predictions.csv'\n")


#IMPUTING AND MISSINGNESS 

# ===== LOAD YOUR IMPUTATION FUNCTION =====
impute.knn <- function(x, k=.05, distmat){
  if(!is.matrix(x))
    stop("kNN does not work on data with mixed featured types. Therefore as a precaution kNN imputation only accept data in matrix form.")
  
  if(k < 1) k <- max(1, round(.05*nrow(x)))
  if(k > nrow(x)-1) stop("k is larger than the maximal number of neighbors.")
  if(!is.matrix(distmat)) distmat <- as.matrix(distmat)
  if(any(nrow(x) != dim(distmat)))
    stop("Distance matrix does not match dataset.")
  
  na.ind <- which(is.na(unname(x)), arr.ind=TRUE)
  
  NN <- apply(distmat, 1, function(z) order(z))
  fills <- apply(na.ind, 1, function(i){
    mean(na.exclude(x[NN[-1, i[1]], i[2]])[1:k])
  })
  x[na.ind] <- fills
  x
}
cat("✓ Imputation function loaded\n\n")

train_ages <- train_age_data$Age[1:ncol(Train)]
test_ages <- test_age_data$Age[1:ncol(Test)]

x_train_all <- t(Train)
x_test_all <- t(Test)

complete_train <- complete.cases(x_train_all, train_ages)
complete_test <- complete.cases(x_test_all, test_ages)

x_train <- x_train_all[complete_train, ]
y_train <- train_ages[complete_train]

x_test <- x_test_all[complete_test, ]
y_test <- test_ages[complete_test]

# ===== FIT YOUR OPTIMAL MODEL =====
cat("Fitting optimal Lasso model...\n")
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)

optimal_lambda <- 0.7  # Your optimal lambda
optimal_n_features <- 12  # Your feature count


# ===== GET SELECTED 12 CpGs =====
optimal_coefs <- coef(lasso_model, s = optimal_lambda)
selected_cpgs <- rownames(optimal_coefs)[optimal_coefs[,1] != 0]
selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]

cat("✓ Model fitted\n")
cat("Selected CpGs (", length(selected_cpgs), "):\n", sep = "")
print(selected_cpgs)

# ===== EXTRACT ONLY YOUR 12 SELECTED CpGs =====
# This is your "modelCpGs" matrix
modelCpGs_train <- Train[selected_cpgs, ]
modelCpGs_test <- Test[selected_cpgs, ]

cat("\nModel CpGs matrix:\n")
cat("  Training:", dim(modelCpGs_train), "(CpGs × samples)\n")
cat("  Test:", dim(modelCpGs_test), "(CpGs × samples)\n\n")


# ===== BASELINE: NO MISSING DATA =====
# FIX: Use the FULL x_test (all 50,000 CpGs), not just selected ones
baseline_pred_test <- predict(lasso_model, 
                              newx = x_test,  # ← All CpGs
                              s = optimal_lambda)[,1]

baseline_mae <- mean(abs(baseline_pred_test - y_test))

cat("=== Baseline (No Missing Data) ===\n")
cat("Test MAE:", round(baseline_mae, 2), "years\n\n")

# ===== TEST IMPUTATION WITH MISSING DATA =====
cat("=== Testing Imputation ===\n\n")

missing_percentages <- seq(5, 95, by = 5)

imputation_results <- data.frame()

for(miss_pct in missing_percentages) {
  
  cat("Testing", miss_pct, "% missing...\n")
  
  # Create copy of FULL test data
  Test_missing <- Test
  
  # Set miss_pct% of values in SELECTED CpGs to NA
  for(cpg in selected_cpgs) {
    n_samples <- ncol(Test_missing)
    n_to_miss <- round(n_samples * miss_pct / 100)
    set.seed(123 + match(cpg, selected_cpgs))  # Reproducible
    miss_samples <- sample(1:n_samples, n_to_miss)
    
    Test_missing[cpg, miss_samples] <- NA
  }
  
  # ===== IMPUTE ONLY THE SELECTED CpGs =====
  # Extract selected CpGs
  modelCpGs_missing <- Test_missing[selected_cpgs, ]
  
  # Transpose (samples as rows)
  met.data <- t(modelCpGs_missing)
  
  # Distance matrix
  all.dist <- dist(met.data)
  
  # Impute
  met.data.impute <- impute.knn(met.data, distmat = all.dist)
  
  # Transpose back
  modelCpGs_imputed <- t(met.data.impute)
  
  # Put imputed values back into full test matrix
  Test_imputed <- Test
  Test_imputed[selected_cpgs, ] <- modelCpGs_imputed
  
  # ===== PREDICT WITH IMPUTED DATA =====
  x_test_imputed <- t(Test_imputed)[complete_test, ]
  
  pred_imputed <- predict(lasso_model, 
                          newx = x_test_imputed,  # All 50k CpGs
                          s = optimal_lambda)[,1]
  
  mae_imputed <- mean(abs(pred_imputed - y_test))
  
  
  # ===== STORE RESULTS =====
  imputation_results <- rbind(imputation_results,
                              data.frame(
                                Percent_Missing = miss_pct,
                                MAE = round(mae_imputed, 2),
                                MAE_Increase = round(mae_imputed - baseline_mae, 2)
                              ))
  
  cat("  MAE:", round(mae_imputed, 2), 
      "years (Δ =", round(mae_imputed - baseline_mae, 2), ")\n")
}
  

# ===== RESULTS TABLE =====
cat("\n")
cat("========================================\n")
cat("IMPUTATION RESULTS\n")
cat("========================================\n")
cat("Baseline (0% missing):", round(baseline_mae, 2), "years\n\n")

print(imputation_results, row.names = FALSE)


# ===== PLOT =====
png("imputation_robustness.png", width = 1000, height = 700, res = 120)

plot(imputation_results$Percent_Missing, 
     imputation_results$MAE,
     type = "b",
     pch = 19,
     col = "red",
     lwd = 3,
     cex = 2,
     xlab = "Percentage of Missing Data in Selected CpGs (%)",
     ylab = "Test MAE (years)",
     main = paste0("Model Robustness to Missing Data\n",
                   "Optimal Model (12 CpGs)"),
     ylim = c(baseline_mae * 0.9, max(imputation_results$MAE) * 1.05))
# Baseline
abline(h = baseline_mae, col = "darkblue", lwd = 2, lty = 2)

# Error bars from baseline to each point
segments(imputation_results$Percent_Missing,
         baseline_mae,
         imputation_results$Percent_Missing,
         imputation_results$MAE,
         col = "gray50",
         lwd = 2)

# Caps on error bars
segments(imputation_results$Percent_Missing - 1.5,
         imputation_results$MAE,
         imputation_results$Percent_Missing + 1.5,
         imputation_results$MAE,
         col = "gray50",
         lwd = 2)
# Add difference labels
text(imputation_results$Percent_Missing, 
     imputation_results$MAE,
     labels = paste0("+", round(imputation_results$MAE - baseline_mae, 2)),
     pos = 4,
     cex = 0.9)

legend("topleft",
       legend = c(
         paste0("Baseline: ", round(baseline_mae, 2), " years"),
         "With kNN imputation",
         "Increase from baseline"
       ),
       col = c("darkblue", "red", "gray50"),
       lwd = c(2, 3, 2),
       lty = c(2, 1, 1),
       pch = c(NA, 19, NA),
       cex = 1,
       bg = "white")
dev.off()


#MISSINGNESS X100 EACH LEVEL

# Load imputation function
impute.knn <- function(x, k=.05, distmat){
  if(!is.matrix(x))
    stop("kNN does not work on data with mixed featured types. Therefore as a precaution kNN imputation only accept data in matrix form.")
  
  if(k < 1) k <- max(1, round(.05*nrow(x)))
  if(k > nrow(x)-1) stop("k is larger than the maximal number of neighbors.")
  if(!is.matrix(distmat)) distmat <- as.matrix(distmat)
  if(any(nrow(x) != dim(distmat)))
    stop("Distance matrix does not match dataset.")
  
  na.ind <- which(is.na(unname(x)), arr.ind=TRUE)
  
  NN <- apply(distmat, 1, function(z) order(z))
  fills <- apply(na.ind, 1, function(i){
    mean(na.exclude(x[NN[-1, i[1]], i[2]])[1:k])
  })
  x[na.ind] <- fills
  x
}

# ===== FIT OPTIMAL MODEL =====
cat("Fitting Lasso model...\n")
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)

optimal_lambda <- 0.7
optimal_n_features <- 12

# Get selected CpGs
optimal_coefs <- coef(lasso_model, s = optimal_lambda)
selected_cpgs <- rownames(optimal_coefs)[optimal_coefs[,1] != 0]
selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]

cat("✓ Model fitted with", length(selected_cpgs), "CpGs\n\n")

# ===== BASELINE PERFORMANCE =====
baseline_pred <- predict(lasso_model, newx = x_test, s = optimal_lambda)[,1]
baseline_mae <- mean(abs(baseline_pred - y_test))
baseline_rmse <- sqrt(mean((baseline_pred - y_test)^2))
baseline_r2 <- summary(lm(baseline_pred ~ y_test))$r.squared

cat("Baseline (0% missing):\n")
cat("  MAE:", round(baseline_mae, 2), "\n")
cat("  RMSE:", round(baseline_rmse, 2), "\n")
cat("  R²:", round(baseline_r2, 3), "\n\n")

# ===== MISSINGNESS TESTING WITH 100 ITERATIONS =====
missing_percentages <- seq(0, 90, by = 10)
n_iterations <- 100  # 100 runs per missingness level

cat("Running", n_iterations, "iterations for each missingness level...\n")
cat("This may take a few minutes...\n\n")

# Store ALL individual results
all_results <- data.frame()

for(miss_pct in missing_percentages) {
  
  cat("Testing", miss_pct, "% missing (", n_iterations, "iterations)...\n")
  
  for(iter in 1:n_iterations) {
    if(miss_pct == 0) {
      # Baseline - no imputation needed
      mae_val <- baseline_mae
      rmse_val <- baseline_rmse
      r2_val <- baseline_r2
      
    } else {
      # Create missing data with unique random seed
      Test_missing <- Test
      
      for(cpg_idx in 1:length(selected_cpgs)) {
        cpg <- selected_cpgs[cpg_idx]
        n_samples <- ncol(Test_missing)
        n_to_miss <- round(n_samples * miss_pct / 100)
        
        set.seed(iter * 1000 + cpg_idx * 10 + miss_pct)
        miss_samples <- sample(1:n_samples, n_to_miss)
        Test_missing[cpg, miss_samples] <- NA
      }
      # Impute
      modelCpGs_missing <- Test_missing[selected_cpgs, ]
      met.data <- t(modelCpGs_missing)
      all.dist <- dist(met.data)
      met.data.impute <- impute.knn(met.data, distmat = all.dist)
      
      # Put back
      Test_imputed <- Test
      Test_imputed[selected_cpgs, ] <- t(met.data.impute)
      
      # Predict
      x_test_imputed <- t(Test_imputed)[complete_test, ]
      pred_imputed <- predict(lasso_model, newx = x_test_imputed, s = optimal_lambda)[,1]
      
      # Calculate metrics
      errors <- pred_imputed - y_test
      mae_val <- mean(abs(errors))
      rmse_val <- sqrt(mean(errors^2))
      r2_val <- summary(lm(pred_imputed ~ y_test))$r.squared
    }
    # Store this iteration's results
    all_results <- rbind(all_results,
                         data.frame(
                           Percent_Missing = miss_pct,
                           Iteration = iter,
                           MAE = mae_val,
                           RMSE = rmse_val,
                           R2 = r2_val
                         ))
  }
  
  cat("  ✓ Complete\n")
}

cat("\n✓ All iterations complete!\n")
cat("Total runs:", nrow(all_results), "\n\n")

# ===== CALCULATE SUMMARY STATISTICS - FIXED VERSION =====
cat("Calculating summary statistics from", n_iterations, "iterations...\n")

summary_table <- data.frame()

for(miss_pct in missing_percentages) {
  
  subset_data <- all_results[all_results$Percent_Missing == miss_pct, ]
  
  summary_table <- rbind(summary_table,
                         data.frame(
                           Percent_Missing = miss_pct,
                           MAE_median = median(subset_data$MAE),
                           MAE_q25 = quantile(subset_data$MAE, 0.25),
                           MAE_q75 = quantile(subset_data$MAE, 0.75),
                           MAE_mean = mean(subset_data$MAE),
                           MAE_sd = sd(subset_data$MAE),
                           RMSE_median = median(subset_data$RMSE),
                           RMSE_q25 = quantile(subset_data$RMSE, 0.25),
                           RMSE_q75 = quantile(subset_data$RMSE, 0.75),
                           RMSE_mean = mean(subset_data$RMSE),
                           RMSE_sd = sd(subset_data$RMSE),
                           R2_median = median(subset_data$R2),
                           R2_q25 = quantile(subset_data$R2, 0.25),
                           R2_q75 = quantile(subset_data$R2, 0.75),
                           R2_mean = mean(subset_data$R2),
                           R2_sd = sd(subset_data$R2)
                         ))
}

cat("✓ Summary statistics calculated\n\n")

# ===== SAVE RESULTS =====
write.csv(all_results, "missingness_all_iterations.csv", row.names = FALSE)
write.csv(summary_table, "missingness_summary_stats.csv", row.names = FALSE)

cat("✓ Results saved\n\n")


# ===== CREATE SEPARATE PLOTS FOR MAE, RMSE, R² =====

# ===== PLOT 1: MAE WITH BOXPLOTS =====
png("missingness_MAE_boxplots.png", width = 1200, height = 800, res = 150)

par(mar = c(5, 5, 4, 2))

# Create boxplot
boxplot(MAE ~ Percent_Missing, 
        data = all_results,
        col = "lightblue",
        border = "steelblue",
        outline = FALSE,  # Don't show outliers as points
        xlab = "Percentage of Missing Data (%)",
        ylab = "Test MAE (years)",
        main = paste0("Model Performance vs Missing Data (MAE)\n",
                      "100 iterations per missingness level"),
        cex.lab = 1.3,
        cex.main = 1.3,
        cex.axis = 1.1,
        las = 1,
        ylim = c(min(all_results$MAE) * 0.95, max(all_results$MAE) * 1.05))
# Add median line connecting boxplots
lines(1:length(missing_percentages), 
      summary_table$MAE_median,
      col = "red",
      lwd = 4,
      lty = 1)
# Add median points
points(1:length(missing_percentages), 
       summary_table$MAE_median,
       pch = 19,
       col = "red",
       cex = 1)

# Baseline reference
abline(h = baseline_mae, col = "darkgreen", lwd = 3, lty = 2)
# Add threshold lines
abline(h = baseline_mae + 0.5, col = "orange", lwd = 2, lty = 3)
abline(h = baseline_mae + 1, col = "darkred", lwd = 2, lty = 3)

# Add shaded zones
rect(par("usr")[1], 0, par("usr")[2], baseline_mae + 0.5, 
     col = rgb(0.5, 0.9, 0.5, 0.1), border = NA)
rect(par("usr")[1], baseline_mae + 0.5, par("usr")[2], baseline_mae + 1, 
     col = rgb(0.9, 0.9, 0.5, 0.1), border = NA)
rect(par("usr")[1], baseline_mae + 1, par("usr")[2], par("usr")[4], 
     col = rgb(0.9, 0.5, 0.5, 0.1), border = NA)
# Redraw boxplot on top
boxplot(MAE ~ Percent_Missing, 
        data = all_results,
        col = "lightblue",
        border = "steelblue",
        outline = FALSE,
        add = TRUE,
        axes = FALSE)

# Redraw median line on top


points(1:length(missing_percentages), 
       summary_table$MAE_median,
       pch = 19,
       col = "red",
       cex = 1)
# Legend
legend("topleft",
       legend = c(
         "Median (100 iterations)",
         "Interquartile range (box)",
         paste0("Baseline: ", round(baseline_mae, 2), " years"),
         "Acceptable (+0.5y)",
         "Warning (+1.0y)"
       ),
       col = c("red", "lightblue", "darkgreen", "orange", "darkred"),
       lwd = c(4, 10, 3, 2, 2),
       lty = c(1, 1, 2, 3, 3),
       pch = c(19, 15, NA, NA, NA),
       pt.cex = c(2, 3, NA, NA, NA),
       cex = 1,
       bg = "white",
       box.lwd = 2)
# Add value labels at median
text(1:length(missing_percentages), 
     summary_table$MAE_median,
     labels = round(summary_table$MAE_median, 2),
     pos = 3,
     cex = 0.9,
     font = 2)

grid(col = "white", lwd = 1.5)

dev.off()

cat("✓ MAE boxplot saved as 'missingness_MAE_boxplots.png'\n")

# ===== PLOT 2: RMSE WITH BOXPLOTS =====
png("missingness_RMSE_boxplots.png", width = 1200, height = 800, res = 150)

par(mar = c(5, 5, 4, 2))

boxplot(RMSE ~ Percent_Missing, 
        data = all_results,
        col = "lightcoral",
        border = "darkred",
        outline = FALSE,
        xlab = "Percentage of Missing Data (%)",
        ylab = "Test RMSE (years)",
        main = paste0("Model Performance vs Missing Data (RMSE)\n",
                      "100 iterations per missingness level"),
        cex.lab = 1.3,
        cex.main = 1.3,
        cex.axis = 1.1,
        las = 1,
        ylim = c(min(all_results$RMSE) * 0.95, max(all_results$RMSE) * 1.05))

# Median line
lines(1:length(missing_percentages), 
      summary_table$RMSE_median,
      col = "darkblue",
      lwd = 4)

points(1:length(missing_percentages), 
       summary_table$RMSE_median,
       pch = 19,
       col = "darkblue",
       cex = 1)

# Baseline
abline(h = baseline_rmse, col = "darkgreen", lwd = 3, lty = 2)

# Threshold lines
abline(h = baseline_rmse + 0.5, col = "orange", lwd = 2, lty = 3)
abline(h = baseline_rmse + 1, col = "darkred", lwd = 2, lty = 3)

# Shaded zones
rect(par("usr")[1], 0, par("usr")[2], baseline_rmse + 0.5, 
     col = rgb(0.5, 0.9, 0.5, 0.1), border = NA)
rect(par("usr")[1], baseline_rmse + 0.5, par("usr")[2], baseline_rmse + 1, 
     col = rgb(0.9, 0.9, 0.5, 0.1), border = NA)
rect(par("usr")[1], baseline_rmse + 1, par("usr")[2], par("usr")[4], 
     col = rgb(0.9, 0.5, 0.5, 0.1), border = NA)

# Redraw boxplot on top
boxplot(RMSE ~ Percent_Missing, 
        data = all_results,
        col = "lightcoral",
        border = "darkred",
        outline = FALSE,
        add = TRUE,
        axes = FALSE)
# Redraw median line
lines(1:length(missing_percentages), 
      summary_table$RMSE_median,
      col = "darkblue",
      lwd = 4)

points(1:length(missing_percentages), 
       summary_table$RMSE_median,
       pch = 19,
       col = "darkblue",
       cex = 1)

legend("topleft",
       legend = c(
         "Median (100 iterations)",
         "Interquartile range",
         paste0("Baseline: ", round(baseline_rmse, 2), " years")
       ),
       col = c("darkblue", "lightcoral", "darkgreen"),
       lwd = c(4, 10, 3),
       lty = c(1, 1, 2),
       pch = c(19, 15, NA),
       pt.cex = c(2, 3, NA),
       cex = 1,
       bg = "white")

# Add median values
text(1:length(missing_percentages), 
     summary_table$RMSE_median,
     labels = round(summary_table$RMSE_median, 2),
     pos = 3,
     cex = 0.9,
     font = 2)

grid(col = "white", lwd = 1.5)
dev.off()

cat("✓ RMSE boxplot saved as 'missingness_RMSE_boxplots.png'\n")

# ===== PLOT 3: R² WITH BOXPLOTS =====
png("missingness_R2_boxplots.png", width = 1200, height = 800, res = 150)

par(mar = c(5, 5, 4, 2))

boxplot(R2 ~ Percent_Missing, 
        data = all_results,
        col = "lightgreen",
        border = "darkgreen",
        outline = FALSE,
        xlab = "Percentage of Missing Data (%)",
        ylab = "Test R² (Coefficient of Determination)",
        main = paste0("Model Fit vs Missing Data (R²)\n",
                      "100 iterations per missingness level"),
        cex.lab = 1.3,
        cex.main = 1.3,
        cex.axis = 1.1,
        las = 1,
        ylim = c(min(all_results$R2) * 0.98, 1))
# Median line
lines(1:length(missing_percentages), 
      summary_table$R2_median,
      col = "hotpink",
      lwd = 4)

points(1:length(missing_percentages), 
       summary_table$R2_median,
       pch = 19,
       col = "hotpink",
       cex = 1)

# Baseline
abline(h = baseline_r2, col = "darkgreen", lwd = 3, lty = 2)

# Reference lines
abline(h = 0.9, col = "orange", lwd = 2, lty = 3)
abline(h = 0.85, col = "darkred", lwd = 2, lty = 3)
legend("bottomleft",
       legend = c(
         "Median (100 iterations)",
         "Interquartile range",
         paste0("Baseline: ", round(baseline_r2, 3)),
         "R² = 0.9 threshold",
         "R² = 0.85 threshold"
       ),
       col = c("hotpink", "lightgreen", "darkgreen", "orange", "darkred"),
       lwd = c(4, 10, 3, 2, 2),
       lty = c(1, 1, 2, 3, 3),
       pch = c(19, 15, NA, NA, NA),
       pt.cex = c(2, 3, NA, NA, NA),
       cex = 1,
       bg = "white")

# Add median values
text(1:length(missing_percentages), 
     summary_table$R2_median,
     labels = round(summary_table$R2_median, 3),
     pos = 3,
     cex = 0.9,
     font = 2)

grid(col = "white", lwd = 1.5)

dev.off()

cat("✓ R² boxplot saved as 'missingness_R2_boxplots.png'\n")

# ===== COMBINED 3-PANEL FIGURE =====
png("missingness_comprehensive_3panel.png", width = 1800, height = 600, res = 150)

par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))


# Panel 1: MAE
boxplot(MAE ~ Percent_Missing, 
        data = all_results,
        col = "lightblue",
        border = "steelblue",
        outline = FALSE,
        xlab = "% Missing",
        ylab = "MAE (years)",
        main = "A. Mean Absolute Error",
        cex.lab = 1.2,
        cex.main = 1.3,
        las = 1)
lines(1:length(missing_percentages), summary_table$MAE_median,
      col = "red", lwd = 4)
points(1:length(missing_percentages), summary_table$MAE_median,
       pch = 19, col = "red", cex = 1)

abline(h = baseline_mae, col = "darkgreen", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Median", "IQR", "Baseline"),
       col = c("red", "lightblue", "darkgreen"),
       lwd = c(4, 10, 2),
       pch = c(19, 15, NA),
       pt.cex = c(2, 3, NA),
       cex = 0.9,
       bg = "white")
# Panel 2: RMSE
boxplot(RMSE ~ Percent_Missing, 
        data = all_results,
        col = "lightcoral",
        border = "darkred",
        outline = FALSE,
        xlab = "% Missing",
        ylab = "RMSE (years)",
        main = "B. Root Mean Square Error",
        cex.lab = 1.2,
        cex.main = 1.3,
        las = 1)

lines(1:length(missing_percentages), summary_table$RMSE_median,
      col = "steelblue", lwd = 4)
points(1:length(missing_percentages), summary_table$RMSE_median,
       pch = 19, col = "steelblue", cex = 1)

abline(h = baseline_rmse, col = "darkgreen", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Median", "IQR", "Baseline"),
       col = c("steelblue", "lightcoral", "darkgreen"),
       lwd = c(4, 10, 2),
       pch = c(19, 15, NA),
       pt.cex = c(2, 3, NA),
       cex = 0.9,
       bg = "white")


# Panel 3: R²
boxplot(R2 ~ Percent_Missing, 
        data = all_results,
        col = "lightgreen",
        border = "darkgreen",
        outline = FALSE,
        xlab = "% Missing",
        ylab = "R²",
        main = "C. Coefficient of Determination",
        cex.lab = 1.2,
        cex.main = 1.3,
        las = 1,
        ylim = c(min(all_results$R2) * 0.98, 1))

lines(1:length(missing_percentages), summary_table$R2_median,
      col = "hotpink", lwd = 4)
points(1:length(missing_percentages), summary_table$R2_median,
       pch = 19, col = "hotpink", cex = 1)

abline(h = baseline_r2, col = "darkgreen", lwd = 2, lty = 2)


legend("bottomleft",
       legend = c("Median", "IQR", "Baseline"),
       col = c("hotpink", "lightgreen", "darkgreen"),
       lwd = c(4, 10, 2, 2),
       pch = c(19, 15, NA, NA),
       pt.cex = c(2, 3, NA, NA),
       cex = 0.9,
       bg = "white")

par(mfrow = c(1, 1))
dev.off()

cat("✓ Combined 3-panel saved as 'missingness_comprehensive_3panel.png'\n")

# ===== DETAILED SUMMARY TABLE =====
cat("\n========================================\n")
cat("MISSINGNESS ANALYSIS SUMMARY\n")
cat("========================================\n\n")

cat("Iterations per missingness level:", n_iterations, "\n")
cat("Total runs:", nrow(all_results), "\n\n")

cat("Performance at Key Missingness Levels:\n")
cat("─────────────────────────────────────────\n")

for(i in 1:nrow(summary_table)) {
  miss <- summary_table$Percent_Missing[i]
  
  cat(sprintf("\n%d%% Missing:\n", miss))
  cat(sprintf("  MAE:  Median=%.2f (IQR: %.2f-%.2f)\n",
              summary_table$MAE_median[i],
              summary_table$MAE_q25[i],
              summary_table$MAE_q75[i]))
  cat(sprintf("  RMSE: Median=%.2f (IQR: %.2f-%.2f)\n",
              summary_table$RMSE_median[i],
              summary_table$RMSE_q25[i],
              summary_table$RMSE_q75[i]))
  cat(sprintf("  R²:   Median=%.3f (IQR: %.3f-%.3f)\n",
              summary_table$R2_median[i],
              summary_table$R2_q25[i],
              summary_table$R2_q75[i]))
}

# ===== FIND ACCEPTABLE MISSINGNESS THRESHOLD =====
cat("\n========================================\n")
cat("ROBUSTNESS THRESHOLDS\n")
cat("========================================\n\n")

# MAE increase <0.5y
acceptable_05 <- summary_table[summary_table$MAE_q75 < baseline_mae + 0.5, ]
if(nrow(acceptable_05) > 0) {
  max_05 <- max(acceptable_05$Percent_Missing)
  cat("MAE increase <0.5y (75th percentile):", max_05, "%\n")
}

# MAE increase <1y
acceptable_1 <- summary_table[summary_table$MAE_q75 < baseline_mae + 1, ]
if(nrow(acceptable_1) > 0) {
  max_1 <- max(acceptable_1$Percent_Missing)
  cat("MAE increase <1y (75th percentile):", max_1, "%\n")
}
# R² stays >0.85
r2_threshold <- summary_table[summary_table$R2_q25 > 0.85, ]
if(nrow(r2_threshold) > 0) {
  max_r2 <- max(r2_threshold$Percent_Missing)
  cat("R² stays >0.85 (25th percentile):", max_r2, "%\n")
}

# ===== PUBLICATION-READY SUMMARY TABLE =====
pub_table <- data.frame(
  Missing = paste0(summary_table$Percent_Missing, "%"),
  MAE_Median = round(summary_table$MAE_median, 2),
  MAE_IQR = paste0(round(summary_table$MAE_q25, 2), "-", round(summary_table$MAE_q75, 2)),
  RMSE_Median = round(summary_table$RMSE_median, 2),
  RMSE_IQR = paste0(round(summary_table$RMSE_q25, 2), "-", round(summary_table$RMSE_q75, 2)),
  R2_Median = round(summary_table$R2_median, 3),
  R2_IQR = paste0(round(summary_table$R2_q25, 3), "-", round(summary_table$R2_q75, 3))
)

cat("\n========================================\n")
cat("PUBLICATION TABLE\n")
cat("========================================\n\n")

print(pub_table, row.names = FALSE)

write.csv(pub_table, "missingness_publication_table.csv", row.names = FALSE)

#===CPG HEATMAP===

# ===== GET YOUR 12 SELECTED CpGs FROM LASSO MODEL =====
# Option 1: If you have them saved
# selected_cpgs <- readRDS("selected_cpgs.RDS")

# Option 2: Extract from model
library(glmnet)

# Load your trained model
x_train <- t(Train)
y_train <- ages_all

complete_idx <- complete.cases(x_train, y_train)
x_train <- x_train[complete_idx, ]
y_train <- y_train[complete_idx]

# Fit or load model
lasso_model <- cv.glmnet(x_train, y_train, alpha = 1, nfolds = 10)

optimal_lambda <- 0.7
# Get selected CpGs
optimal_coefs <- coef(lasso_model, s = optimal_lambda)
selected_cpgs <- rownames(optimal_coefs)[optimal_coefs[,1] != 0]
selected_cpgs <- selected_cpgs[selected_cpgs != "(Intercept)"]

cat("Selected CpGs (", length(selected_cpgs), "):\n", sep = "")
print(selected_cpgs)

# ===== EXTRACT AND ORDER DATA =====
# Get methylation data for 12 CpGs
heatmap_data <- Train[selected_cpgs, ]

# Order samples by age
age_order <- order(ages_all)
data_ordered <- heatmap_data[, age_order]
ages_ordered <- ages_all[age_order]
studies_ordered <- studies_all[age_order]

# Order CpGs by their correlation with age
cpg_cors <- apply(heatmap_data, 1, function(cpg) {
  cor(cpg, ages_all, use = "complete.obs")
})

# Sort CpGs: most negative to most positive correlation
cpg_order <- order(cpg_cors)
data_ordered <- data_ordered[cpg_order, ]
cpgs_ordered <- selected_cpgs[cpg_order]
cors_ordered <- cpg_cors[cpg_order]

cat("\nCpG ordering (by age correlation):\n")
print(data.frame(CpG = cpgs_ordered, Correlation = round(cors_ordered, 3)))

# ===== CREATE HEATMAP =====
png("age_ordered_heatmap.png", width = 2000, height = 1200, res = 150)

# Layout: age bar + study bar + main heatmap + legend
layout(matrix(c(1, 1, 1, 1,
                2, 2, 2, 2,
                3, 3, 3, 4), 
              nrow = 3, byrow = TRUE),
       heights = c(0.3, 0.3, 5))
# ===== TOP BAR: AGE GROUPS =====
par(mar = c(0, 10, 3, 2))

# Create age groups
age_breaks <- c(19, 25, 35, 45, 60)
age_groups <- cut(ages_ordered, breaks = age_breaks, 
                  labels = c("19-25", "26-35", "36-45", "46-60"),
                  include.lowest = TRUE)

# Group colors
group_colors <- c("19-25" = "#d4e7f7",
                  "26-35" = "#f7f2d4", 
                  "36-45" = "#f7dcd4",
                  "46-60" = "#f7d4d4")

age_color_vec <- group_colors[as.character(age_groups)]

# FIX: Correct dimensions for image()
image(1:ncol(data_ordered), 1,
      matrix(1:ncol(data_ordered), ncol = 1),  # ← FIXED: ncol instead of nrow
      col = age_color_vec,
      axes = FALSE,
      xlab = "",
      ylab = "")

# Add age group labels
group_centers <- tapply(1:length(age_groups), age_groups, median)
text(group_centers, 1, names(group_centers), cex = 1.3, font = 2)

box(lwd = 2)
# ===== BAR 2: STUDY (FIXED) =====
par(mar = c(0, 10, 0, 2))

study_colors <- c("GSE185920" = "#b3d9ff",
                  "GSE185445" = "#ffb3d9",
                  "GSE149318" = "#b3ffb3")

study_color_vec <- study_colors[studies_ordered]

# FIXED: Use rect() instead of image()
plot(c(0.5, ncol(data_ordered) + 0.5), c(0, 1),
     type = "n", axes = FALSE, xlab = "", ylab = "")

for(i in 1:ncol(data_ordered)) {
  rect(i - 0.5, 0, i + 0.5, 1,
       col = study_color_vec[i],
       border = NA)
}

box(lwd = 2)
# ===== MAIN HEATMAP =====
par(mar = c(4, 10, 0, 2))

# Use ACTUAL beta values (0-1) not z-scores - matches image better
# Or use z-scores for better contrast
use_raw <- FALSE  # Set to TRUE for raw beta values

if(use_raw) {
  plot_data <- data_ordered
  color_breaks <- seq(0, 1, length.out = 100)
} else {
  plot_data <- t(scale(t(data_ordered)))  # Z-score per CpG
  color_breaks <- seq(-2.5, 2.5, length.out = 100)
}

# Color scheme matching the image (blue-white-red)
heatmap_colors <- colorRampPalette(c("#0033CC", "#6699FF", "#99CCFF", 
                                     "#FFFFFF", 
                                     "#FF9999", "#FF6666", "#CC0000"))(100)
# FIXED: Remove breaks parameter to avoid error
image(1:ncol(data_scaled), 
      1:nrow(data_scaled),
      t(data_scaled),
      col = heatmap_colors,
      axes = FALSE,
      xlab = "",
      ylab = "")

# Add CpG labels on left (LARGE)
axis(2, at = 1:nrow(plot_data), 
     labels = cpgs_ordered,
     las = 2,
     cex.axis = 1.4,
     font = 2,
     tick = FALSE)

# Add "Methylation status" label
mtext("Methylation status", side = 2, line = 8, cex = 1.5, font = 2)
# Add age on bottom
mtext("Age (years) →", side = 1, line = 2.5, cex = 1.5, font = 2)

# Subtle gridlines (optional)
abline(h = seq(0.5, nrow(plot_data) + 0.5, 1), 
       col = "white", lwd = 0.3)

# Box
box(lwd = 3)


# ===== LEGEND: METHYLATION SCALE =====
par(mar = c(4, 1, 0, 4))

# Create vertical color bar
legend_vals <- color_breaks
plot(c(0, 1), range(color_breaks), type = "n", 
     axes = FALSE, xlab = "", ylab = "")

# Draw rectangles
rect(0.2, color_breaks[-length(color_breaks)],
     0.8, color_breaks[-1],
     col = heatmap_colors,
     border = NA)

# Add scale
if(use_raw) {
  axis(4, at = c(0, 0.25, 0.5, 0.75, 1),
       labels = c("0", "0.25", "0.5", "0.75", "1"),
       las = 1,
       cex.axis = 1.2)
  mtext("Methylation level", side = 4, line = 3, cex = 1.3, font = 2)
} else {
  axis(4, at = c(-2, -1, 0, 1, 2),
       labels = c("-2", "-1", "0", "+1", "+2"),
       las = 1,
       cex.axis = 1.2)
  mtext("Z-score", side = 4, line = 3, cex = 1.3, font = 2)
}

box(lwd = 2)

dev.off()

cat("\n✓ Age-ordered heatmap saved as 'age_ordered_heatmap.png'\n")

# ===== PRINT CpG PATTERN SUMMARY =====
cat("\n=== Methylation Patterns ===\n")
cat("CpGs showing HYPERMETHYLATION with age (r > 0):\n")
hyper <- cpg_cor_table[cpg_cor_table$Correlation > 0, ]
print(hyper)

cat("\nCpGs showing HYPOMETHYLATION with age (r < 0):\n")
hypo <- cpg_cor_table[cpg_cor_table$Correlation < 0, ]
print(hypo)

cat("\nExpected pattern: Blue→Red (hypermethylation) or Red→Blue (hypomethylation)\n")