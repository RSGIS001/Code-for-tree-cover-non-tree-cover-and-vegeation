# ====================== Libraries ======================
library(raster)
library(glmnet)
library(foreach)
library(doParallel)
library(matrixStats)

# 固定随机种子（可复现）
set.seed(123)

# 设置工作目录（按需修改）
setwd("E:/Work22/08_Regress/Detrend02")

# 读取栅格栈
veg_stack <- stack("R_EC-Earth3-CC_Detrend.tif")
x1_stack  <- stack("TC_BG_EC-Earth3-CC_alog_Detrend.tif")
x2_stack  <- stack("NTC_BG_EC-Earth3-CC_alog_Detrend.tif")
x3_stack  <- stack("AI_EC-Earth3-CC_Detrend.tif")
x4_stack  <- stack("TAS_EC-Earth3-CC_Detrend.tif")
x5_stack  <- stack("RSS_EC-Earth3-CC_Detrend.tif")
x6_stack  <- stack("WS_EC-Earth3-CC_Detrend.tif")
x7_stack  <- stack("CO2_Detrend_512_256.tif")

veg_mat <- getValues(veg_stack)
x1_mat  <- getValues(x1_stack)
x2_mat  <- getValues(x2_stack)
x3_mat  <- getValues(x3_stack)
x4_mat  <- getValues(x4_stack)
x5_mat  <- getValues(x5_stack)
x6_mat  <- getValues(x6_stack)
x7_mat  <- getValues(x7_stack)

ncell <- nrow(veg_mat)
nt    <- ncol(veg_mat)

# ====================== Z-score (pixel-wise) ======================
zscore_rows <- function(mat){
  mu <- rowMeans2(mat, na.rm = TRUE)
  s  <- rowSds(mat, na.rm = TRUE)
  s[s == 0] <- NA
  sweep(sweep(mat, 1, mu, "-"), 1, s, "/")
}

y_z  <- zscore_rows(veg_mat)
x1_z <- zscore_rows(x1_mat)
x2_z <- zscore_rows(x2_mat)
x3_z <- zscore_rows(x3_mat)
x4_z <- zscore_rows(x4_mat)
x5_z <- zscore_rows(x5_mat)
x6_z <- zscore_rows(x6_mat)
x7_z <- zscore_rows(x7_mat)

# ====================== Output containers ======================
coef_mat   <- matrix(0,  ncell, 7)
vif_mat    <- matrix(NA, ncell, 7)
y_sim_mat  <- matrix(NA, ncell, nt)

r_mat   <- rep(NA, ncell)
p_mat   <- rep(NA, ncell)
r2_mat  <- rep(NA, ncell)
r2p_mat <- rep(NA, ncell)

# 岭回归稳定性指标
beta_stab_TC  <- rep(NA, ncell)
beta_stab_NTC <- rep(NA, ncell)
vif_path_TC   <- rep(NA, ncell)
vif_path_NTC  <- rep(NA, ncell)
deltaR2_TC    <- rep(NA, ncell)
deltaR2_NTC   <- rep(NA, ncell)

# ====================== Valid pixels ======================
valid_idx <- which(rowSums(!is.na(veg_mat)) >= 3)

# ====================== Parallel ======================
cl <- makeCluster(max(1, detectCores(FALSE) - 1))
registerDoParallel(cl)

results <- foreach(i = valid_idx, .packages = c("glmnet", "matrixStats"), 
                   .errorhandling = "pass") %dopar% {
                     
                     tryCatch({
                       yraw <- veg_mat[i, ]
                       yz   <- y_z[i, ]
                       
                       Xz_all <- cbind(
                         x1_z[i,], x2_z[i,], x3_z[i,],
                         x4_z[i,], x5_z[i,], x6_z[i,], x7_z[i,]
                       )
                       
                       ok <- complete.cases(cbind(yz, Xz_all))
                       if(sum(ok) < 5) return(NULL)
                       
                       yz <- yz[ok]
                       Xz <- Xz_all[ok, ]
                       n  <- length(yz)
                       
                       # ===== 岭回归 CV（alpha = 0） =====
                       fit_cv <- cv.glmnet(
                         Xz, yz,
                         alpha = 0,                  # <<< 岭回归
                         nfolds = min(3, n),
                         standardize = FALSE,
                         intercept = TRUE
                       )
                       
                       lambda_min <- fit_cv$lambda.min
                       lambda_seq <- fit_cv$lambda
                       fit_path   <- fit_cv$glmnet.fit
                       
                       beta_path <- as.matrix(coef(fit_path))[-1, ]
                       
                       stab_fun <- function(b){
                         if(length(b) == 0 || all(is.na(b)) || all(b == 0)) return(0)
                         maxb <- max(abs(b), na.rm = TRUE)
                         if(maxb == 0) return(0)
                         mean(abs(b) > 0.1 * maxb, na.rm = TRUE)
                       }
                       
                       stab_TC  <- stab_fun(beta_path[1, ])
                       stab_NTC <- stab_fun(beta_path[2, ])
                       
                       XtX <- crossprod(scale(Xz, center = FALSE, scale = FALSE))
                       vif_path <- matrix(NA, nrow = ncol(Xz), ncol = length(lambda_seq))
                       
                       for(j in seq_along(lambda_seq)) {
                         lam <- lambda_seq[j]
                         if(lam > 0) {
                           inv <- try(solve(XtX + lam * diag(ncol(Xz))), silent = TRUE)
                           if(!inherits(inv, "try-error")) {
                             vif_path[, j] <- diag(XtX) * diag(inv)
                           }
                         }
                       }
                       
                       vif_path_TC  <- mean(vif_path[1, ], na.rm = TRUE)
                       vif_path_NTC <- mean(vif_path[2, ], na.rm = TRUE)
                       
                       R2_fun <- function(Xsub){
                         if(ncol(Xsub) == 0) return(0)
                         fit <- glmnet(Xsub, yz, alpha = 0, lambda = lambda_min,
                                       standardize = FALSE, intercept = TRUE)
                         ypred <- predict(fit, Xsub)
                         if(all(is.na(ypred))) return(0)
                         cor(ypred, yz, use = "complete.obs")^2
                       }
                       
                       R2_full <- R2_fun(Xz)
                       dR2_TC  <- if(ncol(Xz) > 1) R2_full - R2_fun(Xz[, -1, drop = FALSE]) else 0
                       dR2_NTC <- if(ncol(Xz) > 2) R2_full - R2_fun(Xz[, -2, drop = FALSE]) else 0
                       
                       # ===== final 岭回归 fit =====
                       fit <- glmnet(Xz, yz, alpha = 0, lambda = lambda_min,
                                     standardize = FALSE, intercept = TRUE)
                       
                       coef_all   <- as.numeric(coef(fit))
                       intercept_z <- coef_all[1]
                       coef_z      <- coef_all[-1]
                       
                       coef_z[is.na(coef_z)] <- 0
                       
                       sd_y <- sd(yraw, na.rm = TRUE)
                       sd_x <- c(
                         sd(x1_mat[i,], na.rm = TRUE),
                         sd(x2_mat[i,], na.rm = TRUE),
                         sd(x3_mat[i,], na.rm = TRUE),
                         sd(x4_mat[i,], na.rm = TRUE),
                         sd(x5_mat[i,], na.rm = TRUE),
                         sd(x6_mat[i,], na.rm = TRUE),
                         sd(x7_mat[i,], na.rm = TRUE)
                       )
                       
                       coef_raw <- coef_z * sd_y / sd_x
                       coef_raw[!is.finite(coef_raw)] <- 0
                       
                       ysim <- rep(NA, nt)
                       if(sum(ok) > 0) {
                         ypred_z <- predict(fit, Xz)
                         ysim[ok] <- ypred_z * sd_y + mean(yraw, na.rm = TRUE)
                       }
                       
                       vif_min <- rep(NA, ncol(Xz))
                       if(lambda_min > 0) {
                         inv_min <- try(solve(XtX + lambda_min * diag(ncol(Xz))), silent = TRUE)
                         if(!inherits(inv_min, "try-error")) {
                           vif_min <- diag(XtX) * diag(inv_min)
                         }
                       }
                       
                       ok2 <- !is.na(yraw) & !is.na(ysim)
                       r <- r2 <- p_val <- r2p_val <- NA
                       
                       if(sum(ok2) > 2) {
                         cor_test <- try(cor.test(ysim[ok2], yraw[ok2]), silent = TRUE)
                         if(!inherits(cor_test, "try-error")) {
                           r <- cor_test$estimate
                           p_val <- cor_test$p.value
                           r2 <- r^2
                           df1 <- 1
                           df2 <- sum(ok2) - 2
                           if(df2 > 0 && r2 < 1) {
                             F_stat <- (r2 * df2) / (1 - r2)
                             r2p_val <- pf(F_stat, df1, df2, lower.tail = FALSE)
                           }
                         }
                       }
                       
                       list(
                         i = i,
                         coef = coef_raw,
                         vif = vif_min,
                         ysim = ysim,
                         r = r, p = p_val,
                         r2 = r2, r2p = r2p_val,
                         stab_TC = stab_TC, stab_NTC = stab_NTC,
                         vifp_TC = vif_path_TC, vifp_NTC = vif_path_NTC,
                         dR2_TC = dR2_TC, dR2_NTC = dR2_NTC,
                         coef_z = coef_z,
                         intercept_raw = intercept_z * sd_y + mean(yraw, na.rm = TRUE)
                       )
                       
                     }, error = function(e) {
                       list(i = i, error = e$message)
                     })
                   }

stopCluster(cl)

# ====================== Collect results ======================
for(res in results) {
  if(is.null(res)) next
  if(!is.null(res$error)) {
    cat(sprintf("Error at pixel %d: %s\n", res$i, res$error))
    next
  }
  
  i <- res$i
  
  coef_mat[i, ] <- res$coef
  vif_mat[i, ]  <- res$vif
  y_sim_mat[i, ] <- res$ysim
  r_mat[i] <- res$r
  p_mat[i] <- res$p
  r2_mat[i] <- res$r2
  r2p_mat[i] <- res$r2p
  
  beta_stab_TC[i]  <- res$stab_TC
  beta_stab_NTC[i] <- res$stab_NTC
  vif_path_TC[i]   <- res$vifp_TC
  vif_path_NTC[i]  <- res$vifp_NTC
  deltaR2_TC[i]    <- res$dR2_TC
  deltaR2_NTC[i]   <- res$dR2_NTC
}

# ====================== Residual (Z-score) ======================
y_sim_z <- zscore_rows(y_sim_mat)
residual_z <- y_sim_z - y_z

# ====================== 汇总统计 ======================
cat("岭回归结果汇总:\n")
cat(sprintf("总像元数: %d\n", ncell))
cat(sprintf("有效像元数: %d\n", length(valid_idx)))
cat(sprintf("平均R2: %.3f\n", mean(r2_mat, na.rm = TRUE)))

# ====================== Write rasters ======================
template <- veg_stack[[1]]
outdir <- "outputs_ridge"
dir.create(outdir, showWarnings = FALSE)

writeRaster(stack(lapply(1:7, function(k) setValues(raster(template), coef_mat[,k]))),
            file.path(outdir, "coef_stack.tif"), overwrite = TRUE)

for(k in 1:7){
  writeRaster(setValues(raster(template), coef_mat[,k]),
              file.path(outdir, paste0("coef_var", k, ".tif")), overwrite = TRUE)
  writeRaster(setValues(raster(template), vif_mat[,k]),
              file.path(outdir, paste0("vif_var", k, ".tif")), overwrite = TRUE)
}

writeRaster(stack(lapply(1:nt, function(t) setValues(raster(template), y_sim_mat[,t]))),
            file.path(outdir, "y_sim_stack.tif"), overwrite = TRUE)

writeRaster(stack(lapply(1:nt, function(t) setValues(raster(template), residual_z[,t]))),
            file.path(outdir, "residual_z_stack.tif"), overwrite = TRUE)

writeRaster(setValues(raster(template), r_mat),  file.path(outdir, "pearson_r.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), p_mat),  file.path(outdir, "pearson_p.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), r2_mat), file.path(outdir, "R2.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), r2p_mat), file.path(outdir, "R2_pvalue.tif"), overwrite = TRUE)

writeRaster(setValues(raster(template), beta_stab_TC),
            file.path(outdir, "beta_stability_TC.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), beta_stab_NTC),
            file.path(outdir, "beta_stability_NTC.tif"), overwrite = TRUE)

writeRaster(setValues(raster(template), vif_path_TC),
            file.path(outdir, "vif_path_TC.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), vif_path_NTC),
            file.path(outdir, "vif_path_NTC.tif"), overwrite = TRUE)

writeRaster(setValues(raster(template), deltaR2_TC),
            file.path(outdir, "deltaR2_TC.tif"), overwrite = TRUE)
writeRaster(setValues(raster(template), deltaR2_NTC),
            file.path(outdir, "deltaR2_NTC.tif"), overwrite = TRUE)

# 输出统计表
summary_stats <- data.frame(
  Variable = c("TC_BG", "NTC_BG", "AI", "Ta", "CO2", "SR", "WS"),
  Mean_Coefficient = colMeans(coef_mat[valid_idx, ], na.rm = TRUE),
  SD_Coefficient = apply(coef_mat[valid_idx, ], 2, sd, na.rm = TRUE),
  Mean_VIF = colMeans(vif_mat[valid_idx, ], na.rm = TRUE)
)
write.csv(summary_stats, file.path(outdir, "summary_statistics.csv"), row.names = FALSE)

# 系数分布图
png(file.path(outdir, "coefficient_distributions.png"), width = 1200, height = 800)
par(mfrow = c(2, 4))
for(k in 1:7) {
  hist(coef_mat[valid_idx, k], main = paste("Coefficient for Var", k),
       xlab = "Coefficient value", col = "lightblue", breaks = 30)
  abline(v = 0, lwd = 2)
}
dev.off()

cat("\n✅ 岭回归全部完成！\n")
cat(sprintf("输出保存在: %s\n", outdir))
