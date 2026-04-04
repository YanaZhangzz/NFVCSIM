rm(list=ls())
setwd("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code")
library(Matrix)
library(poweRlaw)
library(MultiRNG)
library(MASS)
library(mvtnorm)
library(KernSmooth)
library(VGAM)
library(methods)
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM_functions_rcpp.cpp")
source("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/net_type_functions.R")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/NFVCSIM_alpha_test_func_rcpp.cpp")
Rcpp::sourceCpp("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/power_alpha_analysis.cpp")

# 设置参数
M <- 100
B_bootstrap <- 500
n <- 100
max_points <- 50
alpha_level <- 0.05
effect_sizes <- c(0.2,0.4,0.6,0.8,1)
# effect_sizes = 0

# 初始化参数
D0 = matrix(0, 2, 2)
B0 = matrix(0, 2, 2)
Sige1 = diag(c(0.4, 0.6), nrow = 2)
Sige1[1,2] = 0.1
Sige1[2,1] = 0.1
Sige = Sige1

beta10 = c(1, -0.5, 0.8)
beta20 = c(0.5, 0.8, -1)
beta0 = cbind(beta10, beta20)
g0 = matrix(0, nrow = n, ncol = 4)

# 转换为arma格式
D0_arma <- as.matrix(D0)
B0_arma <- as.matrix(B0)
Sige_arma <- as.matrix(Sige)
g0_arma <- as.matrix(g0)
beta0_arma <- as.matrix(beta0)
load("~/Documents/NFVCSIM/Code_NFVCSIM/NFVCSIM_code/bootstrap_alpha_test_t3_0_10050_new")
bootstrap_test_0


t = seq(0,1,length=max_points)
alpha11 = matrix(t^2,nrow=1)
alpha12 = matrix((1-t)^2,nrow=1)
alpha21 = matrix(5*(t-0.5)^2,nrow=1)
alpha22 = matrix(t^0.5,nrow=1)

alpha11m = matrix(rep(mean(alpha11),max_points),nrow=1)
alpha12m = matrix(rep(mean(alpha12),max_points),nrow=1)
alpha21m = matrix(rep(mean(alpha21),max_points),nrow=1)
alpha22m = matrix(rep(mean(alpha22),max_points),nrow=1)
alpham = c(alpha11m,alpha21m,alpha12m,alpha22m)

# 运行Rcpp版本的Power分析
cat("开始Rcpp版本的Power分析...\n")
system.time({
  results_rcpp <- power_analysis_rcpp(M, B_bootstrap, n, max_points, alpha_level,
                                      effect_sizes, D0_arma, B0_arma, 
                                      Sige_arma, g0_arma, beta0_arma,bootstrap_test_0,
                                      alpham)
})

# save(results_rcpp,file="bootstrap_alpha_test_t3_5025_sizepower_results")
# 处理结果
power_results_rcpp <- as.data.frame(results_rcpp$power_results)
colnames(power_results_rcpp) <- c("EffectSize", "RejectionRate", "StandardError", 
                                  "Type", "CI_lower", "CI_upper", "MeanPValue")
power_results_rcpp$Type <- ifelse(power_results_rcpp$Type == 0, "Size", "Power")

# 显示结果
print(power_results_rcpp)

# 保存结果
save(results_rcpp, file = "power_analysis_rcpp_results.RData")
cat("分析完成！结果已保存。\n")


## 结果分析和可视化
analyze_power_results <- function(power_results) {
  cat("\n========== Power分析最终结果 ==========\n")
  print(power_results)
  
  # 绘制power曲线
  library(ggplot2)
  
  # 提取power部分（排除size）
  power_data <- power_results[power_results$Type == "Power", ]
  
  p <- ggplot(power_data, aes(x = EffectSize, y = RejectionRate)) +
    geom_line(size = 1.2, color = "blue") +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.02, color = "blue", alpha = 0.6) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", size = 1) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "green", size = 1) +
    labs(title = "FVCSIM模型的Power曲线",
         subtitle = paste(M, "次蒙特卡洛 ×", B_bootstrap, "次bootstrap"),
         x = "效应大小 (k)", y = "拒绝率") +
    theme_minimal() +
    ylim(0, 1) +
    annotate("text", x = max(power_data$EffectSize)*0.7, y = 0.1, 
             label = paste("Size =", round(power_results$RejectionRate[1], 4)), 
             color = "red", size = 5) +
    annotate("text", x = max(power_data$EffectSize)*0.7, y = 0.85, 
             label = "80% Power目标", color = "green", size = 5)
  
  print(p)
  
  return(p)
}

# 执行分析并绘图
final_plot <- analyze_power_results(power_results_rcpp)

# 保存完整结果
final_output <- list(
  parameters = list(
    M = M,
    B_bootstrap = B_bootstrap,
    n = n,
    max_points = max_points,
    alpha_level = alpha_level,
    effect_sizes = effect_sizes
  ),
  power_results = power_results,
  computation_info = list(
    total_bootstrap_iterations = M * B_bootstrap,
    total_model_fits = M * (1 + B_bootstrap)  # 原始 + bootstrap
  )
)

save(final_output, file = "complete_power_analysis.RData")

cat("\n========== 计算统计 ==========\n")
cat("总蒙特卡洛次数:", M, "\n")
cat("总bootstrap次数:", M * B_bootstrap, "\n")
cat("总模型拟合次数:", M * (1 + B_bootstrap), "\n")
cat("结果已保存到: complete_power_analysis.RData\n")

# 输出简要报告
cat("\n========== 简要报告 ==========\n")
for(i in 1:nrow(power_results)) {
  row <- power_results[i, ]
  cat(sprintf("效应大小 %.1f: %s = %.3f (%.1f%%)\n", 
              row$EffectSize, row$Type, row$RejectionRate, row$RejectionRate*100))
}