library(ggplot2)

set.seed(42)
n_sim <- 5000  
se_val <- 0.1   
# 扩大范围到 1.0 以观察曲线完全饱和的过程
diff_range <- seq(0, 1.0, length.out = 40) 

results_k2 <- do.call(rbind, lapply(diff_range, function(d) {
  true_betas <- c(d/2, -d/2)
  dist_val <- sqrt(sum(true_betas^2)) 
  
  beta1_hats <- rnorm(n_sim, mean = true_betas[1], sd = se_val)
  beta2_hats <- rnorm(n_sim, mean = true_betas[2], sd = se_val)
  
  # Q 统计量计算
  q_stats <- (beta1_hats - beta2_hats)^2 / (2 * se_val^2)
  power <- mean(pchisq(q_stats, df = 1, lower.tail = FALSE) < 0.05)
  
  data.frame(diff = d, dist = dist_val, power = power)
}))

# 查找 Power 到达 80% 时的临界距离
critical_dist <- results_k2$dist[which.min(abs(results_k2$power - 0.8))]

# 绘图输出
p <- ggplot(results_k2, aes(x = dist, y = power)) +
  geom_line(color = "#e34a33", size = 1.2) + 
  geom_point(color = "#e34a33", alpha = 0.6) +
  # 5% I型错误线
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  # 80% Power 阈值线
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "blue") +
  annotate("text", x = 0.1, y = 0.83, label = "80% Power Threshold", color = "blue", size = 3) +
  theme_minimal(base_size = 14) +
  labs(title = paste("Power Analysis (K=2, SE =", se_val, ")"),
       subtitle = "High noise scenario: Assessment of minimum detectable heterogeneity",
       x = "Euclidean Distance (L2 Norm)",
       y = "Statistical Power")

ggsave("/data2/qiuyue/environment/L2_vs_power_K_2_SE_0.1.pdf", plot = p, width = 7, height = 5)
