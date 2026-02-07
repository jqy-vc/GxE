library(ggplot2)
library(viridis)

# ==========================================
# 1. 核心计算函数 (固定效应 Q 检验)
# ==========================================
calculate_q_test <- function(beta_hat, se) {
  w <- 1 / (se^2)
  beta_fe <- sum(w * beta_hat) / sum(w)
  q_stat <- sum(w * (beta_hat - beta_fe)^2)
  k <- length(beta_hat)
  p_val <- pchisq(q_stat, df = k - 1, lower.tail = FALSE)
  return(p_val)
}

# ==========================================
# 2. 模拟逻辑：非均衡样本量 + 随机异质性
# ==========================================
set.seed(456)
n_sim <- 1000
alpha <- 0.05
k_values <- c(3, 4, 5, 6, 7)
distance_steps <- seq(0, 1, length.out = 50)

all_results <- data.frame()

for (K in k_values) {
  message(paste("正在模拟 K =", K, "(样本量不平衡 + 随机位置)..."))
  
  for (target_dist in distance_steps) {
    rejections <- 0
    
    for (i in 1:n_sim) {
      # --- A. 随机生成不平衡的 SE (标准误) ---
      # 模拟样本量差异：SE 在 0.02 (大样本) 到 0.15 (小样本) 之间随机
      se_vec <- runif(K, min = 0.01, max = 0.2)
      
      # --- B. 随机生成异质性模式 ---
      if (target_dist == 0) {
        true_beta <- rep(0, K)
      } else {
        # 随机决定有几个组存在效应差异 (至少1组，最多K-1组)
        num_diff_groups <- sample(1:(K-1), 1)
        diff_idx <- sample(1:K, num_diff_groups)
        
        raw_beta <- rep(0, K)
        raw_beta[diff_idx] <- rnorm(num_diff_groups)
        
        # 同样进行中心化和距离缩放
        centered_beta <- raw_beta - mean(raw_beta)
        current_dist <- sqrt(sum(centered_beta^2))
        true_beta <- (centered_beta / current_dist) * target_dist
      }
      
      # --- C. 观测值生成与检验 ---
      # 每个组的噪声标准差由它随机生成的 se_vec 决定
      beta_hat <- rnorm(K, mean = true_beta, sd = se_vec)
      p_val <- calculate_q_test(beta_hat, se_vec)
      
      if (p_val < alpha) rejections <- rejections + 1
    }
    
    all_results <- rbind(all_results, data.frame(
      K = as.factor(K),
      distance = target_dist,
      power = rejections / n_sim
    ))
  }
}

# ==========================================
# 3. 可视化
# ==========================================


pdf("/data2/qiuyue/environment/L2_vs_power_multiK.pdf", width=8, height=6)
ggplot(all_results, aes(x = distance, y = power, color = K, group = K)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_smooth(method = "loess", se = FALSE, linetype = "dashed", size = 0.5) +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Power Analysis: Imbalanced Sample Size & Random Patterns",
    subtitle = "Standard errors (SE) are randomly assigned per group (0.01 - 0.2)",
    x = "Total Vector Distance (Heterogeneity Intensity)",
    y = "Statistical Power",
    color = "Number of Groups (K)"
  ) +
  geom_hline(yintercept = 0.05, linetype = "dotted")
dev.off()
