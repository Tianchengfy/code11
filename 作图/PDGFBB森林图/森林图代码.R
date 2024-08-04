library(ggplot2)
setwd("E://代谢物最新//作图//MR森林图炎症因子")

data <- read.csv("森林图示例数据.csv", header = TRUE)


head(data)


data$group <- factor(data$group, levels = rev(unique(data$group)))


p <- ggplot(data) +
  geom_hline(yintercept = 1, linewidth = 0.5) +
  geom_linerange(aes(x = group, y = med, ymin = min, ymax = max, color = p_col),
                 position = position_dodge(width = 0.5), show.legend = FALSE) +
  geom_point(aes(x = group, y = med, color = p_col), size = 1.5,
             position = position_dodge(width = 0.5)) +
  geom_text(aes(x = group, y = max + 0.05, label = p, color = p_col),
            position = position_dodge(width = 0.5), show.legend = FALSE) +
  scale_color_manual(name = "", values = c("MR Egger" = "#00CD00", "Weighted median" = "#0000FF", "Inverse variance weighted" = "#FF0000","Simple mode"="#A020F0","Weighted mode"="#00FFFF")) +
  scale_y_continuous(expand = c(0, 0.15)) +
  xlab("") +
  ylab("Odds Ratio") +
  theme_bw() +
  coord_flip()



# 保存图片为PDF格式
pdf("森林图1.pdf", width = 6.3, height = 10.9)  
print(p)
dev.off()
