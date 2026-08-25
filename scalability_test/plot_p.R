m500 = read.table(file = "scalability_test/result_p/p_500.txt")[1:20,]
m500 = cbind(rep(500, 20), m500)

m1000 = read.table(file = "scalability_test/result_p/p_1000.txt")[1:20,]
m1000 = cbind(rep(1000, 20), m1000)

m2000 = read.table(file = "scalability_test/result_p/p_2000.txt")[1:20,]
m2000 = cbind(rep(2000, 20), m2000)

m3000 = read.table(file = "scalability_test/result_p/p_3000.txt")[1:20,]
m3000 = cbind(rep(3000, 20), m3000)

m4000 = read.table(file = "scalability_test/result_p/p_4000.txt")[1:20,]
m4000 = cbind(rep(4000, 20), m4000)

m5000 = read.table(file = "scalability_test/result_p/p_5000.txt")[1:20,]
m5000 = cbind(rep(5000, 20), m5000)

m6000 = read.table(file = "scalability_test/result_p/p_6000.txt")[1:20,]
m6000 = cbind(rep(6000, 20), m6000)

m7000 = read.table(file = "scalability_test/result_p/p_7000.txt")[1:20,]
m7000 = cbind(rep(7000, 20), m7000)

run_time = as.data.frame(rbind(m500, m1000, m2000, m3000, m4000, m5000))
colnames(run_time) = c("p", "time")
run_time$time = as.numeric(run_time$time)*10


df = aggregate(run_time$time, list(run_time$p), FUN=mean)
colnames(df) = c("p", "time")


library(ggplot2)
p3 = ggplot() + geom_boxplot(data = run_time, aes(x = p, y = time, group = p)) + 
  geom_point(data = df, aes(x = p, y = time), color = "red", size = 3) + 
  geom_line(data = df, aes(x = p, y = time), linetype = 2) + theme_bw() + 
  xlab("Number of genes") + ylab("Running time (seconds)") +
  theme(panel.grid = element_blank()) 


test = data.frame(X= c(0.5, (1:5)), Y = df$time)
model = lm(Y ~ X, data = test)
summary(model)



library(cowplot)
p <- plot_grid(p2, p3, labels=c('(a)', '(b)'))

gg = ggdraw() + 
  draw_plot(p2, x=0, y=0, width = 1/2, height = 1, scale = 1) +
  draw_plot(p3, x=1/2, y=0, width = 1/2, height = 1, scale = 1) 

#plotLegend = plot_complete(500, legend.position = "bottom")
#legend = get_legend(plotLegend[[1]])

#plot_final = plot_grid(gg, legend, nrow = 2,rel_heights = c(9, 1))

#ggsave(gg, filename = "scalability_test/runtime.pdf", width = 8, height = 3)
ggsave(p, filename = "scalability_test/runtime.pdf", width = 8, height = 3)

