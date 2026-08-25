m20 = read.table(file = "scalability_test/result_n/20by20.txt")[1:20,]-2
m20 = cbind(rep(20, 20), m20)

m30 = read.table(file = "scalability_test/result_n/30by30.txt")[1:20,]
m30 = cbind(rep(30, 20), m30)

m40 = read.table(file = "scalability_test/result_n/40by40.txt")[1:20,]
m40 = cbind(rep(40, 20), m40)

m50 = read.table(file = "scalability_test/result_n/50by50.txt")[1:20,]
m50 = cbind(rep(50, 20), m50)

m60 = read.table(file = "scalability_test/result_n/60by60.txt")[1:20,]
m60 = cbind(rep(60, 20), m60)

m70 = read.table(file = "scalability_test/result_n/70by70.txt")[1:12,]
m70 = cbind(rep(70, 20), m70)

m80 = read.table(file = "scalability_test/result_n/80by80.txt")[1:20,]
m80 = cbind(rep(80, 20), m80)

m90 = read.table(file = "scalability_test/result_n/90by90.txt")[1:20,]
m90 = cbind(rep(90, 20), m90)

run_time = as.data.frame(rbind(m20, m30, m40, m50, m60, m70, m80, m90))
colnames(run_time) = c("m", "time")
run_time$time = as.numeric(run_time$time)*10
run_time$n = (run_time$m)^2

df = aggregate(run_time$time, list(run_time$m), FUN=mean)
colnames(df) = c("m", "time")
df$n = (df$m)^2



library(ggplot2)
p1 = ggplot() + geom_boxplot(data = run_time, aes(x = n, y = time, group = n)) + 
  geom_point(data = df, aes(x = n, y = time), color = "red", size = 3) + 
  geom_line(data = df, aes(x = n, y = time), linetype = 2) + theme_bw() + 
  xlab("Square lattice size") + ylab("Running time (seconds)") + 
  theme(panel.grid = element_blank()) + 
  scale_x_discrete(limit = c(20^2, 30^2, 40^2, 50^2, 60^2, 70^2, 80^2, 90^2),
                   labels = c("20X20", "30X30", "40X40", "50X50", "60X60", "70X70", "80X80", "90X90"))





test = data.frame(X= ((2:6))^2, Y = df$time[1:5])
model = lm(Y ~ X, data = test)
summary(model)




run_time = as.data.frame(rbind(m20, m30, m40, m50, m60))
colnames(run_time) = c("m", "time")
run_time$time = as.numeric(run_time$time)*10
run_time$n = (run_time$m)^2

df = aggregate(run_time$time, list(run_time$m), FUN=mean)
colnames(df) = c("m", "time")
df$n = (df$m)^2



library(ggplot2)
p2 = ggplot() + geom_boxplot(data = run_time, aes(x = n, y = time, group = n)) + 
  geom_point(data = df, aes(x = n, y = time), color = "red", size = 3) + 
  geom_line(data = df, aes(x = n, y = time), linetype = 2) + theme_bw() + 
  xlab("Square lattice size") + ylab("Running time (seconds)") + 
  theme(panel.grid = element_blank()) + 
  scale_x_continuous(breaks = c(20^2, 30^2, 40^2, 50^2, 60^2),
                     labels = c("20X20", "30X30", "40X40", "50X50", "60X60"))




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
ggsave(p, filename = "reproduce/img/FigureS11.jpg", width = 8, height = 3)

