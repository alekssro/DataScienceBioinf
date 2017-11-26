library(tidyverse)
source("nature_theme.R")

#
# Function to make boxplots
#
boxplotgroups = function(pd) {
  plot1 = ggplot(pd, aes(x=group, y=bp))+
    geom_boxplot()+
    geom_jitter(width = 0.1, height = 0, aes(color=group))+
    scale_color_manual(guide="none",values = discrete_colors )+
    xlab("")+nature_theme()
  return(plot1)
}

# Function to make a histogram plot of groups
plotgroups = function(pd) {
  plot1 = ggplot(pd, aes(x=bp, fill=group))+
    geom_histogram(bins=100)+
    ylab("")+
    scale_fill_manual(guide="none", values = discrete_colors)+
    facet_wrap( ~ group, ncol=1, scales = "free_y")+nature_theme()
  return(plot1)
}


#### Data from an experiment.  ####

#
# Group 1 is the blood pressure of males with Aortic dissection 
# Group 2 is the blodd pressure of normal males
#
# Link: http://en.wikipedia.org/wiki/Aortic_dissection
#

d  = read.table("Aortic_dissection.tsv")
g1 = d$bp[d$group != "Normal males"]
g2 = d$bp[d$group == "Normal males"]

observed_difference = abs(mean(g1)-mean(g2))

#
# Plot observations as boxplot
#

plot1 = boxplotgroups(d)

ggsave(filename = "observed.data.boxplot.png", plot = plot1, width = 112, height = 63, units = "mm", dpi = 300)


#
# Plot observations as histograms 
#
plot1 = plotgroups(d)+ 
  ggtitle(paste("Observed difference", observed_difference))

ggsave(filename = "observed.data.histogram.png", plot = plot1, width = 112, height = 63, units = "mm", dpi = 300)


#### Parametric t test ####
t.test(g1,g2, var.equal=F)



#### Simulation ####


# Estimate mean assuming H0 is true - that two two groups come from the same population
mean_h0 = mean(d$bp)
sd_h0   = sd(d$bp)

# Make a histogram of the huge population
g3    = data.frame("bp"=rnorm(100000, mean=mean_h0, sd=sd_h0), "group"="Simulated H0")
pd    = rbind(d,g3)

plot1 = plotgroups(pd)+
  ggtitle("With simulated data added")

ggsave(filename = "hypothesis.simulations.png", plot = plot1, width = 80, height = 126, units = "mm", dpi = 300)

#
# Run simulations and get test statistics
#

set.seed(1)
myplots = list()
i=1
myplots[[i]] = plotgroups(d)+
  ggtitle(paste("The observed data", "\nDifference in means=", abs(mean(g1) - mean(g2))))

for (i in 2:4) {
  l1 = length(g1)
  l2 = length(g2)
  ng1   = rnorm(l1, mean=mean_h0, sd=sd_h0)
  ng2   = rnorm(l2, mean=mean_h0, sd=sd_h0)
  pd1 = data.frame(bp=ng1, "group"="Group 1")
  pd2 = data.frame(bp=ng2, "group"="Group 2")
  pd  = rbind(pd1,pd2)
  diff  = abs(mean(ng1)-mean(ng2))
  myplots[[i]] = plotgroups(pd) + 
    ggtitle(paste("Simulation ",i-1, "\nDifference in means=", format(diff,2)))
}

library(gridExtra)
plot1 = arrangeGrob(grobs = myplots, ncol=2)    # Layout 6 plots
ggsave(filename = "simulation.multiplot.png", plot=plot1, width = 223, height = 123, units = "mm", dpi = 300)


#
# Simulate h0
#

set.seed(1)
ts = rep(0,1000000)
l1 = length(g1)
l2 = length(g2)
for (i in 1:length(ts)) {
  ng1   = rnorm(l1, mean=mean_h0, sd=sd_h0)
  ng2   = rnorm(l2, mean=mean_h0, sd=sd_h0)
  ts[i] = abs(mean(ng1)-mean(ng2))
  if(i %% 10000==0) { 
    cat(i,"\n")
    flush.console()
  }
}

#
# Now estimate p value
#
extreme_values      = sum(ts >= observed_difference)
pval_estimated      = (1+extreme_values) / (length(ts) + 1)


#
# Plotting results of simulation
#
pd = data.frame("ts"=ts, "i"=1:length(ts))


getplot = function(pd,rows) {
  obs = sum(observed_difference < pd$ts[1:rows])
  p = (obs + 1) / (rows + 1)
  binwidth = 0.5
  plot1 = ggplot(pd[1:rows,], aes(x = ts, color = ts >= observed_difference, fill=ts >= observed_difference))+ 
    geom_dotplot(binwidth=binwidth)+  
    geom_vline(xintercept = observed_difference, lty=2)+
    ggtitle(paste("The first ",rows, " test values\n",
                  "Number of observations > ", observed_difference, " = ", obs ,"\n",
                  "Estimated p < ", p, sep="")
    ) + 
    coord_cartesian(xlim = c(0, ceiling(observed_difference*1.10)))+
    ylab("")+
    scale_color_manual(guide="none", values = discrete_colors)+
    scale_fill_manual(guide="none", values=discrete_colors)+
    nature_theme()
  return(plot1)
}


# Stacked dots plot
plot1 = getplot(pd,4)
plot2 = getplot(pd,49)
plot3 = getplot(pd,1999)

plot4 = ggplot(pd, aes(x=ts, fill = ts >= observed_difference))+
  geom_histogram(binwidth=0.1)+
  geom_vline(xintercept = observed_difference, lty=2)+
  ylab("")+
  ggtitle(paste("Number of extreme values: ", extreme_values," of ", length(ts), "\n",
                "Simulated H0 p < ", pval_estimated, sep=""))+
  nature_theme()+
  scale_fill_manual(guide="none", values = discrete_colors)
  

f1 = arrangeGrob(grobs = list(plot1,plot2,plot3,plot4), ncol=2)
ggsave(filename = "simulation.pvalues.multiplot.png", plot=f1, width = 223, height = 112, units = "mm", dpi = 300)


#### Permutation ####

# Plot the estimated population using the permuted data

# Make a histogram of the huge population 
g3    = data.frame("bp"= sample(x = d$bp, size = 100000, replace = T) , "group"="Permuted ")
pd    = rbind(d,g3)

plot1 = plotgroups(pd)+
  ggtitle("With permuted data added")

ggsave(filename = "hypothesis.permutation.high.png", plot = plot1, width = 80, height = 126, units = "mm", dpi = 300)


set.seed(1)
myplots = list()
i=1
myplots[[i]] = plotgroups(d) + 
  ggtitle(paste("The observed data", "\nDifference in means=", abs(mean(g1) - mean(g2))))

for (i in 2:4) {
  permuted = sample(d$group)
  ng1      = d$bp[permuted != "Normal males"]
  ng2      = d$bp[permuted == "Normal males"]
  pd1 = data.frame(bp=ng1, "group"="Group 1")
  pd2 = data.frame(bp=ng2, "group"="Group 2")
  pd  = rbind(pd1,pd2)
  diff  = abs(mean(ng1)-mean(ng2))
  myplots[[i]] = plotgroups(pd) + 
  ggtitle(paste("Permutation",i-1, "\nDifference in means=", format(diff,2)))
}


plot1 = arrangeGrob(grobs = myplots, ncol=2)
ggsave(filename = "permutation.multiplot.png", plot=plot1, width = 223, height = 123, units = "mm", dpi = 300)


#
# Permutation of data to estimate distribution of test statistic under H0
#
set.seed(1)
ts = rep(NA, 1000000)

for (i in 1:length(ts)) {
  pmg = sample(d$group)                 # Mix samples between groups
  ng1 = d$bp[pmg != "Normal males"]     # The values from the sick group
  ng2 = d$bp[pmg == "Normal males"]     # The values from the normal group
  ts[i]    = abs(mean(ng1) - mean(ng2)) # Test statistic (H0 is true)
  if(i %% 10000==0) { 
    cat(i,"\n")
    flush.console()
  }
}

# Now estimate p value
extreme_values      = sum(ts >= observed_difference)
pval_estimated      = (1+extreme_values) / (length(ts) + 1)

pd = data.frame("ts"=ts, "i"=1:length(ts))

# Stacked dots plot
plot1 = getplot(pd, 4)
plot2 = getplot(pd, 49)
plot3 = getplot(pd, 1999)

plot4 = ggplot(pd, aes(x=ts, fill = ts >= observed_difference))+
  geom_histogram(binwidth=0.1)+
  geom_vline(xintercept = observed_difference, lty=2)+
  ylab("")+
  ggtitle(paste("Number of extreme values: ", extreme_values," of ", length(ts), "\n",
                "Permuted H0 p < ", pval_estimated, sep=""))+
  nature_theme()+
  scale_fill_manual(guide="none", values = discrete_colors)
  
f1 = arrangeGrob(grobs = list(plot1,plot2,plot3,plot4), ncol=2)

ggsave(filename = "permutation.pvalues.multiplot.png", plot=f1, width = 223, height = 112, units = "mm", dpi = 300)


