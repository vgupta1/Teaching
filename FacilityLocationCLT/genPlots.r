library(ggplot2)
#Load up the facility and demand locations for the plot
demandLocs = read.csv(file="demandLocs.csv")
demandLocs = as.data.frame( t(demandLocs) )
demandLocs = cbind(demandLocs, "Type" = "Demand")
names(demandLocs) <- c("X", "Y", "Type")

facLocs = read.csv(file="facilityLocs.csv")
facLocs = as.data.frame( t(facLocs))
facLocs = cbind(facLocs, "Type"="Facility")
names(facLocs) <- c("X", "Y", "Type")

locs = data.frame(rbind(demandLocs, facLocs))
ggplot(aes(x=X, y=Y, color=Type), data=locs) + 
  geom_point(size=4) + 
  theme_minimal(base_size=18) + 
  theme(legend.position=c(.2, .9), 
        legend.title=element_blank()) + 
  xlab("") + ylab("")

####  Inspect Capacities so that they make sense
capacity = read.csv("capacities.csv")

#Histogram of Naive approach
ggplot(aes(x=-1 * value, group=Method, color=Method, fill=Method), 
       data=subset(profits.melt, Method %in% c("Nominal", "Naive"))) + 
  geom_histogram(alpha=.5, position="identity") + 
  xlab("Profits ($)") +
  theme_minimal(base_size=18) + 
  ylab("") + xlab("Profit") + 
  theme(legend.title=element_blank(), 
        legend.position=c(.2, .9))




#####
#Start with the easy histograms
library(reshape)
profits = read.csv("profits.csv")
profits = cbind(1:50, profits)
names(profits)<- c(c("Run", "Nominal", "Naive"), seq(1.5, 4, .25) )
profits.melt = melt(profits, id.vars="Run", variable_name="Method")

ggplot(aes(x=Method, y=-1*value), data=profits.melt) + 
  geom_boxplot(fill="blue", alpha=.5) + 
  theme_minimal(base_size=18) + 
  ylab("Profit ($)") + xlab("")




ggplot(aes(x=-1 * value, group=Method, color=Method, fill=Method), 
       data=subset(profits.melt, Method %in% c("Nominal", "Naive", "3"))) + 
  geom_histogram(alpha=.5, position="identity") + 
  xlab("Profits ($)") + 
  facet_grid(Method~.) + 
  theme_minimal(base_size=18) + 
  ylab("") + xlab("Profit") + 
  theme(legend.position="none")

