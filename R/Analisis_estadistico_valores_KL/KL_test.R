
# Script to check the statistical significance of the KL values ??????between the groups
# of replicas and those of non-replicas (that is, two samples that are not replicas)

setwd("C:/Users/yolib/Documents/TFM/Linf_T/Datos")

replica <-  read.csv("repl_norm.csv", sep = ",", header= FALSE)
no_replica <-  read.csv("no_replicas_norm.csv", sep = ",", , header= FALSE)

length(unlist(no_replica))
length(unlist(replica))

# I load the data of the replicas and the random data (no replicas)

# The data is unbalanced. 


t.test(as.numeric(unlist(replica)), as.numeric(unlist(no_replica)))
t.test(as.numeric(unlist(no_replica)), as.numeric(unlist(replica)))

wilcox.test(as.numeric(unlist(replica)),as.numeric(unlist(no_replica)), alternative = "two.sided") # we try the wilcox test




require(dplyr)
#random_no_replic <- sample_n(no_replica, size= length(replica[,1]))
random_no_replic <- sample_n(no_replica, size= 1000) # I will take 1000 random points to have the same
# number of samples in both groups.
random_replic <- sample_n(replica, size= 1000)

#length(random_no_replic)
#random_no_replic

hist(x = as.numeric(unlist(random_no_replic)))

qqnorm(as.numeric(unlist(random_no_replic)), xlab = "", ylab = "",
       main = "no_replica", col = "springgreen4")
qqline(as.numeric(unlist(random_no_replic)))


hist(x = as.numeric(unlist(random_replic)))

qqnorm(as.numeric(unlist(random_replic)), xlab = "", ylab = "",
       main = "replica", col = "springgreen4")
qqline(as.numeric(unlist(random_replic)))

require(car)
fligner.test(list(as.numeric(unlist(random_replic))), list(as.numeric(unlist(random_no_replic))))

# The replica data does not follow a normal distribution.
# I shouldn't do a t test, but since the number of samples is very large, it doesn't matter as much.
# Anyway, I'm going to also do a wilcox to compare.

t.test(as.numeric(unlist(random_replic)), as.numeric(unlist(random_no_replic)))

wilcox.test(as.numeric(unlist(random_replic)),as.numeric(unlist(random_no_replic)), alternative = "two.sided")

# With both tests I get that p < 2.2e-16

wilcox.test(as.numeric(unlist(random_replic)),as.numeric(unlist(random_no_replic)),
            exact = FALSE, alternative = "less")

bp <- boxplot(c(random_replic, random_no_replic), at = c(1,2), las = 2, names = c("replicas", "no replicas"))


bp_replic <- boxplot(replica)


replic_no_out <-replica[!(replica$V1 %in% bp_replic$out),]

qqnorm(as.numeric(unlist(replic_no_out)), xlab = "", ylab = "",
       main = "replica", col = "firebrick")
qqline(as.numeric(unlist(replic_no_out)))

