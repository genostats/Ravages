require(oz)
set.seed(1)
x <- random.bed.matrix( rep(0.001,20), 200 )
Beta.M(x, factor(rep(1:2,each=100)), factor(rep("R1", 20)))

oz:::ex_Beta.M( x, factor(rep(1:2,each=100)), factor(rep("R1", 20)), groups = 1)
as.vector(which(rowSums(as.matrix(x)) == 0)) - 1

oz:::ex_Beta.M( x, factor(rep(1:2,each=100)), factor(rep(c("R1","R2"), each = 10)), groups = 1)
as.vector(which(rowSums(as.matrix(x[,1:10])) == 0)) - 1



C.ALPHA(x, factor(rep(1:2,each=100)), factor(rep("R1", 20)))


