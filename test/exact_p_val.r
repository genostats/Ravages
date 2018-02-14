require(oz)
set.seed(1)
x <- random.bed.matrix( rep(0.002,20), 200 )


Beta.M(x, factor(rep(1:2,each=100)), factor(rep("R1", 20)))
Beta.M.exact( x, factor(rep(1:2,each=100)), factor(rep("R1", 20)))

Beta.M(x, factor(rep(1:2,each=100)), factor(rep(c("R1","R2"), each = 10)))
Beta.M.exact( x, factor(rep(1:2,each=100)), factor(rep(c("R1","R2"), each = 10)))
Beta.M.exact( x, factor(rep(1:2,each=100)), factor(rep(c("R1","R2"), each = 10)), regions.to.test = "R1")


