mafs <- function(maf, OR, baseline = 0.01) {
  if(length(OR) == 1 & length(maf) > 1)
    OR <- rep(OR, length(maf))
  if(length(maf) == 1 & length(OR) > 1)
    maf <- rep(maf, length(OR))
  # calcul p_tem : eq du second degré
  alpha <- (OR-1)*(1-baseline)
  beta  <- 1+(OR-1)*(baseline - maf)
  p.tem <- ifelse(OR == 1, maf, (sqrt(beta**2 + 4*alpha*maf) - beta)/(2*alpha) )
  uu <- OR*p.tem/(1-p.tem)
  p.cas <- uu/(1+uu)
  list(maf.controls = p.tem, maf.cases = p.cas)
}
