constructFund <- function(rho, q, priceMat, sharesMat, unique_tickers, unique_dates){
  n = length(unique_tickers)
  d = length(unique_dates)
  c = c(as.vector(rho),rep(0,n))
  A = matrix(0, n^2+n+1, n^2+n)
  A[1,(n^2+1):(n^2+n)] = rep(1,n)
  for (i in 1:n){
    A[(i+1), (n*(i-1)+1):(n*i)] = rep(1,n)
  }
  A[(n+2):(n^2+n+1), 1:n^2] = diag(1,n^2)
  A[(n+2):(n^2+n+1), (n^2+1):(n^2+n)] = matrix(rep(diag(-1,n),n), nrow=n^2, byrow=T)
  b = c(q, rep(1,n), rep(0,n^2))
  dir = c(rep('=',(n+1)),rep('<=',n^2))
  s <- lp('max', c, A, dir, b, all.bin=TRUE)
  share_last = sharesMat[d,]
  price_last = priceMat[d,]
  mkt_cap = share_last*price_last
  stk_similar = matrix(0,n,n)
  for (i in 1:n){
    stk_similar[i,] = mkt_cap[i]*s$solution[(n*(i-1)+1):(n*i)]
  }
  weights = colSums(stk_similar)
  weights_adj = weights/sum(weights)
  return(weights_adj)
}