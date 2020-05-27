similarityMat <- function(priceMat, sharesMat, unique_tickers,unique_dates){
  market_cap = priceMat*sharesMat
  corrMatrix = cor(market_cap, use = "pairwise.complete.obs")
  return(corrMatrix)
}