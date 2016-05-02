#Cross Validation for glasso and modified glasso along with the plots of the cross validation

rm(list=ls())
library(glasso) 
library(infotheo)
library(bnlearn)
library(igraph)

scoreCalculator = function(graphModelWi,X)
{
  AdjMatEigen = eigen(graphModelWi);
  
  LogDet = 0;
  for(value in AdjMatEigen$values){
    LogDet = LogDet+log(value);
  }
  S = t(X)%*%X;
  TraceSTheta = sum(diag(S%*%graphModelWi));
  
  graphScore = LogDet - TraceSTheta;
  return(graphScore);
}

crossValidator = function(X,rho,k)
{
  numberOfRows = nrow(X);
  dataIndices = c(1:numberOfRows);
  testingSize = as.integer(numberOfRows/k);
  iter = k;
  testStartingPoint = 1;
  scoreVector = c();
  while(iter > 0)
  {
    testIndexVector = c();
    for(i in 1:testingSize)
    {
      testIndexVector = append(testIndexVector,i + testStartingPoint -1)
    }
    trainIndexVector = setdiff(dataIndices,testIndexVector);
    trainX = X[trainIndexVector,];
    testX = X[testIndexVector,];
    GraphModelNormal = glasso(t(trainX)%*%trainX, 
                              rho, 
                              zero=NULL, 
                              thr=1.0e-4, 
                              maxit=10,  
                              approx=FALSE, 
                              penalize.diagonal=TRUE, 
                              start="cold", 
                              w.init=NULL,
                              wi.init=NULL, 
                              trace=FALSE);
    likelihoodScore = scoreCalculator(GraphModelNormal$wi,testX);
    testStartingPoint = testStartingPoint + testingSize;
    iter = iter-1;
    if(is.nan(likelihoodScore))
    {
      next
    }
    if(is.complex(likelihoodScore))
    {
      next
    }
    scoreVector = append(scoreVector,likelihoodScore);
  }
  len = length(scoreVector);
  totalSum = sum(scoreVector);
  averageScore = totalSum/len;
  return(averageScore);
}


X = read.csv("WorldMarketsIndices.csv");
X = data.matrix(X);
X = scale(X, center = TRUE, scale = TRUE);
rhoVector = c();
scoreVector = c();
rho = 0.001;
k=10;
bestScore = -Inf;
while (rho <= 0.1)
{
  likelihoodScore = crossValidator(X,rho,k)
  z = paste("rho = ",rho);
  v = paste("likelihood = ", likelihoodScore);
  print(paste(z,v));
  rho = rho + 0.001;
  
  if(is.nan(likelihoodScore))
  {
    next
  }
  if(is.complex(likelihoodScore))
  {
    next
  }
  if(likelihoodScore > bestScore)
  {
    bestScore = likelihoodScore;
    bestRho = rho;
  }
  rhoVector = append(rhoVector,rho);
  scoreVector = append(scoreVector,likelihoodScore);
}