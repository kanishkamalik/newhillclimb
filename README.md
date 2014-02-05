newhillclimb
============
twoinputmodel <- function(B) 
  #takes a budget and creates a 3D landscape for 2 inputs
  #uses 1000 different set of random coefficients between -1 & 1
{
  library(scatterplot3d)
  X3=NULL
  X4=NULL
  P2=NULL
  for (j in 1:1000) #1000 sets will be generated
  {
    c1=c(runif(5, -1, 1)) #vector of five randomly generated coefficients. 
    X1=NULL
    X2=NULL
    P=NULL
    for(i in 1:B)
    {
      X1=c(X1, i) #vector with all possible allocations of input 1
      X2=c(X2, B-i) # vector with all possible allocations of input 2
      X=c(i, B-i, (i^2), ((B-i)^2), (i*(B-i)))
      P=c(P, sum(c1*X)) #profit function
    }
    X3=c(X3, X1) #all the values of the input, X1, get stored in one vector
    X4=c(X4, X2) #all the values of the input, X2, get stored in one vector
    P2=c(P2, P) #all the profit values get stored in one vector. 
    #This allows us to print all the landscapes at once
  }
  scatterplot3d(X3, X4, P2) #3D graph
}

library(mvtnorm)
library(splines)
library(survival)
library(TH.data)

profitfunction <- function(Budget, N) 
  #function takes a budget and a number of inputs from the user
  #gives mean profits of all allocations for a 1000 different landscapes   
{
  library(EMT)
  X=findVectors(N, Budget) 
  #creates every possible allocation that meets the budget constraint
  Profits=NULL 
  #stores the mean profits from each set of coefficients
  for (j in 1:100)  
    #loop for creating a 100 different set of coefficients
  {
    c1=c(runif(N, -1,1))   #randomly generated set of coefficients between -1 & 1
    c2=c(runif(N, -1,1))   #set of coefficients for the squares
    c3=c(runif(choose(N-1,2), -1, 1)) 
    #coefficients of the crossproducts
    #the number of elements of the cross prod is always choose(N-1,2)
    P.alloc=NULL 
    #stores profit from each set of allocations 
    #for the jth set of coefficients
    for (i in 1: nrow(X))
    {
      X1=X[i, ]
      X2=(X1^2)
      X3=NULL
      for (k in 1:(N-1)) 
      { 
        for (l in (k+1):N) 
        { 
          if(k!=l) { X3=c(X3, X1[k]*X1[l]) } 
        } #generates the crossproduct 
      } 
      P.alloc=c(P.alloc, sum(c1* X1)+sum(c2*X2)+sum(c3*X3)) 
      #profits from every possible allocation
    }
    Profits <- c(Profits, (mean(P.alloc))) 
    #vector with mean profits from each set of coefficients
  }
  print(Profits)
}

pf <- function(X1, c1, c2, c3) 
  #generates a profit value for a given allocation and set of coefficients; 
  #this is needed for the next functions
{
  N=length(X1)
  X3=NULL
  for(k in 1:(N-1))
  {
    for(l in (k+1):N)
    {if(k!=1)
    {X3=c(X3, X1[k]*X1[l])}
    }
  }
  profit=sum(c1*X1)+sum(c2*(X1^2))+sum(c3*X3)
  print(profit)
}

# now we test the three different hillclimbing methods
# we use coefficients between -1 & 1 that are randomly generated

steepascent <- function(Budget, N)
  #function that makes climbs to the most profitable neighborhood
  #only one connection between the input variables
{
  library(EMT)
  library(multcomp)
  X=findVectors(N, Budget) 
  Profit=NULL 
  #vector containing local maxima of each set of coefficients
  for (j in 1:1000)
  {
    c1=c(runif(N, -1,1))
    c2=c(runif(N, -1,1))
    c3=c(runif(choose(N-1,2), -1,1))
    for(i in 1:nrow(X))
    {
      XM=X[i,]
      X1=XM
      pf1 <- pf(X1, c1, c2, c3)
      nbd=c(pf1) 
      #stores profit values of the most profitable neighborhood
      #we move; first element is the profit at the current allocation
      for(t in 1:1000)
      {
        mat1 <- multcomp:::contrMat(1:length(X1), "Tukey")	
        #creates a special type of matrix, 
        #when added to the the current allocation, generates every neighborhood
        #1 is added to the first term and subtracted from the second term 
        mat2 <- -1*mat1	#inverts mat1 as 
        #1 is instead subtracted from the first term and added to the second term
        nbd1 <- c(apply(mat1, 1, function(i)
        {X1 + i}))
        nbd2 <- c(apply(mat2, 1, function(i)
        {X1 + i}))
        Prof=NULL 
        #vector with profit values of all neighborhoods of an allocation
        #from this we will select the most profitable point
        for(g in 1:ncol(nbd1)) #iterates through every neighborhood
        {
          pf2 <- pf(nbd1[,g], c1, c2, c3)
          pf3 <- pf(nbd2[,g], c1, c2, c3)
          if (pf2>pf3)
          {Prof=c(Prof, pf2)} 
          else if (pf2<pf3)
          {Prof=c(Prof, pf3)}
        }
        m = seq(1,length(Prof))
        m=i[Prof==max(Prof)] 
        #stores position of the most profitable column vector in the neighborhoods matrix
        if(max(Prof)>max(nbd)) 
          #checks if the neighborhood is more profitable than our current position
        {
          nbd=c(nbd, max(Prof))
          if(pf(nbd1[,m], c1, c2, c3)>=pf(nbd2[,m], c1, c2, c3)) 
          {X1=nbd1[,m]} 
          #this allows us to move from our current position the most profitable neighborhood
          else if(pf(nbd1[,m], c1, c2, c3)<pf(nbd2[,m], c1, c2, c3))
          {X1=nbd2[,m]}
        }
        else if(max(Prof)<max(nbd)) 
          #if the neighborhood is < current position, this is the local maxima
        {Profit=c(Profit, max(nbd))
         break}
      }
    }
  }
  print(mean(Profit)) #prints mean value of the local maxima 
}

medianascent <- function(Budget, N) 
  #climb to the neighborhood whose profit is the median for that set of neighborhoods
{
  library(EMT)
  library(multcomp)
  X=findVectors(N, Budget) 
  Profit=NULL 
  #vector containing local maxima of each set of coefficients
  for (j in 1:1000)
  {
    c1=c(runif(N, -1,1))
    c2=c(runif(N, -1,1))
    c3=c(runif(choose(N-1,2), -1,1))
    for(i in 1:nrow(X))
    {
      XM=X[i,]
      X1=XM
      pf1 <- pf(X1, c1, c2, c3)
      nbd=c(pf1) 
      #stores profit values of the most profitable neighborhood onto which we move
      #first element is the profit at the current allocation
      for(t in 1:1000)
      {
        mat1 <- multcomp:::contrMat(1:length(X1), "Tukey")	
        #creates a matrix which is added to the current allocation
        #generates every neighborhood where 1 is added to the 1st term and subtracted from the 2nd 
        mat2 <- -1*mat1	
        #inverts mat1 as 1 is now subtracted from the first term and added to the 2nd
        nbd1 <- c(apply(mat1, 1, function(i)
        {X1 + i}))
        nbd2 <- c(apply(mat2, 1, function(i)
        {X1 + i}))
        Prof=NULL 
        #vector with profit values of all neighborhoods of an allocation
        #from this we select the most profitable point
        for(g in 1:ncol(nbd1)) 
          #iterates through every neighborhood
        {
          pf2 <- pf(nbd1[,g], c1, c2, c3)
          pf3 <- pf(nbd2[,g], c1, c2, c3)
          if (pf2>pf3)
          {Prof=c(Prof, pf2)} 
          else if (pf2<pf3)
          {Prof=c(Prof, pf3)}
        }
        m = seq(1,length(Prof))
        m=i[Prof==max(Prof)] 
        #stores position of the most profitable column vector in the neighborhoods matrix
        if(median(Prof)>max(nbd)) 
          #checks if the neighborhood is more profitable than our current position
        {
          nbd=c(nbd, median(Prof))
          if(pf(nbd1[,m], c1, c2, c3)>=pf(nbd2[,m], c1, c2, c3)) 
          {X1=nbd1[,m]} 
          #this moves us from our current position to the most profitable neighborhood
          else if(pf(nbd1[,m], c1, c2, c3)<pf(nbd2[,m], c1, c2, c3))
          {X1=nbd2[,m]}
        }
        else if(median(Prof)<max(nbd)) 
          #if the neighborhood < current position, we have reached the local maxima
        {Profit=c(Profit, max(nbd))
         break}
      }
    }
  }
  print(mean(Profit)) #prints mean value of the local maxima 
}

leastascent <- function(Budget, N) 
  #climbs to the least profitable neighborhood
{
  library(EMT)
  library(multcomp)
  X=findVectors(N, Budget) 
  Profit=NULL 
  #vector containing local maxima of each set of coefficients
  for (j in 1:1000)
  {
    c1=c(runif(N, -1,1))
    c2=c(runif(N, -1,1))
    c3=c(runif(choose(N-1,2), -1,1))
    for(i in 1:nrow(X))
    {
      XM=X[i,]
      X1=XM
      pf1 <- pf(X1, c1, c2, c3)
      nbd=c(pf1) 
      #stores profit values of the most profitable neighborhood onto which we move
      #first element is the profit at the current allocation
      for(t in 1:1000)
      {
        mat1 <- multcomp:::contrMat(1:length(X1), "Tukey")	
        mat2 <- -1*mat1	
        nbd1 <- c(apply(mat1, 1, function(i)
        {X1 + i}))
        nbd2 <- c(apply(mat2, 1, function(i)
        {X1 + i}))
        Prof=NULL 
        for(g in 1:ncol(nbd1))
        {
          pf2 <- pf(nbd1[,g], c1, c2, c3)
          pf3 <- pf(nbd2[,g], c1, c2, c3)
          if (pf2>pf3)
          {Prof=c(Prof, pf2)} 
          else if (pf2<pf3)
          {Prof=c(Prof, pf3)}
        }
        m = seq(1,length(Prof))
        m=i[Prof==min(Prof)]
        if(min(Prof)>max(nbd))
        {
          nbd=c(nbd, min(Prof))
          if(pf(nbd1[,m], c1, c2, c3)>=pf(nbd2[,m], c1, c2, c3)) 
          {X1=nbd1[,m]}
          else if(pf(nbd1[,m], c1, c2, c3)<pf(nbd2[,m], c1, c2, c3))
          {X1=nbd2[,m]}
        }
        else if(min(Prof)<max(nbd))
        {Profit=c(Profit, max(nbd))
         break}
      }
    }
  }
  print(mean(Profit))
}

package.skeleton("newhillclimb")
