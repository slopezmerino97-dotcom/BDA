RMSECV = function(X,LV,split){
  # INPUT   X: Data matrix (samples x features)
  #         LV: Total number ofcomponents that will be calculated
  #         split: In how many groups the samples will be divided
  # OUTPUT
  #         RMSECV: Root Mean Squared Error of Cross Validation
  # 
  # AHB 18/11/2022
  
  X = as.matrix(X)     #data matrix
  size = dim(X)
  I = size[1]          # samples
  J = size[2]          # features
  NTest = ceiling(I/split)   # samples per group
  TESTSET = 1:(NTest*split)*NA  # initialize vector with NAs (can be longer than number of samples)
  TESTSET[1:I] = (1:I) # fill vector with 1:number of samples (except positions > number of samples)
  TEST = matrix(TESTSET,nrow = split) #fill matrix with 1:number of samples + NA (as many rows as groups)
  PRESS = 1:LV * 0    # Initiate PRESS vector (LVx1) with zeros
  RMSECV = PRESS      # Initiate RMSECV vector (LVx1) with zeros
  Xhat = X*0          # Initiate Xhat (IxJ) of same dimensions as X, filled with 0s
  XX = X              # copy of data matrix, to iteratively remove components
  
  for(lv in 1:LV){        # for each PC 
    for(s in 1:split){       # for each group
      Train = 1:I            # values 1:number of samples  
      Test = TEST[s,]        # the sample numbers of the s'th split
      Test = Test[!is.na(Test)]    # keep only non-NA samples in s'th split (must be at least one)
      Train = Train [! Train %in% Test]  # keep only samples in training set that are not in test set (all others)
      Xtrain = XX[Train,]    # original data for training set (selected samples)
      Xtest = XX[Test,]      # original data for test set (selected samples)

      # PCA on XTrain using SVD
      svd = svd(Xtrain);    # SVD on training set
      v = svd$v[,1]         # loadings of training set
      
      for (j in 1:J){
        # for each variable j; predict score of testset without using this variable
        # and then predict Xhat test set for this variable
        # This guarantees independence between Xhat and X.
        vj = v[j]           # loading of one variable
        v_j = v[-j]         # loading of all other variables
        #AHB        t_test_j = as.matrix(Xtest[,-j]) %*% as.vector(v_j) / sum(v_j^2)   # Xtest is not a matrix, if split > 1/2 number of samples
        if(length(Test)>1){ # test data predicted for j by all others
          t_test_j = as.matrix(Xtest[,-j]) %*% as.vector(v_j) / sum(v_j^2)  
        }else{
          t_test_j = t(as.matrix(Xtest[-j])) %*% as.vector(v_j) / sum(v_j^2)
        }
        Xhat[Test,j] = Xhat[Test,j] + t_test_j * vj   #fill Xhat with the predicted data times loading of the variable
                      
      } # end J
      
    } # end split 
    
    E = X - Xhat # difference between the real data and the predicted data
    PRESS[lv] = sum(as.vector(E^2))    # sum of squares of residuals
                                       # removed now: if LV > PRESS, it just gets appended, 
                                       # removed now: if LV < PRESS, the 0s stay
    RMSECV[lv] = sqrt(PRESS[lv]/(I*J))   # square root of average sum of squares per data point
    
    # Deflate XX with 1 component
    svd_all = svd(XX)    
    XX = XX - svd_all$d[1] * (svd_all$u[,1]%*%t(svd_all$v[,1]))
  } # end LV
 
  return(RMSECV)  
} # end function
