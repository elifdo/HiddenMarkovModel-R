#install.packages("ppclust")
library(ppclust)
#install.packages("protr")
library(protr)
#install.packages("stats")
library(stats)
#install.packages("cluster")
library(cluster)
#install.packages("caret")
library(caret)
#install.packages("fclust")
library(fclust)
#install.packages("kmed")
library(kmed)
#install.packages("aphid")
library(aphid)
#if (!require("devtools"))
# install.packages("devtools")
#devtools::install_github("shiny", "rstudio")
#install.packages("shiny")
library(shiny)



##getting and preparing data
mydata = read.table("hiv746data.txt")

#cleaved and non cleaved data sequences
data_cle <- as.character(mydata[mydata[,2]==1,1])
data_noncle <- as.character(mydata[mydata[,2]==0,1])

#function that seperates sequence into 8 bit long vectors of aa's
#eg. "AVKLCAVI" to "A","V","K","L","C","A","V","I"
seperate_string <- function(data){
  data_mat <- matrix(, nrow = length(data), ncol = 8)
  for(i in 1:length(data)){
    data_mat[i,] <- unlist(strsplit(data[i],split=""))
  }
  return(data_mat)
}
# seperating data to 90 percent training and 10 percent test data
data_cle_all <- seperate_string(data_cle)
set.seed(100)
indices <- sample(1:dim(data_cle_all)[1],round(dim(data_cle_all)[1]*0.9),replace = FALSE  )
data_cle_mat<- data_cle_all[indices,]
data_cle_test<-data_cle_all[-indices,]
set.seed(110)
data_noncle_all <- seperate_string(data_noncle)
indices <- sample(1:dim(data_noncle_all)[1],round(dim(data_noncle_all)[1]*0.9),replace = FALSE  )
data_noncle_mat<- data_noncle_all[indices,]
data_noncle_test<-data_noncle_all[-indices,]

##Creating 10 folds for cross validation
num_fold <- 10
set.seed(100)
folds_cle <-createFolds(1:dim(data_cle_mat)[1], k = num_fold, list = TRUE, returnTrain = TRUE)
set.seed(100)
folds_noncle <- createFolds(1:dim(data_noncle_mat)[1], k = num_fold, list = TRUE, returnTrain = TRUE)

#Aminoacid data to be used while preparing hidden states
#GETTING AMINOACID DATA
data(AAindex)
#SLICING ONLY THE FEATURE SUBSET
aadata <- AAindex[,7:26]
#NAMES OF THE FEATURES
names<- AAindex[,2]
#to get rid of one duplicate feature
duplicate_indices <- which(duplicated(names))
names <- names[-duplicate_indices]
#keeping names as numbers
numbernames <- 1:544
aadata <- aadata[-duplicate_indices,]
numbernames <- numbernames[-duplicate_indices]
row.names(aadata) <- names
aadata <- as.data.frame(aadata)
#to get rid of rows which contains NAs
row.has.na <- apply(aadata, 1, function(x){any(is.na(x))})
aadata <- aadata[as.vector(  which(row.has.na==0)),]
numbernames <- numbernames[as.vector(  which(row.has.na==0))]
names <- names[as.vector(  which(row.has.na==0))]
row.names(aadata) <- numbernames

######### HMM model ##########################################
## HMM producing function, with given data and hidden states
hmm_producing_fnc<- function(data,S ){
  N <- length(S)      #number of states
  M <- 20      #number of observations
  t <- 8       #sequence length
  R <- dim(data)[1]   #number of data sequences
  O <-c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  #initial Emission probabilities according to data
  E <- matrix(0, nrow = N, ncol = M)      #emission probabilities
  Y <- 0.9  #Probability that staying on the states that the aa belongs to
  Z <- 1-Y  #Probability that not belonging to the states that aa belongs to
  for(i in 1:N){
    aa_not_belonging_to_i <-  setdiff(O, S[[i]])
    aa_belonging_to_i <- S[[i]]
    for(j in 1:M){
      if(O[j]  %in% aa_belonging_to_i ){
        E[i,j] <- Y / length(aa_belonging_to_i)
      }
      else {
        E[i,j] <- Z / length(aa_not_belonging_to_i)
      }
    }
  }
  #Initial start Probabilities according to data
  PI <- rep(0,times=N)    #starting probabilities
  for(i in 1:N){
    for(j in 1:R){
      if(data[j,1] %in%  S[[i]]){
        PI[i] <- PI[i]+1
      }
    }
  }
  PI <- PI/ sum(PI)
  
  #Initial Transition Matrix according to data
  A <- matrix(0, nrow = N, ncol = N)      #transition matrix of the states
  for(i in 1:R){
    for(j in 1:7){
      for(k in 1:N){
        for(l in 1:N){
          if(  (  data[i,j] %in% S[[k]] )  &
               (data[i,j+1])%in%  S[[l]]   ){
            A[k,l] <- A[k,l]+1
          }
        }
      }
    }
  }
  for(i in 1:N){
    if(sum(A[i,]) !=0 ){sum <- sum(A[i,]) 
    for(j in 1:N){
      A[i,j] <- A[i,j]/sum
    }
    }
    else{sum <- 1 
    for(j in 1:N){
      A[i,j] <- sum/N
    }}
    
  }
  #### Training the HMM
  datavector <- c()
  for(i in 1:R){
    datavector <-append(datavector,data[i,])
  }
  #names of the states
  States <- c("begin")
  for(i in 1:N){
    States <- c(States,paste(c("gr",as.character(i)), collapse = ""))
  }
  ### Define the transition probability matrix with addition of begin state
  A <- cbind(rep(0,N), A)
  PI<- append(0,PI)
  A <- rbind(PI,A)
  dimnames(A) <- list(from = States, to = States)
  ### Define the emission probability matrix
  dimnames(E) <- list(states = States[-1], residues = O)
  ### Build the starting HMM object
  hmm_cle_start <- structure(list(A = A, E = E), class = "HMM")
  ### Trained hmm object
  hmm_cle <- train(hmm_cle_start, list(datavector), method = "BaumWelch",
                   deltaLL = 0.001,maxiter=1000)
  return(hmm_cle  )
}




#feature can be tree, kmeans, kmedoids or without(without any feature selection)
#state selection can be fkmmed(fuzzy k medoids),fkm(fuzzy k means),fkmgk(Gustafson and Kessel - like fuzzy k-means),single(from optimally connected hmm paper),multi(same paper)
modelfnc <- function(feature,feature_num, state_selection, state_num ){
  
  ###### feature selection given feature selection method and number of features to be used
  
  ###hierercihal clustering based on euclidean distance of all data
  if(feature=="tree"){
    fit <- hclust(as.dist(1-abs(cor(t(aadata)))), method="complete")
    groups <- cutree(fit, k=feature_num)  
    
    #feature selection, one random feature from each cluster
    selected_features <- c()
    for(i in 1:feature_num){
      gr <- groups[which(groups==i)]
      what <- as.integer(names(gr))
      set.seed(i)
      random_index <- sample(1:length(what), size=1, replace = FALSE, prob = NULL)
      selected_features <- c(selected_features,which(numbernames==what[random_index]))
    }
    selected_features <- sort(selected_features)
    feat_data = t(aadata[selected_features,])  
  }  
  
  ### k means clustering for feature selection
  if(feature=="kmeans"){
    clusters <- kmeans(as.dist(1-abs(cor(t(aadata)))),feature_num)
    groups<- clusters$cluster
    #feature selection, one random feature from each cluster
    selected_features <- c()
    for(i in 1:feature_num){
      gr <- groups[which(groups==i)]
      what <- as.integer(names(gr))
      set.seed(i)
      random_index <- sample(1:length(what), size=1, replace = FALSE, prob = NULL)
      selected_features <- c(selected_features,which(numbernames==what[random_index]))
    }
    selected_features <- sort(selected_features)
    feat_data = t(aadata[selected_features,])
  }
  
  ###without feature selection
  if(feature=="without"){feat_data = t(aadata)  }
  
  #### k medoids, choosing medoids as representatives of the clusters
  if(feature=="kmedoids"){
    result <- fastkmed(as.dist(1-abs(cor(t(aadata)))), ncluster = feature_num, iterate = 100)
    selected_features <- sort(result$medoid)  
    feat_data = t(aadata[selected_features,])
  }
  
  
  
  ######### state production with the reduced feature set, fuzzy clustering with 0.1 treshold###
  
  k=state_num  #number of states
  
  if(state_selection=="fkmmed"){
    clustering <-FKM.med (X=feat_data, k=k, m=1.5,  stand=1,seed=100)
    U <- clustering$U            #membership degrees matrix 
    S <- vector("list", k)
    for(j in 1:k){
      for(i in 1:20){
        if(U[i,j]>0.1  ){S[[j]]<- c(S[[j]], rownames(U)[i])   }
      }
    }
  }
  
  if(state_selection=="fkmgk"){
    clustering2 <-FKM.gk (X=feat_data, k=k, stand=1, seed=100)
    U <- clustering2$U
    S <- vector("list", k)
    for(j in 1:k){
      for(i in 1:20){
        if(U[i,j]>0.1  ){S[[j]]<- c(S[[j]], rownames(U)[i])   }
      }
    }
  }
  
  if(state_selection=="fkm"){
    clustering3 <- FKM (X=feat_data, k=k, stand=1, seed=100)
    U <- clustering3$U
    S <- vector("list", k)
    for(j in 1:k){
      for(i in 1:20){
        if(U[i,j]>0.1  ){S[[j]]<- c(S[[j]], rownames(U)[i])   }
      }
    }
  }
  
  
  
  if(state_selection=="multi"){
    ########## multi property states, optimally connected paper
    S <- list()
    S[[1]] <- c("L","I","V")
    S[[2]] <- c("V","C","G","A")
    S[[3]] <- c("A","S","T","G","C")
    S[[4]] <- c("N","C","S","T","D")
    S[[5]] <- c("E","D","R","K","H")
    S[[6]] <- c("K","C","Y","W","H")
    S[[7]] <- c("Y","F","W","H")
    S[[8]] <- c("M")
    S[[9]] <- c("Q")
    S[[10]] <- c("P")
  }
  
  if(state_selection=="single"){
    ######### Single property states
    S <- list()
    S[[1]] <- c("L","I","V")
    S[[2]] <- c("F","Y","W","H")
    S[[3]] <- c("L","I","V","M","C","F","Y","W","H","K","A","G")
    S[[4]] <- c("Y","W","H","K","R","D","E","C","S","T","N","Q")
    S[[5]] <- c("E","D","R","K","H")
    S[[6]] <- c("K","R","H")
    S[[7]] <- c("A","G","C","S","T")
    S[[8]] <- c("V","C","P","A","G","S","T","D","N")
  }
  
  
  TP <- rep(0,num_fold) #true positives for each fold
  TN <-rep(0,num_fold)  #true negative
  FP <- rep(0,num_fold) #false positive
  FN <- rep(0,num_fold) #false negative
  
  #Producing HMMs, comparing likelihoods, assigning sequences into classes and filling TP,TN,FP and FN
  
  for(i in 1:num_fold){
    hmm_cle <- hmm_producing_fnc(data_cle_mat[unlist(folds_cle[[i]]),],S)  #hmm trained with cleaved data
    hmm_noncle <- hmm_producing_fnc(data_noncle_mat[unlist(folds_noncle[[i]]),],S)  #hmm trained with noncleaved data
    for(j in 1:(dim(data_cle_mat)[1]-length(folds_cle[[i]])) ){
      
      logForwardProbabilities = forward(hmm_cle, data_cle_mat[-unlist(folds_cle[[i]]),] [j,] )
      forwardProbability_cle = exp(logForwardProbabilities[[1]])
      
      logForwardProbabilities = forward(hmm_noncle, data_cle_mat[-unlist(folds_cle[[i]]),][j,])
      forwardProbability_noncle = exp(logForwardProbabilities[[1]])
      
      ##comparing probabilities and assigning to classes
      if(forwardProbability_cle >  forwardProbability_noncle ){
        TP[i] <- TP[i]+1
      }
      else{  FN[i] <-FN[i] +1   }
    }
    for(j in 1:(dim(data_noncle_mat)[1]-length(folds_noncle[[i]]))  ){
      logForwardProbabilities = forward(hmm_cle, data_noncle_mat[-unlist(folds_noncle[[i]]),][j,])
      forwardProbability_cle = exp(logForwardProbabilities[[1]])
      logForwardProbabilities = forward(hmm_noncle, data_noncle_mat[-unlist(folds_noncle[[i]]),][j,])
      forwardProbability_noncle = exp(logForwardProbabilities[[1]])
      
      ##comparing probabilities and assigning to classes
      if(forwardProbability_cle >  forwardProbability_noncle ){
        FP[i] <- FP[i]+1
      }
      else{  TN[i] <-TN[i] +1   }
    }
  }
  
  #### precision
  precisio <- c()
  for(i in 1:num_fold){
    precisio[i] <- TP[i]/(TP[i] +FP[i] )
  }
  pre <-  mean(precisio)
  
  #### recall
  recal <- c()
  for(i in 1:num_fold){
    recal[i] <- TP[i] /(TP[i] +FN[i] )
  }
  rec <- mean(recal)
  
  #### f score
  fscore <- c()
  for(i in 1:num_fold){
    fscore[i] <- 2*((precisio[i] *recal[i])/(precisio[i]+recal[i] )  )
  }
  fsc <- mean(fscore )
  
  #### mcc
  mcc <- c()
  for(i in 1:num_fold){
    mcc[i] <- ((TP[i] *TN[i] )-(FP[i] *FN[i] )  ) /(sqrt((TP[i] +FP[i] )*(TP[i] +FN[i] )*(TN[i] +FP[i] )*(TN[i] +FN[i] )   )  )
  }
  mc <- mean(mcc)
  accuracy <- c()
  for(i in 1:num_fold){
    accuracy[i] <- (TP[i] +TN[i] )/(TP[i] +FP[i] +TN[i] + FN[i] )
  }
  acc <- mean(accuracy)
  
  #### saving results as a text file
  
  mylist <- list(TP=TP,TN=TN,FP=FP,FN=FN,"precision"=pre,"recall"=rec,"fscore"=fsc,"mcc"=mc, "acc"=acc)
  filename <- c()
  filename <- c(filename, paste(c(feature,",",as.character(feature_num),",",state_selection,",",as.character(state_num)), collapse = ""))
  capture.output(mylist, file=filename)
  
  return(list(hmm_cle=hmm_cle,hmm_noncle=hmm_noncle))
}

hmms= modelfnc(feature="kmedoids",feature_num=5, state_selection="multi", state_num=10)

for(i in c(30,35,40,45,50,55,60)){
  for(j in c(5,6,7,8,9,10)){
    modelfnc(feature="kmedoids",feature_num=i, state_selection="fkmmed", state_num=j)
  }
}
for(i in c(45,50,55,60)){
  for(j in c(5,6,7,8,9,10)){
    modelfnc(feature="kmedoids",feature_num=i, state_selection="fkm", state_num=j)
  }
}
for(i in c(30,35,40,45,50,55,60)){
  for(j in c(5,6,7,8,9,10)){
    modelfnc(feature="kmedoids",feature_num=i, state_selection="fkmgk", state_num=j)
  }
}

