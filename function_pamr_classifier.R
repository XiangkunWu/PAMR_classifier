pamr_classifier <- function(train_file, test_file, subtype_file, output_file) {
  
  library(pamr)
  
  # Read in the train and test data
  train.exp <- read.csv(train_file, row.names = 1, check.names = FALSE)
  train.exp <- as.matrix(train.exp)
  test.exp <- read.csv(test_file, row.names = 1, check.names = FALSE)
  test.exp <- as.matrix(test.exp)
  
  # Read in the subtype information
  subtype <- read.csv(subtype_file, row.names = 1, check.names = FALSE)
  
  # Extract common samples and features
  comsam <- intersect(colnames(train.exp), rownames(subtype))
  comfea <- intersect(rownames(train.exp), rownames(test.exp))
  train.exp <- train.exp[comfea, comsam]
  test.exp <- test.exp[comfea, ]
  subtype <- subtype[comsam, , drop = FALSE]
  
  set.seed(11)
  # Train pam classifier
  indata <- t(scale(t(train.exp)))
  mylist <- list(x = indata, 
                 y = as.vector(subtype$CCCRC), geneid = comfea,
                 genenames = comfea, sep = "")
  model <- pamr.train(mylist)
  
  # Perform cross-validation and choose a threshold value
  set.seed(11)
  model.cv <- pamr.cv(model, mylist, nfold = 10)
  threshold_errors <- data.frame(threshold = model.cv$threshold,
                           Overall_error_rate = model.cv$error,
                           signature_size=model.cv$size,
                           stringsAsFactors = F)
  write.table(threshold_errors, file = "threshold_errors.txt", sep="\t", quote=F, col.names=T,row.names = F)
  Delta <- model.cv$threshold[which.min(model.cv$error)]
  
  # Plot centroids and save to a file
  pamr.plotcen(model, mylist, Delta)
  
  # Print a confusion table to evaluate the classifier's performance
  pamr.confusion(model.cv, Delta)
  
  # Plot probability of classification for each sample
  pamr.plotcvprob(model, mylist, Delta)
  
  # Plot the genes that contribute most to each centroid
  pamr.geneplot(model, mylist, Delta)
  
  # Write the centroid information to a file
  write.table(pamr.listgenes(model, mylist, threshold = Delta, genenames = TRUE), file = "centroid.txt",sep="\t", quote=F, col.names=T,row.names = F)
  
  # Use the trained classifier to predict the new data set and write the results to a file
  indata2 <- t(scale(t(test.exp)))
  pamr_test <- pamr.predict(model, indata2, threshold = Delta)
  pamr_CCCRC <- data.frame(ID = colnames(indata2),
                           CCCRC = pamr_test,
                           stringsAsFactors = FALSE)
  write.table(pamr_CCCRC, file = output_file, sep="\t", quote=F, col.names=T,row.names = F)
}
