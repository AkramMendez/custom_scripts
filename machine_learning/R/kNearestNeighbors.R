
library(tidyverse)
library(data.table)
library(caret)



K-Nearest Neighbors Algorithm:
Supervised learning classification algorithm


# Loading reference dataset (COX-2 Activity Data) -------------------------
# COX-2 Activity Data
# From Sutherland, O’Brien, and Weaver (2003): A set of 467 cyclooxygenase-2 (COX-2) inhibitors has been assembled from the published work of a single research group, with in vitro activities against human recombinant enzyme expressed as IC50 values ranging from 1 nM to >100 uM (53 compounds have indeterminate IC50 values).
# 
# A set of 255 descriptors (MOE2D and QikProp) were generated. To classify the data, we used a cutoff of 2^{2.5} to determine activity.
# 
# Using data(cox2) exposes three R objects: cox2Descr is a data frame with the descriptor data, cox2IC50 is a numeric vector of IC50 assay values and cox2Class is a factor vector with the activity results.

data(cox2)
head(cox2Descr)

Splitting data as training and set datasets

set.seed(42)
# Create a subset dataframe to predict cox2 activity based on molecular weight or IC50 of each compound:

cox2DescrSubset <- data.frame(weight=cox2Descr$moe2D_Weight,ic50=cox2IC50, class=cox2Class) 

idxTrain <- createDataPartition(y=cox2DescrSubset$class, p=0.8, list =FALSE)
train_data <- cox2DescrSubset[idxTrain,]
test_data <- cox2DescrSubset[-idxTrain,]

# Verify the distribution of the training data and the partitioned data:
prop.table(table(train_data$class))



prop.table(table(test_data$class))


Pre-processing: Feature scaling
Scaling and centering training and testing sets

#Using base functions:
#train_scaled <- scale(train_data[,c("weight","ic50")],center = TRUE,scale = TRUE)
#test_scaled <- scale(test_data[,c("weight","ic50")],center = TRUE,scale = TRUE)
#train_scaled <- scale(train_data[,c("weight","ic50")],center = TRUE,scale = TRUE)
#test_scaled <- scale(test_data[,c("weight","ic50")],center = TRUE,scale = TRUE)

# Using caret functions:
preProc_cox2 <- preProcess(train_data[,c("weight","ic50")], method = c("center","scale"))
train_scaled <- predict(preProc_cox2,train_data)
test_scaled <- predict(preProc_cox2,test_data)



Hyper-parameter tuning using cross-validation for finding the optimal value of k neighbors to use for classification.

set.seed(42) # set random seed for reproduciblility

# control object for training the kNN model 
ctrl <- trainControl(method = "repeatedcv",
                     summaryFunction = defaultSummary,
                     classProbs = TRUE,
                     number = 10, # Number of folds in the cross-validation
                     repeats = 10) # Number of times to repeat the cross-validation

# Use the train() function to perform the model training/tuning of the k hyperparameter.
knn_cv <- caret::train(class ~ .,
                       data=train_data,
                       method="knn",
                       trControl = ctrl, 
                       metric = "Accuracy", # Performance metric
                       tuneGrid = data.frame(k = seq(1,30,by = 2))) # Grid for searching the optimal hyperparameter k

knn_cv

### Plot kNN model tunning

plot(knn_cv)



## Training KNN Classifier:
After finding the optimal value of k number of neighbours, we will train the kNN classification model using the training dataset:


#library(class)
# Select k according to the hyper-parameter tunning result (k with the largest accuracy)
k <- knn_cv$results$k[which(knn_cv$results$Accuracy==max(knn_cv$results$Accuracy))]

# Using the class:knn() function:
# test_pred <- class::knn( train = train_scaled,
#                          test = test_scaled,
#                          cl = train_data$class, #	factor of true classifications of training set
#                          k=k # Number of neighbours considered (selected optimal k in hyper-parameter tuning)
#                          )

# Using caret::train() function:
best_kNN_model <- knn3(class ~ .,
                       data = train_scaled,
                       k=k)





Model Evaluation
We can compare the predicted Cox2 activity classification labels (Active/Inactive) vs the original class labels from the Cox2 dataset. 
Generating a confusion matrix:

actual_class <- test_data$class
predicted_class <- predict(best_kNN_model,
                           test_scaled,
                           type = "class")
# Using base table() function:
#confusion_matrix_cox2 <- table(actual_class,predicted_class)

# Using caret::confusionMatrix() function:
confusion_matrix_cox2 <- caret::confusionMatrix(actual_class,predicted_class)
confusion_matrix_cox2

Summary of classification performance:

summary_knn <- data.frame(Accuracy = confusion_matrix_cox2$overall["Accuracy"],
           sensitivity= confusion_matrix_cox2$byClass["Sensitivity"],
           specificity= confusion_matrix_cox2$byClass["Specificity"],
           precision = confusion_matrix_cox2$byClass["Precision"],
           recall = confusion_matrix_cox2$byClass["Recall"]
           )
summary_knn


library(pROC)
library(hutils)

predictions_cox2 <- predict(knn_cv,
                           test_data[,c("weight","ic50")],
                           type = "prob")

# ROC curve:
ground_truth_int <- as.numeric(hutils::Switch(as.character(test_data$class),"Active"="1","Inactive"="0", DEFAULT = "0"))
roc_curve_knn <- roc(response=ground_truth_int,# Vector of responses (true class)
                     predictor=predictions_cox2[,"Active"]) # Predicted values for each observation

# Plot ROC curve
plot(roc_curve_knn, main = "ROC Curve for KNN Predictions", col = "steelblue")




