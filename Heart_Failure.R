#Setting the working directory

getwd()
setwd("~/Documents/Heart_Failure_Project")
getwd()

#Installing the required packages

#install.packages('tidyverse')
#install.packages('ggplot2')
#install.packages("ggpubr")
#install.packages("corrplot")
#install.packages("mltools")
#install.packages('fastDummies')
#install.packages('PerformanceAnalytics')
#install.packages('plotROC')

library('fastDummies')
library('mltools')
library('tidyverse')
library('ggplot2')
library("ggpubr")
library("corrplot")
library("plotROC")
library("ISLR")
library("SmartEDA")

#Reading the data
my_data <- read.csv('heart_failure_clinical_records_dataset.csv')

#Looking at the data
head(my_data)
tail(my_data)

#Looking for null values

sum(is.na(my_data))
mean(is.na(my_data))

#Looking at the summary of the data variables

summary(my_data)
sapply(my_data, sd)

#EXPLORATORY DATA ANALYSIS

#Numerical Variables Analysis

numsum <- ExpNumStat(my_data,by="A",gp=NULL,Qnt=seq(0,1,0.1),MesofShape=2,Outlier=TRUE,round=2,Nlim=10)
numsum

#Categorical Variable Analysis

catsum <- ExpCTable(my_data,Target=NULL,margin=1,clim=10,nlim=3,round=2,bin=NULL,per=T)
catsum


#VISUALIZATIONS

#Density plot (Univariate)

plot1 <- ExpNumViz(my_data,target=NULL, nlim=10,Page=c(3,2),sample=6)
plot1[[1]]

#Box plots for all categorical variables

plot2 <- ExpNumViz(my_data,target="DEATH_EVENT",type=1,nlim=3,fname=NULL,col=c("darkgreen", 'blue', 'red'),Page=c(2,2),sample=4)
plot2[[1]]

#Plotting histogram for the continous variables

attach(my_data)
plot(age, creatinine_phosphokinase, main="Age & Creatinine Phosphokinase",
     xlab="Age ", ylab="creatinine phosphokinase", pch=20, col="red", cex = 1)
plot(age, ejection_fraction, main = "Age and Ejection Fraction",
     xlab= "Age", ylab = "ejection_fraction", pch=21, col="blue", cex = 1)
plot(age, platelets, main = "Age and Platelets",
     xlab= "Age", ylab = "platelets", pch=17, col="blue", cex = 1)
plot(age, serum_sodium, main="Age & Serum Sodium",
     xlab="Age ", ylab="serum sodium", pch=16, col="blue", cex=1)
plot(serum_sodium, serum_creatinine, main="Serum Sodium & Serum Creatinine",
     xlab="Serum Sodium", ylab="Serum Creatinine", pch=24, col="blue", cex=1)
plot(serum_sodium, ejection_fraction, main="Serum Sodium & Ejection Fraction ",
     xlab="Serum Sodium", ylab="Ejection Fraction", pch=19, col="blue", cex=1)
plot(serum_creatinine, ejection_fraction, main="Serum Creatinine & Ejection Fraction",
     xlab="Serum Creatinine", ylab="Ejection Fraction", pch=24, col="blue", cex=1)

#The column named 'time' does not seems to make any sense for analysis,
#thus we dropped that column from the data frame

df = subset(my_data, select = -c(time))
df

#Correlation Analysis between variables to see how they are correlated.

cc = cor(df, method = "spearman")

#Plotting the correlation heat map

corrplot(cc, order = "hclust", hclust.method = "average", addrect = 4, tl.cex = 0.7)

chart.Correlation(my_data,
                  method="spearman",
                  histogram=TRUE,
                  pch=16)

#QQ PLot

options(width = 150)
CData = data1
qqp <- ExpOutQQ(CData,nlim=10,fname=NULL,Page=c(3,2),sample=6)
qqp[[1]]

#Compute two-samples Wilcoxon test  

###The Mann–Whitney U test (or Wilcoxon rank–sum test), applied to each feature in relation 
#  to the death event target, detects whether we can reject the null hypothesis that the distribution of the each feature for the groups of samples defined by death event are the same.

# For Age
res <- wilcox.test(age ~ DEATH_EVENT , data = data1,
                   exact = FALSE)
res

#The p-value is 0.0001668, which is less than 0.05 or close to zero indicating that Age is strongly related to death event


#For ejection fraction
res <- wilcox.test(ejection_fraction ~ DEATH_EVENT , data = data1,
                   exact = FALSE)
res

# The p-value is 7.368e-07, which is less than 0.05 or close to zero indicating that ejection_fraction is strongly related to death event

# For Sex

res <- wilcox.test(sex ~  DEATH_EVENT , data = data1,
                   exact = FALSE)
res

# The p-value is 0.9413, which is greater than 0.05 or close to 1 indicating that sex is not related to death event.

# For Serum creatinine

res <- wilcox.test(serum_creatinine ~  DEATH_EVENT, data = data1,
                   exact = FALSE)
res

# The p-value is 1.581e-10, which is less than 0.05 or close to zero indicating  that serum_creatinine is strongly related to death event.


# For Serum Sodium
res <- wilcox.test(serum_sodium ~ DEATH_EVENT, data = data1,
                   exact = FALSE)
res

# The p-value is 0.0002928, which is less than 0.05 or close to zero indicating that serum_sodium is strongly related to death event

#Creating the Rank Table for features using the Wilcoxon Rank Sum Test results

Variables <- c('Sr. Creatinine','Ejection Fraction','Age', 'Sr. Sodium', 'Sex')
P_Value <- c(0.000000000158, 0.0000007368, 0.0001668, 0.0002928, 0.9413)
Rank <- (c(1, 2, 3, 4, 5))
wilcoxrank <- data.frame(Variables, P_Value, Rank)
wilcoxrank

# Scatter plot to visualize the relationship between serum creatinine and ejection fraction, differentiated by the occurrence of a death event.

scatterplot(serum_creatinine ~ ejection_fraction | DEATH_EVENT, data=my_data,
            main="Serum Creatinine & Ejection Fraction ",
            xlab="Ejection Fraction", ylab="Serum Creatinine")


#DATA PREPROCESSING
#One-hot encoding the categorical variables

df <- dummy_cols(df, select_columns = c('anaemia', 'diabetes', 'high_blood_pressure', 'sex', 'smoking'), remove_selected_columns = TRUE)
df

#LOGISTIC REGRESSION

table(my_data$DEATH_EVENT)


# Creating Training Data
input_death <- df[which(df$DEATH_EVENT == 1), ]  # all 1's
input_alive <- df[which(df$DEATH_EVENT == 0), ]  # all 0's
set.seed(100)  # for repeatability of samples
input_death_training_rows <- sample(1:nrow(input_death), 0.7*nrow(input_death))  # 1's for training
input_alive_training_rows <- sample(1:nrow(input_alive), 0.7*nrow(input_alive))  # 0's for training. Pick as many 0's as 1's
training_death <- input_death[input_death_training_rows, ]  
training_alive <- input_alive[input_alive_training_rows, ]
trainingData <- rbind(training_death, training_alive)  # row bind the 1's and 0's 

# Creating Test Data
test_death <- input_death[-input_death_training_rows, ]
test_alive <- input_alive[-input_alive_training_rows, ]
testData <- rbind(test_death, test_alive)  # row bind the 1's and 0's 

model.final = glm(DEATH_EVENT ~ age + creatinine_phosphokinase + ejection_fraction +
                    serum_creatinine + serum_sodium + anaemia_0 + anaemia_1 + diabetes_0 + diabetes_1
                  + high_blood_pressure_0 + high_blood_pressure_1 + sex_0 + sex_1 + smoking_0 + smoking_1, 
                  data=trainingData,
                  family = binomial(link="logit"))

predicted <- plogis(predict(model.final, testData))

summary(model.final)

#MODEL EVALUATION
# Assuming 'predicted' contains your predicted values
roc_obj <- roc(testData$DEATH_EVENT, predicted)

# Plot ROC curve
plot(roc_obj)











