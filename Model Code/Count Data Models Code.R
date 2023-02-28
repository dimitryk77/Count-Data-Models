
rm(list = ls())

library(stargazer)
library(pROC)
library(MASS)
library(Metrics)
library(writexl)
library(pscl)


#Reading in the data

setwd("~/Data Science/MSDS_410/Module_9/Assignment_9/")

my.data <- read.table('medical_care.txt')
colnames(my.data)<-c('ofp','ofnp','opp', 'opnp','emr','hosp','exclhlth','poorhlth','numchron','adldiff',
                       'noreast','midwest','west','age','black','male','married','school','faminc',
                       'employed','privins', 'medicaid')
                       
drop.list <-c('ofnp','opp','opnp','emr','hosp')
my.data <-my.data[,!(names(my.data) %in% drop.list)]


str(my.data)
head(my.data)


#####################################################################
# Train Test Set Creation
#####################################################################
set.seed(789)
my.data$u <- runif(n=dim(my.data)[1],min=0,max=1);

train.df <- subset(my.data, u<0.50)
test.df <- subset(my.data, u>=0.50)

train.df<-subset(train.df,select=-c(u))
test.df<-subset(test.df,select=-c(u))

dim(my.data)[1]
dim(train.df)[1]
dim(test.df)[1]
dim(train.df)[1]+dim(test.df)[1]


######################
###Poisson Regression
######################

upper_poisson.glm <- glm(ofp ~ .,data=train.df, family="poisson");
summary(upper_poisson.glm)
#using backward variable selection with stepAIC
backward_poisson.glm <- stepAIC(object=upper_poisson.glm,direction=c('backward'));
summary(backward_poisson.glm)

#creating model output using the stargazer package
model.path <- 'C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\'
file.name <- 'backward_poisson.html';
stargazer(backward_poisson.glm, type=c('html'),out=paste(model.path,file.name,sep=''),
          title=c('Poisson Regression Model'),
          align=TRUE, digits=2, digits.extra=2, initial.zero=TRUE, intercept.bottom=FALSE)


#####################################
###Poisson Regression with dispersion
#####################################

#stepAIC does not work as it has no loglikelihood
#using result of poisson regression backward variable selection model

poisson_dispersion.glm <- glm(formula(backward_poisson.glm),data=train.df, family="quasipoisson")
summary(poisson_dispersion.glm)


#creating model output using the stargazer package
model.path <- 'C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\'
file.name <- 'poisson_dispersion.html';
stargazer(poisson_dispersion.glm, type=c('html'),out=paste(model.path,file.name,sep=''),
          title=c('Poisson Regression with Dispersion  Model'),
          align=TRUE, digits=2, digits.extra=2, initial.zero=TRUE, intercept.bottom=FALSE)


###############################
###Negative Binomial Regression
###############################

upper_negative_binomial.glm <- glm.nb(ofp ~ .,data=train.df);
summary(upper_negative_binomial.glm)

#using backward variable selection with stepAIC
backward_negative_binomial.glm <- stepAIC(object=upper_negative_binomial.glm,direction=c('backward'));
summary(backward_negative_binomial.glm)


#creating model output using the stargazer package
model.path <- 'C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\'
file.name <- 'backward_negative_binomial.html';
stargazer(backward_negative_binomial.glm, type=c('html'),out=paste(model.path,file.name,sep=''),
          title=c('Negative Binomial Regression Model'),
          align=TRUE, digits=2, digits.extra=2, initial.zero=TRUE, intercept.bottom=FALSE)


######################
###Hurdle Regression
######################


upper_hurdle.glm <- hurdle(ofp ~ ., dist="poisson",data=train.df);
summary(upper_hurdle.glm)

#using backward variable selection with stepAIC
backward_hurdle.glm <- stepAIC(object=upper_hurdle.glm,direction=c('backward'));
summary(backward_hurdle.glm)

#creating model output using the stargazer package
model.path <- 'C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\'
file.name <- 'backward_hurdle.html';
stargazer(backward_hurdle.glm, type=c('html'),out=paste(model.path,file.name,sep=''),
          title=c('Hurdle Regression Model'),
          align=TRUE, digits=2, digits.extra=2, initial.zero=TRUE, intercept.bottom=FALSE)




###########################
###Zero Inflated Regression
###########################

upper_zero_inflated.glm <- zeroinfl(ofp ~ ., dist="poisson",data=train.df);
summary(upper_zero_inflated.glm)

#using backward variable selection with stepAIC
backward_zero_inflated.glm <- stepAIC(object=upper_zero_inflated.glm,direction=c('backward'));
summary(backward_zero_inflated.glm)

#creating model output using the stargazer package
model.path <- 'C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\'
file.name <- 'backward_zero_inflated.html';
stargazer(backward_zero_inflated.glm, type=c('html'),out=paste(model.path,file.name,sep=''),
          title=c('Zero-Inflated Regression Model'),
          align=TRUE, digits=2, digits.extra=2, initial.zero=TRUE, intercept.bottom=FALSE)




########################
###In-sample statistics
########################



gof.ins<- function(model, dataset,  model_name, poisson_dispersion_model=FALSE ) {
  
  print(model_name)
  #creating empty dataframe to store metrics
  df<-data.frame(matrix(NA, nrow = 4, ncol = 2))
  #getting AIC and BIC scores from models
  
  df[1,1]="AIC"
  #if poisson_dispersion_model flag is true, return "NA" as the model does not have
  #an AIC score, otherwise return AIC score
  df[1,2]<-ifelse(poisson_dispersion_model==TRUE,"NA",round(AIC(model),2))
  
  df[2,1]="BIC"
  #if poisson_dispersion_model flag is true, return "NA" as the model does not have
  # BIC score, otherwise return BIC score
  df[2,2]<-ifelse(poisson_dispersion_model==TRUE,"NA",round(BIC(model),2))
  
  #computing fitted values from model
  fitted_values_in_sample<-fitted(model)
  
  #MSE, MAE and MAPE measures are from the Metrics library
  df[3,1]<-"MSE"
  df[3,2]<- format(round(mse(dataset$ofp, fitted_values_in_sample),4), scientific=FALSE)
  df[4,1]<-"MAE"
  df[4,2]<- format(round(mae(dataset$ofp, fitted_values_in_sample),4), scientific=FALSE)

  
  #printing dataframe
  print(df)
  
  #returning dataframe for later use to create a bigger dataframe with output for all models
  df
}


#Creating dataframe with all in-sample results will be called model_metrics_in_sample_df

#poisson model metrics
poisson_df_in_sample<-gof.ins(backward_poisson.glm, train.df, "poisson")
model_metrics_in_sample_df<-poisson_df_in_sample
colnames(model_metrics_in_sample_df)<-c("Metrics","Poisson")


#poisson_dispersion model metrics
poisson_dispersion_df_in_sample<-gof.ins(poisson_dispersion.glm, train.df, "poisson_dispersion",poisson_dispersion_model=TRUE)
model_metrics_in_sample_df["poisson_dispersion"]<-poisson_dispersion_df_in_sample$X2


#negative_binomial model metrics
negative_binomial_df_in_sample<-gof.ins(backward_negative_binomial.glm, train.df, "negative_binomial")
model_metrics_in_sample_df["negative_binomial"]<-negative_binomial_df_in_sample$X2


#hurdle model metrics
hurdle_df_in_sample<-gof.ins(backward_hurdle.glm , train.df, "hurdle")
model_metrics_in_sample_df["hurdle"]<-hurdle_df_in_sample$X2

#zero_inflated model metrics
zero_inflated_df_in_sample<-gof.ins(backward_zero_inflated.glm , train.df, "backwards poisson")
model_metrics_in_sample_df["zero_inflated"]<-zero_inflated_df_in_sample$X2

#taking transpose of the metrics and adding in the model names as the column (R kept erasing it
#when it was the index)
model_metrics_in_sample_df_xls<-as.data.frame(t(model_metrics_in_sample_df))
model_metrics_in_sample_df_xls<- cbind(newColName = rownames(model_metrics_in_sample_df_xls), model_metrics_in_sample_df_xls)

#writing dataframe to excel
write_xlsx(as.data.frame(model_metrics_in_sample_df_xls, row.names=colnames(model_metrics_in_sample_df)),
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\in_sample_metrics.xlsx")


######################
###Out-sample stats
######################


gof.oos<- function(model, dataset,  model_name ) {
  
  print(model_name)
  #creating empty dataframe to store metrics
  df<-data.frame(matrix(NA, nrow = 2, ncol = 2))

  
  #computing predicted values on the test set from model
  predicted_values_out_of_sample<-predict(model, newdata=dataset,type="response")
  
  #MSE, MAE and MAPE measures are from the Metrics library
  df[1,1]<-"MSE"
  df[1,2]<- format(round(mse(dataset$ofp, predicted_values_out_of_sample),4), scientific=FALSE)
  df[2,1]<-"MAE"
  df[2,2]<- format(round(mae(dataset$ofp, predicted_values_out_of_sample),4), scientific=FALSE)

  
  #printing dataframe
  print(df)
  
  #returning dataframe for later use to create a bigger dataframe with output for all models
  df
}


#Creating dataframe with all in-sample results will be called model_metrics_out_of_sample_df


#poisson model metrics
poisson_df_out_of_sample<-gof.oos(backward_poisson.glm, test.df, "poisson")
model_metrics_out_of_sample_df<-poisson_df_out_of_sample
colnames(model_metrics_out_of_sample_df)<-c("Metrics","Poisson")



#poisson_dispersion model metrics
poisson_dispersion_df_out_of_sample<-gof.oos(poisson_dispersion.glm, test.df, "poisson_dispersion")
model_metrics_out_of_sample_df["poisson_dispersion"]<-poisson_dispersion_df_out_of_sample$X2


#negative_binomial model metrics
negative_binomial_df_out_of_sample<-gof.oos(backward_negative_binomial.glm, test.df, "negative_binomial")
model_metrics_out_of_sample_df["negative_binomial"]<-negative_binomial_df_out_of_sample$X2


#hurdle model metrics
hurdle_df_out_of_sample<-gof.oos(backward_hurdle.glm , test.df, "hurdle")
model_metrics_out_of_sample_df["hurdle"]<-hurdle_df_out_of_sample$X2

#zero_inflated model metrics
zero_inflated_df_out_of_sample<-gof.oos(backward_zero_inflated.glm , test.df, "backwards poisson")
model_metrics_out_of_sample_df["zero_inflated"]<-zero_inflated_df_out_of_sample$X2

#taking transpose of the metrics and adding in the model names as the column 
model_metrics_out_of_sample_df_xls<-as.data.frame(t(model_metrics_out_of_sample_df))
model_metrics_out_of_sample_df_xls<- cbind(newColName = rownames(model_metrics_out_of_sample_df_xls), model_metrics_out_of_sample_df_xls)

#writing dataframe to excel
write_xlsx(as.data.frame(model_metrics_out_of_sample_df_xls, row.names=colnames(model_metrics_out_of_sample_df)),
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\out_of_sample_metrics.xlsx")



##############################################
###Accuracy and confusion matrices-in-sample
##############################################

conf_matrix_in_sample<- function(model, dataset,  model_name ) {
  
  print(model_name)
  
  #computing fitted values from model and rounding integers
  fitted_values_is<-round(fitted(model),0)
  
  fitted_values_cat_is<-ifelse(fitted_values_is<6,1,ifelse(fitted_values_is<11,2,3))
  actual_values_cat_is<-ifelse(dataset$ofp<6,1,ifelse(dataset$ofp<11,2,3))
  
  # confusion matrix
  t <- table(actual_values_cat_is, fitted_values_cat_is);
  # Compute row totals;
  r <- apply(t,MARGIN=1,FUN=sum);
  # Normalize confusion matrix to rates;
  print(round(t/r,4)) 
  
  
  print(accuracy(actual_values_cat_is, fitted_values_cat_is))
  table_is<-round(t/r,4)
}


#writing confusion matrices to excel

c_m_backward_poisson_ins <- as.data.frame(unclass(conf_matrix_in_sample(backward_poisson.glm, train.df, "poisson")))
write_xlsx(c_m_backward_poisson_ins,
  "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_poisson_ins.xlsx")

c_m_poisson_dispersion_ins <-as.data.frame(unclass(conf_matrix_in_sample(poisson_dispersion.glm, train.df, "poisson dispersion")))
write_xlsx(c_m_poisson_dispersion_ins,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_poisson_dispersion_ins.xlsx")



c_m_backward_negative_binomial_ins <-as.data.frame(unclass(conf_matrix_in_sample(backward_negative_binomial.glm, train.df, "negative_binomial")))
write_xlsx(c_m_backward_negative_binomial_ins,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_negative_binomial_ins.xlsx")


c_m_backward_hurdle_ins <-as.data.frame(unclass(conf_matrix_in_sample(backward_hurdle.glm, train.df, "hurdle")))
write_xlsx(c_m_backward_hurdle_ins ,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_hurdle_ins.xlsx")


c_m_backward_zero_inflated_ins <-as.data.frame(unclass(conf_matrix_in_sample(backward_zero_inflated.glm, train.df, "zero_inflated")))
write_xlsx(c_m_backward_zero_inflated_ins ,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_zero_inflated_ins.xlsx")

###################################################
###Accuracy and confusion matrices - Out-of-sample
###################################################

conf_matrix_out_of_sample<- function(model, dataset,  model_name ) {
  
  print(model_name)
  
  #computing fitted values from model and rounding integers
  predicted_values_oos<-round(predict(model, newdata=dataset,type="response"),0)
  
  predicted_values_cat_oos<-ifelse(predicted_values_oos<6,1,ifelse(predicted_values_oos<11,2,3))
  actual_values_cat_oos<-ifelse(dataset$ofp<6,1,ifelse(dataset$ofp<11,2,3))
  
  #  confusion matrix
  t <- table(actual_values_cat_oos, predicted_values_cat_oos);
  # Compute row totals;
  r <- apply(t,MARGIN=1,FUN=sum);
  # Normalize confusion matrix to rates;
  print(round(t/r,4)) 

  
  print(accuracy(actual_values_cat_oos, predicted_values_cat_oos))
  table_is<-round(t/r,4)
  
}

#writing confusion matrices to excel


c_m_backward_poisson_oos <- as.data.frame(unclass(conf_matrix_out_of_sample(backward_poisson.glm, test.df, "poisson")))
write_xlsx(c_m_backward_poisson_oos,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_poisson_oos.xlsx")



c_m_poisson_dispersion_oos <-as.data.frame(unclass(conf_matrix_out_of_sample(poisson_dispersion.glm, test.df, "poisson dispersion")))
write_xlsx(c_m_poisson_dispersion_oos,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_poisson_dispersion_oos.xlsx")


c_m_backward_negative_binomial_oos <-as.data.frame(unclass(conf_matrix_out_of_sample(backward_negative_binomial.glm, test.df, "negative_binomial")))
write_xlsx(c_m_backward_negative_binomial_oos,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_negative_binomial_oos.xlsx")


c_m_backward_hurdle_oos <-as.data.frame(unclass(conf_matrix_out_of_sample(backward_hurdle.glm, test.df, "hurdle")))
write_xlsx(c_m_backward_hurdle_oos ,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_hurdle_oos.xlsx")


c_m_backward_zero_inflated_oos <-as.data.frame(unclass(conf_matrix_out_of_sample(backward_zero_inflated.glm, test.df, "zero_inflated")))
write_xlsx(c_m_backward_zero_inflated_oos ,
           "C:\\Users\\dimit\\Documents\\Data Science\\MSDS_410\\Module_9\\Assignment_9\\confusion_matrix\\c_m_backward_zero_inflated_oos.xlsx")





