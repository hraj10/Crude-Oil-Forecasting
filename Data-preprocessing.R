library(tseries)
library(stats)
library(forecast)
library(readxl)


# Dependent variable: Spot Crude Oil Price: West Texas Intermediate (WTI) (WTISPLC)
# Available under https://fred.stlouisfed.org/series/WTISPLC

data<-read.csv("data/MCOILWTICO FINAL.csv")[-(397:402),]
inflation<-read.csv("data/CPIAUCSL.csv")[-(913:917),]
kilian_data<-as.data.frame(read_xlsx("data/kilian.xlsx")[-(661:666),])
stock_change<-as.data.frame(read_xlsx("data/Stock change.xlsx",col_names = FALSE))[-(253:257),]
production<-as.data.frame(read_xlsx("data/Global Production.xlsx",col_names = FALSE)[-(253:257),])

Data_Transform<-function(data,transform=FALSE){
  #input:   data: csv file found in data folder , i.e. "data/MCOILWTICO.csv"
  #         transform={FALSE (default): applies no transformation
  #                    logReturn: log(P_t)-log(P_{t-1}), 
  #                    difflogReturn:  log(P_t)-log(P_{t-1}) - (log(P_{t-1})-log(P_{t-2}))
  #                    }
  #
  #output;  y vector
  y_raw<-data[,-1]
  if (transform=="logReturn") {
    y<-diff(log(y_raw))
  } else if (transform=="difflogReturn"){
    y<-diff(diff(log(y_raw)))
  } else if (transform=="diff"){
    y<-diff(y_raw)
  } else {
    y<-y_raw
  }
  return(y)
}
#Data_Transform(data,"diff")

Data_Select_Period<-function(data, transform=FALSE, end){
  #input:   data: csv file found in data folder , i.e. "data/MCOILWTICO.csv"
  #         end: end data in the form (2023,end month)
  #
  #output:  y: a vector from January 2003 until (2023, end month)
  y_transformed<-Data_Transform(data,transform)
  n<-length(y_transformed)
  start_idx<-12*20  
  y_2003<-y_transformed[(n - start_idx+1):n]
  y<-ts(y_2003,frequency = 12,start=c(2003,01), end = end)
  return(y)
}



# Independent variable 
y_ts<-Data_Select_Period(data,"diff",end=c(2022,12)) #logreturn
# lag 1
data_lag1<-data[-dim(data)[1],]
y_lag1<-Data_Select_Period(data_lag1,"diff",end=c(2022,12))
data_lag2<-data_lag1[-dim(data_lag1)[1],]
y_lag2<-Data_Select_Period(data_lag2,"diff",end=c(2022,12))
y_lag3<-Data_Select_Period(data[-c((dim(data)[1]-2):dim(data)[1]),],"diff",end=c(2022,12))
y_lag6<-Data_Select_Period(data[-c((dim(data)[1]-5):dim(data)[1]),],"diff",end=c(2022,12))
y_lag12<-Data_Select_Period(data[-c((dim(data)[1]-11):dim(data)[1]),],"diff",end=c(2022,12))


# inflation
inflation_lag1<-inflation[-dim(inflation)[1],]
cpi<-Data_Select_Period(inflation,"logReturn",end=c(2022,12))
cpi_lag1<-Data_Select_Period(inflation_lag1,"logReturn",end=c(2022,12))
cpi_lag3<-Data_Select_Period(inflation[-c((dim(inflation)[1]-2):dim(inflation)[1]),],"logReturn",end=c(2022,12))
cpi_lag6<-Data_Select_Period(inflation[-c((dim(inflation)[1]-5):dim(inflation)[1]),],"logReturn",end=c(2022,12))
cpi_lag12<-Data_Select_Period(inflation[-c((dim(inflation)[1]-11):dim(inflation)[1]),],"logReturn",end=c(2022,12))


# kilian index
kilian_idx<-kilian_data[-dim(kilian_data)[1],]
kilian<-Data_Select_Period(kilian_data,"diff", end=c(2022,12))
kilian_lag1<-Data_Select_Period(kilian_idx,"diff", end=c(2022,12))
kilian_lag3<-Data_Select_Period(kilian_data[-c((dim(kilian_data)[1]-2):dim(kilian_data)[1]),],"diff",end=c(2022,12))
kilian_lag6<-Data_Select_Period(kilian_data[-c((dim(kilian_data)[1]-5):dim(kilian_data)[1]),],"diff",end=c(2022,12))
kilian_lag12<-Data_Select_Period(kilian_data[-c((dim(kilian_data)[1]-11):dim(kilian_data)[1]),],"diff",end=c(2022,12))


# Stock change
delta_stock<-Data_Select_Period(stock_change,end=c(2022,12))/1000
stock_change_lag1<-stock_change[-252,]
delta_stock_lag1<-Data_Select_Period(stock_change_lag1,end=c(2022,12))/1000
delta_stock_lag3<-Data_Select_Period(stock_change[-c((dim(stock_change)[1]-2):dim(stock_change)[1]),],end=c(2022,12))/1000
delta_stock_lag6<-Data_Select_Period(stock_change[-c((dim(stock_change)[1]-5):dim(stock_change)[1]),],end=c(2022,12))/1000
delta_stock_lag12<-Data_Select_Period(stock_change[-c((dim(stock_change)[1]-11):dim(stock_change)[1]),],end=c(2022,12))/1000

# Production
prod<-Data_Select_Period(production,"logReturn",end=c(2022,12))
production_lag1<-production[-dim(production)[1],]
prod_lag1<-Data_Select_Period(production_lag1,"logReturn",end=c(2022,12))
prod_lag3<-Data_Select_Period(production[-c((dim(production)[1]-2):dim(production)[1]),],"logReturn",end=c(2022,12))
prod_lag6<-Data_Select_Period(production[-c((dim(production)[1]-5):dim(production)[1]),],"logReturn",end=c(2022,12))
prod_lag12<-Data_Select_Period(production[-c((dim(production)[1]-11):dim(production)[1]),],"logReturn",end=c(2022,12))


#block<-cbind(y_lag1,y_lag2,y_lag12, pi_lag1, kilian_lag1, delta_stock, prod)
# Train test split

Train_Test_Split<-function(X=FALSE,y,split){
  #input:     X matrix with potential covariates (does not have to be included, default=FALSE)
  #           y vector that is partitioned
  #           split proportion of train/test set i.e. c(0.7,0.3)
  #
  #output:    y_train, y_test splitted vectors and X_train, X_test if specified
  
  n<-length(y)
  split_idx<-round(n * split[1])
  #second_split_idx<-round(n * (split[1]+split[2]))
  y_train <- window(y, start = c(2003,01), end = c(2003,split_idx))
  #y_val <- window(y, start = c(2003,(first_split_idx+1)), end = c(2003,second_split_idx))
  y_test <- window(y, start = c(2003,(split_idx+1)))
  if (!is.matrix(X)) {
    return(list(y_train=y_train,y_test=y_test))
  } else {
    X_train <- window(X, start = c(2003,01), end = c(2003,split_idx))
    #X_val <- window(X, start = c(2003,(first_split_idx+1)), end = c(2003,second_split_idx))
    X_test <- window(X, start = c(2003,(split_idx+1)))
    return(list(X_train=X_train,y_train=y_train,
              X_test=X_test,y_test=y_test))
  }
}
#abc<-Train_Test_Split(block,y_ts, c(0.75,0.25))


