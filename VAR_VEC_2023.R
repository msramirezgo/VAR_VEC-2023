# Paquetes 
if (1) {
  library(vars)
  library(urca)
  library(ggplot2)
  library(ggfortify)
  library(gridExtra)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(tseries)
  library(tsDyn)
  library(dynlm)
  library(aTSA)
}
# Parámetros del  modelo
frequency     = c("Quarterly", "Yearly")[2]
if(frequency=="Quarterly"){
  file.name  = "Belize_Quarterly.xlsx"
  # Nombres de columnas que contiene las variables para el VAR.
  variables  = c("Revenue Current",	"Revenue Current - Base 2000",	"Total Revenue and Grants - Base 2000",
                 "Total Revenue and Grants",	"Expenditure Current",	"Expenditure Current - Base 2000",	
                 "Total Expenditure",	"Total Expenditure - Base 2000",	"GDP - Base 2000",	"GDP - Base 2014")[c(2,9)]
  start.date = c(2000,1)
  frequency  = 4
}

#For yearly data
if (frequency=="Yearly"){
  file.name    = "Yearly_Current.xlsx"
  # Nombres de columnas que contiene las variables para el VAR.
  variables    = c( "Total expenditures"
                    ,"Total recurrent expenditure"
                    ,"Personal emoluments"
                    ,"Pensions and ex-gratia"
                    ,"Goods and services"
                    ,"Debt service-interest and other charges"
                    ,"GDP")[c(6,7)]
  start.date = c(2001)
  frequency  = 1
} 

log.all      = c(TRUE,FALSE)[1]                # Verdadero transforma todas las series en logaritmos.
diff.all     = c(TRUE,FALSE)[2]                # Verdadero obtiene la primera diferencia de las series.
max.lags     = 10                              # Número máximo posible para el orden p del VAR. 
lags.pt.test = c(10,15,20,30,50,75)            # Rezagos a usar para las pruebas de autocorrelación serial. (Portmanteau statistic)
n.ahead      = 4                               # Número de pasos adelante para el pronóstico.
eigen.confidence = c("10pct","5pct","1pct")[2] # Nivel de significancia para la prueba de johanssen
EG.procedure = c(TRUE,FALSE)[2]                # Ejecutar la metodología de Engle & Granger.

# Datos -------------------------------------------------------------------
Data = read_xlsx(paste0(getwd(),"/",file.name))
Data = ts(Data[,variables], start = start.date,   frequency = frequency)
#if (variables[1]=="Total expenditures"|variables[1]=="Total recurrent expenditure") {
# Data=Data[3:nrow(Data),] 
#}

# Graficación -------------------------------------------------------------
x11()
par(mfrow=c(length(variables),1))
for (i in variables) {
  plot.ts(Data[,i], col="steelblue", xlab="",ylab="", main=i)
} # Hacer análisis manual de resultados del test. 


# Transformaciones --------------------------------------------------------
if (log.all==TRUE)  {
  for (i in variables) {
    Data[,i]=log(Data[,i])
  }
}
if (diff.all==TRUE) {
  for (i in variables) {
    Data[-1,i]=diff(Data[,i])
  }
  Data=Data[-1,] #Eliminamos la primera fila correspondiente al dato perdido en la diferencia. 
}
x11()
par(mfrow=c(ncol(Data),1))
for (i in variables) {
  plot.ts(Data[,i], col="steelblue", xlab="",ylab="", main=i)
}
# Pruebas de raíz unitaria ------------------------------------------------
type=matrix(NA, nrow=length(variables), ncol=1, dimnames=list(variables, "type"))
Int.Order=c()
for (i in variables) {
  
  trend = ur.df(Data[,i], lags=6, selectlags = "AIC",type="trend")
  drift = ur.df(Data[,i], lags=6, selectlags = "AIC",type="drift")
  none  = ur.df(Data[,i], lags=6, selectlags = "AIC",type="none")
  
  #---- Cuando se rechaze H0 de raíz unitaria ----#
  if (trend@teststat[,"phi3"]>trend@cval["phi3","5pct"] & trend@teststat[,"tau3"]<trend@cval["tau3","5pct"]) {
    type[i,]="trend"
    cat(i,"\n Type: Trend \n tau: Null hyphotesis for unitary root rejected \n")
    Int.Order=0
  } else if (drift@teststat[,"phi1"]>drift@cval["phi1","5pct"] & drift@teststat[,"tau2"]<drift@cval["tau2","5pct"]){
    type[i,]="drift"
    cat(i,"\n Type: Drift \n tau: Null hyphotesis for unitary root rejected \n")
    Int.Order=0
  } else if (none@teststat[,"tau1"]<none@cval["tau1","5pct"]){
    type[i,]="none"
    cat(i,"\n Type: None \n tau: Null hyphotesis for unitary root rejected \n")
    Int.Order[i]=0
  }
  #---- Cuando no se rechaze H0 en cualquier caso ----# 
  if (trend@teststat[,"tau3"]>trend@cval["tau3","5pct"] & 
      drift@teststat[,"tau2"]>drift@cval["tau2","5pct"] &
      none@teststat[,"tau1"] >none@cval["tau1", "5pct"]) {
    Int.Order[i]=1
    cat(i,"\n", "tau: Null hyphotesis for unitary root NOT REJECTED","\n")
  }
}


# Orden p del VAR ---------------------------------------------------------
if (sum(Int.Order)==0) {
  cat("Series are I(0):Estimating a VAR in levels.")
}else if (sum(Int.Order)==length(variables)) cat("Series are I(1): Estimating VAR and  evaluating for cointegration (Johansen) \n")

if (frequency=="Yearly")    lagmax=2
if (frequency=="Quarterly") lagmax=6

# Extraemos los rezagos óptimos para el modelo con tendencia y deriva.
P.tr.cons = VARselect(Data, lag.max=lagmax, type="both", season=NULL)
AIC.tr    = which(t(P.tr.cons$criteria)==min(P.tr.cons$criteria["AIC(n)",]))
# Extraemos los rezagos óptimos para el modelo con deriva.
P.cons    = VARselect(Data, lag.max=lagmax, type="const", season=NULL)
AIC.cons  = which(t(P.cons$criteria)   ==min(P.cons$criteria["AIC(n)",]))
# Extraemos los rezagos óptimos para el modelo sin términos deterministicos.
P.none    = VARselect(Data, lag.max=lagmax, type="none", season=NULL) 
AIC.none  = which(t(P.none$criteria)   ==min(P.none$criteria["AIC(n)",]))


# Estimación del VAR ------------------------------------------------------

# Estimamos los modelos con tendencia y constante, sólo constante, y sólo términos deterministicos.
VAR.both  = VAR(Data, p=AIC.tr  , type="both" )
VAR.const = VAR(Data, p=AIC.cons, type="const")
VAR.none  = VAR(Data, p=AIC.none, type="none" )

# Creamos dataframe con el respectivo AIC de cada modelo para elegir el que mejor se ajuste. 
AIC.VAR   = matrix(c(AIC(VAR.both),AIC(VAR.const), AIC(VAR.none)), nrow=1, ncol=3, dimnames=list("AIC", c("both", "const", "none")))
VAR.type  = colnames(AIC.VAR)[which(AIC.VAR==min(AIC.VAR))]   # Se recomiendo ver la significancia de los parámetros "const" y "trend" de forma manual para
# dar robustes a los resultados o para corregir de ser necesario. 

if (VAR.type=="both")  {
  VAR=VAR.both
  cat("Using VAR trend and constant. \n")}
if (VAR.type=="const") {
  VAR=VAR.const
  cat("Using VAR with constant. \n")}
if (VAR.type=="none")  {
  VAR=VAR.none
  cat("Usig VAR with no deterministic terms. \n")}

# Prueba de Cointegración (Rango de Pi)---------------------------------------

eigen=ca.jo(Data, ecdet = if (VAR.type=="both") {
              "trend"
            }else{
              VAR.type
            }, 
            type = "eigen", K = 2, 
            spec = "longrun", season = NULL)


summary(eigen) 
critical.values=rev(eigen@cval[,eigen.confidence])
for(i in 1:length(eigen@teststat)){
  if (rev(eigen@teststat)[i]>critical.values[i]) {
    cat("Matrix has rank",i,"\n")
    rank=i
    rank.type="reduced"
  }
  if(rev(eigen@teststat)[1]<critical.values[1]&rev(eigen@teststat)[2]<critical.values[2]){
    rank.type = "zero"
    rank      = 0
  } else if(eigen@teststat[1]>critical.values[1]&rev(eigen@teststat)[2]>critical.values[2]){
    rank.type = "complete"
    rank=length(variables)
  }
}
if(rank==0&sum(Int.Order)==length(variables)){
  cat("Rank 0 & series are I(1): Estimate a VAR in differences \n")
  rank.type="zero"
}
if(rank>0 & rank<length(variables) & sum(Int.Order)==length(variables)){
  cat("Reduced Rank: Estimate a VEC(p-1)")
  rank.type="reduced"
}
if(rank==length(variables)&sum(Int.Order)==0){
  cat("Complete Rank & series are I(0): Estimate a VAR in levels \n")
  rank.type="complete"
}
if(rank==length(variables)&sum(Int.Order)==length(variables)){
  cat("Complete Rank & series are I(1): Estimate a VAR in first differences \n")
  rank.type="complete"
}

# VEC ---------------------------------------------------------------------
if (rank.type=="reduced") {
  
  VEC = cajorls(eigen, r=rank) 
  VEC
  # Vector de cointegración normalizado
  coefB(VEC)
  # Coeficientes de velocidad de ajuste
  coefA(VEC)
  # Representación VAR
  VAR.rep = vec2var(eigen, r=rank)

  # Pronóstico de la representación VAR
  cat("Forecasting VECM in it's VAR representation \n")
  predict =predict(VAR.rep, n.ahead=n.ahead)
  forecast    = cbind(as.matrix(predict$fcst[[1]][,"fcst"]),as.matrix(predict$fcst[[2]][,"fcst"]))
  if(log.all==TRUE)levels = exp(as.matrix(rbind(Data, forecast)))
  write.xlsx(levels, file = paste0("VECM_",n.ahead,"_step_ahead_forecast_",variables[1],"_With_",variables[2],".xlsx"))
  x11()
  plot(predict)
  
  
}else if (rank.type=="complete" & sum(Int.Order)==0){
  # Pronóstico con VAR en niveles
  cat("Forecasting VAR in levels")
  predict=predict(VAR, n.ahead=n.ahead)
  forecast    = cbind(as.matrix(predict$fcst[[1]][,"fcst"]),as.matrix(predict$fcst[[2]][,"fcst"]))
  levels      = as.matrix(rbind(Data, forecast))
  for (i in nrow(Data):(nrow(Data)+n.ahead-1)) {
    levels[i+1,]=levels[i,]+levels[i+1,]
  }
  levels = ts(levels, start = start.date,   frequency = frequency)
  if(log.all==TRUE)levels = exp(levels)
  write.xlsx(levels, file = paste0("VAR_",n.ahead,"_step_ahead_forecast_",variables[1],"_With_",variables[2],".xlsx"))
  
  x11()
  plot(predict)
  
  x11()
  for (i in lags.pt.test) {
    plot(serial.test(VAR, lags.pt = i, type = "PT.asymptotic"), title=paste(i," lags"))
  } 
  
  # Navegue por el dispositivo gráfico para ver el 
  # resumen de resultados del Test para cada variable
  # y para cada orden de rezagos. 
  
  
  ##Test Jarque-Bera multivariado
  norm.test=normality.test(VAR.IR) #
  if(norm.test$jb.mul$JB$p.value<0.05)cat("Rejection for normality.")
  if(norm.test$jb.mul$JB$p.value>0.05)cat("No rejection for normality.")
  
}else if (rank.type=="zero" & sum(Int.Order)==length(variables)){
  # Estimación del VAR ------------------------------------------------------
  
  # Estimamos los modelos con tendencia y constante, sólo constante, y sólo términos deterministicos.
  VAR.both  = VAR(diff(Data), p=AIC.tr  , type="both" )
  VAR.const = VAR(diff(Data), p=AIC.cons, type="const")
  VAR.none  = VAR(diff(Data), p=AIC.none, type="none" )
  
  # Creamos dataframe con el respectivo AIC de cada modelo para elegir el que mejor se ajuste. 
  AIC.VAR   = matrix(c(AIC(VAR.both),AIC(VAR.const), AIC(VAR.none)), nrow=1, ncol=3, dimnames=list("AIC", c("both", "const", "none")))
  VAR.type  = colnames(AIC.VAR)[which(AIC.VAR==min(AIC.VAR))]   # Se recomiendo ver la significancia de los parámetros "const" y "trend" de forma manual para
  # dar robustes a los resultados o para corregir de ser necesario. 
  
  if (VAR.type=="both")  {
    VAR=VAR.both
    cat("Using VAR trend and constant. \n")}
  if (VAR.type=="const") {
    VAR=VAR.const
    cat("Using VAR with constant. \n")}
  if (VAR.type=="none")  {
    VAR=VAR.none
    cat("Usig VAR with no deterministic terms. \n")}
  
  cat("Forecasting VAR in first differences \n")
  predict     = predict(VAR, n.ahead=n.ahead)
  forecast    = cbind(as.matrix(predict$fcst[[1]][,"fcst"]),as.matrix(predict$fcst[[2]][,"fcst"]))
  levels      = as.matrix(rbind(Data, forecast))
  for (i in nrow(Data):(nrow(Data)+n.ahead-1)) {
    levels[i+1,]=levels[i,]+levels[i+1,]
  }
  levels = ts(levels, start = start.date,   frequency = frequency)
  if(log.all==TRUE)levels = exp(levels)
  write.xlsx(levels, file = paste0("VAR_",n.ahead,"_step_ahead_forecast_",variables[1],"_With_",variables[2],".xlsx"))
  x11()
  plot(predict)
}


# Engle & Granger ---------------------------------------------------------

if (EG.procedure==TRUE) {
  spread    = Data[,1]-Data[,2]
  Data.plot = cbind(Data, spread)
  x11()
  plot(as.zoo(Data.plot[,c(1,2)]),
       plot.type = "single",
       lty = c(2, 1),
       lwd = 2)
  x11()
  plot.ts(Data.plot[,"spread"],
          col = "steelblue",
          lwd = 2,
          ylab="spread")
  
  type=matrix(NA, nrow=length(variables), ncol=1, dimnames=list(variables, "type"))
  Int.Order=c()
  
  for (i in variables) {
    
    trend = ur.df(Data[,i], lags=6, selectlags = "AIC",type="trend")
    drift = ur.df(Data[,i], lags=6, selectlags = "AIC",type="drift")
    none  = ur.df(Data[,i], lags=6, selectlags = "AIC",type="none")
    
    #---- Cuando se rechaze H0 de raíz unitaria ----#
    if (trend@teststat[,"phi3"]<trend@cval["phi3","5pct"] & trend@teststat[,"tau3"]<trend@cval["tau3","5pct"]) {
      type[i,]="trend"
      cat(i,"\n Type: Trend \n tau: Null hyphotesis for unitary root rejected \n")
      Int.Order=0
    } else if (drift@teststat[,"phi1"]<drift@cval["phi1","5pct"] & drift@teststat[,"tau2"]<drift@cval["tau2","5pct"]){
      type[i,]="drift"
      cat(i,"\n Type: Drift \n tau: Null hyphotesis for unitary root rejected \n")
      Int.Order=0
    } else if (none@teststat[,"tau1"]<none@cval["tau1","5pct"]){
      type[i,]="none"
      cat(i,"\n Type: None \n tau: Null hyphotesis for unitary root rejected \n")
      Int.Order[i]=0
    }
    #---- Cuando no se rechaze H0 en cualquier caso ----# 
    if (trend@teststat[,"tau3"]> trend@cval["tau3","5pct"] & 
        drift@teststat[,"tau2"]> drift@cval["tau2","5pct"] &
        none@teststat[,"tau1"] > none@cval["tau1", "5pct"]) {
      Int.Order[i]=1
      cat(i,"\n", "tau: Null hyphotesis for unitary root NOT REJECTED","\n")
    }
  }
  reg <- dynlm(Data[,1] ~ Data[,2]) 
  summary(reg)
  
  # Obtenemos los residuales
  z_hat <- residuals(reg)
  # Los graficamos
  x11()
  plot.ts(z_hat)
  
  lags=nrow(Data/4)
  par(mfrow=c(1,2))
  acf(z_hat,lag.max=lags,plot=T,lwd=2,xlab='',main='ACF de los residuales') 
  pacf(z_hat,lag.max=lags,plot=T,lwd=2,xlab='',main='PACF de los residuales')
  
  #Hacemos la prueba sobre los residuales
  coint.test(Data[,1], Data[,2], nlag=6)
  
}




#Pruebas
if(0){
  # Pruebas de correlación serial -------------------------------------------
  #Pruebas
  for (i in lags.pt.test) {
    PT.test=serial.test(VAR, lags.pt = i, type = "PT.asymptotic")
    cat(i,"lags", "\n","p-value =" ,PT.test$serial$p.value, if(PT.test$serial$p.value>0.05){
      "No se rechaza el supuesto"
    }else "Se rechaza el supuesto","\n")
  } 
  #Graficación
  x11()
  for (i in lags.pt.test) {
    plot(serial.test(VAR, lags.pt = i, type = "PT.asymptotic"), title=paste(i," lags"))
  } # Navegue por el dispositivo gráfico para ver el 
  # resumen de resultados del Test para cada variable
  # y para cada orden de rezagos. 
  
  
  
}

