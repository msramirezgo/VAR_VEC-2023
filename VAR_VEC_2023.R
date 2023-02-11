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
}

# Parámetros del  modelo
#file.name    = "C:/Users/maico/OneDrive - Universidad Nacional de Colombia/Oliver Pardo/Modelos Multiecuaciones (VAR)/2_datos/Enders.xlsx"
file.name    = "Belize_Quarterly.xlsx"
#For yearly data
#variables    = c("Total Revenue and Grants","Total Expenditure")[1:2] # Nombres de columnas que contiene las variables para el VAR.
#For quarterly data
variables    = c("Revenue Current",	"Revenue Current - Base 2000",	"Total Revenue and Grants - Base 2000",
                 "Total Revenue and Grants",	"Expenditure Current",	"Expenditure Current - Base 2000",	
                 "Total Expenditure",	"Total Expenditure - Base 2000",	"GDP - Base 2000",	"GDP - Base 2014")[c(2,9)]
log.all      = c(TRUE,FALSE)[1]           # Verdadero transforma todas las series en logaritmos.
diff.all     = c(TRUE,FALSE)[2]           # Verdadero obtiene la primera diferencia de las series.
max.lags     = 10                         # Número máximo posible para el orden p del VAR. 
lags.pt.test = c(10,15,20,30,50,75)       # Rezagos a usar para las pruebas de autocorrelación serial. (Portmanteau statistic)
n.ahead      = 15                         # Número de pasos adelante para el pronóstico

# Datos -------------------------------------------------------------------
Data = read_xlsx(file.name)
Data = ts(Data[,variables], start = c(2000,1), frequency = 4)
Data = ts(Data[1:nrow(Data)-1,], start = c(2000,1), frequency = 4)

# Graficación -------------------------------------------------------------
x11()
par(mfrow=c(length(variables),1))
for (i in variables) {
  plot.ts(Data[,i], col="steelblue", xlab="",ylab="", main=i)
  
} # Hacer análisis manual de resultados del test. 


# Transformaciones --------------------------------------------------------
if (log.all==TRUE)    {
  for (i in variables) {
    Data[,i]=log(Data[,i])
  }
}
if (diff.all==TRUE)   {
  for (i in variables) {
    Data[-1,i]=diff(Data[,i])
  }
  Data=Data[-1,] #Eliminamos la primera fila correspondiente al dato perdido en la diferencia. 
}
x11()
par(mfrow=c(ncol(Data),1))
for(i in variables){
plot(Data[,i], col="red", xlab="",ylab="", main=i)
}
# Pruebas de raíz unitaria ------------------------------------------------
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
  if (trend@teststat[,"tau3"]>trend@cval["tau3","5pct"] & 
      drift@teststat[,"tau2"]>drift@cval["tau2","5pct"] &
      none@teststat[,"tau1"] >none@cval["tau1", "5pct"]) {
    Int.Order[i]=1
    cat(i,"\n", "tau: Null hyphotesis for unitary root NOT REJECTED","\n")
  }
}


# Orden p del VAR ---------------------------------------------------------
if (sum(Int.Order)==0) {
  cat("Series are I(0):Estimating a VAR(p-1) in first differences.")
  Data=diff(Data)
}else if (sum(Int.Order)==length(variables)) cat("Series are I(1): Estimating VAR and  evaluating for cointegration (Johansen)")

# Extraemos los rezagos óptimos para el modelo con tendencia y deriva.
P.tr.cons = VARselect(Data, lag.max=6, type="both", season=NULL)
AIC.tr    = which(t(P.tr.cons$criteria)==min(P.tr.cons$criteria["AIC(n)",]))
# Extraemos los rezagos óptimos para el modelo con deriva.
P.cons    = VARselect(Data, lag.max=6, type="const", season=NULL)
AIC.cons  = which(t(P.cons$criteria)   ==min(P.cons$criteria["AIC(n)",]))
# Extraemos los rezagos óptimos para el modelo sin términos deterministicos.
P.none    = VARselect(Data, lag.max=6, type="none", season=NULL) 
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
  cat("Se elige el VAR con tendencia y constante")}
if (VAR.type=="const") {
  VAR=VAR.const
  cat("Se elige el VAR con constante")}
if (VAR.type=="none")  {
  VAR=VAR.none
  cat("Se elige el VAR sin términos deterministicos")}

# Prueba de Cointegración (Rango de Pi)---------------------------------------

eigen=ca.jo(Data, 
            ecdet = if (VAR.type=="both") {
              "trend"
            }else{
              VAR.type
            }, 
            type = "trace", K = if (VAR.type=="both") {
              AIC.tr
            }else if(VAR.type=="const"){AIC.cons} else{
              AIC.none
            }, 
            spec = "longrun", season = NULL)


summary(eigen) 

for(i in 1:ncol(Data)){
  if (eigen@teststat[i]>eigen@cval[i,"5pct"]) {
    cat("Matrix has rank",i,"\n")
    rank=i
    rank.type="reduced"
  } else {
    rank=length(variables)
    rank.type="complete"
  }
}
if(eigen@teststat[1]>eigen@cval[i,"5pct"]&sum(Int.Order)==0){
  cat("Rank 0: Estimate a VAR in differences")
  rank.type="zero"
}
if(rank==length(variables)&sum(Int.Order)==length(variables)){
  cat("Complete Rank: Estimate a VAR in levels")
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
  
  ################################################
  # Validación de supuestos: representación VAR
  ################################################
  
  ## Autocorrelación: PT.asymptotic es para muestra grande y "PT.adjusted" es corrección para muestra pequeña.
  
  for (i in lags.pt.test) {
    PT.test=serial.test(VAR.rep, lags.pt = i, type = if(nrow(Data>=100))"PT.asymptotic" else "PT.adjusted")
    cat(i,"lags", "\n","p-value =" ,PT.test$serial$p.value, if(PT.test$serial$p.value>0.05){
      "No se rechaza el supuesto"
    }else "Se rechaza el supuesto","\n")
  } 
  
  #Graficación
  
  x11()
  for (i in lags.pt.test) {
    plot(serial.test(VAR.IR, lags.pt = i, type = "PT.asymptotic"), title=paste(i," lags"))
  } 
  # Navegue por el dispositivo gráfico para ver el 
  # resumen de resultados del Test para cada variable
  # y para cada orden de rezagos. 
  
  
  ##Test Jarque-Bera multivariado
  norm.test=normality.test(VAR.IR) #
  if(norm.test$jb.mul$JB$p.value<0.05)cat("Rejection for normality.")
  if(norm.test$jb.mul$JB$p.value>0.05)cat("No rejection for normality.")
  
  
  ###########################################
  # Pronóstico de la representación VAR
  ###########################################
  
  # Debido al incumplimiento de normalidad, los intervalos de confianza deben computarse por bootstrapping.
  x11()
  predict=predict(VAR.IR, n.ahead=n.ahead)
  plot(predict)
  
}else
{
  predict=predict(VAR, n.ahead=n.ahead)
  for (i in names(predict$fcst)) {
    write.csv(predict$fcst[i], file=paste0(n.ahead, " step ahead forecast of ",i))
  }
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