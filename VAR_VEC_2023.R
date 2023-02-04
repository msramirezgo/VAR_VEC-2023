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
}

# Parámetros del modelo
#file.name    = "C:/Users/maico/OneDrive - Universidad Nacional de Colombia/Oliver Pardo/Modelos Multiecuaciones (VAR)/2_datos/Enders.xlsx"
file.name    = "BELIZE_REV_EXP_DATA.xlsx"
variables    = c("Total Revenue and Grants","Total Expenditure")[1:2] # Nombres de columnas que contiene las variables para el VAR.
log.all      = c(TRUE,FALSE)[1]           # Verdadero transforma todas las series en logaritmos.
diff.all     = c(TRUE,FALSE)[2]           # Verdadero obtiene la primera diferencia de las series.
max.lags     = 10                         # Número máximo posible para el orden p del VAR. 
lags.pt.test = c(10,15,20,30,50,75)          # Rezagos a usar para las pruebas de autocorrelación serial. (Portmanteau statistic)


# Datos -------------------------------------------------------------------
Data = read_xlsx(file.name)
Data = ts(Data[,variables], start = c(1989,2), frequency = 4)

# Graficación -------------------------------------------------------------
x11()
par(mfrow=c(length(variables),1))
for (i in variables) {
  plot.ts(Data[,i], size=2,ts.colour="red", xlab="",ylab="", main=i)
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
plot.ts(Data[,i],size=2, color="red", xlab="",ylab="", main=i)
}
# Pruebas de raíz unitaria ------------------------------------------------
type=matrix(NA, nrow=length(variables), ncol=1, dimnames=list(variables, "type"))
for (i in variables) {
  
  trend = ur.df(Data[,i], lags=6, selectlags = "AIC",type="trend")
  drift = ur.df(Data[,i], lags=6, selectlags = "AIC",type="drift")
  none  = ur.df(Data[,i], lags=6, selectlags = "AIC",type="none")
  
  #---- Cuando se rechaze H0 de raíz unitaria ----#
  if (trend@teststat[,"phi3"]<trend@cval["phi3","5pct"] & trend@teststat[,"tau3"]<trend@cval["tau3","5pct"]) {
    type[i,]="trend"
    cat(i,"\n","Type: Trend","\n","tau: Null hyphotesis for unitary root rejected", "\n")
  } else if (drift@teststat[,"phi1"]<drift@cval["phi1","5pct"] & drift@teststat[,"tau2"]<drift@cval["tau2","5pct"]){
    type[i,]="drift"
    cat(i,"\n","Type: Drift","\n","tau: Null hyphotesis for unitary root rejected", "\n")
  } else if (none@teststat[,"tau1"]<none@cval["tau1","5pct"]){
    type[i,]="none"
    cat(i,"\n","Type: None","\n","tau: Null hyphotesis for unitary root rejected", "\n")
  }
  
  #---- Cuando no se rechaze H0 en cualquier caso ----# 
  if (trend@teststat[,"tau3"]>trend@cval["tau3","5pct"] & 
      drift@teststat[,"tau2"]>drift@cval["tau2","5pct"] &
      none@teststat[,"tau1"] >none@cval["tau1", "5pct"]) {
    cat(i,"\n", "tau: Null hyphotesis for unitary root NOT REJECTED","\n")
  }
}


# Orden p del VAR ---------------------------------------------------------

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


# Prueba de COintegración -------------------------------------------------

eigen=ca.jo(Data, 
            ecdet = if (VAR.type=="both") {
              "trend"
            }else{
              VAR.type
            }, 
            type = "eigen", K = if (VAR.type=="both") {
              AIC.tr
            }else if(VAR.type=="const"){AIC.cons} else{
              AIC.none
            }, 
            spec = "longrun", season = NULL)


summary(eigen) 

for(i in 1:ncol(Data)){
  if (eigen@teststat[i]<eigen@cval[i,"5pct"]) {
    cat(i, "Cointegration relationship")
    coint=i
  }
}

# VEC ---------------------------------------------------------------------
if (coint!=0) {
  VEC = cajorls(eigen, r=coint) 
  VEC
  # Vector de cointegración normalizado
  coefB(VEC)
  # Coeficientes de velocidad de ajuste
  coefA(VEC)
  VAR.IR = vec2var(eigen, r = coint)
  
  ################################################
  # Validación de supuestos: representación VAR
  ################################################
  
  ## Autocorrelación: PT.asymptotic es para muestra grande y "PT.adjusted" es corrección para muestra pequeña.
  
  for (i in lags.pt.test) {
    PT.test=serial.test(VAR.IR, lags.pt = i, type = "PT.asymptotic")
    cat(i,"lags", "\n","p-value =" ,PT.test$serial$p.value, if(PT.test$serial$p.value>0.05){
      "No se rechaza el supuesto"
    }else "Se rechaza el supuesto","\n")
  } 
  
  #Graficación
  x11()
  for (i in lags.pt.test) {
    plot(serial.test(VAR.IR, lags.pt = i, type = "PT.asymptotic"), title=paste(i," lags"))
  } # Navegue por el dispositivo gráfico para ver el 
  # resumen de resultados del Test para cada variable
  # y para cada orden de rezagos. 
  
  
  #Graficamos los residuales para 20 lags: se grafican los residuales, su distribución, la ACF y PACF de los residuales y
  #la ACF y PACF de los residuales al cuadrado (proxy para heterocedasticidad)
  x11()
  for (i in colnames(Data)) {
    x11()
    plot(P.20, names = i, title=i) #Bien comportados, salvo por los residuales al cuadrado
  } 
  #Homocedasticidad: Test tipo ARCH multivariado
  arch.test(VAR.IR, lags.multi = 24, multivariate.only = TRUE) #Rechazo, no se cumple el supuesto.
  arch.test(VAR.IR, lags.multi = 12, multivariate.only = TRUE) #Rechazo, no se cumple el supuesto
  
  ##Test Jarque-Bera multivariado
  normality.test(VAR.IR) #Rechazo, no se cumple el supuesto. 
  
  
  ###########################################
  # Pronóstico de la representación VAR
  ###########################################
  
  # Debido al incumplimiento de normalidad, los intervalos de confianza deben computarse por bootstrapping.
  x11()
  predict(VAR.IR, n.ahead=12)
  plot(predict(VAR.IR, n.ahead=12))
  
}



