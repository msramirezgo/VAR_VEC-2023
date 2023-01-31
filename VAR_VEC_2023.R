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
file.name = "C:/Users/maico/OneDrive - Universidad Nacional de Colombia/Oliver Pardo/Modelos Multiecuaciones (VAR)/2_datos/Enders.xlsx"
variables = c("IPI","CPI","Unem")[1:3] # Nombres de columnas que contiene las variables para el VAR
log.all   = c(TRUE,FALSE)[1]
diff.all  = c(TRUE,FALSE)[1]

# Datos -------------------------------------------------------------------
Data = read_xlsx(file.name)
Data = ts(Data[,variables])

# Graficación -------------------------------------------------------------
x11()
par(mfrow=c(length(variables),1))
for (i in variables) {
  plot.ts(Data[,i], size=2,ts.colour="red", xlab="",ylab="", main=i)
}


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
par(mfrow=c(3,1))
for(i in variables){
plot.ts(Data[,i])
}
# Pruebas de raíz unitaria ------------------------------------------------
type=matrix(NA, nrow=length(variables), ncol=1, dimnames=list(variables, "type"))
for (i in variables) {
  
  trend = ur.df(Data[,i], lags=6, selectlags = "AIC",type="trend")
  drift = ur.df(Data[,i], lags=6, selectlags = "AIC",type="drift")
  none  = ur.df(Data[,i], lags=6, selectlags = "AIC",type="none")
  
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
  if (trend@teststat[,"tau3"]>trend@cval["tau3","5pct"] & 
      drift@teststat[,"tau2"]>drift@cval["tau2","5pct"] &
      none@teststat[,"tau1"] >none@cval["tau1", "5pct"]) cat(i,"\n", "tau: Null hyphotesis for unitary root NOT REJECTED")
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

VAR(Y, p=AIC.tr  , type="both" )
VAR(Y, p=AIC.cons, type="const")
VAR(Y, p=AIC.none, type="const")




# VEC ---------------------------------------------------------------------




