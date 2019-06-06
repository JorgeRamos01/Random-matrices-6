rm(list=ls())
#leer symbolos
library(quantmod)
library(rvest)
pagina <- read_html("https://en.wikipedia.org/wiki/List_of_S%26P_500_companies")

noticias_texto_principal <- pagina %>%
  html_nodes("#mw-content-text table") %>%
  html_table(fill=T)

etiquetassp<-noticias_texto_principal[[1]][1]

#Obtenemos las series con valores diarios, debido a que el semanal solo toma los de los martes
lista_datos_wiki<-list()
for(i in 1:nrow(etiquetassp)){
  if( !(as.character(etiquetassp[i,1]) %in% c("CBG","BF.B","BRK.B","CSRA","DPS","GGP","LUK","MON","PCLN","SNI","TWX","HCN","WYN")) ){
    lista_datos_wiki[[as.character(etiquetassp[i,1])]] = getSymbols(as.character(etiquetassp[i,1]),auto.assign = FALSE, from = "2008-01-01",periodicity="daily")[,4] 
  }
}

mydate <- as.character(seq.Date(from = as.Date("2008-01-04"), 
                                to = as.Date("2018-11-16"),
                                by = "week"))


#Revisamos el numero de semanas para considerar las de mayor longitud
n<-vector("integer", length = 0)
for (i in names(lista_datos_wiki)){
  print((nrow(lista_datos_wiki[[i]])))
  n[i]<-(nrow(lista_datos_wiki[[i]]))
}
n<-max(n)

#Guardamos los datos de cierre en un data frame
data=matrix(0L, nrow=n)
for (i in names(lista_datos_wiki)){
  if (nrow(lista_datos_wiki[[i]])>=n){
    data=cbind(data, lista_datos_wiki[[i]][,1])
  }
}

data=data[mydate,]  #Generamos el dataframe con solo los viernes
#Removemos lista con NAs
na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
data=data[, na_count==0]

#Conseguimos los datos del indicador global
indicador_global= getSymbols("^GSPC",auto.assign = FALSE, from = "2008-01-07",periodicity="weekly")[,1]
#Matcheamos la lista del indicador global con las fechas que tenemos datos pero considerando para los lunes
indicador_global<-indicador_global$GSPC.Open[time(data)+3]
#Este se usa en caso de tener el ultimo dato repetido
#indicador_global<-indicador_global$GSPC.Open[1:nrow(indicador_global)-1]  

#Realizamos pruebas de estacionariedad
library(tseries)
library(forecast)
library(urca)
data=data[,-1]
estacionarias<-vector(mode="numeric", length = 0) #No todas las series son estacionarias
for (i in 2:ncol(data)){
  estacionarias<-c(estacionarias,adf.test(data[,i])$p.value)
}
sum(estacionarias<0.05) #La suma es menor que el numero de variables por lo que hay variables con raices unitarias

data2<-diff(data)
data2<-data2[-1,]
estacionarias<-vector(mode="numeric", length = 0) #Se genera estacionariedad con la primera diferencia con confianza del 95%
for (i in 1:ncol(data2)){
  estacionarias<-c(estacionarias,adf.test(data2[,i])$p.value)
}
sum(estacionarias<0.05, na.rm=TRUE) #Todas las series son estacionarias

summary(ur.df(as.ts(indicador_global), type="none"))  #Aplicamos prueba de estacionariedad para el indicador global
                                                    #y encontramos la presencia de una raiz unitaria
summary(ur.df(as.ts(diff(indicador_global)[-1,]), type="none")) #Encontramos un p-value peque?o y parece estacionaria
                                                              #en la primera diferencia
#Comprobamos estacionariedad graficando
par(mfrow=c(2,1))
ts.plot(as.ts(indicador_global))
ts.plot(as.ts(diff(indicador_global)[-1,])) #Estacionaria en la primera diferencia

par(mfrow=c(1,1))
#Calculamos PCA a la matriz de datos diferenciados
eigenv<-prcomp(data2[,-1])

#Aplicamos un test de Marcenko-Pastur
Marc.Pastur<-(1+1/sqrt(ncol(data2)/nrow(data2)))^2  #Aparecen 4.41 y redondeamos a 5 componentes
which(cumsum((eigenv$sdev)^2/sum(eigenv$sdev^2))>Marc.Pastur/10)[1]
#Consideramos el criterio de varianza
crit.var<-which(cumsum((eigenv$sdev)^2/sum(eigenv$sdev^2))>0.80)[1] #Con 19 componentes principales se genera el 80% de la varianza

plot(cumsum((eigenv$sdev)^2/sum(eigenv$sdev^2)), ylab="Varianza", xlab="Indice", main="Varianza explicada")
points(crit.var,cumsum((eigenv$sdev)^2/sum(eigenv$sdev^2))[crit.var], col="red", pch=19)

############Aplicando PLS
library(pls)
data3<-matrix(0L, nrow=nrow(data2))
data3=cbind(data3, as.vector(diff(indicador_global$GSPC.Open)[-1,]))

for (i in 1:ncol(data2)){
  data3<-cbind(data3, as.vector(data2[,i]))
}
data3<-data3[,-1]
data3<-as.data.frame(data3)

#Aplicando PLS usando las componentes por M-P hasta el 2017 para predecir el 2018
pcr.marc.past<-pcr(V1 ~ ., ncomp =5, data = data3[1:501,], validation = "CV")

#Aplicando PLS usando las componentes por el criterio de varianza hasta el 2017 para predecir el 2018
pcr.variance<-plsr(V1 ~ ., ncomp =19, data = data3[1:501,], validation = "CV")

#Verificamos los modelos
summary(pcr.marc.past)
summary(pcr.variance)

plot(RMSEP(pcr.marc.past), legendpos = "topright", main="Criterio M-P")
plot(RMSEP(pcr.variance), legendpos = "topright", main="Criterio de varianza")

#Verificamos la prediccion del primer a?o

ts.plot(predict(pcr.marc.past, ncomp=5, newdata=data3[502:nrow(data3),-1]), main="M-P criteria", xlab="Tiempo", ylab="Diff pred")
lines(as.ts(data3[502:nrow(data3),1]), col="red")
legend(20,-50, legend=c("predicho", "real"), col=c("black", "red"), lty=1)

ts.plot(predict(pcr.variance, ncomp=19, newdata=data3[502:nrow(data3),-1]), main="Covariance criteria", xlab="Tiempo", ylab="Diff pred")
lines(as.ts(data3[502:nrow(data3),1]), col="red")
legend(20,-100, legend=c("predicho", "real"), col=c("black", "red"), lty=1)


#Verificamos la efectividad con las diferencias
mp.effect<-sign(predict(pcr.marc.past, ncomp=5, newdata=data3[503:nrow(data3),-1])-data3[502:(nrow(data3)-1),1])

real.effect<-sign(diff(data3[502:nrow(data3),1]))
sum(mp.effect==real.effect)/length(real.effect)

#Verificamos con los datos en niveles

pred<-predict(pcr.marc.past, ncomp=5, newdata=data3[502:nrow(data3),-1])
pred[1]<-pred[1]+indicador_global[501,]
pred.niv<-cumsum(pred)

#Graficamos
ts.plot(pred.niv, main="Criterio de M-P", sub="Datos en niveles")
lines(as.vector(indicador_global$GSPC.Open[502:nrow(data3)]), col="red")
legend(20,2650, legend=c("predicho", "real"), col=c("black", "red"), lty=1)

#Graficamos la efectividad de la prediccion para cada uno de los valores en el 2018
plot(1-(abs(pred.niv-as.vector(indicador_global$GSPC.Open[502:nrow(data3)]))/as.vector(indicador_global$GSPC.Open[502:nrow(data3)])), type="l", ylab="Efectividad", xlab="Indice", main="Efectividad de la prediccion: 2018")

#Promedio de efectividad
mean(1-(abs(pred.niv-as.vector(indicador_global$GSPC.Open[502:nrow(data3)]))/as.vector(indicador_global$GSPC.Open[502:nrow(data3)])))*100

#Revisamos la precision
mp.effect.niv<-sign(pred.niv[2:46]-as.vector(indicador_global$GSPC.Open[502:546]))
real.effect2<-sign(diff(as.vector(indicador_global$GSPC.Open[502:547])))
100*sum(mp.effect.niv==real.effect2)/length(real.effect2)

#Revisamos la efectividad para el lunes 19 de noviembre
(1-abs(pred.niv[46]-as.vector(indicador_global$GSPC.Open)[548])/as.vector(indicador_global$GSPC.Open)[548])*100

