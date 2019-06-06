rm(list=ls())  #Limpiamos la memoria

set.seed(100) #Fijamos una semilla para permitir la reproducibilidad

#Funcion que genera las matrices para cada ensamble y obtiene los valores propios de dichas matrices
generadora<-function(N,M,beta){ 
  if(beta == 1){
    H<- matrix(rnorm(n = N*M,mean = 0,sd = 1),nrow = N,ncol = M )
    H_t<-t(H)
    W<-H%*%H_t
    VP<-eigen(W)$values
    }
  if(beta == 2) {
    Real<- rnorm(N*M,mean = 0,sd = 1)
    Imaginaria<- rnorm(N*M,mean = 0,sd = 1)
    H<- matrix(complex(real = Real,imaginary = Imaginaria),nrow = N,ncol = M )
    H_t<- t(matrix(complex(real = Real,imaginary = -Imaginaria),nrow = N,ncol = M ))
    W<-H%*%H_t
    VP<-eigen(W)$values
  }
  if(beta == 4){
    library(QZ)
    Real<- rnorm(N*M,mean = 0,sd = 1)
    Imaginaria<- rnorm(N*M,mean = 0,sd = 1)
    A<- matrix(complex(real = Real,imaginary = Imaginaria),nrow = N,ncol = M )
    A_c <- matrix(complex(real = Real,imaginary = -Imaginaria),nrow = N,ncol = M )
    
    Real<- rnorm(N*M,mean = 0,sd = 1)
    Imaginaria<- rnorm(N*M,mean = 0,sd = 1)
    B<- matrix(complex(real = Real,imaginary = Imaginaria),nrow = N,ncol = M )
    B_c<- matrix(complex(real = Real,imaginary = -Imaginaria),nrow = N,ncol = M )
    aux1<-cbind(A,B)
    aux2<-cbind(-B_c,A_c)
    H<-rbind(aux1,aux2)
    W<-H%*%H(H)
    VP<-eigen(W)$values
  }
  return(VP)
}



#######Esta parte dibuja el histograma con la N,M y beta de parametros (c = N/M con N<M)
N<-100
M<-200
beta<-1
L<-array(0,dim = c(length(generadora(N,M,beta)),0))

library(RMTstat)
for (i in 1:5000) {
  L<-cbind(L,generadora(N,M,beta))
}

#Realizamos la prueba de Kolmogorov-Smirnov
temp<-rmp(5000, ndf=200, pdim = 100)
ks.test(L/200, temp)

#Realizamos un qqplot para observar graficamente el ajuste de las 2 distribuciones
qqplot(L/200, temp, col="dodgerblue1", main="Qqplot de las 2 distribuciones", xlab = "Percentiles simulada", ylab="Percentiles teorica")

#Dibujamos las distribuciones
plot(density(L/200), main="Densidades: simulada y M-P")
lines(density(temp), col="red")
legend(2, 0.8, legend=c("Simulada-escalada", "M-P"), col=c("black", "red"), lty=1)

#Generamos la grafica que se establece en el libro
#Funcion que genera la distribucion teorica
p_MP<-function(y){
  zeta_menos<-(1-c^(-1/2))^2
  zeta_mas<-(1+c^(-1/2))^2
  1/(2*pi*y)*sqrt((y-zeta_menos)*(zeta_mas-y))
}

#Distribucion teorica con respecto a cada beta
p<-function(x,beta,N){
  1/(beta*N)*p_MP(x/(beta*N))
}

#Graficas con Beta=1,2,4 y proporcion 1/2
N<-100
M<-200
c<-N/M
VP<-generadora(N=N,M=M,beta=1)
temp<-sample(VP, 20)
plot(temp/100,p(temp,beta = 1,N = 100)*100,pch= 20,col="darkmagenta",xlim = c(0,16),ylab = "p(x)", main="Distribucion de M-P con c=1/2 y c=1/8", sub="Tres ensembles")
lines(VP/100, p_MP(VP/100), col="dodgerblue4")

VP<-generadora(N=N,M=M,beta=2)
temp<-sample(VP, 20)
points(temp/(2*100),p(temp,beta = 2,N = 100)*(2*100),pch= 18,col="deepskyblue4")

VP<-generadora(N=N,M=M,beta=4)
temp<-sample(VP, 20)
points(temp/(4*100),p(temp,beta = 4,N = 100)*(4*100),pch= 4,col="gray24")

#Graficas con Beta=1,2,4 y proporcion 1/8
N<-100
M<-800
c<-N/M
VP<-generadora(N=N,M=M,beta=1)
temp<-sample(VP, 20)
points(temp/100,p(temp,beta = 1,N = 100)*100,pch= 20,col="palegreen4")
lines(VP/100, p_MP(VP/100), col="darkorange1")

VP<-generadora(N=N,M=M,beta=2)
temp<-sample(VP, 20)
points(temp/(2*100),p(temp,beta = 2,N = 100)*(2*100),pch= 18,col="cadetblue")

VP<-generadora(N=N,M=M,beta=4)
temp<-sample(VP, 20)
points(temp/(4*100),p(temp,beta = 4,N = 100)*(4*100),pch= 4,col="firebrick4")

legend("topright", inset=0,c("c = 1/2", "beta= 1","beta= 2", "beta=4", "c=1/8", "beta = 1", "beta =2","beta =4"), col=c("dodgerblue4","darkmagenta","deepskyblue4","gray24","darkorange1","palegreen4","cadetblue", "firebrick4"), horiz=F, cex=1,lty=c(1, NA, NA, NA, 1, NA, NA, NA), pch = c(NA,20,18,4,NA,20,18,4),xjust = 0,yjust = 0)
#######



