#variation des taxes

#I- structures

n<-3 #nombre de biens
p<-2 #nombre de facteurs
q<-5 #nombre de classes sociales
m<-matrix(0,n,n)  #consommations intermediaires en volume
a<-matrix(0,n,n)  #coefficients de Leontieff
r<-matrix(0,p,n)	#rémunérations unitaires des facteurs
v<-matrix(0,p,n)  #composantes des valeurs ajoutées
d<-matrix(0,n,q)	#demande finale par classe sociale en volume, un bien par ligne
s<-matrix(0,p,q)  #clés de répartition des revenus par facteur et par classe
C<-rep(0,n)		#consommations minimales vitales selon la fonction d'utilité Stone-Geary
E<-rep(0,n)		#exposants de la fonction d'utilité Stone-Geary
P<-rep(1,n)		#prix
V<-rep(0,n)		#valeur ajoutée unitaire, par secteur
Y<-rep(0,n)		#offre en volume, par secteur
D<-rep(0,n)		#demande finale en volume par secteur
M<-rep(0,q)		#revenu total par classe sociale

to <- matrix(0,n+p,n+q) #matrice des taux des taxes/subventions indirectes
ro <- rep(0,n+q)		#vecteur des proportions des redistributions directes par rapport aux recettes fiscales indirectes totales

#II- situation de référence et paramètres

#matrice des comptes sociaux:
m[1,1]<-10; m[1,2]<-40; m[1,3]<-6;
m[2,1]<-43; m[2,2]<-50; m[2,3]<-10;
m[3,1]<-15; m[3,2]<-40; m[3,3]<-2;

v[1,1]<-42.5; v[1,2]<-19; v[1,3]<-10.2;
v[2,1]<-7.5; v[2,2]<-19; v[2,3]<-40.8;

#structure sociale: clés de répartition des revenus des facteurs
s[1,1]<-0.2;s[1,2]<-0.20;s[1,3]<-0.2;s[1,4]<-0.2;s[1,5]<-0.2
s[2,1]<-0.013;s[2,2]<-0.028;s[2,3]<-0.095;s[2,4]<-0.23;s[2,5]<-0.63

#matrice de Leontieff et va unitaires:
#calcul de l'offre en volume par secteur
for (i in 1:3) {Y[i]<-rep(1,3)%*%m[1:3,i]+rep(1,2)%*%v[1:2,i]}

#calcul de la matrice des coeff de Leontief-volumes
for (i in 1:3) {for (j in 1:3) {a[i,j]<-m[i,j]/Y[j]}}

#calcul des remunérations unitaires et va unitaires
for (i in 1:2) {for (j in 1:3) {r[i,j]<-v[i,j]/Y[j]}}
V<-rep(1,2)%*%r

#taux de taxe et de redistribution
to[(n+2),] <- rep(0.05,(n+q))
ro[1] <- 0.2
ro[2] <- 0.1
ro[4] <- 0.7

#paramètres de la fonction d'utilité Stone-Geary(calculés à partir du calibrage)
C[1]<-25
C[2]<-5.8596449
C[3]<-0.5248405
E[1]<-0.3481126
E[2]<-0.5343837
E[3]<-0.1175019

#matrice des coefficients de Leontief ttc et taxes unitaires
ap <- a*(matrix(1,n,n)+to[1:n,1:n])
toa <- a*to[1:n,1:n]

#taux des redistributions directes aux producteurs et aux consommateurs
roy <- ro[1:n]
roc <- ro[(n+1):(n+q)]

#revenus des consommateurs à taxes nulles
M <- t(s)%*%v%*%rep(1,n)

#consommation finale à taxes nulles
depmin<-t(C)%*%P
d0<-matrix(0,n,q)
for (j in 1:q) {
  for (i in 1:n){
    d0[i,j]<-(C[i]+E[i]*((M[j]/s[1,j])/P[i]-depmin/P[i]))*s[1,j]
  }
}

#III- définition des fonctions nécessaires:

#matrice identité	
Id<-function(nn){
  Idn<-matrix(0,nn,nn)
  for (i in 1:nn) {Idn[i,i]<-1}
  return(Idn)
}
Idn <- Id(n)

#fonction inversion de matrice
inv<-function(mat){
  n<-nrow(mat)
  matmoins1<-matrix(0,n,n)
  for (i in 1:n) {
    matmoins1[,i]<-solve(mat,Idn[,i])
  }
  return(matmoins1)
}

#fonction pour calculer la nouvelle demande et d'autres parametres
dem <- function(d){
  #mise à jour des parametres selon le nouveau Y
  D<-d%*%rep(1,q)
  Y<-inv(Idn-a)%*%D
  alpha <- (a*to[1:n,1:n])%*%Y + (d*to[1:n,(n+1):(n+q)])%*%rep(1,q)
  B <- t(alpha%*%t(roy/Y))
  
  #taxe sur l'utilisation des facteurs et taxe sur les revenus
  Tf <- as.vector(to[n+1,1:n]*r[1,1:n])%*%as.vector(Y) + as.vector(to[n+2,1:n]*r[2,1:n])%*%as.vector(Y)
  Tr <- (as.vector(t(r[1,1:n]*Y)%*%rep(1,n))*s[1,1:q])%*%to[(n+1),1:q] + (as.vector(t(r[2,1:n]*Y)%*%rep(1,n))*s[2,1:q])%*%to[(n+2),1:q]
  
  #prix
  Vp <- matrix(1,1,p)%*%((matrix(1,p,n)+to[(n+1):(n+p),1:n])*r)
  P <- solve(Idn-t(ap)+B,t(Vp)-as.vector(Tf+Tr)*(roy/Y))
  
  #taxe sur la demande et taxe totale
  Td <- t(alpha)%*%P
  T <- Td+Tf+Tr
  
  #mise à jour des remunérations des facteurs
  for (m in 1:p) {for (i in 1:n) {v[m,i]<-r[m,i]*Y[i]}} 
  #revenus des consommateurs apres taxes sur les revenus des consommateurs
  M <- t((matrix(1,p,q)-to[(n+1):(n+p),(n+1):(n+q)])*s)%*%v%*%rep(1,n)
  #revenus des consommateurs apres redistribution
  M <- M + as.vector(T)*roc
  
  #consommation en volume selon les revenus M et les prix P
  depmin<-t(C)%*%P
  d_out<-matrix(0,n,q)
  for (j in 1:q) {
    for (i in 1:n){
      d_out[i,j]<-(C[i]+E[i]*((M[j]/s[1,j])/P[i]-depmin/P[i]))*s[1,j]
    }
  }
  sortie_dem <- list(d_out,P,Y,Tr,Tf,alpha)
  return(sortie_dem)
}

#IV- calcul du gradient numerique
mc <- function(d){
  sortie_dem <- dem(d)
  dp <- sortie_dem[[1]]
  x<-sqrt(sum((d-dp)^2))
  return(x)
}
gradnum <- function(d,delta) {
  x <- matrix(0,n,q)
  y <- mc(d)
  for (i in 1:n) { for (j in 1:q) {
    dpij <- d
    dpij[i,j] <- dpij[i,j] + delta
    ypij <- mc(dpij)
    x[i,j]<-(ypij-y)/delta
  }}
  return(x)
}
d_in <- d0
gradn <- gradnum(d_in,1e-03)	

#V- descente avec gradient numerique
mu <- 1e-4
delta <- 1e-05
nn <- 1
err <- 1
deseq <- 10000
delta_deseq <- -1
d_in <- d0
while (err > 10^-8 & nn <= 1000000 & delta_deseq <= 0) {
  gradn <- gradnum(d_in,delta)
  normgrad<-sqrt(sum(gradn^2))
  err <- mu*normgrad
  deseqnmoins1 <- deseq
  deseq <- mc(d_in) #desequilibre
  delta_deseq <- deseq-deseqnmoins1 #pour arreter la descente au cas où deseq augmente
  print(deseq) 
  print(normgrad)
  print(err)
  print(nn)
  flush.console()
  d_in <- d_in - mu*gradn
  nn <- nn+1
}  
