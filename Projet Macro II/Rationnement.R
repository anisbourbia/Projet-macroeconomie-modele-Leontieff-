#recherche d'optimum avec rationnement de l'industrie

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

#II- situation de référence et paramètres

#matrice des comptes sociaux:
m[1,1]<-10; m[1,2]<-40; m[1,3]<-6;
m[2,1]<-43; m[2,2]<-50; m[2,3]<-10;
m[3,1]<-15; m[3,2]<-40; m[3,3]<-2;

v[1,1]<-42.5; v[1,2]<-19; v[1,3]<-10.2;
v[2,1]<-7.5; v[2,2]<-19; v[2,3]<-40.8;

#structure sociale: clés de répartition des revenus des facteurs
s[1,1]<-0.29;s[1,2]<-0.70;s[1,3]<-0.01
s[2,1]<-0.01;s[2,2]<-0.5;s[2,3]<-0.49

#matrice de Leontieff et va unitaires:
#calcul de l'offre en volume par secteur
for (i in 1:3) {Y[i]<-rep(1,3)%*%m[1:3,i]+rep(1,2)%*%v[1:2,i]}

#calcul de la matrice des coeff de Leontief-volumes
for (i in 1:3) {for (j in 1:3) {a[i,j]<-m[i,j]/Y[j]}}

#calcul des va unitaires
for (i in 1:2) {for (j in 1:3) {r[i,j]<-v[i,j]/Y[j]}}
V<-rep(1,2)%*%r

#fonction d'utilité Stone-Geary
#paramètres (calculés à partir du calibrage)
C[1]<-40
C[2]<-0.8333333
C[3]<-5.6566208
E[1]<-0.2378092
E[2]<-0.6936101
E[3]<-0.0685807

stone_geary<-function(di) {
  u<-((di[1]-C[1])^E[1])*((di[2]-C[2])^E[2])*((di[3]-C[3])^E[3])
  return(u)
}
#fonction d'utilité sociale, en input la demande par bien et par classe
fsoc<-function(d) {
  U<-0
  for (j in 1:q) {U<-U+s[1,j]*stone_geary(d[,j]/s[1,j])}
  return(U)
}

#III- définition des fonctions nécessaires:

#fonction inversion de matrice
inv<-function(mat){
  n<-nrow(mat)
  Idn<-matrix(0,n,n)
  matmoins1<-matrix(0,n,n)
  for (i in 1:n) {Idn[i,i]<-1}
  for (i in 1:n) {
    matmoins1[,i]<-solve(mat,Idn[,i])
  }
  return(matmoins1)
}

#gradient de la fonction d'utilité sociale
gradfsoc<-function(d){
  g<-matrix(0,n,q)
  for (j in 1:q) {
    for (i in 1:n){
      g[i,j]<-E[i]*(1/(d[i,j]/s[1,j]-C[i]))*stone_geary(d[,j]/s[1,j])
    }
  }
  return(g)
}

#fonction de demande sans rationnement, par bien et par classe
dem<-function(M,P) {
  depmin<-t(C)%*%P
  d<-matrix(0,n,q)
  for (k in 1:q) {
    for (i in 1:n){
      d[i,k]<-(C[i]+E[i]*((M[k]/s[1,k])/P[i]-depmin/P[i]))*s[1,k]
    }
  }
  return(d)
}

#fonction de demande avec rationnement à 55,25 de la demande finale de l'industrie
dem_r<-function(d,M,P)  {
  mat_cont<-matrix(0,4,9)
  Id9<-matrix(0,9,9)
  for (i in 1:9) {Id9[i,i]<-1}
  mat_cont[1,1:3]<-P
  mat_cont[2,4:6]<-P
  mat_cont[3,7:9]<-P
  mat_cont[4,]<-c(0,1,0,0,1,0,0,1,0)
  vect_cont<-rbind(M,55.25)
  
  #opérateur de projection sur le sous espace vectoriel (sev) des contraintes saturées
  Proj<-mat_cont%*%t(mat_cont)
  Proj<-inv(Proj)
  Proj<-Id9-t(mat_cont)%*%Proj%*%mat_cont
  
  #recherche d'un point d'amorçage proche de la situation initiale et respectant les contraintes.
  #on modifie successivement les consommations du bien 2 de sorte à satisfaire le rationnement,
  #puis on modifie les consommations des bien 1 puis 3 pour satisfaire les contraintes budgétaires.
  #facteur correctif sur le bien 2
  d0<-d #d0: futur point d'amorçage
  corr2<-as.numeric(vect_cont[4]/(d[2,]%*%rep(1,3)))
  d0[2,]<-corr2*d[2,]
  #report du budget dégagé sur le bien 1
  d0[1,]<-d[1,]+d[2,]*(1-corr2)*P[2]/P[1]
  #ajustement du budget du cons 1 par le bien 3
  d0[3,1]<-d0[3,1]-as.numeric((d0[,1]%*%P-vect_cont[1])/P[3])
  #ajustement du budget du cons 2 par le bien 3
  d0[3,2]<-d0[3,2]-as.numeric((d0[,2]%*%P-vect_cont[2])/P[3])
  #ajustement du budget du cons 3 par le bien 3
  d0[3,3]<-d0[3,3]-as.numeric((d0[,3]%*%P-vect_cont[3])/P[3])
  
  #vectorisation de l'argument de gradfsoc (pour pouvoir utiliser les opérations matricielles)
  vect_gradfsoc<-function(vdx){
    dx<-matrix(vdx,n,q)
    res_gradfsoc<-gradfsoc(dx)
    res_gradfsoc<-as.vector(res_gradfsoc)
    return(res_gradfsoc)
  }
  #descente du gradient projeté
  mu <- 0.01
  y <- as.vector(d0)
  nn <- 1
  err <- 1
  while (err > 10^-8 & nn <= 1000000) { 
    x<-y
    u <- Proj%*%vect_gradfsoc(x)
    #normu<-sqrt(sum(u^2))
    y <- x+mu*u
    err <- sqrt(sum((x-y)^2))
    nn <- nn+1
  }
  my<-matrix(y,3,3)
  return(my)
}

#fonction de calcul des revenus des classes à partir de la demande finale à prix unitaires constants
rev<-function(d,r) {
  #calcul de la demande finale totale par bien
  D<-d%*%rep(1,q)
  
  #mise à jour de l'offre totale par bien en volume
  Id3<-matrix(0,3,3)
  for (i in 1:3) {Id3[i,i]<-1}
  
  Y<-inv(Id3-a)%*%D
  
  #composantes des valeurs ajoutées
  for(j in 1:p) { for (i in 1:n) {v[j,i]<-r[j,i]*Y[i]}}
  
  #partage des revenus
  M<-t(s)%*%v%*%rep(1,n)
  
  return(M)
}

#IV- calculs d'initialisation

#partage des revenus dans la situation de référence
M<-t(s)%*%v%*%rep(1,n)

#calcul de la demande finale par bien et par classe sans rationnement
d<-dem(M,P)

#calcul de la demande finale par bien et par classe avec rationnement de l'industrie à 55,25
#dr<-dem_r(d,M,P)

#V- calcul de l'équilibre

#demande->revenus->demande...
F<-function(vd){
  d<-matrix(vd,3,3)
  M<-rev(d,r)
  dr<-dem_r(d,M,P)
  vdr<-as.vector(dr)
  return(vdr)
}
#descente du gradient numérique
vd <- as.vector(d)
vdr<- F(vd)
delta<-1e-05 #pour le gradient numérique
mu <- 0.01   #pas de la descente
y<-vdr
nn <- 1
err <- 1
GradF<-matrix(0,9,9)
u<-rep(0,9) #u: gradient de la fonction norme(F(x)-x)^2
Id9<-matrix(0,9,9)
for (i in 1:9) {Id9[i,i]<-1}
while (err > 10^-10 & nn <= 600) { 
  x<-y
  #calcul du gradient numérique
  Fx<-F(x)
  for (i in 1:9) { GradF[,i]<-(F(x+delta*Id9[i,])-Fx)/delta }
  for (i in 1:9) { u[i] <- 2*(GradF[,i]-Id9[,i])%*%(Fx-x) }
  normu<-sqrt(sum(u^2))
  print('norme de u:')
  print(normu)
  y <- x-mu*u
  err <- sqrt(sum((x-y)^2))
  print('y:')
  print(y)
  print(nn)
  print(err)
  flush.console()
  nn <- nn+1
}


#rm(list = ls(all = TRUE))#pour tout effacer