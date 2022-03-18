#calibrage de l'utilité stone-geary
#inputs: d = matrice des consommations individuelles en valeur par bien x classe
#m = revenus individuels par classe
#C[1] = consommation alimentaire minimale vitale

#I - Déclarations

n <- 3		  #n doit être >= 2
q <- 5
d <- matrix (0,n,q) #consommation individuelle des classes en valeur
E <- rep(0,n)       #elasticités
C <- rep(0,n)	  #consommations minimales vitales
m <- rep(0,q)	  #budget individuelle par classe

#II- Données de références

C[1] <- 0 #on doit choisir une comp. de C car sans le système est indéterminé
d <- matrix(c(5.8,7.4,10.4,13.6,24.8,2.4,4,7.6,12.8,38.2,0.375,0.75,1.25,4.5,5.125),n,q)
m <- c(8.575	,12.150,	19.250,	30.900,	68.125)






#III- Fonctions

m_c <- function(x) { #somme des carrés des différences à minimiser
  #x est la concaténation de C[2:n] et de E
  mc <- 0
  #deconcatenation:
  C[2:n] <- x[1:(n-1)]
  E <- x[n:(2*n-1)]
  for (i in 1:n) {
    for (j in 1:q){
      mc <- mc+(d[i,j]-C[i]-E[i]*(m[j]-sum(C)))^2
    }
  }
  return(mc)
}

grad_mc <- function(x) { #gradient de m_c
  #deconcatenation:
  C[2:n] <- x[1:(n-1)]
  E <- x[n:(2*n-1)]
  A <- matrix(0,n,q)
  for (i in 1:n) {
    for (j in 1:q){
      A[i,j] <- d[i,j]-C[i]-E[i]*(m[j]-sum(C))
    }
  }
  B <- rep(0,q)
  B <- m - sum(C)
  grad <- rep(0,(2*n))
  for (i in 1:n) {
    for (ip in 1:n) {grad[i] <- grad[i] + 2*sum(A[ip,])*E[ip]}
    grad[i] <- grad[i] - 2*sum(A[i,])
    grad[n+i] <- grad[n+i] - 2*A[i,]%*%B
  }
  grad <- grad[2:(2*n)] #on supprime la première composante car C1 est connu
  return(grad)
}

dem <- function(m,C,E) { #calcule la matrice des consommations estimées à partir des revenus et des coeffs C & E
  d_es <- matrix(0,n,q)
  for (i in 1:n) {
    for (j in 1:q){
      d_es[i,j] <- d_es[i,j] + C[i] + E[i]*(m[j]-sum(C))
    }
  }
  return(d_es)
}

m_c1 <- function(C1) {
  #écart entre cons reelle et estimée pour mieux choisir C1
  #on utilise: C1 et x la sortie de la boucle while; m le vecteur des revenus; d la matrice des consommations réelles
  #il s'avère que ça sert à rien car le choix de C1 n'augmente ni ne diminue l'erreur
  #il ne fait que decaler les autres Ci
  #fonction à supprimer pour les élèves
  C <- rep(0,n)
  C[1] <- C1
  C[2:n] <- x[1:(n-1)]
  E <- x[n:(2*n-1)]
  d_es <- dem(m,C,E) #calcul des consommations simulées
  #calcul des Ci à partir de C1 et d_es
  sigmac <- m[1]-(d_es[1,1]-C[1])/E[1] #somme des Ci
  for (i in 2:n) {
    C[i] <- d_es[i,1]-E[i]*(m[1]-sigmac)
  }
  mc <- 0
  for (i in 1:n) {
    for (j in 1:q){
      mc <- mc+(d[i,j]-C[i]-E[i]*(m[j]-sum(C)))^2
    }
  }
  return(mc)
}

#IV- Point d'amorçage (calculé sur les individus j = 1 et j = 2)

#Point d'amorçage
#calcul des Ei
Ei <- matrix(0,q,q)
meanE <- rep(0,n)
for (i in 1:n) {
  compteur <- 0     #pour calculer la moyenne des Ei
  sumEi <- 0        #pour calculer la moyenne des Ei
  sumsquaresEi <- 0 #pour calculer la variance des Ei
  for (j1 in 1:q) {
    for (j2 in 1:q) {
      if (j1 != j2) {
        Ei[j1,j2] <- (d[i,j2]-d[i,j1])/(m[j2]-m[j1])
        compteur <- compteur + 1
        sumEi <- sumEi + Ei[j1,j2]
        sumsquaresEi <- sumsquaresEi + (Ei[j1,j2])^2
      }
    }
  }
  meanE[i] <- sumEi / compteur
  varEi <- sumsquaresEi / compteur - (meanE[i])^2
  relative_devEi <- sqrt(varEi)/meanE[i]
  print(Ei) #pour voir la stabilité des Ei
  print(relative_devEi)
}
E <- meanE
#pour avoir sum(Ei)=1
for (i in 1:n) {
  E[i] <- E[i]/sum(E)
}
#calcul des Ci
sigmac <- m[1]-(d[1,1]-C[1])/E[1] #somme des Ci
for (i in 2:n) {
  C[i] <- d[i,1]-E[i]*(m[1]-sigmac)
}

#V- Descente du gradient
x <- c(C[2:n],E) #initialisation au point d'amorçage
mu <- 1e-1	     #on peut accélerer à mu = 1e-1 par moment
mc <- m_c(x)
nn <- 1
dmc <- 1
f<-1e4 #pour le rescaling de certaines composantes du gradient
while (dmc > 10^-11 & nn <= 1000000) {
  u <- grad_mc(x)
  #projection du gradient sur le sev de la contrainte sum(E)=1
  moyE <- mean(u[n:(2*n-1)])
  u[n:(2*n-1)] <- u[n:(2*n-1)] - moyE
  u[n:(2*n-1)] <- u[n:(2*n-1)]/f #rescaling des 3 dernières composantes du gradient
  x[n:(2*n-1)] <- x[n:(2*n-1)]*f
  normu<-sqrt(sum(u^2))
  x <- x-mu*u
  x[n:(2*n-1)] <- x[n:(2*n-1)]/f
  dmc <- mc-m_c(x)
  mc <- mc-dmc
  err <- mu*sqrt(sum(u^2))
  print(mc)
  print(dmc)
  print(normu)
  print(err)
  print(nn)
  flush.console()
  nn <- nn+1
}

