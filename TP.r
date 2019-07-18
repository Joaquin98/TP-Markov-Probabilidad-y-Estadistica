library(markovchain)
library(ggplot2)
library(MASS)
#Distribución inicial
pi <- c(2/5,1/5,2/5)
#Estados
states <- c("A","B","C")
#Matriz de transición
transition <- matrix(c(0,1/3,2/3,1/4,3/4,0,2/5,0,3/5),byrow=TRUE,nrow=3)
#Creo la cadena de markov. Chequea automáticamente las propiedades necesarias
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio5")

#Probabilidad de ir de un estado a otro
transitionProbability(mc,"A","B")
mc[1,2]

#Grafo asociado
plot(mc)

#¿Es irreducible?
is.irreducible(mc)
#Período de la cadena
period(mc)
#Estados transitorios
transientStates(mc)
#Distribución estacionaria
steadyStates(mc)
#Estados absorbentes
absorbingStates(mc)

#Matriz de transición en un paso
transition1 <- mc
#Matriz de transición en dos pasos
transition2 <- mc * mc
#Matriz de transición en tres pasos
transition3 <- mc ^ 3

#Simula una trayectoria de 100 estados partiendo desde C
rmarkovchain(mc,n=100)

# Ejercicio 1

# Distribución inicial
pi <- c(0,0,1,0,0,0)
#Estados
states <- c("0","1","2","3","4","5")
p = 0.5
#Matriz de transición
transition <- matrix(c(1,0,0,0,0,0,
                       1-p,0,p,0,0,0,
                       0,1-p,0,p,0,0,
                       0,0,1-p,0,p,0,
                       0,0,0,1-p,0,p,
                       0,0,0,0,0,1),byrow=TRUE,nrow=6)

#Creo la cadena de markov. Chequea automáticamente las propiedades necesarias
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 1")

#Grafo asociado
plot(mc)

transitionProbability(mc^100,"2","0")

# gamblersruin.R
# Example 1.11

# gamble(k, n, p)
#   k: Gambler's initial state
#   n: Gambler plays until either $n or Ruin
#   p: Probability of winning $1 at each play
#   Function returns 1 if gambler is eventually ruined
#                    returns 0 if gambler eventually wins $n

gamble <- function(k,n,p) {
  while (k > 0 & k < n) 
    k <- k + sample(c(-1,1),1,prob=c(1-p,p))
  k == 0
}   
mean(replicate(1000, gamble(2, 5, 1/2)))
mean(replicate(1000, gamble(2, 10, 1/2)))
mean(replicate(1000, gamble(2, 50, 1/2)))
mean(replicate(100, gamble(250, 500, 1/2)))


canonicForm(mc)

Q = matrix(c( 0.0, 0.5, 0.0, 0.0,
              0.5, 0.0, 0.5, 0.0,
              0.0, 0.5, 0.0, 0.5,
              0.0, 0.0, 0.5, 0.0),byrow=TRUE,nrow=4)

B = matrix(c( 0.5, 0.0, 
              0.0, 0.0,
              0.0, 0.0, 
              0.0, 0.5 ),byrow=TRUE,nrow=4)

S = ginv(diag(4)-Q)
G = S %*% B
G


r = rep(0,1000)
for(i in 1:1000){
  x = rmarkovchain(mc,n=100,t0="2")
  for(j in 1:100){
    if(x[j] == "0"){
      r[i] = 1
      break  
    } 
    if(x[j] == "5") break
  }
}
mean(r)


generarMatrizJuego <- function(n,p){
  n = n+1
  P = matrix(rep(0,n*n),byrow=TRUE,nrow=n)
  P[1,1] = 1
  P[n,n] = 1
  for(i in 2:(n-1)){
    P[i,i-1] = 1-p
    P[i,i+1] = p
  }
  P
}

states <- c("0","1","2","3","4","5","6","7","8","9","10")
transition <- generarMatrizJuego(10,0.5)
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 1")
transitionProbability(mc^100,"2","0")

canonicForm(mc)

Q = matrix(c( 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0),byrow=TRUE,nrow=9)

B = matrix(c( 0.5, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.0,
              0.0, 0.5),byrow=TRUE,nrow=9)

S = ginv(diag(9)-Q)
G = S %*% B
G


# Ejercicio 2

# Distribución inicial
pi <- c(0,0,0,0,1)
#Estados
states <- c("80","135","139","145","NA")
#Matriz de transición
transition <- matrix(c(0,0,0,0,1,
                       0,8/13,3/13,1/13,1/13,
                       1/16,3/16,3/8,1/4,1/8,
                       0,1/11,4/11,5/11,1/11,
                       0,1/8,1/2,1/8,1/4),byrow=TRUE,nrow=5)

#Creo la cadena de markov. Chequea automáticamente las propiedades necesarias
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 2")

#Grafo asociado
plot(mc)

# Probabilidad de estados luego de 2 semanas.
semana2 = pi * (mc^2)
semana2
max(semana2) # 0.3868007 -> Puerto 139
min(semana2) # 0.03125 -> Puerto 80

#Distribución estacionaria
steadyStates(mc)

# Ejercicio 3

# Distribución inicial
pi <- c(0,0,0,0,1)
#Estados
states <- c("s","t","m","a","r")
states <- c("0","1","2","3","4")
#Matriz de transición
transition <- matrix(c(0.84,0.11,0.01,0.04,0,
                       0.03,0.8,0.04,0.1,0.03,
                       0.01,0.15,0.7,0.07,0.07,
                       0.03,0.19,0.02,0.75,0.01,
                       0.03,0.09,0.05,0,0.83),byrow=TRUE,nrow=5)

#Creo la cadena de markov. Chequea automáticamente las propiedades necesarias
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 3")

#Grafo asociado
plot(mc)

is.irreducible(mc)
period(mc)
steadyStates(mc)


n = 1
sim = matrix(rep(0,len = n*100),ncol=n)
for(i in (1:n))
sim[,i] = rmarkovchain(mc,n=100)
matplot(sim, type='l', lty=1, col=3:5, ylim=c(0,4), ylab='Actividades', xlab='Tiempo')


# Ejercicio 4

states <- c("a","b","c","d","e","f")
transition <- matrix(c(0,1,0,0,0,0,
                       1/4,0,1/4,1/4,1/4,0,
                       0,1/4,0,1/4,1/4,1/4,
                       0,1/4,1/4,0,1/4,1/4,
                       0,1/3,1/3,1/3,0,0,
                       0,0,1/2,1/2,0,0 ),byrow=TRUE,nrow=6)
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 4")



is.irreducible(mc)
period(mc)
transientStates(mc)
steadyStates(mc)

  cont = 0
  for(i in 1:10){
    x = rmarkovchain(mc,n=50)
    if(tail(x, n = 1) == "c") cont = cont + 1
  }
  cont


plot(mc)


# Ejercicio 5

states <- c("0","1","2","3","4")
transition <- matrix(c(0,1,0,0,0,
                       0,1/4,3/4,0,0,
                       0,0,2/4,2/4,0,
                       0,0,0,3/4,1/4,
                       0,0,0,0,1 ),byrow=TRUE,nrow=5)
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 5")

transitionProbability(mc^6,"0","4")


F(6,"0","4",states) + F(5,"0","4",states) + F(4,"0","4",states)

r = rep(0,1000)
for(i in 1:1000){
  x = rmarkovchain(mc,n=6,t0="0")
  for(j in 1:6){
    if(x[j] == "4"){
      r[i] = 1
      break  
    } 
  }
}
mean(r)

# Ejercicio 6

Binpy <- function(n,p,y)
  choose(n,y) * p^y * (1-p)^(n-y)

P = matrix( rep( 0, len=36), nrow = 6)
for(i in 1:6)
  for(j in 1:6)
    if(i<=j) P[i,j] = Binpy(5-i+1,1/6,j-i)
P

pi <- c(1,0,0,0,0,0)
states <- c("0","1","2","3","4","5")
transition <- P
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 4")

tercerPaso = mc^3

tercerPaso[1,6]

mc^100

is.irreducible(mc)
steadyStates(mc)

F(3,"0","5",states) + F(2,"0","5",states) + F(1,"0","5",states)

r = rep(0,1000)
for(i in 1:1000){
  x = rmarkovchain(mc,n=3,t0="0")
  for(j in 1:3){
    if(x[j] == "5"){
      r[i] = 1
      break  
    } 
  }
}
mean(r)

# Ejercicio 7

pi <- c(1/4,1/4,1/4,1/4)
states <- c("0","1","2","3")
transition <- matrix(c(0,1,0,0,
                       3/4,0,1/4,0,
                       0,3/4,0,1/4,
                       0,0,1,0),byrow=TRUE,nrow=4)
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 7")

transitionProbability(mc^3,"2","1")

transitionProbability(mc^2,"1","3") * ((pi*(mc^3))[2])

period(mc)

P(3,"2","1")

# Ejercicio 8
pi <- c(0,1,0,0,0,0)
states <- c("A","1","2","3","4","R")
transition <- matrix(c(1,0,0,0,0,0,
                       6/100,3/100,91/100,0,0,0,
                       6/100,0,3/100,91/100,0,0,
                       4/100,0,0,3/100,93/100,0,
                       4/100,0,0,0,3/100,93/100,
                       0,0,0,0,0,1),byrow=TRUE,nrow=6)
mc <- new("markovchain",states = states,transitionMatrix = transition,name="Ejercicio 8")

# Tipo de grafo
is.irreducible(mc)
absorbingStates(mc)

# Distribución límite
steadyStates(mc)
markovchain ::: sort(mc)

plot(mc)

# Promedio de años
sum = 0
l = 10000
for(i in 1:l){
  x = rmarkovchain(mc,n=20,t0=1)
  for(j in 1:length(x)){
    if(x[j] == "A" || x[j] == "R"){
      sum = sum + j - 1
      break
    }
  }
}
res = 1 + sum/l

res

Q = matrix(c( 3/100 , 91/100 , 0 , 0 ,
              0 , 3/100 , 91/100 , 0,
              0 , 0 , 3/100 , 93/100,
              0 , 0 , 0 , 3/100 ),byrow=TRUE,nrow=4)

N = ginv(diag(4)-Q)
sum = 0
for(i in 1:4)
  sum = sum + N[1,i] 
sum

transientStates(mc)

P = matrix(c(1,0,0,0,0,0,
         6/100,3/100,91/100,0,0,0,
         6/100,0,3/100,91/100,0,0,
         4/100,0,0,3/100,93/100,0,
         4/100,0,0,0,3/100,93/100,
         0,0,0,0,0,1),byrow=TRUE,nrow=6)

Q = matrix(c( 3/100 , 91/100 , 0 , 0 ,
              0 , 3/100 , 91/100 , 0,
              0 , 0 , 3/100 , 93/100,
              0 , 0 , 0 , 3/100 ),byrow=TRUE,nrow=4)

N = ginv(diag(4)-Q)
sum(N[1,])

F(10,"1","4",states)

# Ejercicio 9

Pij <- function(i,j,k)
  choose(k,j) * (i/k)^j * (1-(i/k))^(k-j)

P = matrix(rep( 1, len=36), nrow = 6)

for(i in 1:6)
  for(j in 1:6)
    if(!(i==1 && j==1) && !(i==6 && j==6)) 
      P[i,j] = Pij(i-1,j-1,5)

states <- c("0","1","2","3","4","5")
transition <- P
mc <- new("markovchain",states = states,transitionMatrix = transition)
mc
rmarkovchain(mc,n=10)

simulacion <- function(x0){
  r = rep(0,1000)
  for(i in 1:1000){
    x = rmarkovchain(mc,n=100,t0 = x0)
    for(j in 1:100){
      if(x[j] == "5"){
        r[i] = 1
        break  
      } 
    }
  }
  mean(r)
}

r = rep(0,1000)
for(i in 1:1000){
  x = rmarkovchain(mc,n=100)
  for(j in 1:100){
    if(x[j] == "5"){
      r[i] = 1
      break  
    } 
  }
}
mean(r)

canonicForm(mc)

Q = matrix(c( 0.4096, 0.2048, 0.0512, 0.0064,
              0.2592, 0.3456, 0.2304, 0.0768,
              0.0768, 0.2304, 0.3456, 0.2592,
              0.0064, 0.0512, 0.2048, 0.4096),byrow=TRUE,nrow=4)

B = matrix(c( 0.32768, 0.00032, 
              0.07776, 0.01024, 
              0.01024, 0.07776, 
              0.00032, 0.32768 ),byrow=TRUE,nrow=4)

S = ginv(diag(4)-Q)
G = S %*% B
G

# Ejercicio 10

x = rpois(2*60, 100/60)
x[1] = 150 + x[i]
for(i in 2:120)
  x[i] = x[i-1] + x[i]
x
matplot(x, type='l', lty=1, col=3:5, ylim=c(150,400), ylab='Cantindad de Votantes', xlab='Tiempo(Minutos)')

r = 0
for(i in 0:199)
  r = r + dpois(i, lambda=200)
1-r

  
  P <- function(n,i,j){
    transitionProbability(mc^n,i,j)
  }
  
  F <- function(n,i,j,E){
      if(n==1) return (P(1,i,j))
      else {
        sum = 0
        for(b in 1:length(E)){
          if(E[b] != j) sum = sum + P(1,i,E[b]) * F(n-1,E[b],j,E)
        }
        return (sum)
      }
  }
  


























pi <- c(0,1,0,0,0,0)
states <- c("1","2","3")
Q <- matrix(c(0.565,0.085,0,
                       0,0.7875,0.0125,
                       0,0,0.95),byrow=TRUE,nrow=4)

S = ginv(diag(3)-transition)
S



# PRACTICA

# Ej 14

# A
states = c("0","1","2","3")
P = matrix(c(0.3,0.5,0.2,0,
             0.6,0,0.4,0,
             0,0.2,0.4,0.4,
             0,0,0.2,0.8),byrow=TRUE,nrow=4)
mc <- new("markovchain",states = states,transitionMatrix = P)
is.irreducible(mc)
period(mc)
steadyStates(mc)
mc^100

# B
states = c("0","1","2","3","4")
P = matrix(c(0,0,0.5,0,0.5,
             0,0,1,0,0,
             0,0,0,1,0,
             0,0,0.6,0,0.4,
             0.8,0.2,0,0,0),byrow=TRUE,nrow=5)
mc <- new("markovchain",states = states,transitionMatrix = P)
is.irreducible(mc)
period(mc)
steadyStates(mc)
mc^100
mc^101

# C
states = c("0","1","2","3")
P = matrix(c(0,0.2,0.8,0,
             0.2,0,0,0.8,
             0.2,0,0,0.8,
             0,0,1,0),byrow=TRUE,nrow=4)
mc <- new("markovchain",states = states,transitionMatrix = P)
is.irreducible(mc)
period(mc)
steadyStates(mc)
mc^100
mc^101

# Ej 15

states = c("1","2","3","S")
P = matrix(c(0.565   ,0.085  ,0       ,0.35,
             0       ,0.7875 ,0.0125   ,0.2,
             0       ,0      ,0.95    ,0.05,
             0       ,0      ,0       ,1),byrow=TRUE,nrow=4)
mc <- new("markovchain",states = states,transitionMatrix = P)
is.irreducible(mc)
period(mc)
steadyStates(mc)

canonicForm(mc)


Q = matrix(c( 0.565, 0.0850, 0.0000,
              0.000, 0.7875, 0.0125,
              0.000, 0.0000, 0.9500),byrow=TRUE,nrow=3)
R = ginv(diag(3) - Q)
R
sum(R[1,])
R[1,3] / R[3,3]


rexp(16, 200)

############################################################################
# Simulate a homogeneous Poisson process - uses the 'poisson' package in R
############################################################################

# Install poisson package
library("poisson")


####################################################
# Understand the following functions for simulation
####################################################
hpp.event.times
hpp.sim


# Recreate the code above - simulate event times
PP.sim <- function (rate, num.events, num.sims = 1, t0 = 0) 
{
  if (num.sims == 1) {
    x = t0 + cumsum(rexp(n = num.events, rate = rate))
    return(c(t0,x))
  }
  else {
    xtemp = matrix(rexp(n = num.events * num.sims, rate = rate), num.events)
    x = t0 + apply(xtemp, 2, cumsum)
    return(rbind(rep(t0, num.sims), x))
  }
}


# Run the function above, num.sims default is 1
rate = 100
x = PP.sim(rate, num.events = 200)
dif = x
x
sum(x)
for(i in 2:201)
  x[i] = x[i-1] + x[i]
x
plot(x)
plot(dif)
curve(exp)
hist(dif)
# Compare to hpp.sim
hpp.sim(rate, num.events=10, num.sims = 5)



sim_pp2 <- function(t, rate) {
  
  path <- matrix(0, nrow = 1, ncol = 2)
  
  jumps_number <- rpois(1, lambda = rate * t)
  jumps_time <- runif(n = jumps_number, min = 0, max = t) %>% sort()
  
  for(j in seq_along(jumps_time)) {
    jump <- matrix(c(jumps_time[j], path[nrow(path), 2],
                     jumps_time[j], path[nrow(path), 2]  + 1),
                   nrow = 2, ncol = 2, byrow = TRUE)
    path <- rbind(path, jump)
  }
  
  path <- rbind(path,
                c(t, path[nrow(path), 2]))
  
  list(path, jumps_time)
  
}


library("ggplot2")
library("magrittr")

set.seed(1)

path1 <- sim_pp1(1000, 1)
mean(diff(path1[[2]])); var(diff(path1[[2]]))
# [1] 1.029312
# [1] 0.9722406

path2 <- sim_pp2(2, 100)
path2
sum(path2)
mean(diff(path2[[2]])); var(diff(path2[[2]]))
# [1] 1.006302
# [1] 1.066079

data.frame(it = diff(path2[[2]])) %>%
  ggplot() +
  geom_histogram(aes(it, y = ..density..)) +
  stat_function(fun = dexp) +
  theme_bw() + 
  theme(text = element_text(size = 24))


x = rep(0,1000)
for(i in 1:1000){
  datos = rexp(200,100)
  if(sum(datos)<=2) x[i] = 1
}
mean(x)

datos = rexp(200,100)
sum(datos)
dif = datos
for(i in 2:200) datos[i] = datos[i-1] + datos[i]
plot(datos,151:350,xlab="Tiempo(Horas)",ylab = "Cantidad de Votantes")

hist(dif,freq = FALSE,xlab = "Intervalos de Tiempo(Horas)", ylab = "Densidad",main = "Distribución de Probabilidad del Tiempo entre Arrivos")
curve(dexp(x,100),add = TRUE)
?plot
?rexp
