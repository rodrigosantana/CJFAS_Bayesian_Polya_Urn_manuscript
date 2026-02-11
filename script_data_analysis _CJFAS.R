########################################################################
## Description: Analyses to go with the CJFAS article on Bayesian
##              Finite Population Sampling in Artisanal Fishery
## 
## Maintainer: Datenkraft Data Science / University of Vale do Itajaí
## Author: Paul Gerhard Kinas and Rodrigo Sant'Ana
## Version: 0.0.1
## 
## URL: 
## Doc URL: 
## 
## Database info: 
## 
### Commentary: 
## 
### Code:
########################################################################

rm(list=ls())
#@
#@ Data were generated with 
#@ 'dkft_PopulacaoSimulada_DoisEstagios_v2.R'
load(file="dkft_CensoVirtual_perCapita.Rda")
ls()
str(censo_lista)
names(censo_lista)
#@ Data from DGD master dissertation
#load(file="PopVir_list.Rda")
#ls()
#str(PopVir_list)
#names(PopVir_list)
#@-----
#@
dd <- censo_lista$censo_p_cap
( ppf <- censo_lista$ParPopFin )
dim(dd)
( N <- dim(dd)[1] )

##@@@@@@@@@@@
##@ some Finite Population (FP) summary statistics
head(dd)
( ctab <- table(dd$cluster) )

##@@@----------------
##@ Figure 'PopPairs'
par(mfrow=c(2,3))
##@ all 2x2 plots
plot(dd$u1.1,dd$u2.1, xlab= "Sp1.w1",ylab="Sp2.w1", main="(a)")
plot(dd$u1.1,dd$u3.1, xlab= "Sp1.w1",ylab="Sp3.w1", main="(b)")
plot(dd$u1.1,dd$u4.1, xlab= "Sp1.w1",ylab="Sp4.w1", main="(c)")
#
plot(dd$u2.1,dd$u3.1, xlab= "Sp2.w1",ylab="Sp3.w1", main="(d)")
plot(dd$u2.1,dd$u4.1, xlab= "Sp2.w1",ylab="Sp4.w1", main="(e)")
plot(dd$u3.1,dd$u4.1, xlab= "Sp3.w1",ylab="Sp4.w1", main="(f)")
##@@@@@
par(mfrow=c(1,1))
##@@@--------------


##@--
#@ Total for Category 1 by Cluster and Week
##@ a. start stacking all weeks
dst1 <- data.frame( weeks = factor(rep(c(1,2,3,4),c(N,N,N,N))),
                    clu = factor(rep(dd$cluster,4)),
                    ff = c(dd$f.1,dd$f.2,dd$f.3,dd$f.4),
                    u1 = c(dd$u1.1,dd$u1.2,dd$u1.3,dd$u1.4)) 
##@ b. define factor combinations and calculate 'sum' by cell
fatores <- list(dst1$weeks, dst1$clu)
round(tapply(dst1$ff, fatores, sum ), 1)
round(tapply(dst1$ff, fatores, mean ),2)
round(tapply(dst1$ff, fatores, mad ),2)
##@
##@ Could do the same for Categories 2 to 4

##@ For Category 2
dst2 <- data.frame( weeks = factor(rep(c(1,2,3,4),c(N,N,N,N))),
                    clu = factor(rep(dd$cluster,4)),
                    ff = c(dd$f.1,dd$f.2,dd$f.3,dd$f.4),
                    u2 = c(dd$u2.1,dd$u2.2,dd$u2.3,dd$u2.4)) 
##@ b. define factor combinations and calculate 'sum' by cell
fatores <- list(dst2$weeks, dst2$clu)
round(tapply(dst2$ff, fatores, sum ),1)
round(tapply(dst2$ff, fatores, mean ),2)
round(tapply(dst2$ff, fatores, mad ),2)

##@ For Category 3
dst3 <- data.frame( weeks = factor(rep(c(1,2,3,4),c(N,N,N,N))),
                    clu = factor(rep(dd$cluster,4)),
                    ff = c(dd$f.1,dd$f.2,dd$f.3,dd$f.4),
                    u3 = c(dd$u3.1,dd$u3.2,dd$u3.3,dd$u3.4)) 
##@ b. define factor combinations and calculate 'sum' by cell
fatores <- list(dst3$weeks, dst3$clu)
round(tapply(dst3$ff, fatores, sum ),1)
round(tapply(dst3$ff, fatores, mean ),2)
round(tapply(dst3$ff, fatores, mad ),2)

##@ For Category 4
dst4 <- data.frame( weeks = factor(rep(c(1,2,3,4),c(N,N,N,N))),
                    clu = factor(rep(dd$cluster,4)),
                    ff = c(dd$f.1,dd$f.2,dd$f.3,dd$f.4),
                    u4 = c(dd$u4.1,dd$u4.2,dd$u4.3,dd$u4.4)) 
##@ b. define factor combinations and calculate 'sum' by cell
fatores <- list(dst4$weeks, dst4$clu)
round(tapply(dst4$ff, fatores, sum ),1)
round(tapply(dst4$ff, fatores, mean ),2)
round(tapply(dst4$ff, fatores, mad ),2)

###@@
##@ Totals per week (all clusters) and by Category 
xx1 <- rowSums( tapply(dst1$u1, fatores, sum ) )
#round( rowSums(xx1),2)
xx2 <- rowSums( tapply(dst2$u2, fatores, sum ) )
xx3 <- rowSums( tapply(dst3$u3, fatores, sum ) )
xx4 <- rowSums( tapply(dst4$u4, fatores, sum ) )

totals <- data.frame(rbind(xx1,xx2,xx3,xx4))
totals
names(totals) <- c("Sp1", "Sp2","Sp3","Sp4")
rownames(totals) <- c("w1","w2","w3","w4")
##@ ready for table in Paper:
Tab1Sp <- round(totals,1)
save(Tab1Sp, file="Tab1Sp.Rda")
##@
ddm <- as.matrix(dd)
ddTab <- data.frame(rbind( ddm[,c(1,7,11:14)],
                           ddm[,c(1,8,15:18)],
                           ddm[,c(1,9,19:22)],
                           ddm[,c(1,10,23:26)]))
ddTab <- cbind(week = rep(c("w1","w2","w3","w4"),c(350,350,350,350)), ddTab)
names(ddTab) <- c("week", "Com", "Eff","Sp1", "Sp21", "Sp3", "Sp4" )
ddTab$week <- as.factor(ddTab$week)
ddTab$Com <- as.factor(ddTab$Com)
levels(ddTab$Com) <- c( "a","b","c","d","e","f","g","h","i","j","k","l")
str(ddTab)
Tab1b <- round(tapply(ddTab[,3:7],ddTab[,c(2,1)], sum ),1)
ComSize <- as.numeric(table(ddTab$Com)/4)
##@ Tabela com as Somas Sp1 a Sp4 por Com e por semana
( Tab1sum <- cbind(Tab1b,ComSize) )
## Tabela com Medias totais por Com e por semana
( Tab1mean <- data.frame(cbind(round(Tab1b[,1:4]/ComSize,1), Size = ComSize)) )
save(Tab1mean, file="Tab1mean.Rda")


#save(dd, file="dd.Rda")
#load(dd, file="dd.Rda")


##@=================================================================
##@---
##@ Obs. For N = 350,the sampling procedure in use today in ELP is 
##@      equivalent to 20 fishers per week totaling 80 per months
##@
##@    Starting july 2025 we are changing to something closer to 
##@      24 fisher per week totaling 96 per months 
##@       
##@---

##@ A function to produce de sample of clusters with probability
##@ proportional to the number of combinations of 'm' fishers that
##@ can be sampled within a cluster
f_cluster_sample <- function(C,m,mc){
  ##@------------------------
  ##@ function 'f_cluster_sample'
  # C: a vector with cluster sizes
  # m: a fixed sample size within clusters
  # mc: number of clusters to be sampled
  ##---------------------------
  a <- lchoose(C,m)
  aM <- max(a)
  b <- exp(a - aM)
  A <- aM + log(sum(b))
  p_cl <- exp(a - A)
  am <-sample(sample(1:length(C),mc, replace=F, prob = p_cl),mc)
  am
}  # end function 'f_cluster_sample'


##@--
weeks <- 4
##@ using 20 per week
nc <- 2       # nr of clusters selected each week
nf <- 10      # nr de fishers selected within each cluster
##@-
##@ using 24 per week
#nc <- 2
#nf <- 12  ## can´t use since the smallest cluster has 10 fishers
#@ OR
#nc <- 3
#nf <- 8
##@--

#ctab <- table(dd$cluster)
cluster_index <- 1:length(ctab)
cluster_size <- as.numeric(ctab)
cluster_names <- names(ctab)

##@ Sampling
##@ select the sample o clusters
set.seed(18062025)
am_cluster <- matrix(f_cluster_sample(cluster_size,nf,nc*weeks),ncol=nc)

##@ select the sample of fishers WITHIN selected clusters and store
##@ it in 'am_pesc_ic' (ic por "in cluster")
#@
line_nr <- 1:dim(dd)[1]
am_pesc_ic <- array(rep(NA,weeks*nc*nf),c(nc,nf,weeks))

for(w in 1:weeks){
  for(i in 1:nc){
    zz <- cluster_names[am_cluster[w,i]]
    line_sel <- line_nr[dd$cluster== zz]
    am_pesc_ic[i,,w] <- sort(sample(line_sel,size=nf, replace=F))
  }
}
am_pesc_ic  #here is the sample of line by cluster and week
##@


##@--------------------------------------------------------------
##@ The next step is to extract de information obtained each week
##@ with sampling

##@ decompose 'dd' into four weeks
#names(dd)
Y_week1 <- dd[, c("f.1","u1.1", "u2.1", "u3.1", "u4.1")]
Y_week2 <- dd[, c("f.2","u1.2", "u2.2", "u3.2", "u4.2")]
Y_week3 <- dd[, c("f.3","u1.3", "u2.3", "u3.3", "u4.3")]
Y_week4 <- dd[, c("f.4","u1.4", "u2.4", "u3.4", "u4.4")]
##@ population size
( N <- nrow(dd) )
#@ number of variables in the sample (for estimation of Finite Pop Parameters)
( kY <- ncol(Y_week1) )    

##@--------------
# Tamanho da amostra MENSAL (1/4 para cada semana)
n.mes <- 80
n <- 20
##@============================
#@ Selecionar Amostra da Populacao para as varias semanas e armazenar
#@ sob forma de matriz
#@ select lines per week
lw1 <- as.vector(t(am_pesc_ic[,,1]))
lw2 <- as.vector(t(am_pesc_ic[,,2]))
lw3 <- as.vector(t(am_pesc_ic[,,3]))
lw4 <- as.vector(t(am_pesc_ic[,,4]))
#@
Ys1 <- as.matrix(Y_week1[lw1,])
Ys2 <- as.matrix(Y_week2[lw2,])
Ys3 <- as.matrix(Y_week3[lw3,])
Ys4 <- as.matrix(Y_week4[lw4,])


##@---
##@ Estimate Population Totals with POLYA posterior
##@ Call library first 
#@ incluir pacote para usar a funcao de amosttrador de Polya
library(polyapost)
#@ And create TWO useful summary functions
fun.sumario.plus <- function(x){
  c(quantile(x,c(0.025,0.5,0.975)),mean=mean(x),sd=sd(x),mad.sd=mad(x),mad.cv=mad(x)/median(x))
}

##--- Function ---------------------------------
bayeshist <- function(x,col.hist="light grey",main="main text",xlab="text x",
                      col.ci = "black",cred=0.95,nclass=20,prob=TRUE,
                      border = "white")
  # bayesian one-dimensional posterior distrinution (histogram)
  # with credibility intervals (default cred=0.95)
{
  hist(x,nclass=nclass,col=col.hist,main=main,xlab=xlab,prob=prob, 
       border=border)
  aa <- (1-cred)/2
  bb <- quantile(x, prob=c(aa,1-aa))
  segments(bb[1],0,bb[2],0,lwd=5,col=col.ci)
  points(mean(x),0, pch=21, bg="white",cex=2.0, col=col.ci)
}
##-- End Function

##@--
##@ Set Seed for Polya
set.seed(18082025)
# Nr of simulations
M <- 1500

#@@
# Simulação Polya: Criar matrizes para as saidas 
Q_polya1 <- matrix(NA, nrow = M, ncol = kY)
Q_polya2 <- matrix(NA, nrow = M, ncol = kY)
Q_polya3 <- matrix(NA, nrow = M, ncol = kY)
Q_polya4 <- matrix(NA, nrow = M, ncol = kY)
Y_sim <- matrix(NA, nrow = N, ncol = kY)
#@
for(i in 1:M){
  # semana 1
  idx_polya <- polyap(1:n,N-n)
  Y_sim <- Ys1[idx_polya,]
  Q_polya1[i, ] <- colSums(Y_sim)
  # semana 2
  idx_polya <- polyap(1:n,N-n)
  Y_sim <- Ys2[idx_polya,]
  Q_polya2[i,] <- colSums(Y_sim)
  # semana 3
  idx_polya <- polyap(1:n,N-n)
  Y_sim <- Ys3[idx_polya,]
  Q_polya3[i,]  <- colSums(Y_sim)
  # semana 4
  idx_polya <- polyap(1:n,N-n)
  Y_sim <- Ys4[idx_polya,]
  Q_polya4[i,] <- colSums(Y_sim)
}  # fim do Loop de Polya

# Resultados em forma de data.frame
QP1 <- as.data.frame(Q_polya1)
names(QP1) <- names(Y_week1)
QP2 <- as.data.frame(Q_polya2)
names(QP2) <- names(Y_week2)
QP3 <- as.data.frame(Q_polya3)
names(QP3) <- names(Y_week3)
QP4 <- as.data.frame(Q_polya4)
names(QP4) <- names(Y_week4)


##@========================
##@ TOTAIS MENSAIS (4 semanas somadas)
##@---
##@ Os parametros REAIS

( ppf <- censo_lista$ParPopFin )
especies <- c("Sp1","Sp2","Sp3","Sp4")

##@ Estimativas por Posterior de Polya
# Juntar todos SOMAR as 4 semanas
QP <- QP1 + QP2 + QP3 + QP4
names(QP) <- c("Eff",especies)

# Verificar a estrutura do data frame completo
#str(QP)
head(QP)

##@ sumarios estatísticos da posterior de polya
Tab3 <- data.frame(rbind(round(apply(QP,2,fun.sumario.plus),2),
              round(ppf[1:5],2)))
row.names(Tab3)[c(2,8)] <- c("median","FPP")
Tab3
#save(Tab3, file="Tab3.Rda")

## Comparar as medianas posteriores com os parametros reais:
round(apply(QP,2, median),2)
#@ parametros reais
ppf[1:5]
names(ppf)

##@---
#@ Histogramas dos totais MENSAIS:

##@ Esforco:
bayeshist(QP[,1], main=paste("Total -", names(QP)[1]),
     xlab="Fishing Days", nclass = 20)
abline(v=ppf[1], 
       col="red", lwd=2)

#@ capturas por especie
par(mfrow=c(2,2))
for(i in 2:5){
  bayeshist(QP[,i]/1000, main=paste("Total -", names(QP)[i]),
       xlab="Tonnes", nclass=20)
  abline(v=ppf[i]/1000, 
         col="red", lwd=2)
}
par(mfrow=c(1,1))

#@ Captura total (todas especies juntas)
Y_tot <- apply(QP[,-1],1,sum)
round(fun.sumario.plus(Y_tot/1000),2)
##@ Valor Real
ppf[6]/1000
#@
bayeshist(Y_tot/1000, main= "Overall Total",
     xlab="Tonnes", nclass = 20)
abline(v=ppf[6]/1000, 
       col="red", lwd=2)




##@-----------------------------------------------
##@ Incluir as estimativas de Estatisticas menos usuais
##@ (embora úteis) e de dificil estimacao por metodos 
##@ fSample-based.
##@
##@  (a) CPUE para "Sp2"
##@
##@  (b) Captura Total excedendo algum Threshold 'G = 62 ton'

##@ Resolvendo (a) CPUE
cpue.Sp2 <- QP[,3]/QP[,1] 
round(fun.sumario.plus(cpue.Sp2),2)
C.Sp2 <- ppf[3]/ppf[1] 

##@ Resolvendo (b) Threshold
##@ Here "True PPF" is 52.5 ton
G <- 62000

x.G <- sum(Y_tot > G)
n.G <- length(Y_tot)
## Using a Non-Informative Beta Prior
## Posterior is Beta(x.G +1, n>G - x.G +1)
curve(dbeta(x,x.G+1,n.G-x.G+1),0.05,.15,lwd=1.6)
pr_sim <- rbeta(M,x.G+1,n.G-x.G+1)
qbeta(c(.025,0.5,0.975),x.G+1,n.G-x.G+1)

##@ Figura
par(mfrow=c(1,2))
bayeshist(cpue.Sp2, main="CPUE Sp2",
          xlab="tons/day", nclass=20)
abline(v=C.Sp2, col="red", lwd=2)
##
bayeshist(pr_sim, main="Pr(TCatch > G)",
          xlab="p", nclass=20)

curve(dbeta(x,x.G+1,n.G-x.G+1),0.05,0.15,lwd=1.6,add=T)
par(mfrow=c(1,1))


##@==============================================
##@ Fim do Arquivo
##@==============================================


########################################################################
## 
##                  Creative Commons License 4.0
##                       (CC BY-NC-SA 4.0)
## 
##  This is a humam-readable summary of (and not a substitute for) the
##  license (https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode)
## 
##  You are free to:
## 
##  Share - copy and redistribute the material in any medium or format.
## 
##  The licensor cannot revoke these freedoms as long as you follow the
##  license terms.
## 
##  Under the following terms:
## 
##  Attribution - You must give appropriate credit, provide a link to
##  license, and indicate if changes were made. You may do so in any
##  reasonable manner, but not in any way that suggests the licensor
##  endorses you or your use.
## 
##  NonCommercial - You may not use the material for commercial
##  purposes.
## 
##  ShareAlike - If you remix, transform, or build upon the material,
##  you must distributive your contributions under the same license
##  as the  original.
## 
##  No additional restrictions — You may not apply legal terms or
##  technological measures that legally restrict others from doing
##  anything the license permits.
## 
########################################################################


