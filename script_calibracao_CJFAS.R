########################################################################
## Description: R code developed for implementing a Polya Urn applied to
##              fisheries data.
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
### Commentary: The R code presents frequentist performance measures for 
##              finite population estimators obtained from two-stage (2E)
##              samples of fishermen extracted from the virtual population
##              constructed for this purpose.
##              The following are measures of performance that can be 
##              calculated with the code: 
##                 1) accuracy of the estimators 'mean' and 'median'; 
##                 2) range and coverage of the percentile and/or HDI 
##                    intervals;
##                 3) mean value of the coefficient of variation (CV) 
##                    and its robust version (RCV)
##                 4) robust version (RCV)
## 
### Code:
########################################################################

##@ Precisa carregar a biblioteca geradora das amostras de polya
library(polyapost)
library(HDInterval)

##@-------------------------------------------
##@ PARTE 1:
##@ Importacao da Populacao Completa e a definicao de termos globais para 
##@ as etapas posteriores.
##@------------------------------------------

##@  1. Importar a populacao em uma  lista denotada 'censo_lista' a partir 
##@  da qual serao extraidas as amostras para efetuar as inferencias 
##@  via posterior de Polya.
rm(list=ls())
load(file="dkft_CensoVirtual_perCapita.Rda")

##@ Carregar os parametros de Pop Finita
( ParPopFin <- censo_lista[[2]] )
censo_pc <- censo_lista[[1]]
names(censo_pc)

##@ Para cada uma das semanas as colunas a serem selecionadas de 'censo' estao 
##@ listadas nos quatro vetores a seguir
##@ Obs.: Vamos extrair apenas as colunas dos cluster, dos esforcos e das 
##@       capturas das 4 categorias de pescado (por semmana)
col_week1 <- c(1,7,11:14)
col_week2 <- c(1,8,15:18)
col_week3 <- c(1,9,19:22)
col_week4 <- c(1,10,23:26)

##@ O tamanho da populacao = nro de linhas do data.frame 'censo'
( N <- dim(censo_pc)[1]  ) 

##@---------
##@ Para o Caso de 2E:
##@ Inicia-se especificando os tamanhos amostrais de cada estagio
nc <- 2     # nro de clusters por semana
nf <- 10     # nro de unidades dentro de cada cluster
##@ Supondo que todos os cluster N[k] >= nf tem-se 
#@  por semana:
( nw <- nf*nc )
#@ e por mes:
( n <- nw*4 )

##@----------
##@ O número total de Clusters no data.frame 'censo' 
##@ No caso em pauta, 'nro_c = 12'
( nro_c <- length(table(censo_pc$cluster)) ) 

##@ O número total de unidades secundarias no data.frame 'censo' 
##@ No caso, 'N_c <- as.numeric(table(censo$cluster))'
##@ (Obs. O comprimento deste vetor deve ser igual a 'nro_c')
( N_c <- as.numeric(table(censo_pc$cluster)) )

##@ O numero de variaveis correspondentes ao numero de colunas
##@ a ser criado em uma amostra semanal de unidades'am_acs'
#(nvar <- length(col_week1) + 1 )
nvar <- 5  # (ff, u1,u2,u3, u4)
##@---------

##@-------
##@ Algumas FUNCOES PROPRIAS:
##@==================
# Bayesian summary function to be used with apply()
fun.sumario.plus <- function(x){
  c(quantile(x,c(0.025,0.5,0.975)),mean=mean(x),sd=sd(x),
    rcv=100*mad(x)/median(x),cv=100*sd(x)/mean(x))
}

### Incluir funcao para amostrar comunidades com probabilidade diferentes
##@ A function to produz de sample of clusters
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

##- teste
#f_cluster_sample(N_c,nf,nc)

weeks <- 4
am_cluster <- matrix(f_cluster_sample(N_c,nf,nc*weeks),ncol=nc)
line_nr <- 1:N   #1:dim(censo)[1]
am_pesc_in_cluster <- array(rep(NA,weeks*nc*nf),c(nc,nf,weeks))
for(w in 1:weeks){
  for(i in 1:nc){
    #    zz <-com_rgd[am_cluster_rgd[w,i]]
    zz <- am_cluster[w,i]
    #    line_sel <- line_nr[drgd$Comunidade==zz]
    line_sel <- line_nr[censo_pc$cluster==zz]
    # precisa excluir if() e mudar amostra para replace=F quando tiver
    # cluster_size > nf para toda comunidade
    #if(length(line_sel) == 1) am_pesc_in_cluster[i,,w] <- rep(line_sel, nf)
    #else
      #      am_pesc_in_cluster[i,,w] <- sort(sample(line_sel,size=nf, replace=T))
    am_pesc_in_cluster[i,,w] <- sort(sample(line_sel,size=nf, replace=F))
  }
}
#am_pesc <- am_pesc_in_cluster
##@

##@ FUNCAO PARA POLYA DA PRODUCAO 
##@ Segue Abaixo a funcao 'polya_prod_mes' para criar as posteriores de Polya 
##@ para as producoes mensais de todas as categorias de pescado e do esforco
##@ ATENCAO: Esta funcao vale somente para 4 Semanas Completas.
##@ N caso de 5 semanas usamos 'polya_prod_mes5' (Anexo Abaixo) 

##@ INICIO DA FUNCAO
polya_prod_mes <- function(N_mun, Nsim, YYYmXw1, YYYmXw2, YYYmXw3, YYYmXw4,aa,bb)
  ## 'aa' e 'bb' denotam  o inicio e fim das COLUNAS de YYYmXwS a ser usado
{
  NN <- N_mun 
  ##@ Fazer tudo com MATRIZES para ver se vai mais rapido que com dataframe
  am_u1 <- as.matrix(YYYmXw1[,c(aa:bb)])
  am_u2 <- as.matrix(YYYmXw2[,c(aa:bb)])
  am_u3 <- as.matrix(YYYmXw3[,c(aa:bb)])
  am_u4 <- as.matrix(YYYmXw4[,c(aa:bb)])
  #@ definir o tamanho da amostra na semana 'ns'
  #@ de cada municipio.
  ns <- c(dim(am_u1)[1],dim(am_u2)[1],dim(am_u3)[1],dim(am_u4)[1])
  #@@---
  #@ criar matriz para quardar as saidas das simulacoes
  #@ de cada semana
  po_u1 <-  matrix( rep(NA,Nsim*dim(am_u1)[2]),ncol=dim(am_u1)[2])
  po_u2 <- po_u1
  po_u3 <- po_u1
  po_u4 <- po_u1
  Qpolya <- matrix(rep(NA,NN*dim(am_u1)[2]),ncol=dim(am_u1)[2] )
  
  #@@ Inicio do loop de Polya
  for(j in 1:Nsim){
    #@ w1
    idx_polya <- polyap(1:ns[1],NN-ns[1])
    for(k in 1:NN){
      Qpolya[k,] <- am_u1[idx_polya[k],]
    }
    po_u1[j,] <- apply(Qpolya,2,sum)
    #@ w2
    idx_polya <- polyap(1:ns[2],NN-ns[2])
    for(k in 1:NN){
      Qpolya[k,] <- am_u2[idx_polya[k],]
    }
    po_u2[j,] <- apply(Qpolya,2,sum)
    #@ w3
    idx_polya <- polyap(1:ns[3],NN-ns[3])
    for(k in 1:NN){
      Qpolya[k,] <- am_u3[idx_polya[k],]
    }
    po_u3[j,] <- apply(Qpolya,2,sum)
    #@ w4
    idx_polya <- polyap(1:ns[4],NN-ns[4])
    for(k in 1:NN){
      Qpolya[k,] <- am_u4[idx_polya[k],]
    }
    po_u4[j,] <- apply(Qpolya,2,sum)
  }
  ##@ Fim do Loop de Polya
  po_YYY_mX <- data.frame(po_u1 + po_u2 + po_u3 + po_u4)
  names(po_YYY_mX) <- names(YYYmXw1)[aa:bb]
  return(po_YYY_mX)
} ##@ end of Function
##@  FIM DA FUNCAO
##@====================


##@-------------------------------------------
##@ PARTE 2:
##@ Fazer LOOP sobre extracao de amostras e calculo das estimativas de 
##@ populacao finita que queremos avaliar com este loop
##@------------------------------------------

##@ 
populacao <- 1:N

##@
weeks <- 4

#@ numero total de LOOPS é dado por 'nloop'
nloop <- 300 

#@ numero total simulacoes para a posterior de Polya DENTRO de cada LOOP é 
#@ dado por 'sim'
sim <- 1000

##@ Preparar matrizes para armazenar saidas do LOOP
n_in_fs <- 7   # nro de sumarios estatisticos definidos em 'fun.sumario.plus'

##@ a_calib(nloop, n_in_fs, nvar)
#a_calib <- array(rep(NA,nloop*n_in_fs*(nvar-1)),dim = c(nloop,n_in_fs,(nvar-1)))
a_calib <- array(rep(NA,nloop*n_in_fs*nvar),dim = c(nloop,n_in_fs,nvar))

##@ criar matrizes para guardar os resultados finais das comparacoes
##@ dos procedimento de Polya replicados com relacao aos parametros 
##@ populacionais reais
##@ Coluna 1: media                    Coluna 2: acuracia absoluta
##@ Coluna 3: amplitude icr percentil  Coluna 4: cobertura percentual
##@ Coluna 5: amplitude icr hdi        Coluna 6: cobertura percentua
#Fr_xxp <- matrix(rep(NA,nloop*6),ncol=6)
#Fr_nnp <- matrix(rep(NA,nloop*6),ncol=6)
FrEsf <- matrix(rep(NA,nloop*6),ncol=6)
#FrUG <- matrix(rep(NA,nloop*6),ncol=6)
FrU1 <- matrix(rep(NA,nloop*6),ncol=6)
FrU2 <- matrix(rep(NA,nloop*6),ncol=6)
FrU3 <- matrix(rep(NA,nloop*6),ncol=6)
FrU4 <- matrix(rep(NA,nloop*6),ncol=6)


###@===================================
###@ Aqui inicia o LOOP para Calibracao
 for(m in 1: nloop){

   ##@-------
   ##@ OPCAO 2: Amostra em DOIS ESTAgios (2E)
   ##@ Sorteio de clusters (estagio 1):

   #weeks <- 4
   am_cluster <- matrix(f_cluster_sample(N_c,nf,nc*weeks),ncol=nc)
   line_nr <- 1:N   #1:dim(censo)[1]
   am_pesc_in_cluster <- array(rep(NA,weeks*nc*nf),c(nc,nf,weeks))
   for(w in 1:weeks){
     for(i in 1:nc){
       #    zz <-com_rgd[am_cluster_rgd[w,i]]
       zz <- am_cluster[w,i]
       #    line_sel <- line_nr[drgd$Comunidade==zz]
       line_sel <- line_nr[censo_pc$cluster==zz]
       # precisa excluir if() e mudar amostra para replace=F quando tiver
       # cluster_size > nf para toda comunidade
       #if(length(line_sel) == 1) am_pesc_in_cluster[i,,w] <- rep(line_sel, nf)
       #else
       #      am_pesc_in_cluster[i,,w] <- sort(sample(line_sel,size=nf, replace=T))
       am_pesc_in_cluster[i,,w] <- sort(sample(line_sel,size=nf, replace=F))
     }
   }
   #am_pesc_in_cluster
   ##@
   
      ##@  Precisa agora criar as matrizes de produções para os pescadores
   ##@  amostrados a cada semana. Ou seja, as 4 matrizes semanais
   ##@  'am_w1' ate 'am_w4', 
   ##@ As matrizes semanais com os CONTEUDOS das entrevistas amostradas.
   ##@ (Obs. Em dois estagios deve inclui o cluster
   for(d in 1:nc){
     if(d == 1) 
       am_w1 <- censo_pc[am_pesc_in_cluster[1,,1],col_week1]
     else
       am_w1 <- rbind(am_w1,censo_pc[am_pesc_in_cluster[d,,1],col_week1] )
   }
   am_w1 <- cbind(id = as.numeric(row.names(am_w1)),am_w1)
   #@
   for(d in 1:nc){
     if(d == 1) 
       am_w2 <- censo_pc[am_pesc_in_cluster[1,,2],col_week2]
     else
       am_w2 <- rbind(am_w2,censo_pc[am_pesc_in_cluster[d,,2],col_week2] )
   }
   am_w2 <- cbind(id = as.numeric(row.names(am_w2)),am_w2)
   #@
   for(d in 1:nc){
     if(d == 1) 
       am_w3 <- censo_pc[am_pesc_in_cluster[1,,3],col_week3]
     else
       am_w3 <- rbind(am_w3,censo_pc[am_pesc_in_cluster[d,,3],col_week3] )
   }
   am_w3 <- cbind(id = as.numeric(row.names(am_w3)),am_w3)
   #@
   for(d in 1:nc){
     if(d == 1) 
       am_w4 <- censo_pc[am_pesc_in_cluster[1,,4],col_week4]
     else
       am_w4 <- rbind(am_w4,censo_pc[am_pesc_in_cluster[d,,4],col_week4] )
   }
   am_w4 <- cbind(id = as.numeric(row.names(am_w4)),am_w4)
   
   #@ criar lista para armazenar
   am_2e <- list(am_w1=am_w1, am_w2=am_w2, am_w3=am_w3, am_w4=am_w4)

#@ construir lista para armazenar
#am_acs <- list(am_w1=am_w1, am_w12=am_w2, am_w3=am_w3, am_w4=am_w4)

##@--------------------------------------------- 
##@ Preparar os dados para efetuar Polya:

##@ Esforco e Captura (4 especies) ja esta per capita 
#@ ... semana 1
dbq1 <- am_2e[[1]]
#dbq1[,5:9] <- dbq1[,5:9]/dbq1[,4]
#@ ... semana 2
dbq2 <- am_2e[[2]]
#dbq2[,5:9] <- dbq2[,5:9]/dbq2[,4]
#@ ... semana 3
dbq3 <- am_2e[[3]]
#dbq3[,5:9] <- dbq3[,5:9]/dbq3[,4]
#@ ... semana 4
dbq4 <- am_2e[[4]]
#dbq4[,5:9] <- dbq4[,5:9]/dbq4[,4]


##@ selecionar os clusters para compor as amostras semanais
#( clu1 <- as.numeric(names(table(dbq1$cluster))) )
#( clu2 <- as.numeric(names(table(dbq2$cluster))) )
#( clu3 <- as.numeric(names(table(dbq3$cluster))) )
#( clu4 <- as.numeric(names(table(dbq4$cluster))) )


###@Fazer o NOVO loop de Polya para as 4 semanas:
##@ Usar a funcao 'polya_prod_mes'

dt_polya <- polya_prod_mes(N, sim, dbq1, dbq2, dbq3, dbq4, 3, 7)
names(dt_polya) <- c("ff","uu1","uu2","uu3","uu4")
apply(dt_polya,2,fun.sumario.plus)
ParPopFin[-6]
#apply(dt_polya,2,hdi)


##@ Aqui sao calculados as estimativas posteriores de Polya de interesse e que
##@ estao apresentados na funcao 'fun.sumario.plus'
a_calib[m,,1:nvar] <- apply(dt_polya,2,fun.sumario.plus)
## sumario para a captura total
#a_calib[m,,nvar-1] <- fun.sumario(cap_tot_polya)

##@------------------------
# 2e. Fazer os sumarios sobre as K replicas de Polya em relacao
#  aos parametros populacionais "reais"
# (i) Esforco Total (Esf)
Tf.est <- dt_polya[,1]
EsfTotal <- ParPopFin[1]
FrEsf[m,1] <- 1 - abs(median(Tf.est) - EsfTotal)/EsfTotal
FrEsf[m,2] <- 1 - abs(mean(Tf.est) - EsfTotal)/EsfTotal
qq <- as.numeric(quantile(Tf.est,prob=c(0.025,0.975)))
FrEsf[m,3] <- qq[2]-qq[1]
FrEsf[m,4] <- (qq[1] <= EsfTotal)*(EsfTotal <= qq[2])
qq <- as.numeric(hdi(Tf.est, 0.95))
FrEsf[m,5] <- qq[2]-qq[1]
FrEsf[m,6] <- (qq[1] <= EsfTotal)*(EsfTotal <= qq[2])
#--
# (ii) Captura Total (UG)
#Ty.est <- cap_tot_polya
#UGTotal <- ParPopFin[6]
#FrUG[m,1] <- 1 - abs(median(Ty.est) - UGTotal)/UGTotal
#FrUG[m,2] <- 1 - abs(mean(Ty.est) - UGTotal)/UGTotal
#qq <- as.numeric(quantile(Ty.est,prob=c(0.025,0.975)))
#FrUG[m,3] <- qq[2]-qq[1]
#FrUG[m,4] <- (qq[1] <= UGTotal)*(UGTotal <= qq[2])
#qq <- as.numeric(hdi(Ty.est, 0.95))
#FrUG[m,5] <- qq[2]-qq[1]
#FrUG[m,6] <- (qq[1] <= UGTotal)*(UGTotal <= qq[2])
#--
# (iii) Captura 1 (U1)
Ty1.est <- dt_polya[,2]
U1Total <- ParPopFin[2]
FrU1[m,1] <- 1 - abs(median(Ty1.est) - U1Total)/U1Total
FrU1[m,2] <- 1 - abs(mean(Ty1.est) - U1Total)/U1Total
qq <- as.numeric(quantile(Ty1.est,prob=c(0.025,0.975)))
FrU1[m,3] <- qq[2]-qq[1]
FrU1[m,4] <- (qq[1] <= U1Total)*(U1Total <= qq[2])
qq <- as.numeric(hdi(Ty1.est, 0.95))
FrU1[m,5] <- qq[2]-qq[1]
FrU1[m,6] <- (qq[1] <= U1Total)*(U1Total <= qq[2])
#--
# (iv) Captura 2 (U2)
Ty2.est <- dt_polya[,3]
U2Total <- ParPopFin[3]
FrU2[m,1] <- 1 - abs(median(Ty2.est) - U2Total)/U2Total
FrU2[m,2] <- 1 - abs(mean(Ty2.est) - U2Total)/U2Total
qq <- as.numeric(quantile(Ty2.est,prob=c(0.025,0.975)))
FrU2[m,3] <- qq[2]-qq[1]
FrU2[m,4] <- (qq[1] <= U2Total)*(U2Total <= qq[2])
qq <- as.numeric(hdi(Ty2.est, 0.95))
FrU2[m,5] <- qq[2]-qq[1]
FrU2[m,6] <- (qq[1] <= U2Total)*(U2Total <= qq[2])
#--
# (v) Captura 3 (U3)
Ty3.est <- dt_polya[,4]
U3Total <- ParPopFin[4]
FrU3[m,1] <- 1 - abs(median(Ty3.est) - U3Total)/U3Total
FrU3[m,2] <- 1 - abs(mean(Ty3.est) - U3Total)/U3Total
qq <- as.numeric(quantile(Ty3.est,prob=c(0.025,0.975)))
FrU3[m,3] <- qq[2]-qq[1]
FrU3[m,4] <- (qq[1] <= U3Total)*(U3Total <= qq[2])
qq <- as.numeric(hdi(Ty3.est, 0.95))
FrU3[m,5] <- qq[2]-qq[1]
FrU3[m,6] <- (qq[1] <= U3Total)*(U3Total <= qq[2])
#--
# (vi) Captura 4 (U4)
Ty4.est <- dt_polya[,5]
U4Total <- ParPopFin[5]
FrU4[m,1] <- 1 - abs(median(Ty4.est) - U4Total)/U4Total
FrU4[m,2] <- 1 - abs(mean(Ty4.est) - U4Total)/U4Total
qq <- as.numeric(quantile(Ty4.est,prob=c(0.025,0.975)))
FrU4[m,3] <- qq[2]-qq[1]
FrU4[m,4] <- (qq[1] <= U4Total)*(U4Total <= qq[2])
qq <- as.numeric(hdi(Ty4.est, 0.95))
FrU4[m,5] <- qq[2]-qq[1]
FrU4[m,6] <- (qq[1] <= U4Total)*(U4Total <= qq[2])
#--
} ###@ Aqui termina o LOOP para Calibracao com 'nloop'


###@-------------------------------------------------
##@ NOTA 1: 'a_calib' eh array com dim=c(nloop,n_in_fs,nvar)
##@ As colunas de n_in_fs sao: 2.5%, 50%, 97.5%, mean, sd, rcv, cv.
##@
##@ As colunas de nvar sao: ff, uu1, uu2, uu3, uu4
##@ --
##@ NOTA 2: A matriz 'FrEsf' tem dim=c(nloop, 6).
##@ As colunas contem: (1) Acuracia percentual da mediana posterior
##@                    (2) Acuracia percentual da media posterior
##@                    (3) Amplitude do ICr nominal de 95%
##@                    (4) Cobertura do ICr
##@                    (5) Amplitude do HDI nominal de 95%
##@                    (6) Cobertura do HDI
##@ A mesma estrutura vale para 'FrU1', 'FrU2', 'FrU3', 'FrU4'.
##@=================================================

x1 <- apply(FrEsf,2,mean)
outEsf <- data.frame( true = EsfTotal,
                      acc.median = x1[1]*100,
                      acc.mean = x1[2]*100,
                      length.icr = x1[3],
                      cov.icr = x1[4]*100,
                      length.hdi = x1[5],
                      cov.hdi = x1[6]*100,
                      row.names="FrEsf")
round(outEsf,2)
# Captura Total (UG)
#x1 <- apply(FrUG,2,mean)
#outUG <- data.frame(  true = UGTotal,
#                      acc.median = x1[1]*100,
#                      acc.mean = x1[2]*100,
#                      length.icr = x1[3],
#                      cov.icr = x1[4]*100,
#                      length.hdi = x1[5],
#                      cov.hdi = x1[6]*100,
#                      row.names="FrUG")
#round(outUG,2)
# Captura 1 (U1)
x1 <- apply(FrU1,2,mean)
outU1 <- data.frame( true = U1Total,
                     acc.median = x1[1]*100,
                     acc.mean = x1[2]*100,
                     length.icr = x1[3],
                     cov.icr = x1[4]*100,
                     length.hdi = x1[5],
                     cov.hdi = x1[6]*100,
                     row.names="FrU1")
#round(outU1,2)
# Captura 2 (U2)
x1 <- apply(FrU2,2,mean)
outU2 <- data.frame( true = U2Total,
                     acc.median = x1[1]*100,
                     acc.mean = x1[2]*100,
                     length.icr = x1[3],
                     cov.icr = x1[4]*100,
                     length.hdi = x1[5],
                     cov.hdi = x1[6]*100,
                     row.names="FrU2")
#round(outU2,2)
# Captura 3 (U3)
x1 <- apply(FrU3,2,mean)
outU3 <- data.frame( true = U3Total,
                     acc.median = x1[1]*100,
                     acc.mean = x1[2]*100,
                     length.icr = x1[3],
                     cov.icr = x1[4]*100,
                     length.hdi = x1[5],
                     cov.hdi = x1[6]*100,
                     row.names="FrU3")
#round(outU3,2)
# Captura 4 (U4)
x1 <- apply(FrU4,2,mean)
outU4 <- data.frame( true = U4Total,
                     acc.median = x1[1]*100,
                     acc.mean = x1[2]*100,
                     length.icr = x1[3],
                     cov.icr = x1[4]*100,
                     length.hdi = x1[5],
                     cov.hdi = x1[6]*100,
                     row.names="FrU4")
#round(outU4,2)

# Criar uma tabela-Resumo única  (Pode ser uma LISTA)
#desempenho <- rbind(outEsf,outUG,outU1,outU2,outU3,outU4)
desempenho <- rbind(outEsf,outU1,outU2,outU3,outU4)

save(desempenho, file = "desempenho_2e.Rda")
load(file = "desempenho_2e.Rda")
#load(file="desempenho_acs.Rda")
round(desempenho,1)
Tab4 <- round(desempenho[,-1],1)
row.names(Tab4) <- c("Eff","Sp1", "Sp2", "Sp3", "Sp4")
names(Tab4)  <- c("ARA(g)", "ARA(t)", "CrI.95", "covCrI", "HDI.95", "covHDI" )
Tab4
save(Tab4, file="Tab4.Rda")


##@==============================================
##@ Estimativas posteriores de Parametros de Pop Finitas para as 'nloop'
##@ replicações de todo o procedimento desde a seleção da amostra ate
##@ o calculo dos sumarios posteriores
##@
#@ para fracao de proprietarios 'xx_p'
#cal_xxp <- data.frame(a_calib[,,1])
#colnames(cal_xxp) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para nro de pessoas pescando juntas 'nnp'
#cal_nnp <- data.frame(a_calib[,,2])
#colnames(cal_nnp) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ oara esforco 'ff'
#cal_ff <- data.frame(a_calib[,,3])
cal_ff <- data.frame(a_calib[,,1])
colnames(cal_ff) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para captura da categoria 1 em 'uu1'
#cal_u1 <- data.frame(a_calib[,,4])
cal_u1 <- data.frame(a_calib[,,2])
colnames(cal_u1) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para captura da categoria 2 em 'uu2'
#cal_u2 <- data.frame(a_calib[,,5])
cal_u2 <- data.frame(a_calib[,,3])
colnames(cal_u2) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para captura da categoria 3 em 'uu3'
#cal_u3 <- data.frame(a_calib[,,6])
cal_u3 <- data.frame(a_calib[,,4])
colnames(cal_u3) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para captura da categoria 4 em 'uu4'
#cal_u4 <- data.frame(a_calib[,,7])
cal_u4 <- data.frame(a_calib[,,5])
colnames(cal_u4) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

#@ para captura total em 'utotal'
#cal_ug <- data.frame(a_calib[,,8])
#colnames(cal_ug) <- c("Q025", "median", "Q975", "mean", "sd", "rcv", "cv")

##@ Juntar tudo em uma lista para analises complementares como nos exemplos
##@ destacados abaixo

#l_calib <- list(xxp = cal_xxp, nnp = cal_nnp, eff = cal_ff, u1 = cal_u1,
#                u2 = cal_u2, u3 = cal_u3, u4 = cal_u4, ug = cal_ug,
#                true  = ParPopFin)
l_calib <- list(eff = cal_ff, u1 = cal_u1,
                u2 = cal_u2, u3 = cal_u3, u4 = cal_u4,
                true  = ParPopFin)

#save(l_calib, file="lista_calib_2e.Rda")

##@------------------------------------------------------------
##@ Possiveis analises complementares:
##@ Se comeca por aquim precisa carrega
##@ os arquivos pertinentes.
load(file="lista_calib_2e.Rda")
ParPopFin <- l_calib$true

##@ Minha funcao de g histigrama posterior bayesiano##--- Function ---------------------------------
bayeshist_md <- function(x,col.hist="light grey",main="main text",xlab="text x",
                      col.ci = "black",cred=0.95,nclass=20,prob=TRUE, median=TRUE,
                      border = "white")
  # bayesian one-dimensional posterior distrinution (histogram)
  # with credibility intervals (default cred=0.95)
{
  hist(x,nclass=nclass,col=col.hist,main=main,xlab=xlab,prob=prob, 
       border=border)
  aa <- (1-cred)/2
  bb <- quantile(x, prob=c(aa,1-aa))
  segments(bb[1],0,bb[2],0,lwd=5,col=col.ci)
  if(median == TRUE)
  points(median(x),0, pch=21, bg="white",cex=2.0, col=col.ci)
  else
  points(mean(x),0, pch=21, bg="white",cex=2.0, col=col.ci)
  
}
##-- End Function



##@ Graficos:
#@ Fazer um histograma das medianas posteriores de esforco
cal_ff <- l_calib$eff
#hist(cal_ff[,2], main="Median Effort",xlab="ff",nclass=12)
bayeshist_md(cal_ff[,2], main="Median Effort",xlab="days", median=FALSE)
abline(v = ParPopFin[1],col="blue", lwd=2)

par(mfrow=c(2,2))
#1
cal_u1 <- l_calib$u1
#hist(cal_u1[,2], main="U1",xlab="kg",nclass=12)
bayeshist_md(cal_u1[,2]/1000, main="U1",xlab="ton", median=FALSE)
abline(v = ParPopFin[2]/1000,col="blue", lwd=2)
#2
cal_u2 <- l_calib$u2
#hist(cal_u2[,2], main="U2",xlab="kg",nclass=12)
bayeshist_md(cal_u2[,2]/1000, main="U2",xlab="ton", median=FALSE)
abline(v = ParPopFin[3]/1000,col="blue", lwd=2)
#3
cal_u3 <- l_calib$u3
#hist(cal_u3[,2], main="U3",xlab="kg",nclass=12)
bayeshist_md(cal_u3[,2]/1000, main="U3",xlab="ton", median=FALSE)
abline(v = ParPopFin[4]/1000,col="blue", lwd=2)
#4
cal_u4 <- l_calib$u4
#hist(cal_u4[,2], main="U4",xlab="kg",nclass=12)
bayeshist_md(cal_u4[,2]/1000, main="U2",xlab="ton", median=FALSE)
abline(v = ParPopFin[5]/1000,col="blue", lwd=2)
par(mfrow=c(1,1))


##@ Estatisticas {Confuso: Melhorar esses sumarios !!!}
# calibration summary function to be used with apply()
fun.calibration <- function(x){
  c(median=median(x),mad=mad(x),mean=mean(x),sd=sd(x))
}
apply(cal_ff,2,fun.calibration)

##@ ou, somente para a mediana da captura total ug (medido em toneladas)
#fun.calibration(cal_ug[,2]/1000)
#@ e o valor real para comparacao
#ParPopFin[6]/1000
ParPopFin[1]
#@
#hist(cal_ug[,2], main="Median Total Catch",xlab="ug",nclass=12)
#abline(v = ParPopFin[6],col="blue", lwd=2)

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


