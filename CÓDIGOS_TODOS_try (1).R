rm(list=ls(all=TRUE))
############GLARMA MA1#####################



# Vetores para armazenar as estimativas


beta0.glm=NULL
beta1.glm=NULL       

beta0.garma=NULL
beta1.garma=NULL
teta.garma=NULL

beta0.glarma=NULL
beta1.glarma=NULL
teta.glarma=NULL

beta0.glarma2=NULL
beta1.glarma2=NULL
teta.glarma2=NULL



######################################
# Programa principal
######################################

# Especificações iniciais
#########################################################################
N = 100
                        # tamanho da serie
                                              
K=50

                        # Burn-in
BetaV=c(1,0.3)          # Parametros reais das covariaveis
                                
PhiV=c()             # Parametros reais dos termos autorregressivos
        
TetaV=c(-0.95)              # Parametros reais dos termos medias moveis
                
c=0.1                   # Constante para o log(Y), caso Y=0

lamb2=1		        # parametro para o residuo de Pearson no GLARMA
lamb=0.5
Ind=0			# Indicadora para Like e LikeG (se Ind=0, retorna a	#
	
MC=10           
#########################################################################

# Declaração de variáveis
X=NULL
                        # Covariavel sem burn-in
nb=length(BetaV)        # Número de covariaveis
np=length(PhiV)         # Numero de componentes AR
nt=length(TetaV)       # Numero de componentes MA
ParametrosV=NULL
npar=nb+np+nt
EST=rep(0,nb+np+nt) # Inicialização dos parâmetros para a estimação por MV
for (i in 1:nb) {
   ParametrosV[i]=BetaV[i]}
if (np>0) {
   for (i in (nb+1):(nb+np)) {
      ParametrosV[i]=PhiV[i-nb]}
}
if (nt>0) {
   for (i in (nb+np+1):npar) {
      ParametrosV[i]=TetaV[i-nb-np]}
}

erroGL5 = NULL
erroGL1 = NULL
okGL5 = NULL
okGL1 = NULL

############################################################################
        set.seed(197)
 
for(r in 1:MC){ # Inicio do Monte Carlo
# Geração da serie
                                                                                                                        
#        set.seed(97)
                                                                
                                                                                          
       # Geracao da covariavel com burn-in
                                

       XE = c(1:(N+K))/(N+K)
                                                                

      #XE = rnorm(N+K)
                                                                

       #XE=arima.sim(list(ar=c(0.4)),n=N+K)
                                        
                                                                                                 

        # Eliminacao dos primeiros K valores da covariavel
                
         for(i in 1:N){   
                                                                        
         X[i]=XE[i+K]
                                                                        

        }                                                                                        

 #      set.seed(97)        
                                                                           
        Y = glarma.serie(N, K, nb, np, nt, ParametrosV, XE, lamb)


######################################################################

        # Geração da Verossimilhança GARMA
                                        
                                                                          
  #      set.seed(97)
                                                                
                                                                                                
        EMV.Garma = LikeGMA(ParametrosV, Y, nb, nt, X, c, Ind)

######################################################################
###################################################################

    # Geração da Verossimilhança GLARMA
                                        
                                                                                    
   #     set.seed(97)
                                           
                                                                                                
        EMV.Glarma = LikeMA(ParametrosV,Y, nb, nt, lamb, X, Ind)

######################################################################
######################################################################

    # Geração da Verossimilhança GLARMA com lamb=1.0
                                        
                                                                                    
        #set.seed(97)
                                                                
                                                                                                
        EMV.Glarma2 = LikeMA(ParametrosV,Y, nb, nt, lamb2, X, Ind)

     
######################################################################
######################################################################


 #                    Ajuste dos modelos                           
	
	# 2) GLM
		M1=glm(Y~X,family=poisson)
		summary(M1[r])
            beta0.glm[r]= M1$coef[1]
		beta1.glm[r]=M1$coef[2]


#*****************************************************************
	# 3) GARMA
		Ind=0	# Maximiza a verossimilhança
		MGA=nlm(LikeGMA, EST, Y, nb, nt, X, c,Ind)
		Ind=1	# Calcula Mi_t
		MiGA=LikeGMA(MGA$est, Y, nb, nt, X, c, Ind)
	
            beta0.garma[r]= MGA$est[1]
		beta1.garma[r]=MGA$est[2]
		#phi.garma[r]=MGA$est[3]
           teta.garma[r]=MGA$est[3]

#*****************************************************************
	# 4) GLARMA
       outGL5 <- try(GLARMA(Y,X,nb,nt,lamb)											)
        if(class(outGL5) != "try-error"){
  		beta0.glarma[r]= outGL5$beta0
		beta1.glarma[r]= outGL5$beta1
#		phi1.glarma[r]= outGL5$phi
		teta.glarma[r]= outGL5$teta
	      okGL5 <- c(okGL5,r)
        }
        else erroGL5 <- c(erroGL5,r)

#*****************************************************************
	# 5) GLARMA com lamb= 1.0
       outGL1 <- try(GLARMA(Y,X,nb,nt,lamb2)											)
        if(class(outGL1) != "try-error"){
  		beta0.glarma2[r]= outGL1$beta0
		beta1.glarma2[r]= outGL1$beta1
#		phi1.glarma2[r]= outGL1$phi
		teta.glarma2[r]= outGL1$teta
	      okGL1 <- c(okGL1,r)
        }
        else erroGL1 <- c(erroGL1,r)

print(r)        

} # Fim do Monte Carlo

### Estimativa GLM###

if(okGL1 & okGL5) okGL

MBeta0GLM=mean(beta0.glm[okGL5])
EQMBeta0GLM=mean(beta0.glm-BetaV[1])^2

glm0=c(MBeta0GLM,EQMBeta0GLM)


MBeta1GLM=mean(beta1.glm)
EQMBeta1GLM=mean(beta1.glm-BetaV[2])^2

glm1=c(MBeta1GLM,EQMBeta1GLM)


### Estimativa GARMA###

MBeta0GARMA=mean(beta0.garma)
EQMBeta0GARMA=mean(beta0.garma-BetaV[1])^2

garma0=c(MBeta0GARMA,EQMBeta0GARMA)


MBeta1GARMA=mean(beta1.garma)
EQMBeta1GARMA=mean(beta1.garma-BetaV[2])^2

garma1=c(MBeta1GARMA,EQMBeta1GARMA)



MTetaGARMA=mean(teta.garma)
EQMTetaGARMA=mean(teta.garma-TetaV)^2

garmateta=c(MTetaGARMA,EQMTetaGARMA)


### Estimativa GLARMA###

MBeta0GLARMA=mean(beta0.glarma,na.rm=T)
EQMBeta0GLARMA=mean((beta0.glarma-BetaV[1])^2,na.rm=T)

glarma0=c(MBeta0GLARMA,EQMBeta0GLARMA)


MBeta1GLARMA=mean(beta1.glarma)
EQMBeta1GLARMA=mean(beta1.glarma-BetaV[2])^2

glarma1=c(MBeta1GLARMA,EQMBeta1GLARMA)


MTetaGLARMA=mean(teta.glarma)
EQMTetaGLARMA=mean(teta.glarma-TetaV)^2

glarmateta=c(MTetaGLARMA,EQMTetaGLARMA)


### Estimativa GLARMA2###

MBeta0GLARMA2=mean(beta0.glarma2)
EQMBeta0GLARMA2=mean(beta0.glarma2-BetaV[1])^2

glarma02=c(MBeta0GLARMA2,EQMBeta0GLARMA2)


MBeta1GLARMA2=mean(beta1.glarma2)
EQMBeta1GLARMA2=mean(beta1.glarma2-BetaV[2])^2

glarma12=c(MBeta1GLARMA2,EQMBeta1GLARMA2)

MTetaGLARMA2=mean(teta.glarma2)
EQMTetaGLARMA2=mean(teta.glarma2-TetaV)^2

glarmateta2=c(MTetaGLARMA2,EQMTetaGLARMA2)



garma0
garma1
garmateta


glarma0
glarma1
glarmateta

glarma02
glarma12
glarmateta2


glm0
glm1

#######Guardando os vetores

b0glama1p=c(beta0.garma,beta0.glarma,beta0.glarma2,beta0.glm)
b1glama1p=c(beta1.garma,beta1.glarma,beta1.glarma2,beta1.glm)
tetaglama1p=c(teta.garma,teta.glarma,teta.glarma2)


########################## GLARMA AR1################################


######################################
# Programa principal
######################################

# Especificações iniciais
#########################################################################
N = 100
                        # tamanho da serie
                                              
K=50

                        # Burn-in
BetaV=c(1,0.3)          # Parametros reais das covariaveis
                                
PhiV=c(0.95)             # Parametros reais dos termos autorregressivos
        
TetaV=c()              # Parametros reais dos termos medias moveis
                
c=0.1                   # Constante para o log(Y), caso Y=0

lamb2=1		        # parametro para o residuo de Pearson no GLARMA
lamb=0.5
Ind=0			# Indicadora para Like e LikeG (se Ind=0, retorna a	#
	
MC=10
#########################################################################

# Declaração de variáveis
X=NULL
                        # Covariavel sem burn-in
nb=length(BetaV)        # Número de covariaveis
np=length(PhiV)         # Numero de componentes AR
nt=length(TetaV)       # Numero de componentes MA
ParametrosV=NULL
npar=nb+np+nt
EST=rep(0,nb+np+nt) # Inicialização dos parâmetros para a estimação por MV
for (i in 1:nb) {
   ParametrosV[i]=BetaV[i]}
if (np>0) {
   for (i in (nb+1):(nb+np)) {
      ParametrosV[i]=PhiV[i-nb]}
}
if (nt>0) {
   for (i in (nb+np+1):npar) {
      ParametrosV[i]=TetaV[i-nb-np]}
}

############################################################################
set.seed(92)
for(r in 1:MC){ # Inicio do Monte Carlo
# Geração da serie
                                                                                                                        
        #set.seed(97)
                                                                
                                                                                          
       # Geracao da covariavel com burn-in
                                

       XE = c(1:(N+K))/(N+K)
                                                                

      #XE = rnorm(N+K)
                                                                

       #XE=arima.sim(list(ar=c(0.4)),n=N+K)
                                        
                                                                                                 

        # Eliminacao dos primeiros K valores da covariavel
                
         for(i in 1:N){   
                                                                        
         X[i]=XE[i+K]
                                                                        

        }                                                                                        

                                                                                   
        Y = glarma.serie(N, K, nb, np, nt, ParametrosV, XE, lamb)


######################################################################

        # Geração da Verossimilhança GARMA
                                        
                                                                          
       # set.seed(97)
                                                                
                                                                                                
        EMV.Garma = LikeGAR(ParametrosV, Y, nb, np, X, c, Ind)

######################################################################
###################################################################

    # Geração da Verossimilhança GLARMA
                                        
                                                                                    
#        set.seed(97)
                                           
                                                                                                
        EMV.Glarma = LikeAR(ParametrosV,Y, nb, np, lamb, X, Ind)

######################################################################
######################################################################

    # Geração da Verossimilhança GLARMA com lamb=1.0
                                        
                                                                                    
        #set.seed(97)
                                                                
                                                                                                
        EMV.Glarma2 = LikeAR(ParametrosV,Y, nb, np, lamb2, X, Ind)

     
######################################################################
######################################################################


 #                    Ajuste dos modelos                           
	
	# 2) GLM
		M1=glm(Y~X,family=poisson)
		summary(M1[r])
            beta0.glm[r]= M1$coef[1]
		beta1.glm[r]=M1$coef[2]


#*****************************************************************
	# 3) GARMA
		Ind=0	# Maximiza a verossimilhança
		MGA=nlm(LikeGAR, EST, Y, nb, np, X, c,Ind)
		Ind=1	# Calcula Mi_t
		MiGA=LikeGAR(MGA$est, Y, nb, np, X, c, Ind)
	
            beta0.garma[r]= MGA$est[1]
		beta1.garma[r]=MGA$est[2]
		phi.garma[r]=MGA$est[3]
      #      teta.garma[r]=MGA$est[2]

#*****************************************************************
	# 4) GLARMA
		Ind=0	# Maximiza a verossimilhança
		MGL=nlm(LikeAR, EST, Y, nb, np, lamb, X, Ind)
            Ind=1	# Calcula Mi_t
		MiGL=LikeAR(MGL$est, Y, nb, np, lamb, X, Ind)


		beta0.glarma[r]= MGL$est[1]
		beta1.glarma[r]=MGL$est[2]
        	phi.glarma[r]=MGL$est[3]
            #teta.glarma[r]=MGL$est[2]

#*****************************************************************
	# 5) GLARMA com lamb= 1.0
		Ind=0	# Maximiza a verossimilhança
		MGL2=nlm(LikeAR, EST, Y, nb, np, lamb2, X, Ind)
            Ind=1	# Calcula Mi_t
		MiGL2=LikeAR(MGL2$est, Y, nb, np, lamb2, X, Ind)


		beta0.glarma2[r]= MGL2$est[1]
		beta1.glarma2[r]=MGL2$est[2]
		phi.glarma2[r]=MGL2$est[3]
            #teta.glarma2[r]=MGL2$est[2]
print(r)        

} # Fim do Monte Carlo

### Estimativa GLM###

MBeta0GLM=mean(beta0.glm)
EQMBeta0GLM=mean(beta0.glm-BetaV[1])^2

glm0=c(MBeta0GLM,EQMBeta0GLM)


MBeta1GLM=mean(beta1.glm)
EQMBeta1GLM=mean(beta1.glm-BetaV[2])^2

glm1=c(MBeta1GLM,EQMBeta1GLM)


### Estimativa GARMA###

MBeta0GARMA=mean(beta0.garma)
EQMBeta0GARMA=mean(beta0.garma-BetaV[1])^2

garma0=c(MBeta0GARMA,EQMBeta0GARMA)


MBeta1GARMA=mean(beta1.garma)
EQMBeta1GARMA=mean(beta1.garma-BetaV[2])^2

garma1=c(MBeta1GARMA,EQMBeta1GARMA)


MPhiGARMA=mean(phi.garma)
EQMPhiGARMA=mean(phi.garma-PhiV)^2

garmaphi=c(MPhiGARMA,EQMPhiGARMA)

#MTetaGARMA=mean(teta.garma)
#EQMTetaGARMA=mean(teta.garma-TetaV)^2

#garmateta=c(MTetaGARMA,EQMTetaGARMA)


### Estimativa GLARMA###

MBeta0GLARMA=mean(beta0.glarma)
EQMBeta0GLARMA=mean(beta0.glarma-BetaV[1])^2

glarma0=c(MBeta0GLARMA,EQMBeta0GLARMA)


MBeta1GLARMA=mean(beta1.glarma)
EQMBeta1GLARMA=mean(beta1.glarma-BetaV[2])^2

glarma1=c(MBeta1GLARMA,EQMBeta1GLARMA)


MPhiGLARMA=mean(phi.glarma)
EQMPhiGLARMA=mean(phi.glarma-PhiV)^2

glarmaphi=c(MPhiGLARMA,EQMPhiGLARMA)

#MTetaGLARMA=mean(teta.glarma)
#EQMTetaGLARMA=mean(teta.glarma-TetaV)^2

#glarmateta=c(MTetaGLARMA,EQMTetaGLARMA)


### Estimativa GLARMA2###

MBeta0GLARMA2=mean(beta0.glarma2)
EQMBeta0GLARMA2=mean(beta0.glarma2-BetaV[1])^2

glarma02=c(MBeta0GLARMA2,EQMBeta0GLARMA2)


MBeta1GLARMA2=mean(beta1.glarma2)
EQMBeta1GLARMA2=mean(beta1.glarma2-BetaV[2])^2

glarma12=c(MBeta1GLARMA2,EQMBeta1GLARMA2)


MPhiGLARMA2=mean(phi.glarma2)
EQMPhiGLARMA2=mean(phi.glarma2-PhiV)^2

glarmaphi2=c(MPhiGLARMA2,EQMPhiGLARMA2)

#MTetaGLARMA2=mean(teta.glarma2)
#EQMTetaGLARMA2=mean(teta.glarma2-TetaV)^2

#glarmateta2=c(MTetaGLARMA2,EQMTetaGLARMA2)



garma0
garma1
garmaphi
#garmateta


glarma0
glarma1
glarmaphi
#glarmateta

glarma02
glarma12
glarmaphi2
#glarmateta2


glm0
glm1

#######Guardando os vetores

b0glaar1=c(beta0.garma,beta0.glarma,beta0.glarma2,beta0.glm)
b1glaar1=c(beta1.garma,beta1.glarma,beta1.glarma2,beta1.glm)
phiglaar1=c(phi.garma,phi.glarma,phi.glarma2)


###############GARMA MA1#######################################

######################################
# Programa principal
######################################

# Especificações iniciais
#########################################################################
N = 100
                        # tamanho da serie
                                                
K=50
                        # Burn-in
BetaV=c(1,0.3)          # Parametros reais das covariaveis
                                
PhiV=c()             # Parametros reais dos termos autorregressivos
        
TetaV=c(0.9)              # Parametros reais dos termos medias moveis
                
c=0.1                   # Constante para o log(Y), caso Y=0

lamb2=1		        # parametro para o residuo de Pearson no GLARMA
lamb=0.5
Ind=0			# Indicadora para Like e LikeG (se Ind=0, retorna a	#
	
MC=10
#########################################################################

# Declaração de variáveis
X=NULL
                        # Covariavel sem burn-in
nb=length(BetaV)        # Número de covariaveis
np=length(PhiV)         # Numero de componentes AR
nt=length(TetaV)        # Numero de componentes MA
ParametrosV=NULL
npar=nb+np+nt
EST=rep(0,nb+np+nt) # Inicialização dos parâmetros para a estimação por MV
for (i in 1:nb) {
   ParametrosV[i]=BetaV[i]}
if (np>0) {
   for (i in (nb+1):(nb+np)) {
      ParametrosV[i]=PhiV[i-nb]}
}
if (nt>0) {
   for (i in (nb+np+1):npar) {
      ParametrosV[i]=TetaV[i-nb-np]}
}


################################################################### 
# Vetores para armazenar as estimativas


beta0.glm=NULL
beta1.glm=NULL       
beta0.garma=NULL
beta1.garma=NULL
phi.garma=NULL
beta0.glarma=NULL
beta1.glarma=NULL
phi.glarma=NULL
beta0.glarma2=NULL
beta1.glarma2=NULL
phi.glarma2=NULL
teta.glarma=NULL
teta.garma=NULL
teta.glarma2=NULL



###################################################################
#  set.seed(97)

for(r in 1:MC){ # Inicio do Monte Carlo

        # Geração da serie
                                        
                                                                                      

       # set.seed(97)
                                                                
                                                                                          
       # Geracao da covariavel com burn-in
                                

       XE = c(1:(N+K))/(N+K)
                                                                
 
      #XE = rnorm(N+K)
                                                                

     # XE=arima.sim(list(ar=c(0.4)),n=N+K)
                                        
                                                                                                 

        # Eliminacao dos primeiros K valores da covariavel
                
         for(i in 1:N){   
                                                                        
         X[i]=XE[i+K]
                                                                        

        }                                                                                        

                                                                                                

        Y = garma.serie(N, K, nb, np, nt, ParametrosV, XE, c)

      
 ###################################################################

    # Geração da Verossimilhança GLARMA com lamb=0.5
                                        
                                                                                    
        #set.seed(97)
                                                                
                                                                                                
        EMV.Glarma = LikeMA(ParametrosV,Y, nb, nt, lamb, X, Ind)

     
######################################################################

    # Geração da Verossimilhança GLARMA com lamb=1.0
                                        
                                                                                    
        #set.seed(97)
                                                                
                                                                                                
        EMV.Glarma2 = LikeMA(ParametrosV,Y, nb, nt, lamb2, X, Ind)

     
######################################################################

        # Geração da Verossimilhança GARMA
                                        
                                                                          
       # set.seed(97)
                                                                
                                                                                                
        EMV.Garma = LikeGMA(ParametrosV, Y, nb, nt, X, c, Ind)

######################################################################


 #                    Ajuste dos modelos                           
	
	# 2) GLM
		M1=glm(Y~X,family=poisson)
		summary(M1[r])
		
            beta0.glm[r]= M1$coef[1]
		beta1.glm[r]=M1$coef[2]


#*****************************************************************
	# 3) GARMA
		Ind=0	# Maximiza a verossimilhança
		MGA=nlm(LikeGMA, EST, Y, nb, nt, X, c,Ind)
		Ind=1	# Calcula Mi_t
		MiGA=LikeGMA(MGA$est, Y, nb, nt, X, c, Ind)
	
            beta0.garma[r]= MGA$est[1]
		beta1.garma[r]=MGA$est[2]
		#phi.garma[r]=MGA$est[2]
            teta.garma[r]=MGA$est[3]

#*****************************************************************
	# 4) GLARMA
		Ind=0	# Maximiza a verossimilhança
		MGL=nlm(LikeMA, EST, Y, nb, nt, lamb, X, Ind)
                Ind=1	# Calcula Mi_t
		MiGL=LikeMA(MGL$est, Y, nb, nt, lamb, X, Ind)


		beta0.glarma[r]= MGL$est[1]
		beta1.glarma[r]=MGL$est[2]
#        	phi.glarma[r]=MGL$est[2]
            teta.glarma[r]=MGL$est[3]

#*****************************************************************
	# 5) GLARMA com lamb= 0.5
		Ind=0	# Maximiza a verossimilhança
		MGL2=nlm(LikeMA, EST, Y, nb, nt, lamb2, X, Ind)
                Ind=1	# Calcula Mi_t
		MiGL2=LikeMA(MGL2$est, Y, nb, nt, lamb2, X, Ind)


		beta0.glarma2[r]= MGL2$est[1]
		beta1.glarma2[r]=MGL2$est[2]
#		phi.glarma2[r]=MGL2$est[4]
            teta.glarma2[r]=MGL2$est[3]
        
print(r)
} # Fim do Monte Carlo

### Estimativa GLM###

MBeta0GLM=mean(beta0.glm)
EQMBeta0GLM=mean(beta0.glm-BetaV[1])^2

glm0=c(MBeta0GLM,EQMBeta0GLM)


MBeta1GLM=mean(beta1.glm)
EQMBeta1GLM=mean(beta1.glm-BetaV[2])^2

glm1=c(MBeta1GLM,EQMBeta1GLM)


### Estimativa GARMA###

MBeta0GARMA=mean(beta0.garma)
EQMBeta0GARMA=mean(beta0.garma-BetaV[1])^2

garma0=c(MBeta0GARMA,EQMBeta0GARMA)


MBeta1GARMA=mean(beta1.garma)
EQMBeta1GARMA=mean(beta1.garma-BetaV[2])^2

garma1=c(MBeta1GARMA,EQMBeta1GARMA)


#MPhiGARMA=mean(phi.garma)
#EQMPhiGARMA=mean(phi.garma-PhiV)^2

#garmaphi=c(MPhiGARMA,EQMPhiGARMA)

MTetaGARMA=mean(teta.garma)
EQMTetaGARMA=mean(teta.garma-TetaV)^2

garmateta=c(MTetaGARMA,EQMTetaGARMA)



### Estimativa GLARMA###

MBeta0GLARMA=mean(beta0.glarma)
EQMBeta0GLARMA=mean(beta0.glarma-BetaV[1])^2

glarma0=c(MBeta0GLARMA,EQMBeta0GLARMA)


MBeta1GLARMA=mean(beta1.glarma)
EQMBeta1GLARMA=mean(beta1.glarma-BetaV[2])^2

glarma1=c(MBeta1GLARMA,EQMBeta1GLARMA)


#MPhiGLARMA=mean(phi.glarma)
#EQMPhiGLARMA=mean(phi.glarma-PhiV)^2

#glarmaphi=c(MPhiGLARMA,EQMPhiGLARMA)

MTetaGLARMA=mean(teta.glarma)
EQMTetaGLARMA=mean(teta.glarma-TetaV)^2

glarmateta=c(MTetaGLARMA,EQMTetaGLARMA)


### Estimativa GLARMA2###

MBeta0GLARMA2=mean(beta0.glarma2)
EQMBeta0GLARMA2=mean(beta0.glarma2-BetaV[1])^2

glarma02=c(MBeta0GLARMA2,EQMBeta0GLARMA2)


MBeta1GLARMA2=mean(beta1.glarma2)
EQMBeta1GLARMA2=mean(beta1.glarma2-BetaV[2])^2

glarma12=c(MBeta1GLARMA2,EQMBeta1GLARMA2)


#MPhiGLARMA2=mean(phi.glarma2)
#EQMPhiGLARMA2=mean(phi.glarma2-PhiV)^2

#glarmaphi2=c(MPhiGLARMA2,EQMPhiGLARMA2)

MTetaGLARMA2=mean(teta.glarma2)
EQMTetaGLARMA2=mean(teta.glarma2-TetaV)^2

glarmateta2=c(MTetaGLARMA2,EQMTetaGLARMA2)


garma0
garma1
#garmaphi
garmateta


glarma0
glarma1
#glarmaphi
glarmateta

glarma02
glarma12
#glarmaphi2
glarmateta2


glm0
glm1


################################# GARMA AR1#####################################


#################################################################################

# Vetores para armazenar as estimativas


beta0.glm=NULL
beta1.glm=NULL       
beta0.garma=NULL
beta1.garma=NULL
phi.garma=NULL

beta0.glarma=NULL
beta1.glarma=NULL
phi.glarma=NULL

beta0.glarma2=NULL
beta1.glarma2=NULL
phi.glarma2=NULL


# Especificações iniciais
#########################################################################
N = 100
                        # tamanho da serie
                                              
K=50

                        # Burn-in
BetaV=c(1,0.3)          # Parametros reais das covariaveis
                                
PhiV=c(-0.5)             # Parametros reais dos termos autorregressivos
        
TetaV=c()              # Parametros reais dos termos medias moveis
                
c=0.1                   # Constante para o log(Y), caso Y=0

lamb2=1		        # parametro para o residuo de Pearson no GLARMA
lamb=0.5
Ind=0			# Indicadora para Like e LikeG (se Ind=0, retorna a	#
	
MC=1000           
#########################################################################

# Declaração de variáveis
X=NULL
                        # Covariavel sem burn-in
nb=length(BetaV)        # Número de covariaveis
np=length(PhiV)         # Numero de componentes AR
nt=length(TetaV)       # Numero de componentes MA
ParametrosV=NULL
npar=nb+np+nt
EST=rep(0,nb+np+nt) # Inicialização dos parâmetros para a estimação por MV
for (i in 1:nb) {
   ParametrosV[i]=BetaV[i]}
if (np>0) {
   for (i in (nb+1):(nb+np)) {
      ParametrosV[i]=PhiV[i-nb]}
}
if (nt>0) {
   for (i in (nb+np+1):npar) {
      ParametrosV[i]=TetaV[i-nb-np]}
}

############################################################################

for(r in 1:MC){ # Inicio do Monte Carlo
# Geração da serie
                                        
                                                                                      

        #set.seed(97)
                                                                
                                                                                          
       # Geracao da covariavel com burn-in
                                

       XE = c(1:(N+K))/(N+K)
                                                                

      #XE = rnorm(N+K)
                                                                

       #XE=arima.sim(list(ar=c(0.4)),n=N+K)
                                        
                                                                                                 

        # Eliminacao dos primeiros K valores da covariavel
                
         for(i in 1:N){   
                                                                        
         X[i]=XE[i+K]
                                                                        

        }                                                                                        

                                                                                   
        Y = garma.serie(N, K, nb, np, nt, ParametrosV, XE, c)


######################################################################

        # Geração da Verossimilhança GARMA
                                        
                                                                          
       # set.seed(97)
                                                                
                                                                                                
        EMV.Garma = LikeGAR(ParametrosV, Y, nb, np, X, c, Ind)

######################################################################
###################################################################

    # Geração da Verossimilhança GLARMA
                                        
                                                                                    
#        set.seed(97)
                                           
                                                                                                
        EMV.Glarma = LikeAR(ParametrosV,Y, nb, np, lamb, X, Ind)

######################################################################
######################################################################

    # Geração da Verossimilhança GLARMA com lamb=1.0
                                        
                                                                                    
        #set.seed(97)
                                                                
                                                                                                
        EMV.Glarma2 = LikeAR(ParametrosV,Y, nb, np, lamb2, X, Ind)

     
######################################################################
######################################################################


 #                    Ajuste dos modelos                           
	
	# 2) GLM
		M1=glm(Y~X,family=poisson)
		summary(M1[r])
            beta0.glm[r]= M1$coef[1]
		beta1.glm[r]=M1$coef[2]


#*****************************************************************
	# 3) GARMA
		Ind=0	# Maximiza a verossimilhança
		MGA=nlm(LikeGAR, EST, Y, nb, np, X, c,Ind)
		Ind=1	# Calcula Mi_t
		MiGA=LikeGAR(MGA$est, Y, nb, np, X, c, Ind)
	
            beta0.garma[r]= MGA$est[1]
		beta1.garma[r]=MGA$est[2]
		phi.garma[r]=MGA$est[3]
           

#*****************************************************************
	# 4) GLARMA
		Ind=0	# Maximiza a verossimilhança
		MGL=nlm(LikeAR, EST, Y, nb, np, lamb, X, Ind)
            Ind=1	# Calcula Mi_t
		MiGL=LikeAR(MGL$est, Y, nb, np, lamb, X, Ind)


		beta0.glarma[r]= MGL$est[1]
		beta1.glarma[r]=MGL$est[2]
        	phi.glarma[r]=MGL$est[3]
           

#*****************************************************************
	# 5) GLARMA com lamb= 0.5
		Ind=0	# Maximiza a verossimilhança
		MGL2=nlm(LikeAR, EST, Y, nb, np, lamb2, X, Ind)
            Ind=1	# Calcula Mi_t
		MiGL2=LikeAR(MGL2$est, Y, nb, np, lamb2, X, Ind)


		beta0.glarma2[r]= MGL2$est[1]
		beta1.glarma2[r]=MGL2$est[2]
		phi.glarma2[r]=MGL2$est[3]
           
print(r)        

} # Fim do Monte Carlo

### Estimativa GLM###

MBeta0GLM=mean(beta0.glm)
EQMBeta0GLM=mean(beta0.glm-BetaV[1])^2

glm0=c(MBeta0GLM,EQMBeta0GLM)


MBeta1GLM=mean(beta1.glm)
EQMBeta1GLM=mean(beta1.glm-BetaV[2])^2

glm1=c(MBeta1GLM,EQMBeta1GLM)


### Estimativa GARMA###

MBeta0GARMA=mean(beta0.garma)
EQMBeta0GARMA=mean(beta0.garma-BetaV[1])^2

garma0=c(MBeta0GARMA,EQMBeta0GARMA)

MBeta1GARMA=mean(beta1.garma)
EQMBeta1GARMA=mean(beta1.garma-BetaV[2])^2

garma1=c(MBeta1GARMA,EQMBeta1GARMA)


MPhiGARMA=mean(phi.garma)
EQMPhiGARMA=mean(phi.garma-PhiV[1])^2

garmaphi=c(MPhiGARMA,EQMPhiGARMA)

MPhiPhiGARMA=mean(phi.phi.garma)
EQMPhiPhiGARMA=mean(phi.phi.garma-PhiV[2])^2

garmaphiphi=c(MPhiPhiGARMA,EQMPhiPhiGARMA)


### Estimativa GLARMA###

MBeta0GLARMA=mean(beta0.glarma)
EQMBeta0GLARMA=mean(beta0.glarma-BetaV[1])^2

glarma0=c(MBeta0GLARMA,EQMBeta0GLARMA)


MBeta1GLARMA=mean(beta1.glarma)
EQMBeta1GLARMA=mean(beta1.glarma-BetaV[2])^2

glarma1=c(MBeta1GLARMA,EQMBeta1GLARMA)


MPhiGLARMA=mean(phi.glarma)
EQMPhiGLARMA=mean(phi.glarma-PhiV[1])^2

glarmaphi=c(MPhiGLARMA,EQMPhiGLARMA)

MPhiPhiGLARMA=mean(phi.phi.glarma)
EQMPhiPhiGLARMA=mean(phi.phi.glarma-PhiV[2])^2

glarmaphiphi=c(MPhiPhiGLARMA,EQMPhiPhiGLARMA)


### Estimativa GLARMA2###

MBeta0GLARMA2=mean(beta0.glarma2)
EQMBeta0GLARMA2=mean(beta0.glarma2-BetaV[1])^2

glarma02=c(MBeta0GLARMA2,EQMBeta0GLARMA2)


MBeta1GLARMA2=mean(beta1.glarma2)
EQMBeta1GLARMA2=mean(beta1.glarma2-BetaV[2])^2

glarma12=c(MBeta1GLARMA2,EQMBeta1GLARMA2)


MPhiGLARMA2=mean(phi.glarma2)
EQMPhiGLARMA2=mean(phi.glarma2-PhiV[1])^2

glarmaphi2=c(MPhiGLARMA2,EQMPhiGLARMA2)

MPhiPhiGLARMA2=mean(phi.phi.glarma2)
EQMPhiPhiGLARMA2=mean(phi.phi.glarma2-PhiV[2])^2

glarmaphiphi2=c(MPhiPhiGLARMA2,EQMPhiPhiGLARMA2)



garma0
garma1
garmaphi
garmaphiphi


glarma0
glarma1
glarmaphi
glarmaphiphi

glarma02
glarma12
glarmaphi2
glarmaphiphi2


glm0
glm1

