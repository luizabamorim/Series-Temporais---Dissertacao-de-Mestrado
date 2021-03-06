
####FUN��O GARMA

garma.serie <- function(n, k, nb, np, nt, Parametros, x, c){

 
        T=k+n
                                # Tamanho da serie + Burn-in

        mi    = rep(0, T)
                
         w     = rep(0, T)

        erro  = rep(0, T) 
         z     = rep(0, T)
                
         serie = rep(0, T)
                # Serie com burn-in
         y     = rep(0, n)
                # Serie sem burn-in

         # In�cializa��o dos par�metros

        if (nb>0) {
                  beta <- Parametros[1:nb]
         }

        if (np>0) {
              fi <- Parametros[(1+nb):(nb+np)]

        }
         if (nt>0) {

             teta <- Parametros[(nb+np+1):length(Parametros)]

        }

         # In�cio do procedimento recursivo

  
        if (np==0) {
           
        if (nb==1) {
                          w[1] = beta[1]
                  }
  
                if (nb>1) {
                          w[1] = beta[1]+beta[2]*x[1]
                  }

        }
           if (np>0) {
            if (nb==1) {

           
        for (j in 1:np) {
                          w[1] = w[1]-beta[1]*fi[j]
                  }

                w[1] = w[1]+beta[1]
            }
 
           if (nb>1) {
            
        for (j in 1:np) {
                          w[1] = w[1]-fi[j]*(beta[1]+beta[2]*x[1])
                  }

           w[1] = w[1]+beta[1]+beta[2]*x[1]

           }
         }

 
        mi[1] = exp(w[1])
         serie[1] = rpois(1, mi[1])

 
        for(t in 2:T){
            if(serie[t-1]>0){

                erro[t-1]<-log(serie[t-1])-w[t-1]

           }
            if(serie[t-1]==0){

                erro[t-1]<-log(c)-w[t-1]
            }

           if (np>0) {
                  for (j in 1:min(c(np,(t-1)))) {
                     if (nb==0) {
                             if(serie[t-j]>0){
                            z[t] <- z[t] + fi[j]*log(serie[t-j])
                             }

                   
        if(serie[t-j]==0){
                            z[t]<-z[t] + fi[j]*log(c)
          
           
        }
                    }
                     if (nb==1) {
                             if(serie[t-j]>0){
                            z[t] <- z[t] + fi[j]*(log(serie[t-j])-beta[1])

                   
        }
                             if(serie[t-j]==0){
                            z[t]<-z[t] + fi[j]*(log(c)-beta[1])
          
           
        }
                    }
                     if (nb>1) {
                             if(serie[t-j]>0){
                            z[t] <- z[t] + fi[j]*(log(serie[t-j])-beta[1]-beta[2]*x[t-j])

                   
        }
                             if(serie[t-j]==0){
                            z[t]<-z[t] + fi[j]*(log(c)-beta[1]-beta[2]*x[t-j])

                             }

                   }
                  }

           }
            if (nt>0) {
                  for (j in 1:min(c(nt,(t-1)))) {
                    z[t] <- z[t] + teta[j]*(erro[t-j])
                  }

           }
        
             if (nb==0) {# S�rie somente com parte AR ou MA

                w[t] = z[t]
            }
  
           if (nb==1) {# S�rie com intercepto + parte AR ou MA

                w[t] = beta[1] + z[t]
            }
  
           if (nb>1) {# S�rie com Intercepto + covariavel + parte AR ou MA

                w[t] =  beta[1]+beta[2]*x[t] + z[t]

           }
            mi[t] = exp(w[t])

           serie[t] = rpois(1, mi[t])
         } # Fim da constru��o da s�rie

 
        #         Eliminacao dos primeiros k valores da serie com outlier

        for(i in 1:n){   
         y[i]=serie[i+k]
         }

 
         return(y)

} # Fim da Fun��o garma.serie


###FUN��O GLARMA

glarma.serie <- function(n, k, nb, np, nt, Parametros, x, lamb){

         T=k+n
                                # Tamanho da serie + Burn-in

         mi    = rep(0, T)      # Vetor para receber os valores de mi  
                
         w     = rep(0, T)     

         erro  = rep(0, T) 

         z     = rep(0, T)
                
         serie = rep(0, T)      # Serie com burn-in
               
         y     = rep(0, n)      # Serie sem burn-in
                

        # In�cializa��o dos par�metros

        if (nb>0) {                            # Tamanho do vetor beta
                  beta <- Parametros[1:nb]
         }

        if (np>0) {                            # Tamanho do vetor phi
              phi <- Parametros[(1+nb):(nb+np)]

        }
         if (nt>0) {                            # Tamanho do vetor teta

             teta <- Parametros[(nb+np+1):length(Parametros)]

        }


         # In�cio do procedimento recursivo



         if (nb==1) {
       
             w[1]=beta[1]
                  }

         if (nb>1) {

             w[1]=beta[1]+beta[2]*x[1]
                                       }


          mi[1]=exp(w[1]) 
        

        
         serie[1]=rpois(1,mi[1])

          erro[1]=(serie[1]-mi[1])/(mi[1]^lamb)

     for(t in 2:T){
     
        if (np>0) {

      
             for (j in 1:min(c(np,(t-1)))){

         z[t]<-z[t] + phi[j]*(z[t-j]+erro[t-j])
                                                }}

         
         if (nt>0){
              
             for (k in 1:min(c(nt,(t-1)))){

             z[t]= z[t]+teta[k]*erro[t-k]

                                            }}
   
          

        if (nb==0) {# S�rie somente com parte AR ou MA

                w[t] = z[t]
            }
  
           if (nb==1) {# S�rie com intercepto + parte AR ou MA

                w[t] = beta[1] + z[t]
            }
  
           if (nb>1) {# S�rie com Intercepto + covariavel + parte AR ou MA

                w[t] =  beta[1]+beta[2]*x[t] + z[t]

           }

            mi[t] = exp(w[t])

           serie[t] = rpois(1, mi[t])

          erro[t]=(serie[t]-mi[t])/(mi[t]^lamb)

         } # fim da constru��o da s�rie

 
        #         Eliminacao dos primeiros k valores da serie com outlier

        for(i in 1:n){   
         y[i]=serie[i+k]
         }

# valores=cbind(erro,z,w,mi,y)
         return(y)

} # fim da Fun��o glarma.serie
#valores

#### FUN��O PARA VEROSSIMILHAN�A AR

   
###############################################
# Fun��o que calcula a verossimilhan�a - GLARMA
###############################################

LikeAR <- function(Parametros, y, nb, np, lamb, x, Ind) {
			
		n     <- length(y)

		if (nb>0) {
			beta <- Parametros[1:nb]
		}
		if (np>0) {
		     fi <- Parametros[(1+nb):length(Parametros)]
		}

		# C�lculo de z, w, mi e erro dados os par�metros fornecidos
				
		z    <- rep(0, n)
		w    <- rep(0, n)
		mi   <- rep(0, n)
		erro <- rep(0, n)

	  	if (nb==1) {
			w[1] = beta[1]
		}
  		if (nb>1) {
			w[1] = beta[1]+beta[2]*x[1]
		}
		mi[1]   <- exp(w[1])

	  # Usando residuos de Pearson
		erro[1] <- (y[1] - mi[1])/(mi[1]^lamb)
		L <- y[1]*w[1] - exp(w[1])
		
		for (t in 2:n) {
	  	   if (np>0) {
			for (k in 1:min(c(np,(t-1)))) {
				z[t]<-z[t] + fi[k]*(z[t-k]+erro[t-k])
	 		}
		   }
	
 		   if (nb==0) {
			w[t] = z[t]
		   }
 		   if (nb==1) {
			w[t] = beta[1] + z[t]
		   }
  		   if (nb>1) {
			w[t] =  beta[1]+beta[2]*x[t] + z[t]
		   }
		   mi[t]   <- exp(w[t])
	  # Usando residuos de Pearson
		   erro[t] <- (y[t] - mi[t])/(mi[t]^lamb)
		   L <- L + y[t]*w[t] - exp(w[t])
		}
      L=L*(-1)

if(Ind==0){return(L)} # Retorna a verossimilhan�a

if(Ind==1){ # Retorna o n�vel da s�rie
return(exp(w))
}

} # Fim da fun��o



##############################################
# Fun��o que calcula a verossimilhan�a - GARMA
##############################################

LikeGAR <- function(Parametros, y, nb, np, x, c, Ind) {

		n     <- length(y)

		if (nb>0) {
			beta <- Parametros[1:nb]
		}
		if (np>0) {
		     fi <- Parametros[(1+nb):length(Parametros)]
		}

		# C�lculo de z, w e erro dados os par�metros fornecidos
				
		z    <- rep(0, n)
		w    <- rep(0, n)
		erro <- rep(0, n)

  		if (np==0) {
	  	   if (nb==1) {
			w[1] = beta[1]
		   }
  		   if (nb>1) {
			w[1] = beta[1]+beta[2]*x[1]
		   }
		}
  		if (np>0) {
	   	   if (nb==1) {
	   		for (j in 1:np) {
			   w[1] = w[1]-beta[1]*fi[j]
			}
			w[1] = w[1]+beta[1]
	   	   }
 	   	   if (nb>1) {
	   		for (j in 1:np) {
			   w[1] = w[1]-fi[j]*(beta[1]+beta[2]*x[1])
			}
	   		w[1] = w[1]+beta[1]+beta[2]*x[1]
	   	   }
		}

		L <- y[1]*w[1] - exp(w[1])

		for (t in 2:n) {
		   if(y[t-1]>0){
			erro[t-1]<-log(y[t-1])-w[t-1]
		   }
	   	   if(y[t-1]==0){
  			erro[t-1]<-log(c)-w[t-1]
	   	   }

		   if (np>0) {
			for (j in 1:min(c(np,(t-1)))) {
  			   if (nb==0) {
			   	if(y[t-j]>0){
				   z[t] <- z[t] + fi[j]*(log(y[t-j]))
			   	}
			   	if(y[t-j]==0){
				   z[t]<-z[t] + fi[j]*(log(c))
	 		   	}
			   }
  			   if (nb==1) {
			   	if(y[t-j]>0){
				   z[t] <- z[t] + fi[j]*(log(y[t-j])-beta[1])
			   	}
			   	if(y[t-j]==0){
				   z[t]<-z[t] + fi[j]*(log(c)-beta[1])
	 		   	}
			   }
  			   if (nb>1) {
			   	if(y[t-j]>0){
				   z[t] <- z[t] + fi[j]*(log(y[t-j])-beta[1]-beta[2]*x[t-j])
			   	}
			   	if(y[t-j]==0){
				   z[t]<-z[t] + fi[j]*(log(c)-beta[1]-beta[2]*x[t-j])
		 	   	}
			   }
			}
		   }

 		   if (nb==0) {
			w[t] = z[t]
		   }
 		   if (nb==1) {
			w[t] = beta[1] + z[t]
		   }
  		   if (nb>1) {
			w[t] =  beta[1]+beta[2]*x[t] + z[t]
		   }

		   L <- L + y[t]*w[t] - exp(w[t])
		}
      L=L*(-1)


if(Ind==0){return(L)} # Retorna a verossimilhan�a

if(Ind==1){ # Retorna o n�vel da s�rie
return(exp(w))
}

} # Fim da fun��o



### FUN��O DE VEROSSIMILHAN�A PARA S�RIE COM PARTE M�DIA M�VEL

###############################################
# Fun��o que calcula a verossimilhan�a - GLARMA MA
###############################################

LikeMA <- function(Parametros, y, nb, nt, lamb, x, Ind) {
			
		n     <- length(y)

		if (nb>0) {
			beta <- Parametros[1:nb]
		}
		if (nt>0) {
		     theta <- Parametros[(1+nb):length(Parametros)]
		}

		# C�lculo de z, w, mi e erro dados os par�metros fornecidos
				
		z    <- rep(0, n)
		w    <- rep(0, n)
		mi   <- rep(0, n)
		erro <- rep(0, n)

	  	if (nb==1) {
			w[1] = beta[1]
		}
  		if (nb>1) {
			w[1] = beta[1]+beta[2]*x[1]
		}
		mi[1]   <- exp(w[1])

	  # Usando residuos de Pearson
		erro[1] <- (y[1] - mi[1])/(mi[1]^lamb)
		L <- y[1]*w[1] - exp(w[1])
		
		for (t in 2:n) {
	  	   if (nt>0) {
			for (k in 1:min(c(nt,(t-1)))) {
				z[t]<-z[t] + theta[k]*erro[t-k]
	 		}
		   }
	
 		   if (nb==0) {
			w[t] = z[t]
		   }
 		   if (nb==1) {
			w[t] = beta[1] + z[t]
		   }
  		   if (nb>1) {
			w[t] =  beta[1]+beta[2]*x[t] + z[t]
		   }
		   mi[t]   <- exp(w[t])
	  # Usando residuos de Pearson
		   erro[t] <- (y[t] - mi[t])/(mi[t]^lamb)
		   L <- L + y[t]*w[t] - exp(w[t])
		}
      L=L*(-1)

if(Ind==0){return(L)} # Retorna a verossimilhan�a

if(Ind==1){ # Retorna o n�vel da s�rie
return(exp(w))
}

} # Fim da fun��o

##############################################
# Fun��o que calcula a verossimilhan�a - GARMA MA
##############################################

LikeGMA <- function(Parametros, y, nb, nt, x, c, Ind) {

		n     <- length(y)

		if (nb>0) {
			beta <- Parametros[1:nb]
		}
		if (nt>0) {
		     theta <- Parametros[(1+nb):length(Parametros)]
		}

		# C�lculo de z, w e erro dados os par�metros fornecidos
				
		z    <- rep(0, n)
		w    <- rep(0, n)
		erro <- rep(0, n)

  	
	  	   if (nb==1) {
			w[1] = beta[1]
		   }
  		   if (nb>1) {
			w[1] = beta[1]+beta[2]*x[1]
		   }
		
  		
		L <- y[1]*w[1] - exp(w[1])

		for (t in 2:n) {
		   if(y[t-1]>0){
			erro[t-1]<-log(y[t-1])-w[t-1]
		   }
	   	   if(y[t-1]==0){
  			erro[t-1]<-log(c)-w[t-1]
	   	   }

		   if (nt>0) {
			for (j in 1:min(c(nt,(t-1)))) {
  			   
			      z[t] <- z[t] + theta[j]*erro[t-j]
			   	
	 		   	}
			  } 
 		   if (nb==0) {
			w[t] = z[t]
		   }
 		   if (nb==1) {
			w[t] = beta[1] + z[t]
		   }
  		   if (nb>1) {
			w[t] = beta[1]+beta[2]*x[t]+z[t]
		   }

		   L <- L + y[t]*w[t] - exp(w[t])
		}
      L=L*(-1)


if(Ind==0){return(L)} # Retorna a verossimilhan�a

if(Ind==1){ # Retorna o n�vel da s�rie
return(exp(w))
}

} # Fim da fun��o

GLARMA <- function(Y,X,nb,nt,lamb){

		Ind=0	# Maximiza a verossimilhan�a
	MGA=nlm(LikeMA, EST, Y, nb, nt, lamb, X, Ind)


           beta0= MGA$est[1]
		beta1=MGA$est[2]
#		phi=MGA$est[3]
           teta=MGA$est[3]
   return(list(beta0=beta0,beta1=beta1,teta=teta))
}



