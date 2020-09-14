######################################################################
############  Modelo espacial efectos mixtos   ######################
######################################################################

#Paquetes necesarios 
library(tidyverse)
library(sp)
library(spdep)
library(tmap)
library(ggplot2)
library(dplyr)
library(haven)
library(numDeriv)
library(rstan)

# Planteamiento General ---------------------------------------------------
        full.data <- read.csv(file.choose())
        # Mapa descargado de https://gadm.org/download_country_v3.html
        #setwd("D://computer_accsess_data_base//REV//Spatial//Spatial")
        peru <- readRDS(file.choose())
        names(peru)
        
        # Regiones contiguas
        peru.nb <- poly2nb(peru,snap=0.0002)

# Enfemeros ---------------------------------------------------------------


          data_enf <- full.data %>%
            filter(prof == "Enfermero/a") 
          
          data_enf$timew<-scale(data_enf$timew)
          data_enf$tweek<-scale(data_enf$tweek)
          
          View(data_enf)
          levels(data_enf$INSTITUCION)
          # Estructura de datos -----------------------------------------------------
          
          #+ tweek : modelo considerando year ! pero no tweek para buscar convergencia 
          modelo_datos_enf <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+timew , data = data_enf) #ver1
          modelo_datos_enf <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) , data = data_enf) #ver2
          modelo_datos_enf <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) + region , data = data_enf) #ver2
                modelo_datos_enf_es <- model.matrix(~ lf+ relevel(INSTITUCION,ref=2) + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) + region , data = data_enf) #esssalud
                modelo_datos_enf_ffa <- model.matrix(~ lf+ relevel(INSTITUCION,ref=3) + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) + region , data = data_enf) #esssalud
                modelo_datos_enf_minsa <- model.matrix(~ lf+ relevel(INSTITUCION,ref=4) + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) + region , data = data_enf) #esssalud
                
          View(modelo_datos_enf_es)
          View(modelo_datos_enf_ffa)
          View(modelo_datos_enf_minsa)
          
          
          
          modelo_datos_enf <-modelo_datos_enf[,-1] #eliminacion del intercepto 
          model_fin_1      <- subset(modelo_datos_enf, modelo_datos_enf[,1]==Inf) # modelo_datos[,2]# ojo cuidado con el indice de columna ( por defecto model matrix te crea intercepto y eso lo elimino)
          model_fin_2      <- subset(modelo_datos_enf, modelo_datos_enf[,1]<Inf) # modelo_datos[,2]# es la 2da columna con la variable (ingreso ) por la cual dividimos la DB
          
          model_fin_1<-model_fin_1[,-1]
          model_fin_2<-model_fin_2[,-1]
          
          
          
          # Revision de la base de datos --------------------------------------------
          
          
          dim(modelo_datos_enf)
          dim(model_fin_1)
          dim(model_fin_2)
          dim(model_fin_1)[1]+dim(model_fin_2)[1]
          length(X1)
          length(X2)
          view(model_fin_1)
          

          # Parametros STAN ---------------------------------------------------------

          
          N <- length(peru$GID_1) #numero de regiones 
          I <- diag(N)
          W <- nb2mat(peru.nb) # Matriz de regiones contiguas( matriz de pesos espaciales)
          e <- rep(1,N)
          k <- dim(model_fin_1)[2] #por que quiero el numero de columnas 
          ##############################################################################
          
          #check missing data and others 
          
          
          
          #cohercecion de la base de datos()
          #enfermeros_list_mix --> interacciones tranversales , la 2 por provincia
          enfermeros_list_mix_2  <- within(list(), {
            N<-N
            I<-I
            W<-W
            e<-e
            N1 <- nrow(data_enf[data_enf$lf==Inf,]) # split database (lf)
            N2 <- nrow(data_enf[data_enf$lf<Inf,])
            X1 <- (data_enf$sexo[data_enf$lf==Inf]=="Mujer") #dummy genero
            X2 <- (data_enf$sexo[data_enf$lf<Inf]=="Mujer")
            z<-model_fin_1
            u<-model_fin_2 # k y z deben tener si o si la misma dimension 
            K <- dim(model_fin_1)[2] # es necesario introducir la dimension xq dentro del stancode lo llamo,num of intec = n-1 per variable where n is the number of levels
            #X3 <- (data_enf$timew[data_enf$lf==Inf])           #Time working 
            #X4 <- (data_enf$timew[data_enf$lf<Inf]) # (key sensitive)
            #X5 <- (data_enf$espec[data_enf$lf==Inf]=="Sí")
            #X6 <- (data_enf$espec[data_enf$lf<Inf]=="Sí") #2110 y 4376
            #S<-K+3+1 # 3 x 1 = (niveles-1)* factor , +1 por interaccion espec sexo
            #lo de arribas antes necesario cuando eran interacciones tranversales 
            d1  <- c(data_enf$departamento[data_enf$lf==Inf])
            d2  <- c(data_enf$departamento[data_enf$lf<Inf])
            L1 <- data_enf$li[data_enf$lf==Inf]/1000 # ingreso prob.
            L2 <- data_enf$li[data_enf$lf<Inf]/1000
            U2 <- data_enf$lf[data_enf$lf<Inf]/1000
          })
          
          
          #########
          #STAN####
          #########
          scode_end_interacciones_ver1(intento) <- "
          data {
          int N;
          matrix<lower=0>[N,N] W;
          vector<upper=1>[N] e;
          matrix<lower=0,upper=1>[N,N] I;
          int N1;
          int N2;
          vector[N1] L1;  // lower bounds for right censored observations
          vector[N2] L2;  // lower bounds for interval censored obervations
          vector[N2] U2;  // upper bounds for interval censored obervations
          
          int<lower=0> K;   // number of predictors 
          int<lower=0> S;   // number of predictors + interactions S mayus ojo 

          matrix[N1, K] z;   // predictor matrix  lower bounds
          matrix[N2, K] u;   // predictor matrix  upper bounds
          
          // Covariables
          vector[N1] X1; // X for right censored observations
          vector[N2] X2; // X for interval censored observations
          // vector[N1] X3;  X for right censored observations time working
          //   vector[N2] X4;  X for interval censored observations time working
          
          
          int<lower = 1, upper = 26> d1[N1];
          int<lower = 1, upper = 26> d2[N2];
          }
          parameters {
          vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
          vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
          vector[26] alpha; // Spatially Varying intercept
          vector[26] beta;  // Spatially Varying gap
          vector[S] omega;   // coefficients for predictors, ojo S= k + interacciones 
          
          real<lower=0> phi;
          real mbeta;
          real malpha;
          real<lower = 0> sigma_b;
          real<lower = 0> sigma_a;
          real<lower=-1,upper=1> lambda_b;
          real<lower=-1,upper=1> lambda_a;
          }
          transformed parameters {
          vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
          vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
          vector[N1] mu1;
          vector[N2] mu2;
          for(i in 1:N1){
          mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+omega[4]*z[i, 4]+ omega[5]*z[i, 5] + omega[6]*z[i, 6]+omega[7]*z[i, 7]+omega[8]*z[i, 8]);
          }
          for(j in 1:N2){
          mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+omega[4]*u[j, 4]+ omega[5]*u[j, 5] + omega[6]*u[j, 6]+omega[7]*u[j, 7]+omega[8]*u[j, 8]);
          }  
          model {
          for(i in 1:N1){
          y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
          }
          for(j in 1:N2){
          y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
          }
          beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
          alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
          }
          "
          
          ##############
#####################          ##############
          ##########
          
          #nota: 
          #func:Modelo funciono adecuamente, sobre este esquema se derivan todas las estimaciones variando solo el num de regresores e interacciones
          
          scode_end_interacciones_Ver3_func <- "
          data {
          int N;
          matrix<lower=0>[N,N] W;
          vector<upper=1>[N] e;
          matrix<lower=0,upper=1>[N,N] I;
          int N1;
          int N2;
          vector[N1] L1;  // lower bounds for right censored observations
          vector[N2] L2;  // lower bounds for interval censored obervations
          vector[N2] U2;  // upper bounds for interval censored obervations
          
          int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
          int<lower=0> S;   //predictors + interactios 
          matrix[N1, K] z;   // predictor matrix  lower bounds
          matrix[N2, K] u;   // predictor matrix  upper bounds
          
          // Covariables
          vector[N1] X1; // X for right censored observations
          vector[N2] X2; // X for interval censored observations
          // vector[N1] X3;  X for right censored observations time working
          //   vector[N2] X4;  X for interval censored observations time working
          
          
          int<lower = 1, upper = 26> d1[N1];
          int<lower = 1, upper = 26> d2[N2];
          }
          parameters {
          vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
          vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
          vector[26] alpha; // Spatially Varying intercept
          vector[26] beta;  // Spatially Varying gap
          vector[S] omega;   // coefficients for predictors, ojo S= k + interacciones  
          
          real<lower=0> phi;
          real mbeta;
          real malpha;
          real<lower = 0> sigma_b;
          real<lower = 0> sigma_a;
          real<lower=-1,upper=1> lambda_b;
          real<lower=-1,upper=1> lambda_a;
          }
          transformed parameters {
          vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
          vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
          vector[N1] mu1;
          vector[N2] mu2;
          for(i in 1:N1){
        mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+omega[4]*X1[i]*z[i, 1] + omega[5]*X1[i]*z[i, 2]+omega[6]*X1[i]*z[i, 3]+omega[7]*z[i, 4]+omega[8]*z[i, 5]+omega[9]*z[i, 6]+omega[10]*z[i, 7]+omega[11]*z[i, 8]+omega[12]*z[i, 9]+omega[13]*z[i, 10]+omega[14]*z[i, 11]+omega[15]*z[i, 12]+omega[16]*z[i, 13]+omega[17]*z[i, 14]+omega[18]*z[i, 15]+omega[19]*z[i, 16]+omega[20]*z[i, 17]+omega[21]*X1[i]*z[i, 4]+omega[22]*z[i, 18]+omega[23]*z[i, 19]+omega[24]*z[i, 20]);
          }
          for(j in 1:N2){
          mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+omega[4]*X2[j]*u[j, 1] + omega[5]*X2[j]*u[j, 2]+omega[6]*X2[j]*u[j, 3]+omega[7]*u[j, 4]+omega[8]*u[j, 5]+omega[9]*u[j, 6]+omega[10]*u[j, 7]+omega[11]*u[j, 8]+omega[12]*u[j, 9]+omega[13]*u[j, 10]+omega[14]*u[j, 11]+omega[15]*u[j, 12]+omega[16]*u[j, 13]+omega[17]*u[j, 14]+omega[18]*u[j, 15]+omega[19]*u[j, 16]+omega[20]*u[j, 17]+omega[21]*X2[j]*u[j, 4]+omega[22]*u[j, 18]+omega[23]*u[j, 19]+omega[24]*u[j, 20]);
          }
          }  
          model {
          for(i in 1:N1){
          y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
          }
          for(j in 1:N2){
          y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
          }
          beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
          alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
          }
          "
          
          
 #######################         
          #interacciones provinciales(enfermeros) : region , dependientes 
          scode_end_interacciones_Ver4_func <- "
          data {
          int N;
          matrix<lower=0>[N,N] W;
          vector<upper=1>[N] e;
          matrix<lower=0,upper=1>[N,N] I;
          int N1;
          int N2;
          vector[N1] L1;  // lower bounds for right censored observations
          vector[N2] L2;  // lower bounds for interval censored obervations
          vector[N2] U2;  // upper bounds for interval censored obervations
          
          int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
             // antes predictors + interactios int<lower=0> S, pero las interacciones provinciales requieren su propio parametro
          matrix[N1, K] z;   // predictor matrix  lower bounds
          matrix[N2, K] u;   // predictor matrix  upper bounds
          
          // Covariables
          vector[N1] X1; // X for right censored observations
          vector[N2] X2; // X for interval censored observations
          // vector[N1] X3;  X for right censored observations time working
          //   vector[N2] X4;  X for interval censored observations time working
          
          
          int<lower = 1, upper = 26> d1[N1];
          int<lower = 1, upper = 26> d2[N2];
          }
          parameters {
          vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
          vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
          vector[26] alpha; // Spatially Varying intercept
          vector[26] beta;  // Spatially Varying gap
          vector[26] niuno;
          vector[26] nidos;
          vector[26] nitres;
          vector[26] nicuatro;
          
          
          vector[K] omega;   // coefficients for predictors, ojo S= k + interacciones  
          
          real<lower=0> phi;
          real mbeta;
          real malpha;
          real<lower = 0> sigma_b;
          real<lower = 0> sigma_a;
          real<lower=-1,upper=1> lambda_b;
          real<lower=-1,upper=1> lambda_a;
          }
          transformed parameters {
          vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
          vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
          vector[N1] mu1;
          vector[N2] mu2;
          for(i in 1:N1){
        mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+niuno[d1[i]]*X1[i]*z[i, 1] + nidos[d1[i]]*X1[i]*z[i, 2]+nitres[d1[i]]*X1[i]*z[i, 3]+omega[4]*z[i, 4]+omega[5]*z[i, 5]+omega[6]*z[i, 6]+omega[7]*z[i, 7]+omega[8]*z[i, 8]+omega[9]*z[i, 9]+omega[10]*z[i, 10]+omega[11]*z[i, 11]+omega[12]*z[i, 12]+omega[13]*z[i, 13]+omega[14]*z[i, 14]+omega[15]*z[i, 15]+omega[16]*z[i, 16]+omega[17]*z[i, 17]+nicuatro[d1[i]]*X1[i]*z[i, 5]+omega[18]*z[i, 18]+omega[19]*z[i, 19]+omega[20]*z[i, 20]);
          }
          for(j in 1:N2){
          mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+niuno[d2[j]]*X2[j]*u[j, 1] + nidos[d2[j]]*X2[j]*u[j, 2]+nitres[d2[j]]*X2[j]*u[j, 3]+omega[4]*u[j, 4]+omega[5]*u[j, 5]+omega[6]*u[j, 6]+omega[7]*u[j, 7]+omega[8]*u[j, 8]+omega[9]*u[j, 9]+omega[10]*u[j, 10]+omega[11]*u[j, 11]+omega[12]*u[j, 12]+omega[13]*u[j, 13]+omega[14]*u[j, 14]+omega[15]*u[j, 15]+omega[16]*u[j, 16]+omega[17]*u[j, 17]+nicuatro[d2[j]]*X2[j]*u[j, 5]+omega[18]*u[j, 18]+omega[19]*u[j, 19]+omega[20]*u[j, 20]);
          }
          }  
          model {
          for(i in 1:N1){
          y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
          }
          for(j in 1:N2){
          y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
          }
          beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
          alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
          }
          "
          
   
          
          inits_fin_enf <- function() {
            list(beta=rep(0,26),mbeta=0,sigma_b=1,lambda_b=0, 
                 alpha=rep(2,26),malpha=2,sigma_a=1,lambda_a=0,
                 phi = 0.6,
                 y_raw1 = rep(0.5,enfermeros_list_mix_2$N1), # init values > 0
                 y_raw2 = rep(0.5,enfermeros_list_mix_2$N2)  # 0 < init values < 1
            )
          }
          
          
          fit_complete_enf <- stan(model_code = scode_end_interacciones_Ver2, data = enfermeros_list_mix, 
                                chains = 3, cores = 3,init=inits_fin_enf)##1er run
          
          fit_complete_enf_intec <- stan(model_code = scode_end_interacciones_Ver2_func, data = enfermeros_list_mix, 
                                   chains = 3, cores = 3,init=inits_fin_enf)##2do run
          
          fit_complete_enf_intec_2 <- stan(model_code = scode_end_interacciones_Ver3_func, data = enfermeros_list_mix, 
                                         chains = 3, cores = 3,init=inits_fin_enf)##3ero run
          
          fit_complete_enf_intec_3 <- stan(model_code = scode_end_interacciones_Ver4_func, data = enfermeros_list_mix_2, 
                                           chains = 3, cores = 3,init=inits_fin_enf)##3ero run
          
          
          n#1er run
          print(fit_complete_enf,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
          #2do run
          print(fit_complete_enf_intec,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
          #3er run
          print(fit_complete_enf_intec_2,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
          
          #4to run
          print(fit_complete_enf_intec_3,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
          
          
          saveRDS(fit_complete_enf,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\modelo_efectos_mixtos_inicial.rds")
          saveRDS(fit_complete_enf_intec,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\modelo_efectos_mixtos_2.rds")
          saveRDS(fit_complete_enf_intec_2,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\modelo_efectos_mixtos_3.rds")
          saveRDS(fit_complete_enf_intec_3,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\modelo_efectos_mixtos_4.rds") #efectos interaccion depdientes region por depa
          
          fit_complete_enf_intec<-readRDS(file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\modelo_efectos_mixtos_2.rds")
          
          #Grafico  modelo interacciones modelo enfermeros
          library(bayesplot)
          
          color_scheme_set("blue") #comando aguanta solo uno por uno 
          fake_rhat_values <- c(1.00, 1.00,1.00,1.00, 1.01, 1.00, 1.01, 1.00,1.00,1.01,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.00,1.00)
          mcmc_intervals(fit_complete_enf_intec_2, 
                         pars=c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                "beta[11]","beta[12]","beta[13]","beta[14]","beta[15]",
                                "beta[16]","beta[17]","beta[18]","beta[19]","beta[20]",
                                "beta[21]","beta[22]","beta[23]","beta[24]","beta[25]","beta[26]"), rhat = fake_rhat_values)+
            ggplot2::labs(
              title = "Intervalos de Credibilidad",
              subtitle = "Al 95 %"
            )+  coord_flip()+ #Esto ultimo falta editar
            xaxis_text(angle = -70, vjust = 1, hjust = 0)+
            theme(
              plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))#,axis.text.x = element_text(angle = 90))
          
          
          
          
          
          
          #4to modelo interacciones  efecto espacial a nivel provincial enfermeroa
          library(ggmcmc)
          resul.gcc.l<-bind_rows(
                  plab("beta", list(depa= peru$NAME_1)) %>%
                          mutate(Coefficient = "coeficiente"))
          
          resultados.intect.prov.espacial.enf<-fit_complete_enf_intec_3
          resul.gcc.intect.prov.espacial.enf<-ggs(resultados.intect.prov.espacial.enf,par_labels=resul.gcc.l,family = "beta\\[(\\d)+\\]")
          
          ggs_caterpillar(resul.gcc.intect.prov.espacial.enf, horizontal=FALSE)+
                  labs(x = expression(beta),y="Departamentos")+
                  theme_tq()+theme(axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0.5))+
                  geom_vline(xintercept=c(0), linetype="dotted")+
                  ggtitle("Brecha Salarial")
          
          
          
###################################################################################################################################################################          
          

# Medicos -----------------------------------------------------------------
              
               full.data <- read.csv(file.choose())
          
          
              data_med <- full.data %>%
                filter(prof == "Médico") 
              
              data_med$timew<-scale(data_med$timew)
              data_med$tweek<-scale(data_med$tweek)
              
              
              
              # Estructura de datos -----------------------------------------------------
              
              
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION  + espec + labor.asis +labor.noasis+ timew+tweek+ labor.docen+tipo.contrato, data = data_med) #2do modelo corrido
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION , data = data_med) #no compila si consideras exactamente la cantidad de parametros necesaria, creo que fue por  el vector theta que lo emplee mal 
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION + timew , data = data_med) #1er modelo corrido 
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION  + espec + labor.asis +labor.noasis+ timew + labor.docen, data = data_med) #3er modelo corrido
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + year, data = data_med)
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year) , data = data_med) #mixtos por depa 
              modelo_datos_med <- model.matrix(~ lf+ INSTITUCION  + espec +depen+ labor.asis +labor.noasis+ timew + tweek + labor.docen + tipo.contrato + factor(year)+ region , data = data_med) # region mas mixtos por depa 
              
              
              modelo_datos_med <-modelo_datos_med[,-1] #eliminacion del intercepto 
              model_fin_1      <- subset(modelo_datos_med, modelo_datos_med[,1]==Inf) # modelo_datos[,2]# ojo cuidado con el indice de columna ( por defecto model matrix te crea intercepto y eso lo elimino)
              model_fin_2      <- subset(modelo_datos_med, modelo_datos_med[,1]<Inf) # modelo_datos[,2]# es la 2da columna con la variable (ingreso ) por la cual dividimos la DB
              
              model_fin_1<-model_fin_1[,-1]
              model_fin_2<-model_fin_2[,-1]
              
              
              
              # Revision de la base de datos --------------------------------------------
              
              
              dim(modelo_datos_med)
              View(model_fin_1)
              dim(model_fin_1)
              dim(model_fin_2)
              dim(model_fin_1)[1]+dim(model_fin_2)[1]
              length(X1)
              length(X2)
              
              
              
              # Parametros STAN ---------------------------------------------------------
              
              
              N <- length(peru$GID_1) #numero de regiones 
              I <- diag(N)
              W <- nb2mat(peru.nb) # Matriz de regiones contiguas( matriz de pesos espaciales)
              e <- rep(1,N)
              k <- dim(model_fin_1)[2] #por que quiero el numero de columnas 
              ##############################################################################
              
              #check missing data and others 
              
              
              
              #cohercecion de la base de datos()
              medicos_list_mix  <- within(list(), {
                N<-N
                I<-I
                W<-W
                e<-e
                N1 <- nrow(data_med[data_med$lf==Inf,]) # split database (lf)
                N2 <- nrow(data_med[data_med$lf<Inf,])
                X1 <- (data_med$sexo[data_med$lf==Inf]=="Mujer") #dummy genero
                X2 <- (data_med$sexo[data_med$lf<Inf]=="Mujer")
                z<-model_fin_1
                u<-model_fin_2
                K <- dim(model_fin_1)[2] # es necesario introducir la dimension xq dentro del stancode lo llamo 
                S<-K+3+1 #+2 para los indices por depa 
                #X3 <- (data_med$timew[data_med$lf==Inf])           #Time working 
                #X4 <- (data_med$timew[data_med$lf<Inf]) # (key sensitive)
                #X5 <- (data_med$espec[data_med$lf==Inf]=="Sí")
                #X6 <- (data_med$espec[data_med$lf<Inf]=="Sí") #2110 y 4376
                
                
                d1  <- c(data_med$departamento[data_med$lf==Inf])
                d2  <- c(data_med$departamento[data_med$lf<Inf])
                L1 <- data_med$li[data_med$lf==Inf]/1000 # ingreso prob.
                L2 <- data_med$li[data_med$lf<Inf]/1000
                U2 <- data_med$lf[data_med$lf<Inf]/1000
              })
              
              
              medicos_list_mix_2  <- within(list(), {
                      N<-N
                      I<-I
                      W<-W
                      e<-e
                      N1 <- nrow(data_med[data_med$lf==Inf,]) # split database (lf)
                      N2 <- nrow(data_med[data_med$lf<Inf,])
                      X1 <- (data_med$sexo[data_med$lf==Inf]=="Mujer") #dummy genero
                      X2 <- (data_med$sexo[data_med$lf<Inf]=="Mujer")
                      z<-model_fin_1
                      u<-model_fin_2
                      K <- dim(model_fin_1)[2] # es necesario introducir la dimension xq dentro del stancode lo llamo 
                      #S<-K+3+1 #+2 para los indices por depa 
                      #X3 <- (data_med$timew[data_med$lf==Inf])           #Time working 
                      #X4 <- (data_med$timew[data_med$lf<Inf]) # (key sensitive)
                      #X5 <- (data_med$espec[data_med$lf==Inf]=="Sí")
                      #X6 <- (data_med$espec[data_med$lf<Inf]=="Sí") #2110 y 4376
                      
                      
                      d1  <- c(data_med$departamento[data_med$lf==Inf])
                      d2  <- c(data_med$departamento[data_med$lf<Inf])
                      L1 <- data_med$li[data_med$lf==Inf]/1000 # ingreso prob.
                      L2 <- data_med$li[data_med$lf<Inf]/1000
                      U2 <- data_med$lf[data_med$lf<Inf]/1000
              })
              
              medicos_list_mix_3_D  <- within(list(), {
                      ncas<-nrow(data_med)
                      n<-N
                      N<-N
                      I<-I
                      W<-W
                      e<-e
                      N1 <- nrow(data_med[data_med$lf==Inf,]) # split database (lf)
                      N2 <- nrow(data_med[data_med$lf<Inf,])
                      X1 <- (data_med$sexo[data_med$lf==Inf]=="Mujer") #dummy genero
                      X2 <- (data_med$sexo[data_med$lf<Inf]=="Mujer")
                      z<-model_fin_1
                      u<-model_fin_2
                      K <- dim(model_fin_1)[2] # es necesario introducir la dimension xq dentro del stancode lo llamo 
                      #S<-K+3+1 #+2 para los indices por depa 
                      #X3 <- (data_med$timew[data_med$lf==Inf])           #Time working 
                      #X4 <- (data_med$timew[data_med$lf<Inf]) # (key sensitive)
                      #X5 <- (data_med$espec[data_med$lf==Inf]=="Sí")
                      #X6 <- (data_med$espec[data_med$lf<Inf]=="Sí") #2110 y 4376
                      
                      
                      d1  <- c(data_med$departamento[data_med$lf==Inf])
                      d2  <- c(data_med$departamento[data_med$lf<Inf])
                      L1 <- data_med$li[data_med$lf==Inf]/1000 # ingreso prob.
                      L2 <- data_med$li[data_med$lf<Inf]/1000
                      U2 <- data_med$lf[data_med$lf<Inf]/1000
              })
              
              
              #########
              #STAN####
              #########
              ##Espec , en el modelo 4 es dependiente  + region aca solo espec sin region
              scode_end_interacciones_med_3 <- "
              data {
              int N;
              matrix<lower=0>[N,N] W;
              vector<upper=1>[N] e;
              matrix<lower=0,upper=1>[N,N] I;
              int N1;
              int N2;
              vector[N1] L1;  // lower bounds for right censored observations
              vector[N2] L2;  // lower bounds for interval censored obervations
              vector[N2] U2;  // upper bounds for interval censored obervations
              
              int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
              int<lower=0> S;   //predictors + interactios
              matrix[N1, K] z;   // predictor matrix  lower bounds
              matrix[N2, K] u;   // predictor matrix  upper bounds
              
              // Covariables
              vector[N1] X1; // X for right censored observations
              vector[N2] X2; // X for interval censored observations
              // vector[N1] X3;  X for right censored observations time working
              //   vector[N2] X4;  X for interval censored observations time working
              
              
              int<lower = 1, upper = 26> d1[N1];
              int<lower = 1, upper = 26> d2[N2];
              }
              parameters {
              vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
              vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
              vector[26] alpha; // Spatially Varying intercept
              vector[26] beta;  // Spatially Varying gap
              vector[26] niuno;
              vector[26] nidos;
              vector[26] nitres; 
              
              vector[S] omega;   // coefficients for predictors 
              
              real<lower=0> phi;
              real mbeta;
              real malpha;
              real<lower = 0> sigma_b;
              real<lower = 0> sigma_a;
              real<lower=-1,upper=1> lambda_b;
              real<lower=-1,upper=1> lambda_a;
              }
              transformed parameters {
              vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
              vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
              vector[N1] mu1;
              vector[N2] mu2;
              for(i in 1:N1){
              mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+niuno[d1[i]]*X1[i]*z[i, 1] + nidos[d1[i]]*X1[i]*z[i, 2]+nitres[d1[i]]*X1[i]*z[i, 3]+omega[7]*z[i, 4]+omega[8]*z[i, 5]+omega[9]*z[i, 6]+omega[10]*z[i, 7]+omega[11]*z[i, 8]+omega[12]*z[i, 9]+omega[13]*z[i, 10]+omega[14]*z[i, 11]+omega[15]*z[i, 12]+omega[16]*z[i, 13]+omega[17]*z[i, 14]+omega[18]*z[i, 15]+omega[19]*z[i, 16]+omega[20]*z[i, 17]+omega[21]*X1[i]*z[i, 4]);
              }
              for(j in 1:N2){
              mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+niuno[d2[j]]*X2[j]*u[j, 1] + nidos[d2[j]]*X2[j]*u[j, 2]+nitres[d2[j]]*X2[j]*u[j, 3]+omega[7]*u[j, 4]+omega[8]*u[j, 5]+omega[9]*u[j, 6]+omega[10]*u[j, 7]+omega[11]*u[j, 8]+omega[12]*u[j, 9]+omega[13]*u[j, 10]+omega[14]*u[j, 11]+omega[15]*u[j, 12]+omega[16]*u[j, 13]+omega[17]*u[j, 14]+omega[18]*u[j, 15]+omega[19]*u[j, 16]+omega[20]*u[j, 17]+omega[21]*X2[j]*u[j, 4]);
              }
              }  
              model {
              for(i in 1:N1){
              y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
              }
              for(j in 1:N2){
              y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
              }
              beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
              alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
              }
              "
              
              scode_end_interacciones_med_4 <- "
              data {
              int N;
              matrix<lower=0>[N,N] W;
              vector<upper=1>[N] e;
              matrix<lower=0,upper=1>[N,N] I;
              int N1;
              int N2;
              vector[N1] L1;  // lower bounds for right censored observations
              vector[N2] L2;  // lower bounds for interval censored obervations
              vector[N2] U2;  // upper bounds for interval censored obervations
              
              int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
              int<lower=0> S;   //predictors + interactios
              matrix[N1, K] z;   // predictor matrix  lower bounds
              matrix[N2, K] u;   // predictor matrix  upper bounds
              
              // Covariables
              vector[N1] X1; // X for right censored observations
              vector[N2] X2; // X for interval censored observations
              // vector[N1] X3;  X for right censored observations time working
              //   vector[N2] X4;  X for interval censored observations time working
              
              
              int<lower = 1, upper = 26> d1[N1];
              int<lower = 1, upper = 26> d2[N2];
              }
              parameters {
              vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
              vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
              vector[26] alpha; // Spatially Varying intercept
              vector[26] beta;  // Spatially Varying gap
              vector[26] niuno;
              vector[26] nidos;
              vector[26] nitres; 

              vector[S] omega;   // coefficients for predictors 
              
              real<lower=0> phi;
              real mbeta;
              real malpha;
              real<lower = 0> sigma_b;
              real<lower = 0> sigma_a;
              real<lower=-1,upper=1> lambda_b;
              real<lower=-1,upper=1> lambda_a;
              }
              transformed parameters {
              vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
              vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
              vector[N1] mu1;
              vector[N2] mu2;
              for(i in 1:N1){
              mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+niuno[d1[i]]*X1[i]*z[i, 1] + nidos[d1[i]]*X1[i]*z[i, 2]+nitres[d1[i]]*X1[i]*z[i, 3]+omega[7]*z[i, 4]+omega[8]*z[i, 5]+omega[9]*z[i, 6]+omega[10]*z[i, 7]+omega[11]*z[i, 8]+omega[12]*z[i, 9]+omega[13]*z[i, 10]+omega[14]*z[i, 11]+omega[15]*z[i, 12]+omega[16]*z[i, 13]+omega[17]*z[i, 14]+omega[18]*z[i, 15]+omega[19]*z[i, 16]+omega[20]*z[i, 17]+omega[21]*X1[i]*z[i, 5]+omega[4]*z[i, 18]+omega[5]*z[i, 19]+omega[6]*z[i, 20]);
              }
              for(j in 1:N2){
              mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+niuno[d2[j]]*X2[j]*u[j, 1] + nidos[d2[j]]*X2[j]*u[j, 2]+nitres[d2[j]]*X2[j]*u[j, 3]+omega[7]*u[j, 4]+omega[8]*u[j, 5]+omega[9]*u[j, 6]+omega[10]*u[j, 7]+omega[11]*u[j, 8]+omega[12]*u[j, 9]+omega[13]*u[j, 10]+omega[14]*u[j, 11]+omega[15]*u[j, 12]+omega[16]*u[j, 13]+omega[17]*u[j, 14]+omega[18]*u[j, 15]+omega[19]*u[j, 16]+omega[20]*u[j, 17]+omega[21]*X2[j]*u[j, 5]+omega[4]*u[j, 18]+omega[5]*u[j, 19]+omega[6]*u[j, 20]);
              }
              }  
              model {
              for(i in 1:N1){
              y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
              }
              for(j in 1:N2){
              y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
              }
              beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
              alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
              }
              "
              ###
              scode_end_interacciones_med_4_ver2 <- "
              data {
              int N;
              matrix<lower=0>[N,N] W;
              vector<upper=1>[N] e;
              matrix<lower=0,upper=1>[N,N] I;
              int N1;
              int N2;
              vector[N1] L1;  // lower bounds for right censored observations
              vector[N2] L2;  // lower bounds for interval censored obervations
              vector[N2] U2;  // upper bounds for interval censored obervations
              
              int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
                 //predictors + interactios   int<lower=0> S;
              matrix[N1, K] z;   // predictor matrix  lower bounds
              matrix[N2, K] u;   // predictor matrix  upper bounds
              
              // Covariables
              vector[N1] X1; // X for right censored observations
              vector[N2] X2; // X for interval censored observations
              // vector[N1] X3;  X for right censored observations time working
              //   vector[N2] X4;  X for interval censored observations time working
              
              
              int<lower = 1, upper = 26> d1[N1];
              int<lower = 1, upper = 26> d2[N2];
              }
              parameters {
              vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
              vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
              vector[26] alpha; // Spatially Varying intercept
              vector[26] beta;  // Spatially Varying gap
              vector[26] niuno;
              vector[26] nidos;
              vector[26] nitres;
              vector[26] nicuatro;

              vector[K] omega;   // coefficients for predictors 
              
              real<lower=0> phi;
              real mbeta;
              real malpha;
              real<lower = 0> sigma_b;
              real<lower = 0> sigma_a;
              real<lower=-1,upper=1> lambda_b;
              real<lower=-1,upper=1> lambda_a;
              }
              transformed parameters {
              vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
              vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
              vector[N1] mu1;
              vector[N2] mu2;
              for(i in 1:N1){
              mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+niuno[d1[i]]*X1[i]*z[i, 1] + nidos[d1[i]]*X1[i]*z[i, 2]+nitres[d1[i]]*X1[i]*z[i, 3]+omega[4]*z[i, 4]+omega[5]*z[i, 5]+omega[6]*z[i, 6]+omega[7]*z[i, 7]+omega[8]*z[i, 8]+omega[9]*z[i, 9]+omega[10]*z[i, 10]+omega[11]*z[i, 11]+omega[12]*z[i, 12]+omega[13]*z[i, 13]+omega[14]*z[i, 14]+omega[15]*z[i, 15]+omega[16]*z[i, 16]+omega[17]*z[i, 17]+nicuatro[d1[i]]*X1[i]*z[i, 5]+omega[18]*z[i, 18]+omega[19]*z[i, 19]+omega[20]*z[i, 20]);
              }
              for(j in 1:N2){
              mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+niuno[d2[j]]*X2[j]*u[j, 1] + nidos[d2[j]]*X2[j]*u[j, 2]+nitres[d2[j]]*X2[j]*u[j, 3]+omega[4]*u[j, 4]+omega[5]*u[j, 5]+omega[6]*u[j, 6]+omega[10]*u[j, 7]+omega[8]*u[j, 8]+omega[9]*u[j, 9]+omega[10]*u[j, 10]+omega[11]*u[j, 11]+omega[12]*u[j, 12]+omega[13]*u[j, 13]+omega[14]*u[j, 14]+omega[15]*u[j, 15]+omega[16]*u[j, 16]+omega[17]*u[j, 17]+nicuatro[d2[j]]*X2[j]*u[j, 5]+omega[18]*u[j, 18]+omega[19]*u[j, 19]+omega[20]*u[j, 20]);
              }
              }  
              model {
              for(i in 1:N1){
              y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
              }
              for(j in 1:N2){
              y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
              }
              beta ~ multi_normal_prec(e * mbeta,  tcrossprod(I - lambda_b * W)/(sigma_b*sigma_b));
              alpha ~ multi_normal_prec(e * malpha,  tcrossprod(I - lambda_a * W)/(sigma_a*sigma_a));
              }
              "
              
              
              scode_end_interacciones_med_5_D_lik <- "
              data {
              int N;
              int n;
              int ncas;
              matrix<lower=0>[N,N] W;
              vector<upper=1>[N] e;
              matrix<lower=0,upper=1>[N,N] I;
              int N1;
              int N2;
              vector[N1] L1;  // lower bounds for right censored observations
              vector[N2] L2;  // lower bounds for interval censored obervations
              vector[N2] U2;  // upper bounds for interval censored obervations
              
              int<lower=0> K;   // number of predictors others apart from gender ,si reconoce mayusculas y minusculas 
                 //predictors + interactios   int<lower=0> S;
              matrix[N1, K] z;   // predictor matrix  lower bounds
              matrix[N2, K] u;   // predictor matrix  upper bounds
              
              // Covariables
              vector[N1] X1; // X for right censored observations
              vector[N2] X2; // X for interval censored observations
              // vector[N1] X3;  X for right censored observations time working
              //   vector[N2] X4;  X for interval censored observations time working
              
              
              int<lower = 1, upper = 26> d1[N1];
              int<lower = 1, upper = 26> d2[N2];
              }
              
              transformed data{
                  vector[n] zeros;
                  matrix<lower = 0>[n, n] D;
                  {
                    vector[n] W_rowsums;
                    for (i in 1:n) {
                      W_rowsums[i] = sum(W[i, ]);
                    }
                    D = diag_matrix(W_rowsums);
                  }
                  zeros = rep_vector(0, n);

              }
              
              
              
              parameters {
              vector<lower=0>[N1] y_raw1;              // latent variable for right censored observations
              vector<lower=0, upper=1>[N2] y_raw2;     // latent variable for interval censored observations
              vector[26] alpha; // Spatially Varying intercept
              vector[26] beta;  // Spatially Varying gap
              vector[26] niuno;
              vector[26] nidos;
              vector[26] nitres;
              vector[26] nicuatro;

              vector[K] omega;   // coefficients for predictors 
              
              real<lower=0> phi;
              real mbeta;
              real malpha;
              real<lower = 0> sigma_b;
              real<lower = 0> sigma_a;
              real<lower=-1,upper=1> lambda_b;
              real<lower=-1,upper=1> lambda_a;
              }
              transformed parameters {
              vector[N1] y1 = L1 + y_raw1;                   // latent variable with lower bound L1
              vector[N2] y2 = L2 + (U2 - L2) .* y_raw2;      // latent variable with lower bound L2 and upper bound U2
              vector[N1] mu1;
              vector[N2] mu2;
              for(i in 1:N1){
              mu1[i] = exp(alpha[d1[i]]+beta[d1[i]]*X1[i]+ omega[1]*z[i, 1] +omega[2]*z[i, 2]+omega[3]*z[i, 3]+niuno[d1[i]]*X1[i]*z[i, 1] + nidos[d1[i]]*X1[i]*z[i, 2]+nitres[d1[i]]*X1[i]*z[i, 3]+omega[4]*z[i, 4]+omega[5]*z[i, 5]+omega[6]*z[i, 6]+omega[7]*z[i, 7]+omega[8]*z[i, 8]+omega[9]*z[i, 9]+omega[10]*z[i, 10]+omega[11]*z[i, 11]+omega[12]*z[i, 12]+omega[13]*z[i, 13]+omega[14]*z[i, 14]+omega[15]*z[i, 15]+omega[16]*z[i, 16]+omega[17]*z[i, 17]+nicuatro[d1[i]]*X1[i]*z[i, 5]);
              }
              for(j in 1:N2){
              mu2[j] = exp(alpha[d2[j]]+beta[d2[j]]*X2[j]+ omega[1]*u[j, 1] +omega[2]*u[j, 2]+omega[3]*u[j, 3]+niuno[d2[j]]*X2[j]*u[j, 1] + nidos[d2[j]]*X2[j]*u[j, 2]+nitres[d2[j]]*X2[j]*u[j, 3]+omega[4]*u[j, 4]+omega[5]*u[j, 5]+omega[6]*u[j, 6]+omega[10]*u[j, 7]+omega[8]*u[j, 8]+omega[9]*u[j, 9]+omega[10]*u[j, 10]+omega[11]*u[j, 11]+omega[12]*u[j, 12]+omega[13]*u[j, 13]+omega[14]*u[j, 14]+omega[15]*u[j, 15]+omega[16]*u[j, 16]+omega[17]*u[j, 17]+nicuatro[d2[j]]*X2[j]*u[j, 5]);
              }
              }  
              model {
              for(i in 1:N1){
              y1[i] ~ gamma(mu1[i]*mu1[i] * phi , mu1[i] * phi);
              }
              for(j in 1:N2){
              y2[j] ~ gamma(mu2[j]*mu2[j] * phi , mu2[j] * phi);
              }
              beta ~ multi_normal_prec(e * mbeta,  (I - lambda_b * inverse(D)* W) *D/(sigma_b*sigma_b) * (I - lambda_b * inverse(D)* W)'  );
              alpha ~ multi_normal_prec(e * malpha,  (I - lambda_a * inverse(D) * W)*D/(sigma_a*sigma_a)*(I - lambda_a * inverse(D)* W)' );
              }
              
              generated quantities {
                  vector[N1] log_lik_1;
                  vector[N2] log_lik_2;
                  vector[ncas] log_lik;
                  
                  for (s in 1:N1) {
                    log_lik_1[s] = gamma_lpdf(y1[s] | mu1[s]*mu1[s] * phi, mu1[s] * phi);
                  }  
                  for (i in 1:N2) {
                    log_lik_2[i] = gamma_lpdf(y2[i] | mu2[i]*mu2[i] * phi, mu2[i] * phi);
                  } 
                  
                  log_lik = append_row(log_lik_1, log_lik_2);
                }
              
              
              
              "
              
              
              
              
              
               
                        inits_fin_med <- function() {
                          list(beta=rep(0,26),mbeta=0,sigma_b=1,lambda_b=0, 
                               alpha=rep(2,26),malpha=2,sigma_a=1,lambda_a=0,
                               phi = 0.6,
                               y_raw1 = rep(0.5,medicos_list_mix$N1), # init values > 0
                               y_raw2 = rep(0.5,medicos_list_mix$N2)  # 0 < init values < 1
                          )
                        }
                        
              #Creado solo para la 4 ver2
                      inits_fin_med_2 <- function() {
                              list(beta=rep(0,26),mbeta=0,sigma_b=1,lambda_b=0, 
                                   alpha=rep(2,26),malpha=2,sigma_a=1,lambda_a=0,
                                   phi = 0.6,
                                   y_raw1 = rep(0.5,medicos_list_mix_2$N1), # init values > 0
                                   y_raw2 = rep(0.5,medicos_list_mix_2$N2)  # 0 < init values < 1
                              )
                      }
                      
                      inits_fin_med_3 <- function() {
                              list(beta=rep(0,26),mbeta=0,sigma_b=1,lambda_b=0, 
                                   alpha=rep(2,26),malpha=2,sigma_a=1,lambda_a=0,
                                   phi = 0.6,
                                   y_raw1 = rep(0.5,medicos_list_mix_3_D$N1), # init values > 0
                                   y_raw2 = rep(0.5,medicos_list_mix_3_D$N2)  # 0 < init values < 1
                              )
                      }
                      
              
              
              
                        #fit_complete2
                        fit_complete_intect_med <- stan(model_code = scode_end_interacciones_med, data = medicos_list_mix, 
                                                chains = 3, cores = 3,init=inits_fin_med)
              
                        fit_complete_intect_med_2 <- stan(model_code = scode_end_interacciones_med_2, data = medicos_list_mix, 
                                              chains = 3, cores = 3,init=inits_fin_med)
                        fit_complete_intect_med_3 <- stan(model_code = scode_end_interacciones_med_3, data = medicos_list_mix, 
                                                          chains = 3, cores = 3,init=inits_fin_med)
                        
                        fit_complete_intect_med_4 <- stan(model_code = scode_end_interacciones_med_4, data = medicos_list_mix, 
                                                          chains = 3, cores = 3,init=inits_fin_med)
                        
                        fit_complete_intect_med_4_ver2 <- stan(model_code = scode_end_interacciones_med_4_ver2, data = medicos_list_mix_2, 
                                                          chains = 3, cores = 3,init=inits_fin_med_2)
                        
                        fit_complete_intect_med_5_D <- stan(model_code = scode_end_interacciones_med_5_D_lik, data = medicos_list_mix_3_D, 
                                                               chains = 3, cores = 3,init=inits_fin_med_3)
                        
                        
                        
                        saveRDS(fit_complete_intect_med,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_1.rds")
                        saveRDS(fit_complete_intect_med_2,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_2.rds")
                        saveRDS(fit_complete_intect_med_3,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_3_intec_por_region.rds")
                        saveRDS(fit_complete_intect_med_4,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_4_intec_por_region.rds") #me huevee sobreeescribite el 3 sobre el 4 hay que hacer nuevo
                        saveRDS(fit_complete_intect_med_4_ver2,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_4_intec_por_region_ver2.rds") #me huevee sobreeescribite el 3 sobre el 4 hay que hacer nuevo
                        saveRDS(fit_complete_intect_med_5_D,file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_5_intect_D_lik.rds") #me huevee sobreeescribite el 3 sobre el 4 hay que hacer nuevo
                        
                        print(fit_complete_intect_med,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
                        print(fit_complete_intect_med_2,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
                        print(fit_complete_intect_med_3,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
                        print(fit_complete_intect_med_4_ver2,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
                        print(fit_complete_intect_med_5_D,pars=c("beta","alpha","sigma_b","lambda_b","lambda_a","mbeta"))
                        
                        
                        fit_complete_intect_med  <-  readRDS(file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_1.rds")
                        fit_complete_intect_med_2<-  readRDS(file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_2.rds")
                        fit_complete_intect_med_3<-  readRDS(file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_3_intec_por_region.rds")
                        fit_complete_intect_med_4<-  readRDS(file="C:\\Users\\JESÚS\\Desktop\\investigacion mb\\runs\\fit_complete_intect_med_4_intec_por_region.rds")
                        
                        print(fit_complete_intect_med,pars=c("lambda_b","lambda_a"))
                        print(fit_complete_intect_med_2,pars=c("lambda_b","lambda_a"))
                        
                        
                        library("loo")
                        
                        # Extract log-likelihood and compute LOO
                        log_lik1 <- extract_log_lik(fit_complete_intect_med_5_D)
                        loo1 <- loo(log_lik1) 
                        
                        
                        
                        
                        #Grafico  modelo interacciones modelo 1 
                        library(bayesplot)
                        
                        color_scheme_set("blue") #comando aguanta solo uno por uno 
                        fake_rhat_values <- c(1.00, 1.00,1.00,1.00, 1.01, 1.00, 1.01, 1.00,1.00,1.01,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.00,1.00)
                        mcmc_intervals(fit_complete_intect_med, 
                                       pars=c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                              "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                              "beta[11]","beta[12]","beta[13]","beta[14]","beta[15]",
                                              "beta[16]","beta[17]","beta[18]","beta[19]","beta[20]",
                                              "beta[21]","beta[22]","beta[23]","beta[24]","beta[25]","beta[26]"), rhat = fake_rhat_values)+
                          ggplot2::labs(
                            title = "Intervalos de Credibilidad",
                            subtitle = "Al 95 %"
                          )+  coord_flip()+ #Esto ultimo falta editar
                          xaxis_text(angle = -70, vjust = 1, hjust = 0)+
                          theme(
                            plot.title = element_text(hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5))#,axis.text.x = element_text(angle = 90))
                        
                        
                        
                        
                        #Grafico  modelo interacciones modelo 2
                        library(bayesplot)
                        
                        color_scheme_set("blue") #comando aguanta solo uno por uno 
                        fake_rhat_values <- c(1.00, 1.00,1.00,1.00, 1.01, 1.00, 1.01, 1.00,1.00,1.01,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.00,1.00)
                        mcmc_intervals(fit_complete_intect_med_2, 
                                       pars=c("beta[1]","beta[2]","beta[3]","beta[4]","beta[5]",
                                              "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                              "beta[11]","beta[12]","beta[13]","beta[14]","beta[15]",
                                              "beta[16]","beta[17]","beta[18]","beta[19]","beta[20]",
                                              "beta[21]","beta[22]","beta[23]","beta[24]","beta[25]","beta[26]"), rhat = fake_rhat_values)+
                          ggplot2::labs(
                            title = "Intervalos de Credibilidad",
                            subtitle = "Al 95 %"
                          )+  coord_flip()+ #Esto ultimo falta editar
                          xaxis_text(angle = -70, vjust = 1, hjust = 0)+
                          theme(
                            plot.title = element_text(hjust = 0.5),
                            plot.subtitle = element_text(hjust = 0.5))#,axis.text.x = element_text(angle = 90))
                        
                        
                        
                        ##########
                        
                        
                        color_scheme_set("green") #comando aguanta solo uno por uno 
                        fake_rhat_values <- c(1.00, 1.00,1.00,1.00, 1.01, 1.00, 1.01, 1.00,1.00,1.01,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.01,1.00,1.00,1.00,1.00,1.00)
                        mcmc_intervals(fit_complete_intect_med_3, 
                                       pars=c("beta[16]","beta[14]","beta[13]","beta[1]","beta[3]","beta[4]","beta[5]",
                                              "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
                                              "beta[11]","beta[15]",
                                              "beta[17]","beta[18]","beta[19]","beta[20]",
                                              "beta[21]","beta[22]","beta[23]","beta[25]","beta[26]","beta[12]","beta[24]","beta[2]"), rhat = fake_rhat_values)+
                                ggplot2::labs(
                                        title = "Intervalos de Credibilidad",
                                        subtitle = "Al 95 %"
                                )+  coord_flip()+ #Esto ultimo falta editar
                                xaxis_text(angle = -70, vjust = 1, hjust = 0)+
                                theme(
                                        plot.title = element_text(hjust = 0.5),
                                        plot.subtitle = element_text(hjust = 0.5))#,axis.text.x = element_text(angle = 90))
                        