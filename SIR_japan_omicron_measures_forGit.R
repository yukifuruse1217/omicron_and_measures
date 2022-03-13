require(deSolve) # for the "ode" function
library(dplyr)



#Equation#######################




sir_1 <- function(beta_eTOc_ini, beta_eTOy2_ini, beta_eTOy4_ini, beta_eTOe_ini,
                  beta_y4TOc_ini, beta_y4TOy2_ini, beta_y4TOy4_ini, beta_y4TOe_ini,
                  beta_y2TOc_ini, beta_y2TOy2_ini, beta_y2TOy4_ini, beta_y2TOe_ini,
                  beta_cTOc_ini, beta_cTOy2_ini, beta_cTOy4_ini, beta_cTOe_ini,
                  

                  alpha_ini, ee_ini, gamma,
                  
                  beta_coeff_adult, beta_coeff_child,
                  alpha_coeff, ee_coeff, vac_move_a_ini, vac_move_c_ini,
                  
                  ne, n4, n2, nc,
                  
                  Se0, Ie0, De0, Re0, Ve0,
                  Sy40, Iy40, Dy40, Ry40, Vy40,
                  Sy20, Iy20, Dy20, Ry20, Vy20,
                  Sc0, Ic0, Dc0, Rc0, Vc0,
                  

                  times) {
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {

    if (time < 30) {
       beta_eTOc <- beta_eTOc_ini / nc
       beta_eTOy2 <- beta_eTOy2_ini / n2
       beta_eTOy4 <- beta_eTOy4_ini / n4
       beta_eTOe <- beta_eTOe_ini / ne
       
       beta_y4TOc <- beta_y4TOc_ini / nc
       beta_y4TOy2 <- beta_y4TOy2_ini / n2
       beta_y4TOy4 <- beta_y4TOy4_ini / n4
       beta_y4TOe <- beta_y4TOe_ini / ne
       
       beta_y2TOc <- beta_y2TOc_ini / nc
       beta_y2TOy2 <- beta_y2TOy2_ini / n2
       beta_y2TOy4 <- beta_y2TOy4_ini / n4
       beta_y2TOe <- beta_y2TOe_ini / ne
       
       beta_cTOc <- beta_cTOc_ini / nc
       beta_cTOy2 <- beta_cTOy2_ini / n2
       beta_cTOy4 <- beta_cTOy4_ini / n4
       beta_cTOe <- beta_cTOe_ini / ne
       

       alpha <- alpha_ini
       ee    <- ee_ini
    } else {
      beta_eTOc <- beta_eTOc_ini / nc
      beta_eTOy2 <- beta_eTOy2_ini / n2 * (1 - (1 - beta_coeff_adult) / 2)
      beta_eTOy4 <- beta_eTOy4_ini / n4 * (1 - (1 - beta_coeff_adult) / 2)
      beta_eTOe <- beta_eTOe_ini / ne * (1 - (1 - beta_coeff_adult) / 2)
      
      beta_y4TOc <- beta_y4TOc_ini / nc
      beta_y4TOy2 <- beta_y4TOy2_ini / n2 * beta_coeff_adult
      beta_y4TOy4 <- beta_y4TOy4_ini / n4 * beta_coeff_adult
      beta_y4TOe <- beta_y4TOe_ini / ne * (1 - (1 - beta_coeff_adult) / 2)
      
      beta_y2TOc <- beta_y2TOc_ini / nc
      beta_y2TOy2 <- beta_y2TOy2_ini / n2 * beta_coeff_adult
      beta_y2TOy4 <- beta_y2TOy4_ini / n4 * beta_coeff_adult
      beta_y2TOe <- beta_y2TOe_ini / ne * (1 - (1 - beta_coeff_adult) / 2)
      
      beta_cTOc <- beta_cTOc_ini / nc * beta_coeff_child
      beta_cTOy2 <- beta_cTOy2_ini / n2
      beta_cTOy4 <- beta_cTOy4_ini / n4
      beta_cTOe <- beta_cTOe_ini / ne
       
       
      alpha <- alpha_ini * alpha_coeff
      ee <- ee_ini + ee_coeff
     }

     if (time >= 30 && time < 31) {
          dSe  <- -beta_cTOe  * (Ic + Dc) * Se  -beta_y2TOe  * (Iy2 + Dy2) * Se  -beta_y4TOe  * (Iy4 + Dy4) * Se  -beta_eTOe  * (Ie + De) * Se      - Se  * vac_move_a_ini
          dSy4 <- -beta_cTOy4 * (Ic + Dc) * Sy4 -beta_y2TOy4 * (Iy2 + Dy2) * Sy4 -beta_y4TOy4 * (Iy4 + Dy4) * Sy4 -beta_eTOy4 * (Ie + De) * Sy4     - Sy4 * vac_move_a_ini
          dSy2 <- -beta_cTOy2 * (Ic + Dc) * Sy2 -beta_y2TOy2 * (Iy2 + Dy2) * Sy2 -beta_y4TOy2 * (Iy4 + Dy4) * Sy2 -beta_eTOy2 * (Ie + De) * Sy2     - Sy2 * vac_move_a_ini
          dSc  <- -beta_cTOc  * (Ic + Dc) * Sc  -beta_y2TOc  * (Iy2 + Dy2) * Sc  -beta_y4TOc  * (Iy4 + Dy4) * Sc  -beta_eTOc  * (Ie + De) * Sc      - Sc  * vac_move_c_ini

          dVe <-  +Se * vac_move_a_ini
          dVy4 <- +Sy4 * vac_move_a_ini
          dVy2 <- +Sy2 * vac_move_a_ini
          dVc <-  +Sc * vac_move_c_ini
     
     } else {
       dSe  <- -beta_cTOe  * (Ic + Dc) * Se  -beta_y2TOe  * (Iy2 + Dy2) * Se  -beta_y4TOe  * (Iy4 + Dy4) * Se  -beta_eTOe  * (Ie + De) * Se
       dSy4 <- -beta_cTOy4 * (Ic + Dc) * Sy4 -beta_y2TOy4 * (Iy2 + Dy2) * Sy4 -beta_y4TOy4 * (Iy4 + Dy4) * Sy4 -beta_eTOy4 * (Ie + De) * Sy4
       dSy2 <- -beta_cTOy2 * (Ic + Dc) * Sy2 -beta_y2TOy2 * (Iy2 + Dy2) * Sy2 -beta_y4TOy2 * (Iy4 + Dy4) * Sy2 -beta_eTOy2 * (Ie + De) * Sy2
       dSc  <- -beta_cTOc  * (Ic + Dc) * Sc  -beta_y2TOc  * (Iy2 + Dy2) * Sc  -beta_y4TOc  * (Iy4 + Dy4) * Sc  -beta_eTOc  * (Ie + De) * Sc
       
       dVe <- 0
       dVy4 <- 0
       dVy2 <- 0
       dVc <- 0
     }
      dIe  <- (1 -ee) * (+beta_cTOe  * (Ic + Dc) * Se  +beta_y2TOe  * (Iy2 + Dy2) * Se  +beta_y2TOe  * (Iy4 + Dy4) * Se  +beta_eTOe  * (Ie + De) * Se)   -gamma * Ie
      dIy4 <- (1 -ee) * (+beta_cTOy4 * (Ic + Dc) * Sy4 +beta_y2TOy4 * (Iy2 + Dy2) * Sy4 +beta_y4TOy4 * (Iy4 + Dy4) * Sy4 +beta_eTOy4 * (Ie + De) * Sy4)  -gamma * Iy4
      dIy2 <- (1 -ee) * (+beta_cTOy2 * (Ic + Dc) * Sy2 +beta_y2TOy2 * (Iy2 + Dy2) * Sy2 +beta_y4TOy2 * (Iy4 + Dy4) * Sy2 +beta_eTOy2 * (Ie + De) * Sy2)  -gamma * Iy2
      dIc  <- (1 -ee) * (+beta_cTOc  * (Ic + Dc) * Sc  +beta_y2TOc  * (Iy2 + Dy2) * Sc  +beta_y4TOc  * (Iy4 + Dy4) * Sc  +beta_eTOc  * (Ie + De) * Sc)   -gamma * Ic
      
      dDe  <- ee * (+beta_cTOe  * (Ic + Dc) * Se  +beta_y2TOe  * (Iy2 + Dy2) * Se  +beta_y2TOe  * (Iy4 + Dy4) * Se  +beta_eTOe  * (Ie + De) * Se)   -alpha * De
      dDy4 <- ee * (+beta_cTOy4 * (Ic + Dc) * Sy4 +beta_y2TOy4 * (Iy2 + Dy2) * Sy4 +beta_y4TOy4 * (Iy4 + Dy4) * Sy4 +beta_eTOy4 * (Ie + De) * Sy4)  -alpha * Dy4
      dDy2 <- ee * (+beta_cTOy2 * (Ic + Dc) * Sy2 +beta_y2TOy2 * (Iy2 + Dy2) * Sy2 +beta_y4TOy2 * (Iy4 + Dy4) * Sy2 +beta_eTOy2 * (Ie + De) * Sy2)  -alpha * Dy2
      dDc  <- ee * (+beta_cTOc  * (Ic + Dc) * Sc  +beta_y2TOc  * (Iy2 + Dy2) * Sc  +beta_y4TOc  * (Iy4 + Dy4) * Sc  +beta_eTOc  * (Ie + De) * Sc)   -alpha * Dc
      
      dRe  <- +gamma * Ie  + alpha * De
      dRy4 <- +gamma * Iy4 + alpha * Dy4
      dRy2 <- +gamma * Iy2 + alpha * Dy2
      dRc  <- +gamma * Ic  + alpha * Dc
      
      
     return(list(c(dSe, dIe, dDe, dRe, dVe,
                   dSy4, dIy4, dDy4, dRy4, dVy4,
                   dSy2, dIy2, dDy2, dRy2, dVy2,
                   dSc, dIc, dDc, dRc, dVc)))
     })
  }
  
  # the parameters values:
  parameters_values <- c(beta_eTOc_ini=beta_eTOc_ini, beta_eTOy2_ini=beta_eTOy2_ini, beta_eTOy4_ini=beta_eTOy4_ini, beta_eTOe_ini=beta_eTOe_ini,
                         beta_y4TOc_ini=beta_y4TOc_ini, beta_y4TOy2_ini=beta_y4TOy2_ini, beta_y4TOy4_ini=beta_y4TOy4_ini, beta_y4TOe_ini=beta_y4TOe_ini,
                         beta_y2TOc_ini=beta_y2TOc_ini, beta_y2TOy2_ini=beta_y2TOy2_ini, beta_y2TOy4_ini=beta_y2TOy4_ini, beta_y2TOe_ini=beta_y2TOe_ini,
                         beta_cTOc_ini=beta_cTOc_ini, beta_cTOy2_ini=beta_cTOy2_ini, beta_cTOy4_ini=beta_cTOy4_ini, beta_cTOe_ini=beta_cTOe_ini,
                         

                         alpha_ini=alpha_ini, ee_ini=ee_ini, gamma=gamma,
                         
                         beta_coeff_adult=beta_coeff_adult, beta_coeff_child=beta_coeff_child,
                         alpha_coeff=alpha_coeff, ee_coeff=ee_coeff, vac_move_a_ini=vac_move_a_ini, vac_move_c_ini=vac_move_c_ini,
                         
                         ne=ne, n4=n4, n2=n2, nc=nc)
  
  # the initial values of variables:
  initial_values <- c(Se=Se0, Ie=Ie0, De=De0, Re=Re0, Ve=Ve0,
                      Sy4=Sy40, Iy4=Iy40, Dy4=Dy40, Ry4=Ry40, Vy4=Vy40,
                      Sy2=Sy20, Iy2=Iy20, Dy2=Dy20, Ry2=Ry20, Vy2=Vy20,
                      Sc=Sc0, Ic=Ic0, Dc=Dc0, Rc=Rc0, Vc=Vc0)
  
  # solving
  out <- ode(method="rk4",initial_values, times, sir_equations, parameters_values)
  #out <- lsode(initial_values, times, sir_equations, parameters_values,atol=1e-20,maxsteps=10000)
  
  #method = c("lsoda", "lsode", "lsodes", "lsodar", "vode", "daspk",
  #           "euler", "rk4", "ode23", "ode45", "radau", 
  #           "bdf", "bdf_d", "adams", "impAdams", "impAdams_d", "iteration")
  # rk4

    
  # returning the output:
  simulationresult <- as.data.frame(out)
  return(simulationresult)
}



#Simulations
 




timeL <- 60

gamma_delta <- 1/5
gamma_omicron <- 1/3



basic_beta_delta <- 0.138    *gamma_delta

basic_beta_child <- 0.209    *gamma_delta

basic_beta_short <- 0.115    *gamma_omicron

basic_beta_omicron <- 0.175  *gamma_omicron





list_beta_eTOc_ini  <- c(basic_beta_delta,     basic_beta_child*1.2,   basic_beta_short,   basic_beta_omicron*1.2)
list_beta_eTOy2_ini <- c(basic_beta_delta,     basic_beta_child,       basic_beta_short,   basic_beta_omicron)
list_beta_eTOy4_ini <- c(basic_beta_delta,     basic_beta_child,       basic_beta_short,   basic_beta_omicron)
list_beta_eTOe_ini  <- c(basic_beta_delta*4,   basic_beta_child*4,     basic_beta_short*4, basic_beta_omicron*4)
                                                  
list_beta_y4TOc_ini  <- c(basic_beta_delta*2,  basic_beta_child*1.2,   basic_beta_short*2,  basic_beta_omicron*1.2)
list_beta_y4TOy2_ini <- c(basic_beta_delta*2,  basic_beta_child,       basic_beta_short*2,  basic_beta_omicron)
list_beta_y4TOy4_ini <- c(basic_beta_delta*8,  basic_beta_child*4,     basic_beta_short*8, basic_beta_omicron*4)
list_beta_y4TOe_ini  <- c(basic_beta_delta*2,  basic_beta_child,       basic_beta_short*4,  basic_beta_omicron)

list_beta_y2TOc_ini  <- c(basic_beta_delta*2,  basic_beta_child*1.2,   basic_beta_short*2,  basic_beta_omicron*1.2)
list_beta_y2TOy2_ini <- c(basic_beta_delta*8,  basic_beta_child*4,     basic_beta_short*8, basic_beta_omicron*4)
list_beta_y2TOy4_ini <- c(basic_beta_delta*2,  basic_beta_child,       basic_beta_short*2,  basic_beta_omicron)
list_beta_y2TOe_ini  <- c(basic_beta_delta*2,  basic_beta_child,       basic_beta_short*2,  basic_beta_omicron)

list_beta_cTOc_ini   <- c(basic_beta_delta*4,  basic_beta_child*4.8,   basic_beta_short*4, basic_beta_omicron*4.8)
list_beta_cTOy2_ini  <- c(basic_beta_delta,    basic_beta_child,       basic_beta_short,   basic_beta_omicron)
list_beta_cTOy4_ini  <- c(basic_beta_delta,    basic_beta_child,       basic_beta_short,   basic_beta_omicron)
list_beta_cTOe_ini   <- c(basic_beta_delta,    basic_beta_child,       basic_beta_short,   basic_beta_omicron)







list_alpha_ini <- c(1/2,1/2,1/2,1/2)
list_ee_ini <- c(0,0,0,0)
list_gamma <- c(gamma_delta,gamma_delta,gamma_omicron,gamma_omicron)





list_beta_coeff_adult <- c(1,1,1,1, 0.51, 0.44,  1, 0.44, 0.44, 0.51, 0.51)
list_beta_coeff_child <- c(1,1,1,1, 0.51,    1,  1,    1,    1, 0.51, 0.51)

list_alpha_coeff <- c(1,1,1,1,1,1,1,1,1,1,1)
list_ee_coeff <-    c(0,0,0,0,0,0,0.73,0.73,0.73,0.73,0.73)

list_vac_move_a_ini <- c(0.48, 0.58, 0.48/2, 0.58/2, 0, 0, 0, 0.58, 0.58/2, 0.48, 0.48/2)
list_vac_move_c_ini <- c(0.48,    0, 0.48/2,      0, 0, 0, 0,    0,      0, 0.48, 0.48/2)

sim <-0

for (i in 1:4){ #1delta, 2child, 3short, 4omicron
  for (ii in 1:11){
    #1:allvaccine, 2:adult vaccine, 3:half-allvaccine, 4:half-adultvaccine,
    #5:all-homo-Beta, 6:adult-homo-Beta,
    #7:clusterEE
    #8:adultvaccine+adult-homo-Beta+clusterEE
    #9:half-adultvaccine+adult-homo-Beta+clusterEE
    #10:allvaccine+all-homo-Beta+clusterEE
    #11:half-allvaccine+all-homo-Beta+clusterEE
    
    sim <- sim + 1
    
    set_beta_eTOc_ini <- list_beta_eTOc_ini[i]
    set_beta_eTOy2_ini <- list_beta_eTOy2_ini[i]
    set_beta_eTOy4_ini <- list_beta_eTOy4_ini[i]
    set_beta_eTOe_ini <- list_beta_eTOe_ini[i]
    
    set_beta_y4TOc_ini <- list_beta_y4TOc_ini[i]
    set_beta_y4TOy2_ini <- list_beta_y4TOy2_ini[i]
    set_beta_y4TOy4_ini <- list_beta_y4TOy4_ini[i]
    set_beta_y4TOe_ini <- list_beta_y4TOe_ini[i]
    
    set_beta_y2TOc_ini <- list_beta_y2TOc_ini[i]
    set_beta_y2TOy2_ini <- list_beta_y2TOy2_ini[i]
    set_beta_y2TOy4_ini <- list_beta_y2TOy4_ini[i]
    set_beta_y2TOe_ini <- list_beta_y2TOe_ini[i]
    
    set_beta_cTOc_ini <- list_beta_cTOc_ini[i]
    set_beta_cTOy2_ini <- list_beta_cTOy2_ini[i]
    set_beta_cTOy4_ini <- list_beta_cTOy4_ini[i]
    set_beta_cTOe_ini <- list_beta_cTOe_ini[i]

    
    
    
    set_alpha_ini <- list_alpha_ini[i] 
    set_ee_ini <- list_ee_ini[i]
    set_gamma <- list_gamma[i]

    

    
    
    
    
    set_beta_coeff_adult <- list_beta_coeff_adult[ii]
    set_beta_coeff_child <- list_beta_coeff_child[ii]
    
    
    set_alpha_coeff <- list_alpha_coeff[ii]
    set_ee_coeff <- list_ee_coeff[ii]
    
    set_vac_move_a_ini <- list_vac_move_a_ini[ii]
    set_vac_move_c_ini <- list_vac_move_c_ini[ii]



    
    
    
    result <- sir_1(beta_eTOc_ini=set_beta_eTOc_ini, beta_eTOy2_ini=set_beta_eTOy2_ini, beta_eTOy4_ini=set_beta_eTOy4_ini, beta_eTOe_ini=set_beta_eTOe_ini,
                      beta_y4TOc_ini=set_beta_y4TOc_ini, beta_y4TOy2_ini=set_beta_y4TOy2_ini, beta_y4TOy4_ini=set_beta_y4TOy4_ini, beta_y4TOe_ini=set_beta_y4TOe_ini,
                      beta_y2TOc_ini=set_beta_y2TOc_ini, beta_y2TOy2_ini=set_beta_y2TOy2_ini, beta_y2TOy4_ini=set_beta_y2TOy4_ini, beta_y2TOe_ini=set_beta_y2TOe_ini,
                      beta_cTOc_ini=set_beta_cTOc_ini, beta_cTOy2_ini=set_beta_cTOy2_ini, beta_cTOy4_ini=set_beta_cTOy4_ini, beta_cTOe_ini=set_beta_cTOe_ini,
                      
                      
                      alpha_ini=set_alpha_ini, ee_ini=set_ee_ini, gamma=set_gamma,
                      
                      beta_coeff_adult=set_beta_coeff_adult, beta_coeff_child=set_beta_coeff_child,
                      alpha_coeff=set_alpha_coeff, ee_coeff=set_ee_coeff, vac_move_a_ini=set_vac_move_a_ini, vac_move_c_ini=set_vac_move_c_ini,
                      
                      ne=3700*(10^4), n4=3350*(10^4), n2=3250*(10^4), nc=2300*(10^4),
                      
                      Se0=3700*(10^4), Ie0=370*2.08, De0=0, Re0=0, Ve0=0,
                      Sy40=3350*(10^4), Iy40=335*2.08, Dy40=0, Ry40=0, Vy40=0,
                      Sy20=3250*(10^4), Iy20=325*2.08, Dy20=0, Ry20=0, Vy20=0,
                      Sc0=2300*(10^4), Ic0=230*2.08, Dc0=0, Rc0=0, Vc0=0,
                      
                      
                      times = seq(0, timeL, by=1))    
    
    
    

    
    
    
    var_i <- rep(i, timeL+1)
    var_ii <- rep(ii, timeL+1)
    var_sim <- rep(sim, timeL+1)
    
    result <- cbind(result,var_i)
    result <- cbind(result,var_ii)
    result <- cbind(result,var_sim)
    
    temp <- c(0, result$Ie + result$De + result$Re)
    temp <- temp[-(timeL+2)]
    newInfE <- result$Ie + result$De + result$Re - temp
    result <- cbind(result,newInfE)

    temp <- c(0, result$Iy4 + result$Dy4 + result$Ry4)
    temp <- temp[-(timeL+2)]
    newInfY4 <- result$Iy4 + result$Dy4 + result$Ry4 - temp
    result <- cbind(result,newInfY4)

    temp <- c(0, result$Iy2 + result$Dy2 + result$Ry2)
    temp <- temp[-(timeL+2)]
    newInfY2 <- result$Iy2 + result$Dy2 + result$Ry2 - temp
    result <- cbind(result,newInfY2)
    
    temp <- c(0, result$Ic + result$Dc + result$Rc)
    temp <- temp[-(timeL+2)]
    newInfC <- result$Ic + result$Dc + result$Rc - temp
    result <- cbind(result,newInfC)
    
    newInf <- newInfE + newInfY4 + newInfY2 + newInfC
    result <- cbind(result,newInf)
    
    propChild <- newInfC / newInf
    result <- cbind(result,propChild)

        
    if (sim == 1) {
      comb_result <- result
    } else {
      comb_result <- rbind(comb_result,result)
    }
    
    
}}

total <- comb_result$Se + comb_result$Ie + comb_result$De + comb_result$Re + comb_result$Ve +
  comb_result$Sy4 + comb_result$Iy4 + comb_result$Dy4 + comb_result$Ry4 + comb_result$Vy4 +
  comb_result$Sy2 + comb_result$Iy2 + comb_result$Dy2 + comb_result$Ry2 + comb_result$Vy2 +
  comb_result$Sc + comb_result$Ic + comb_result$Dc + comb_result$Rc + comb_result$Vc


comb_result <- cbind(comb_result, total)


