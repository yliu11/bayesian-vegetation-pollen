# Bayesian statistical model estimating pollen productivity and dispersal 
# using pollen and vegetation data from 33 sites in 3 different regions.

vp_model <- function(){ 
  for(s in 1:S){ # for s-th pollen site (lake)
    # Likelihood for pollen counts Y:    
    Y[s,1:T] ~ dmulti(P[s,1:T],N[s]) # Eqn [1]
    # Predicted pollen counts:
    Y_rep[s,1:T] ~ dmulti(P[s,1:T],N[s])
    # Calculating observed and predicted pollen proportions:
    Ep_p[s,1:T] <- Y[s,1:T]/sum(Y[s,])
    Ep_rep[s,1:T] <- Y_rep[s,1:T]/sum(Y_rep[s,])
    
    for(t in 1:T){
      # The relative pollen load:
      P[s,t] <- phi[t,Reg[s]] * (gamma[t]*VL[s,t] + (1-gamma[t])*VR[s,t]) # Eqn [2]
      
      # Local distance weighting (vegetation within the 1-km radius from lake shore):
      # Add-up VL_calc from each "ring" to obtain distance weighted local vegetation compositions (relative abundances):
      VL[s,t] <- sum(VL_calc[s,,t]) # This and the calculation of VL_calc (two lines down) correspond to Eqn [7]
      for(i in 1:I){
        # Distance weighted relative abundance of surveyed vegetation VJCK (data reported in Jackson 2019):
        VL_calc[s,i,t] <- VJCK[((s-1)*I+i),t] * wl[s,i,t] # Part of Eqn [7], see above
        # Normalizing wl_raw[s,i,t] to ensure wl[s,i,t] sum to one across all i:
        wl[s,i,t] <- wl_raw[s,i,t]/cl[s,t] # This line and the line below correspond to Eqn [9]
        # Unnormalized local distance-weights:
        wl_raw[s,i,t] <- exp(b[t]*(R[s]^theta-(DL[i]+R[s])^theta))-exp(b[t]*(R[s]^theta-(DL[i+1]+R[s])^theta))
      }
      # Calculating the normalizing constant cl[s,t] for each site-taxon by summing wl[s,j,t] across all distances:
      cl[s,t] <- sum(wl_raw[s,,t]) # Eqn [11]

      
      # Distance weighting for regional vegetation (located 1-298 km from the lake shore):
      # Add-up VR_calc from each "ring" to obtain distance weighted regional vegetation compositions (relative abundances):
      VR[s,t] <- sum(VR_calc[s,,t])   # This and the calculation of VR_calc (two lines down) below correspond to Eqn [8]
      for(j in 1:J){
        # Distance weighted relative abundance of FIA vegetation:
        VR_calc[s,j,t] <- VFIA[((s-1)*J+j),t] * wr[s,j,t] # Part of Eqn [8], also see above
        # Normalizing wr_raw[s,j,t] to ensure wr[s,j,t] sum to one across all j:
        wr[s,j,t] <- wr_raw[s,j,t]/cr[s,t] # This line and the line below correspond to Eqn [10]
        # Unnormalized regional distance-weights:
        wr_raw[s,j,t] <- exp(b[t]*(R[s]^theta-(DR[j]+R[s])^theta))-exp(b[t]*(R[s]^theta-(DR[j+1]+R[s])^theta))
      }
      # Calculating the normalizing constant cr[s,t] for each site-taxon by summing wr[s,j,t] across all distances:
      cr[s,t] <- sum(wr_raw[s,,t]) # Eqn [12]
      
      # The Localness Index (LCI):
      LCI[s,t] <- gamma[t]/gammaTil[s,t] # Eqn [16]
      # Cummulative contribution (unnormalized [gammaR_raw] and normalized [gammaR]) at 1 km predicted by regional pollen dispersal only:
      gammaTil[s,t] <- gammaTil_raw[s,t]/(gammaTil_raw[s,t] + cr[s,t]) # This line and the line below correspond to Eqn [17]
      gammaTil_raw[s,t] <- exp(b[t]*(R[s]^theta-(R[s])^theta))-exp(b[t]*(R[s]^theta-(1000+R[s])^theta))
    }     
  }
  
  # Pollen productivity:
  for(r in 1:R){
    # Regional-level productivities vary around taxon-level productivities:
    phi[1:T,r] ~ ddirch(alpha_phi[1:T]) # This line and the line below correspond to Eqn [3]
  }
  # Prior for regional-level producitivities:
  alpha_phi[1:T] <- alpha * phi_star[1:T] # Part of Eqn [3], see above
  # Prior for taxon-level productivities:
  phi_star[1:T] ~ ddirch(rep(1,T)) # Eqn [5]
  # Flat prior for alpha the scaler parameter:
  alpha ~ dunif(20,1000) # Eqn [4]
  
  for(t in 1:T){
    # Priors for relative contribution (proportion) of local vegetation:
    gamma[t] ~ dunif(0,1) # Eqn [6]
    
    # Monitoring predicted local pollen dispersal (Fl):
    for(i in 1:I){
      Fl[t,i] <- Fl_raw[t,i]/Fl_raw[t,I] # This line and the line below correspond to Eqn [13]
      Fl_raw[t,i] <-  1 - exp(b[t]*(30^theta-(DL[i+1]+30)^theta))
    }
    
    # Monitoring predicted regional pollen dispersal with (Flr) and without (Fr) explicitely considering local vegetation:
    for(j in 1:J){
      Fr[t,j] <- Fr_raw[t,j]/Fr_raw[t,J] # This line and above correspond to Eqn [14]
      Fr_raw[t,j] <-  1 - exp(b[t]*(30^theta-(DR[j+1]+30)^theta))
      Flr[t,j] <- gamma[t] + (1-gamma[t])*Fr[t,j] # Eqn [15]
    }
    
    # Priors for taxon-specific pollen dispersal parameters:
    b[t] ~ dexp(0.1) # Eqn [18]
  } 
  
  # Fixed value for theta (value representing unstable atmosphere):
  theta <- 0.1 # n paratmer = 0.2 in Sutton's equation
}
