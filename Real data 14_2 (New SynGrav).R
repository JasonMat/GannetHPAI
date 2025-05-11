##################################
###       1. Data import       ###
##################################

# Import of real network data (coordinates and size)
colSummaries<-read.csv("Data/Gannet_colonies_overview.csv")
di<-as.matrix(read.csv("Data/Gannet_colony_distances.csv"))
n<-nrow(di) # Number of colonies
# Colony sizes
P0<-as.numeric(round(read.csv("Data/Gannet_posterior_colony_sizes.csv", header=FALSE)))
Access<-t(as.matrix(read.csv("Data/Gannet_colony_accessibility.csv")))
Carcass<-t(as.matrix(read.csv("Data/Gannet_carcass_zeros.csv")))
Carcass<-Carcass[,10:53]/2
dats<-c("28-Mar", "4-Apr", "11-Apr", "18-Apr", "25-Apr", "2-May", "9-May", "16-May", "23-May", "30-May", "6-Jun", "13-Jun", "20-Jun", "27-Jun", "4-Jul", "11-Jul", "18-Jul", "25-Jul", "1-Aug", "8-Aug", "15-Aug", "22-Aug", "29-Aug", "5-Sep", "12-Sep", "19-Sep", "26-Sep", "3-Oct", "10-Oct", "17-Oct", "24-Oct", "31-Oct", "7-Nov", "14-Nov", "21-Nov", "28-Nov", "5-Dec", "12-Dec", "19-Dec", "26-Dec", "2-Jan", "9-Jan", "16-Jan")

##################################
###      2. Model fitting      ###
##################################

library(runjags)
library(coda)
library(stringr)


###### Model statement ######
AIModDet<-  "model{

  # Initialisation of states
  for(i in 1:n)
  {
    P[i,1]<-P0[i]
    S[i,1]<-P0[i]-I[i,1]
    I[i,1]~dpois(a0)
    R[i,1]<-0
    D[i,1]<-0
    U[i,1]<-0
  }
 

  # Connectivity matrices
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      g[i,j]<-exp(-(deltaG*di[i,j]))
      h[i,j]<-exp(-(deltaH*di[i,j]))
    }
    
    H[i]<-sum(h[1:n,i]) #Colsums for synchrony
  }
  
  
  
  ### MAIN LOOP ###
  for(t in 1:(Tmax-1))
  {
    
    fG[1:n,t]<- (t(I[1:n,t])%*%g[1:n,1:n])
    fH[1:n,t]<- (t(U[1:n,t])%*%h[1:n,1:n])
    cumSt[t]<-t-min(3,t-1)
    for(i in 1:n) #destination
    {
      
      Ief[i,t]<- (fG[i,t]+aH*fH[i,t]/H[i])/P[i,t] # Contacts with infected individuals
      GIef[i,t]<-fG[i,t]-I[i,t]
      HIef[i,t]<-aH*fH[i,t]/H[i]
      hod[i,t]~dgamma(theta,theta)
      contagionsD[i,t]~dpois(hod[i,t]*(P[i,t]*c0*Ief[i,t]*(S[i,t]/P[i,t])^c1))
      contagions[i,t]<-min(S[i,t], contagionsD[i,t]) # Double trap???
      
      deathsAI[i,t]<-aiM*I[i,t]

      recoveries[i,t]<-rrr*I[i,t]*(1-aiM)
      
      Sr[i,t]<-max(0, S[i,t]-contagions[i,t]) # Double trap???
      deathsNAS[i,t]<-naM*Sr[i,t]
      deathsNAR[i,t]<-naM*R[i,t]
      R[i,t+1]<-max(0,R[i,t]+recoveries[i,t]-deathsNAR[i,t])
      I[i,t+1]<-max(0, I[i,t]+contagions[i,t]-deathsAI[i,t]-recoveries[i,t])
      S[i,t+1]<-max(0,S[i,t]-contagions[i,t]-deathsNAS[i,t])
      U[i,t+1]<-min(I[i,t+1],1) # Is disease present in the colony?
      
      deaths[i,t]<-deathsAI[i,t]+deathsNAS[i,t]+deathsNAR[i,t]
      D[i,t+1]<-D[i,t]+deaths[i,t]
      P[i,t+1]<-I[i,t]+S[i,t]+R[i,t]

      cumDeaths[i,t]<-sum(deaths[i,cumSt[t]:t])
      Carcass[i,t]~dnorm(cumDeaths[i,t]*esc,1/(cumDeaths[i,t]*esc*(1-esc)+0.001))
      Carcass2[i,t]~dnorm(cumDeaths[i,t]*esc,1/(cumDeaths[i,t]*esc*(1-esc)+0.001))
    }
    Infected[t]<-sum(I[1:n,t])
    GItot[t]<-sum(GIef[1:n,t])
    HItot[t]<-sum(HIef[1:n,t])
  }
  EndInfc<-sum(I[1:n,Tmax])
  EndRecv<-sum(R[1:n,Tmax])
  EndDed<-sum(D[1:n,Tmax])
  EndSsc<-sum(S[1:n,Tmax])
  
  for(i in 1:n) #destination
  {
  NInfc[i]<-I[i,Tmax]
  NRecv[i]<-R[i,Tmax]
  NDed[i]<-D[i,Tmax]
  NSsc[i]<-S[i,Tmax]
  }

  
  for(ts in 1:5)
  {
    uN[1:n,ts]<-U[1:n,timestamps[ts]]
  }

  for(i in 1:nMort) # Final deaths - based on empty nests 
  {
    muM[i]<-min(0.9999,max(0.0001,D[morts[i],Tmax]/P0[morts[i]]))
    obsMort[i]~dbinom(muM[i],testMort[i])

  }
  
  for(i in 1:nRec) # Final recoveries - based on black-eye inspection
  {
    muR[i]<-min(0.9999,max(0.0001,R[recs[i],Tmax]/(R[recs[i],Tmax]+S[recs[i],Tmax])))
    obsRec[i]~dbinom(muR[i],testRec[i])
    
  }
  
  
 ### Priors ###
  # Connectivity
  deltaHd~dnorm(-11,1)
  deltaH<-exp(deltaHd)
  
  deltaGd~dnorm(-11,1)
  deltaG<-exp(deltaGd)
  
  #deltaGd~dgamma(1,10)
  #deltaG<-(1+deltaGd)*deltaH # Synchorny by definition is a more connected system, preforming closer to perfect mixing.



  # Epidemiology
  thetaDum~dbeta(1,100)
  theta<-5+100*(1-thetaDum)
  a0~dgamma(1,10)
  # c0d~dbeta(tig,tig)
  # c0<-0.00001+c0d
  c0cap~dgamma(5,5)
  c0<-(1+c0cap)*rrr*aiM


  aH~dgamma(0.1,0.1)
  c1~dgamma(100,100)

  tig<-10 
  rd~dbeta(tig,tig)
  rrr<-0.1+0.8*rd
  aiMd~dbeta(tig,tig)
  aiM<-0.3+0.7*aiMd
  naM<-0.0001

  # Detectability
  escd~dbeta(1,10) # Detectability prior
  esc<-0.3*escd # Detectability upper limit of 30%

  #data# Tmax,n,P0, Carcass, di , obsMort, obsRec, nMort, nRec, morts, recs, testMort,testRec, timestamps
  #monitor#  c0,c1,deltaG,aH,deltaH,aiM,rrr,esc,theta, muM,muR, 
  #monitor#  Carcass2, Infected, deaths, uN, EndInfc, EndRecv, EndDed, EndSsc, NInfc, NRecv, NDed, NSsc,GItot, HItot,
  #inits# aiMd
  #
}"

####################################################################################################################

# Optional investigation of priors for background and IID infection probability
# a0min<-5
# a0max<-15
# a0sig<-runif(100,10,100)
# mu<- runif(100,0,1)
# sd<- rbeta(100,01,100)
# al<-mu*(mu*(1-mu)/sd^2-1)
# be<-(1-mu)*(mu*(1-mu)/sd^2-1)
# ll0<- -(a0min+(a0max-a0min)*mu)
# ll1<- -(a0min+(a0max-a0min)*rbeta(100,al,be))
# pp1<-exp(ll1)/(1+exp(ll1))
# pp0<-exp(ll0)/(1+exp(ll0))
# plot(ll1,pp1*70000) # Combined background and IID
# points(ll0, pp0*70000, col="red") # Background only


######################################################## Model Running ###
# Final data preparations
Tmax<-ncol(Carcass)

#Provides information for final outcome of outbreak for data-rich colonies.
#Needs direct input from the field
nMort<-1
morts<-c(3)
obsMort<-c(24000)
testMort<-c(55000)

nRec<-2
recs<-c(3,4)
obsRec<-c(24,2)
testRec<-c(87,8)

timestamps<-c(1,10,20,30,40)

a0d<-list(chain1=.6,chain2=.10, chain3=.5, chain4=.20)
aGd<-list(chain1=.10,chain2=.20, chain3=.50, chain4=.60)
aHd<-list(chain1=.10,chain2=.20, chain3=.50, chain4=.60)
aiMd<-list(chain1=0.9,chain2=0.8, chain3=0.7, chain4=0.6)
escd<-list(chain1=0.1,chain2=0.2, chain3=0.1, chain4=0.2)


results <- run.jags(AIModDet, n.chains=4, adapt=50, burnin=50, sample=100, thin=50, method="parallel")
save(results, file="Output/AISpread14_2.rda")

plot(results,plot.type=c("trace","histogram"), 
     vars=c("c0","aH","deltaH","deltaG","aiM","rrr","c1","esc", "theta"))


load("Output/AISpread14_2.rda")
for(blocks in 1:100)
{
  results <- extend.jags(results, adapt=50, sample=1000, thin=50, combine=FALSE)
  plot(results,plot.type=c("trace","histogram"), 
       vars=c("c0","aH","aG","deltaH","deltaG","aiM","rrr","c1","esc","theta"))
  save(results, file="Output/AISpread14_2.rda")
  
  
  su<-summary(results)
  
  ############   4. Plotting of population trajectories   ############
  
  par(mfrow=c(4,2))
  nms<-rownames(su)
  for(i in 1:n)
  {
    kp<-su[str_detect(nms,paste("Carcass2\\[",i,",",sep="")),]
    T<-1:(Tmax-1)
    plot(T,kp[,2], type="l", ylim=c(0,max(kp[,3],na.omit(Carcass[i,]))), xlab="Time", ylab="Carcasses",main=colSummaries[i,2])
    polygon(c(T,rev(T)),c(kp[,1],rev(kp[,3])), border=NA, col="salmon")
    lines(T,kp[,2])
    lines(T,Carcass[i,T], lty=1, col="grey", lwd=3)
    points(T,Carcass[i,T], col="grey", lwd=3)
    kp<-su[str_detect(nms,paste("deaths\\[",i,",",sep="")),]
    lines(T, kp[,2]/P0[i]*100, col="red", lty=2, lwd=2)
  }
  
  
  kp<-su[str_detect(nms,"Infected"),]
  T<-1:(Tmax-1)
  plot(T,kp[,2], type="l", ylim=c(0,max(kp[,3],na.omit(Carcass[i,]))), xlab="Time", ylab="Incoming",main="Force of infection")
  polygon(c(T,rev(T)),c(kp[,1],rev(kp[,3])), border=NA, col="salmon")
  lines(T,kp[,2])
  
  
  
  
  kp1<-su[str_detect(nms,"muM"),2]
  kp2<-su[str_detect(nms,"muR"),2]
  kp3<-su[str_detect(nms,"muS"),2]
  plot(kp1, ylim=c(0,1), xlim=c(0, max(length(kp1),length(kp2))), col="red", main="Final outcomes")
  text(1.3,0.7, "Dead", col="red")
  points(kp2, col="green")
  text(1.3,0.6, "Recovered", col="green")
  points(kp3, col="blue")
  text(1.3,0.5, "Unaffected", col="blue")
  
  par(mfrow=c(1,1))
  
}

plot(results,plot.type=c("trace","histogram"), vars=c("esc"))


