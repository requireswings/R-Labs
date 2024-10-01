library(rhandsontable)
library(shiny)
library(jagsUI)
library(MASS) #needed for mvrnorm()
library(fields) #needed for image.plot()
library(plotrix) #needed for draw.ellipse
library(BB) #needed for dfsane()
library(car) #needed for ellipse()
library(rmarkdown) #needed for generating downloadable reports

#LOAD VALUES AND FUNCTIONS####


#season hatch/fate calculation functions####
h.ex.calc<-function(parms){
  
  lin.f.ex<-exp(parms$alpha.f)
  lin.a.ex<-exp(parms$alpha.a+ parms$beta.a.ex)
  lin.p.ex<-exp(parms$alpha.p  + parms$beta.p.ex)
  survp.ex<-1/(lin.f.ex+lin.a.ex+lin.p.ex+1)
  pred.p.ex<-((lin.p.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  a<-((lin.a.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  flood.p.ex<-((lin.f.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  
  h<-survp.ex^34   #probability of period survival exclosed nest
  o<-flood.p.ex+pred.p.ex 
  m=parms$m
  r2=parms$r2
  r3=parms$r3
  h.ex<-h + (o*r2*h + a*(1-m)*r2*h) + o*r2*(o*r3*h + a*(1-m)*r3*h) + a*(1-m)*r2*(o*r3*h+a*(1-m)*r3*h)
  output<-list(h.ex=h.ex, a=a)
  return(output)
}

h.un.calc<-function(parms){
  lin.f<-exp(parms$alpha.f)
  lin.a<-exp(parms$alpha.a)
  lin.p<-exp(parms$alpha.p)
  survp<-1/(lin.f+lin.a+lin.p+1)   #daily survival probability
  
  #period fate probabilities for each nest
  pred.p<-((lin.p/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #predation unexclosed
  a<-((lin.a/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #abandonment unexclosed
  flood.p<-((lin.f/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #flooding probability
  o<-flood.p+pred.p
  h<-survp^34  #probability of period survival unexclosed nest
  m=parms$m
  r2=parms$r2
  r3=parms$r3
  h.un <- h + (o*r2*h + a*(1-m)*r2*h) + o*r2*(o*r3*h + a*(1-m)*r3*h) + a*(1-m)*r2*(o*r3*h+a*(1-m)*r3*h) #revised; mistake somewhere in original
  output<-list(h.un=h.un, a=a)
  return(output)
}
#logodds####
logodds=function(x){log(x/(1-x))}

#projection function####
lambda.calc <- function(parms, sd.parms, eta.a=0, eta.p=0, n.ex=0, n.iter=1000){ #take out etas - not used
  n.iter=n.iter
  eta.a=eta.a
  eta.p=eta.p
  lambda <- rep(NA, n.iter)
  aban.counts.ex <- rep(NA, n.iter)
  aban.counts.un <- rep(NA, n.iter)
  dead.females.ex <- rep(NA, n.iter)
  
  #run iterations
  for (i in 1:n.iter){
    #draw stochastic values
    parms.iter <- list(
      survmean=plogis(mvrnorm(1,c(qlogis(parms$Phij),qlogis(parms$Phia)),sd.parms$vcovSurv)),
      f = plogis(rnorm(1,qlogis(parms$f),sd.parms$sd.f)),
      ys = rnorm(1,logodds(parms$ys),sd.parms$ys), #probability of second-year bird nesting
      yt = parms$yt, #probability of third-year bird nesting fixed at 0.99
      r2 = rnorm(1, logodds(parms$r2), sd.parms$r2),
      r3 = rnorm(1, logodds(parms$r3), sd.parms$r3),
      m = plogis(rnorm(1,qlogis(parms$m), sd.parms$m)),
      alpha.f = rnorm(1, parms$alpha.f, sd.parms$alpha.f),
      aban.coeff = mvrnorm(1, c(parms$alpha.a, parms$beta.a.ex), sd.parms$vcovAban),
      pred.coeff = mvrnorm(1, c(parms$alpha.p, parms$beta.p.ex), sd.parms$vcovPred)
    )
    parms.iter$alpha.a <- parms.iter$aban.coeff[1] + eta.a
    parms.iter$beta.a.ex <- parms.iter$aban.coeff[2]
    parms.iter$alpha.p <- parms.iter$pred.coeff[1] + eta.p
    parms.iter$beta.p.ex <- parms.iter$pred.coeff[2]
    
    #baseline abandonment rate (unexclosed)
    a.p.un<-h.un.calc(parms=parms.iter)$a
    
    #hatch and exclosure-related abandonment probabilities
    h.un<-h.un.calc(parms=parms.iter)$h.un
    h.ex<-h.ex.calc(parms=parms.iter)$h.ex
    a.p.ex<-h.ex.calc(parms=parms.iter)$a
    
    #fecundity
    Fa<-parms.iter$yt*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))
    Fs<-parms.iter$ys*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))
    
    #breeding-season mortality with abandonment
    m<-parms.iter$m; r2<-parms.iter$r2; r3<-parms.iter$r3
    m.a.un=parms.iter$yt*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m)) 
    m.s.un=parms.iter$ys*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m))
    m.a.ex=parms.iter$yt*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))  
    m.s.ex=parms.iter$ys*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))
    
    aban.counts.ex[i]<-a.p.ex*n.ex*sum(parms$Na0,parms$Ns0)
    aban.counts.un[i]<-a.p.un*(1-n.ex)*sum(parms$Na0,parms$Ns0)
    dead.females.ex[i]<- m.s.ex*parms$Ns0 + m.a.ex*parms$Na0
    
    #annual survival rates	
    Phij.w <- parms.iter$survmean[1]^(10/12) #winter post-fledging juvenile survival
    Phia.w <- parms.iter$survmean[2]^(10/12) #winter adult survival
    Phij.b <- parms.iter$survmean[1]^(2/12)  #second-year breeding season survival
    Phia.b <- parms.iter$survmean[2]^(2/12)  #ASY breeding survival
    Phij.b <- Phij.b/(1-m.s.un)*(1-m.s.ex*n.ex)  #add in probability of surviving exclosure-related abandonment
    Phia.b <- Phia.b/(1-m.a.un)*(1-m.a.ex*n.ex)
    Phij.ann <- Phij.b*Phij.w 
    Phia.ann <- Phia.b*Phia.w
    
    #matrix calculations
    f<-parms.iter$f
    Lesmat<-matrix(c(Phij.ann*Fs*f,Fa*Phia.ann,Phij.ann*f,Phia.ann),2,2,byrow=TRUE)  #found mistake in original; entry [2,1] had adult survival, not juv
    lambda[i]<-eigen(Lesmat)$values[1]
  } #n.iter
  output <- list(lambda=lambda, aban.counts.ex=aban.counts.ex, aban.counts.un=aban.counts.un, dead.females.ex=dead.females.ex)
  return(output)
}

#mean parameters####

mean.parms <- list(
  yt=0.99,  #probability of third-year bird nesting
  ys=0.68,  #probability of second-year bird nesting
  Phij=0.52,  
  Phia=0.74,
  f=0.4,
  m=0.34,
  E = 0.94,
  r2 = 0.7,
  r3 = 0.7,
  beta.a.ex = 1.284,
  beta.p.ex = -2.532,
  alpha.a = -7.447,
  alpha.p = -4.257,
  alpha.f = -6.040
)	

mean.alpha.a = -7.447
mean.alpha.p = -4.257
CV=0.1 #coefficient of variation for most stochastic parameters

#standard deviations for stochasticity	
sd.parms.init <- list(
  vcovSurv=matrix(c((logodds(mean.parms$Phij)*CV)^2,logodds(mean.parms$Phij)*CV*logodds(mean.parms$Phia)*CV, 
                    logodds(mean.parms$Phij)*CV*logodds(mean.parms$Phia)*CV, (logodds(mean.parms$Phia)*CV))^2,2,2),  #had sd in here earlier, not variance
  sd.f = abs(logodds(mean.parms$f))*0.06,
  yt = logodds(mean.parms$yt)*CV,
  ys = logodds(mean.parms$ys)*CV,
  r2 = abs(logodds(mean.parms$r2))*CV,
  r3 = abs(logodds(mean.parms$r3))*CV,
  m = abs(logodds(mean.parms$m))*1.04,
  beta.a.ex = 0.348,
  beta.p.ex = 0.272,
  alpha.a = 0.386,
  alpha.p = 0.116,
  alpha.f = 0.222,
  vcovAban = matrix(c(0.386^2, -0.092, -0.092, 0.348^2),2,2),
  vcovPred = matrix(c(0.116^2, -0.007, -0.007, 0.272^2),2,2)
)

#Example data for download####
exampleDat <- read.csv("sampleData.csv", header=T)

#data needed for 3D plots and confidence intervals####
#for 3D plot surface:
z <- as.matrix(read.csv("zvalues.csv", header=F))
predictors <- read.csv("predictors.csv",header=T)
x <- predictors$x
y<-predictors$y

#for extracting z values for site point and CIs:
df2<-read.csv("predictionFrame.csv", header=T)
lambda.loess2 = loess(lambda.diff~a.prob.ex*p.prob.un, data = df2) #prediction function

#make colors for 3D plot
nrz<-length(x)-1
ncz<-length(x)-1

jet.colors <-  colorRampPalette(c("darkred","red","orange","yellow","green"))
nbcol<-64
color<-jet.colors(nbcol)
zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]
facetcol<-cut(zfacet,nbcol)

#SERVER####
shinyServer(
 function(input, output) {

#upload and edit nest data  
  nestdata <- reactiveValues()
  
  observeEvent(input$file,{
    if (!is.null(input$file)){
    nestdata$nestdata <- read.csv(input$file$datapath, stringsAsFactors = F)
    nestdata$nestdata$Date <- as.Date(nestdata$nestdata$Date, format="%m/%d/%Y")
    nestdata$nestdata$Interval <- as.integer(rep(NA, length(nestdata$nestdata[,1])))
    }
  })
  
  observeEvent(input$hot, {
    if(!is.null(input$hot)){
      nestdata$nestdata <- hot_to_r(input$hot)
    }
    
  })
  
  observeEvent(input$clean,{
    if(input$clean && !is.null(nestdata$nestdata)){
      nestdata$nestdata <- nestdata$nestdata[order(nestdata$nestdata$Nest.ID, nestdata$nestdata$Date),]
          
          Interval<-rep(NA, length(nestdata$nestdata$Nest.ID))
          postHatch <- rep(0, length(nestdata$nestdata$Nest.ID))
          
          for(i in 2:length(nestdata$nestdata$Nest.ID)){
            if (nestdata$nestdata$Nest.ID[i]==nestdata$nestdata$Nest.ID[i-1] && !is.na(nestdata$nestdata$Date[i]) &&
                !is.na(nestdata$nestdata$Date[i-1])) {Interval[i]=(as.numeric(nestdata$nestdata$Date[i]-nestdata$nestdata$Date[i-1]))}
            else {Interval[i]=NA}
            
            if(nestdata$nestdata$Nest.ID[i]==nestdata$nestdata$Nest.ID[i-1] && nestdata$nestdata$Status[i-1]==6)
              postHatch[i] <- 1
          }
          nestdata$nestdata$Interval <- Interval
          nestdata$nestdata <- nestdata$nestdata[postHatch==0,] #remove post-hatch rows
          nestdata$nestdata <- subset(nestdata$nestdata, select=c("Nest.ID","Date","Status","Exclosed","Interval"))
      
    }
  })
  

#NEST FATE ANALYSIS####   
 
  survival <- eventReactive(input$analyze.go,{
    validate(
      need(!is.null(nestdata$nestdata), "Please upload data"),
      need(length(which(is.na(nestdata$Date)))<1, "Warning: missing dates detected. Double-check date formatting.")
    )
    withProgress(message = "Calculating, this may take a few minutes...",{
     nestdata <- nestdata$nestdata
     nest.count <- length(unique(nestdata$Nest.ID))
     
     #in case it hasn't been done earlier
     nestdata<- nestdata[order(nestdata$Nest.ID, nestdata$Date),]
     postHatch <- rep(0, length(nestdata$Nest.ID))
     
     for(i in 2:length(nestdata$Nest.ID)){
       if (nestdata$Nest.ID[i]==nestdata$Nest.ID[i-1] && !is.na(nestdata$Date[i]) &&
           !is.na(nestdata$Date[i-1])) {nestdata$Interval[i]=(as.numeric(nestdata$Date[i]-nestdata$Date[i-1]))}
       if(nestdata$Nest.ID[i]==nestdata$Nest.ID[i-1] && nestdata$Status[i-1]==6)
         postHatch[i] <- 1
     }
 
     nestdata <- nestdata[postHatch==0,] #remove post-hatch rows
     nestdata<-nestdata[!is.na(nestdata$Interval[]),] #remove first nest check
     nestdata<-nestdata[nestdata$Interval[]>=1,] #removes any 0 or negative (=error) nest check intervals
     
     n <- length(nestdata[,1])
     
     Surv=rep(0,n)
     Aban=rep(0,n)
     Flood=rep(0,n)
     Pred=rep(0,n)
     
     for (i in 1:n){
       Surv[i][nestdata[i,"Status"]==1|nestdata[i,"Status"]==6]<-1
       Aban[i][nestdata[i,"Status"]==4]<-1
       Pred[i][nestdata[i,"Status"]==2]<-1
       Flood[i][nestdata[i,"Status"]==3]<-1
     }			
     Fate<-cbind(Surv, Aban, Pred, Flood)
     
     for (i in 1:length(Fate[,1])){
       if (sum(Fate[i,])==0) {Fate[i,]<-NA}
     }
     #replace Y/N in exclosed column to 0 and 1 if necessary
     if(!is.numeric(nestdata$Exclosed)){
       nestdata$Exclosed<-gsub(".*(n).*","0",nestdata$Exclosed, ignore.case=T)
       nestdata$Exclosed<-gsub(".*(y).*", "1", nestdata$Exclosed, ignore.case=T)
       nestdata$Exclosed<-as.numeric(nestdata$Exclosed)
     }
     
     #data validation
     validate(
       need(!is.na(max(nestdata$Exclosed)), "Cannot Proceed: Missing data in exclosure column. Please fill in missing data.")

     )
     
     
     #bundle data and prep for jags
     win.data<-list(ex = nestdata[,"Exclosed"],n=n,interval=nestdata[,"Interval"],Fate=Fate,
                    mean.alpha.a = mean.alpha.a, mean.alpha.p=mean.alpha.p)
     inits<-function(){list(alpha.p=rnorm(1,0,1), alpha.a=rnorm(1,0,1), alpha.f=rnorm(1,0,1), beta.a.ex=rnorm(1,0,1),
                            beta.p.ex=rnorm(1,0,1))}
     
     params<-c("alpha.p","alpha.f","alpha.a", "beta.a.ex","beta.p.ex", "Hatchp.ex", "Hatchp.un",
               "Survp.ex","Survp.un","Abanp.ex","Abanp.un", "Predp.ex","Predp.un", "Flood.p")
     
     nc<-3
     
     out<-autojags(win.data, inits, params, "single_site_model.txt", n.chains=nc, parallel=TRUE)
     
     #extract estimates for later use
     #Period fate probabilities:
     Fate.est <- list(Hatch.Ex = out$mean$Hatchp.ex,
                      Hatch.Un = out$mean$Hatchp.un,
                      Aban.Ex = out$mean$Abanp.ex,
                      Aban.Un = out$mean$Abanp.un,
                      Pred.Ex = out$mean$Predp.ex,
                      Pred.Un = out$mean$Predp.un,
                      Flood = out$mean$Flood.p)
     
     #standard deviations for fate probabilities
     Fate.SD <- list(Hatch.Ex = out$sd$Hatchp.ex,
                     Hatch.Un = out$sd$Hatchp.un,
                     Aban.Ex = out$sd$Abanp.ex,
                     Aban.Un = out$sd$Abanp.un,
                     Pred.Ex = out$sd$Predp.ex,
                     Pred.Un = out$sd$Predp.un,
                     Flood = out$sd$Flood.p)
     
     #Model Estimates:
     alpha.p <- out$mean$alpha.p
     alpha.a <- out$mean$alpha.a
     alpha.f <- out$mean$alpha.f
     beta.a.ex <- out$mean$beta.a.ex
     beta.p.ex <- out$mean$beta.p.ex

     #covariances
     vcovAban <- cov(cbind(out$sims.list$alpha.a, out$sims.list$beta.a.ex))
     vcovPred <- cov(cbind(out$sims.list$alpha.p, out$sims.list$beta.p.ex))
     
     alpha.f.sd <- out$sd$alpha.f
     
     #prepare output
     list(alpha.p=alpha.p, alpha.a=alpha.a, alpha.f=alpha.f, beta.a.ex=beta.a.ex,  beta.p.ex=beta.p.ex, vcovAban=vcovAban,
          vcovPred=vcovPred, Fate.est=Fate.est,Fate.SD=Fate.SD, out=out, alpha.f.sd=alpha.f.sd) 
    })
  }) #survival()
  
output$SmallSampleWarn <- renderText({ 
  if (input.summary()$nest.count<=10){
    paste("Warning: Data Entered Do Not Meet Minimum Recommended Sample Size of 10 Nests")
  }
  })  
  
#update params to use in lambda calculations	####
  parms <- reactive({
    list(
      Na0 = (input$Pairs)/2,
      Ns0 = (input$Pairs)/2,
      yt=0.99,  #probability of third-year bird nesting
      ys=0.68,  #probability of second-year bird nesting
      Phij=0.52,  
      Phia=0.74,
      f=0.4,
      m=0.34,
      E = 0.94,
      r2 = 0.7,
      r3 = 0.7,
      beta.a.ex = survival()$beta.a.ex,
      beta.p.ex = survival()$beta.p.ex,
      alpha.a = survival()$alpha.a,
      alpha.p = survival()$alpha.p,
      alpha.f = survival()$alpha.f
    ) })
  
  sd.parms <- reactive({
    list(
      vcovSurv=matrix(c((logodds(parms()$Phij)*CV)^2,logodds(parms()$Phij)*CV*logodds(parms()$Phia)*CV, 
                        logodds(parms()$Phij)*CV*logodds(parms()$Phia)*CV, (logodds(parms()$Phia)*CV))^2,2,2),
      sd.f = abs(logodds(parms()$f))*0.06,
      yt = logodds(parms()$yt)*CV,
      ys = logodds(parms()$ys)*CV,
      r2 = abs(logodds(parms()$r2))*CV,
      r3 = abs(logodds(parms()$r3))*CV,
      m = abs(logodds(parms()$m))*1.04,  
      alpha.f = survival()$out$sd$alpha.f,
      vcovAban = survival()$vcovAban,
      vcovPred = survival()$vcovPred
    ) })

  ###################################################################################### lambda
  
#LAMBDA####
  
  trajectory<-eventReactive(input$button, {
   
    withProgress(message = "Calculating Projections ..", {

      traj.noEx <- lambda.calc(parms=parms(), sd.parms=sd.parms(), n.ex=0, n.iter=10000) #lambda without exclosures
      traj.Ex <- lambda.calc(parms=parms(), sd.parms=sd.parms(), n.ex=1, n.iter=10000) #lambda with exclosures
      
      #calculate probabilities of decline and growth
      decline.ex <- mean(traj.Ex$lambda<1)
      decline.un <- mean(traj.noEx$lambda<1)
      
      decline.steep.ex <- mean(traj.Ex$lambda<0.95)
      decline.steep.un <- mean(traj.noEx$lambda<0.95)
      
      growth.ex <- mean(traj.Ex$lambda>1)
      growth.un <- mean(traj.noEx$lambda>1)
      
      growth.steep.ex <-mean(traj.Ex$lambda>1.05)
      growth.steep.un <-mean(traj.noEx$lambda>1.05)
      
      #prepare output

      list(lambda.un=traj.noEx$lambda, lambda.ex=traj.Ex$lambda, aban.counts.ex = traj.Ex$aban.counts.ex,
           decline.ex=decline.ex, decline.un=decline.un, decline.steep.ex=decline.steep.ex, decline.steep.un=decline.steep.un,
           growth.ex=growth.ex, growth.un=growth.un, growth.steep.ex=growth.steep.ex, growth.steep.un=growth.steep.un)
      
    })#withProgress
    
  })#trajectory
  
  lambda.summary <- reactiveValues(lambda.ex=NA, lambda.un=NA)
  
 observeEvent(input$button, {
    if(input$button) {
      lambda.summary$lambda.un <- trajectory()$lambda.un; lambda.summary$lambda.ex <- trajectory()$lambda.ex
    }
  })
 
  ################################################################################################### lamdba
 
#THRESHOLDS####
  n.vals <- 200 #number of aban values to draw for use in simulations
  ni<-200 #number of iterations per draw
  thresholds<-eventReactive(input$threshold.data, {
    withProgress(message = "Calculating Thresholds...", {
      #draw abandonment values
      eta.a.vec <- runif(n.vals,-2,4)
      eta.a.vec <- eta.a.vec[sort.list(eta.a.vec)] #sort from low to high; makes plotting later easier
      #define stuff
      lambda.mat.ex <- lambda.mat.un <- abans.ex <- matrix(NA, nrow=ni, ncol=n.vals)
      
      #run projection models
      for(i in 1:n.vals){ 
        growth.ex <- lambda.calc(parms=parms(), sd.parms=sd.parms(), eta.a=eta.a.vec[i], n.ex=1, 
                                 n.iter=ni)
        lambda.mat.ex[,i] <- growth.ex$lambda
        abans.ex[,i] <- growth.ex$aban.counts.ex
        incProgress(1/n.vals) #progress bar for simulations 
      }
      #summarize data  
      ref.line <- lambda.calc(parms=parms(), sd.parms=sd.parms(),   
                              n.ex=0, n.iter=ni) #unexclosed reference
      aban.ex.means <- colMeans(abans.ex)
      prob.decline <- prob.decline.ref <- rep(NA, n.vals)
      prob.decline.ref <- mean(ref.line$lambda[]<1)
      for (i in 1:n.vals){
        prob.decline[i] <- length(which(lambda.mat.ex[,i]<1))/ni 
      }
      prob.curv <- smooth.spline(aban.ex.means, prob.decline)
      diffs <- abs(prob.curv$y-prob.decline.ref)
      aban.tolerance <- round(prob.curv$x[min(which(diffs[]==min(diffs)))],0)
      #package data for output  
      list(prob.decline=prob.decline, prob.decline.ref=prob.decline.ref, aban.ex.means=aban.ex.means, prob.curv=prob.curv,
           aban.tolerance=aban.tolerance)
    }) #thresholds withProgress
  }) #thresholds 
#########################################################################################################

thresh.summary <- reactiveValues(prob.decline=NA, prob.decline.ref=NA, aban.ex.means=NA, prob.curv=NA,aban.tolerance=NA)  

observeEvent(input$threshold.data, {
  if(input$threshold.data) {
    thresh.summary$prob.decline <- thresholds()$prob.decline; thresh.summary$prob.decline.ref <- thresholds()$prob.decline.ref;
    thresh.summary$aban.ex.means <- thresholds()$aban.ex.means; thresh.summary$prob.curv <- thresholds()$prob.curv;
    thresh.summary$aban.tolerance <- thresholds()$aban.tolerance
  }
})  
    
#SCENARIO MODELING####
  
  risk.select<-reactive({  #retrieve selected values for abandonment and predation
    list(
      a=sliderInput$abanRisk,
      p=sliderInput$predRisk
    )
  })
  
#scenario modeling function  
  scenario <- eventReactive(input$scenario,
                            {withProgress(message = "Calculating ..",{
                              
                              #back-transform inputs to get alpha.p and beta.a.ex values
                              back.calc <- function(x){
                                alpha.a <- mean.parms$alpha.a
                                alpha.p <- x[2]
                                alpha.f <- mean.parms$alpha.f
                                beta.a.ex <- x[1]
                                beta.p.ex <- mean.parms$beta.p.ex
                                a <- x[3]
                                p <- x[4]
                                y <- rep(NA,2)
                                #linear predictors:
                                linp.ex <- exp(alpha.p + beta.p.ex)
                                linp <- exp(alpha.p)
                                linf <- exp(alpha.f)
                                lina.ex <- exp(alpha.a + beta.a.ex)
                                lina <- exp(alpha.a)
                                #denominator
                                den.ex <- 1 + linp.ex + linf + lina.ex
                                den <- 1 + linp + linf + lina
                                surv.ex <- 1/(den.ex) #daily survival probability for exclosed nests
                                surv <- 1/den #daily survival probability for unexclosed nests
                                #solve for beta.a.ex given a, with y[1] set to 0:
                                y[1] <- ((lina.ex/den.ex)/(1-surv.ex))*(1-surv.ex^34)-a 
                                #take equation p = (stuff) and substract p (make y[2] 0)
                                y[2] <- ((linp/den)/(1-surv))*(1-surv^34)-p
                                y[3] <- a-a
                                y[4] <- p-p
                                y
                              }
                              
                              #supply 1 and -4 as initial parameter value guesses
                              parm.est <- dfsane(c(1,-4,input$abanRisk/100,input$predRisk/100),back.calc)$par
                              
                              parms.scen <- reactive({
                                list(
                                  Na0 = (input$ScenPairs)/2,
                                  Ns0 = (input$ScenPairs)/2,
                                  yt=0.99,  #probability of third-year bird nesting
                                  ys=0.68,  #probability of second-year bird nesting
                                  Phij=0.52,  
                                  Phia=0.74,
                                  f=input$fledge/100,
                                  m=input$mortality/200,
                                  E = 0.94,
                                  r2 = 0.7,
                                  r3 = 0.7,
                                  beta.a.ex = parm.est[1],
                                  beta.p.ex = mean.parms$beta.p.ex,
                                  alpha.a = mean.parms$alpha.a,
                                  alpha.p = parm.est[2],
                                  alpha.f = mean.parms$alpha.f
                                ) })
                              
                              sd.parms.scen <- reactive({
                                list(
                                  vcovSurv=matrix(c((logodds(parms.scen()$Phij)*CV)^2,logodds(parms.scen()$Phij)*CV*logodds(parms.scen()$Phia)*CV, 
                                                    logodds(parms.scen()$Phij)*CV*logodds(parms.scen()$Phia)*CV, (logodds(parms.scen()$Phia)*CV))^2,2,2),  #had sd in here earlier, not variance
                                  sd.f = abs(logodds(parms.scen()$f))*0.06,
                                  yt = logodds(parms.scen()$yt)*CV,
                                  ys = logodds(parms.scen()$ys)*CV,
                                  r2 = abs(logodds(parms.scen()$r2))*CV,
                                  r3 = abs(logodds(parms.scen()$r3))*CV,
                                  m = min(abs(logodds(parms.scen()$m))*1.04,abs(logodds(parms.scen()$m+0.01))*1.04),  #use 0.01 if value of 0 entered
                                  alpha.f = sd.parms.init$alpha.f,
                                  vcovAban = sd.parms.init$vcovAban,
                                  vcovPred = sd.parms.init$vcovPred
                                ) })
                              
                              traj.noEx <- lambda.calc(parms=parms.scen(), sd.parms=sd.parms.scen(), n.ex=0, n.iter=1000)
                              traj.Ex <- lambda.calc(parms=parms.scen(), sd.parms=sd.parms.scen(), n.ex=1, n.iter=1000)
                              
                              #calculate probabilities of decline and growth
                              decline.ex <- mean(traj.Ex$lambda<1)
                              decline.un <- mean(traj.noEx$lambda<1)
                              
                              decline.steep.ex <- mean(traj.Ex$lambda<0.95)
                              decline.steep.un <- mean(traj.noEx$lambda<0.95)
                              
                              growth.ex <- mean(traj.Ex$lambda>1)
                              growth.un <- mean(traj.noEx$lambda>1)
                              
                              growth.steep.ex <-mean(traj.Ex$lambda>1.05)
                              growth.steep.un <-mean(traj.noEx$lambda>1.05)
                              
                              list(lambda.un=traj.noEx$lambda, lambda.ex=traj.Ex$lambda, decline.ex=decline.ex, 
                                   decline.un=decline.un, decline.steep.ex=decline.steep.ex, decline.steep.un=decline.steep.un,
                                   growth.ex=growth.ex, growth.un=growth.un, growth.steep.ex=growth.steep.ex, 
                                   growth.steep.un=growth.steep.un, parm.est=parm.est)
                            })
                        }#with progress
  )#event reactive
  

  #THRESHOLDS FOR SCENARIO####

  scen.thresholds<-eventReactive(input$Scen.threshold.data, {
    withProgress(message = "Calculating Thresholds...", {
      
      parms.scen.t <- reactive({
        list(
          Na0 = (input$ScenPairs)/2,
          Ns0 = (input$ScenPairs)/2,
          yt=0.99,  #probability of third-year bird nesting
          ys=0.68,  #probability of second-year bird nesting
          Phij=0.52,  
          Phia=0.74,
          f=input$fledge/100,
          m=input$mortality/200,
          E = 0.94,
          r2 = 0.7,
          r3 = 0.7,
          beta.a.ex = scenario()$parm.est[1],
          beta.p.ex = mean.parms$beta.p.ex,
          alpha.a = mean.parms$alpha.a,
          alpha.p = scenario()$parm.est[2],
          alpha.f = mean.parms$alpha.f
        ) })
      
      sd.parms.scen.t <- reactive({
        list(
          vcovSurv=matrix(c((logodds(parms.scen.t()$Phij)*CV)^2,logodds(parms.scen.t()$Phij)*CV*logodds(parms.scen.t()$Phia)*CV, 
                            logodds(parms.scen.t()$Phij)*CV*logodds(parms.scen.t()$Phia)*CV, (logodds(parms.scen.t()$Phia)*CV))^2,2,2),  #had sd in here earlier, not variance
          sd.f = abs(logodds(parms.scen.t()$f))*0.06,
          yt = logodds(parms.scen.t()$yt)*CV,
          ys = logodds(parms.scen.t()$ys)*CV,
          r2 = abs(logodds(parms.scen.t()$r2))*CV,
          r3 = abs(logodds(parms.scen.t()$r3))*CV,
          m = min(abs(logodds(parms.scen.t()$m))*1.04,abs(logodds(parms.scen.t()$m+0.01))*1.04),  #use 0.01 if value of 0 entered
          alpha.f = sd.parms.init$alpha.f,
          vcovAban = sd.parms.init$vcovAban,
          vcovPred = sd.parms.init$vcovPred
        ) })
      
      
      n.vals <- 200 #number of aban values to draw for use in simulations
      ni<-200 #number of iterations per draw
      #draw abandonment values
      eta.a.vec <- runif(n.vals,-2,4)
      eta.a.vec <- eta.a.vec[sort.list(eta.a.vec)] #sort from low to high; makes plotting later easier
      #define stuff
      lambda.mat.ex <- lambda.mat.un <- abans.ex <- matrix(NA, nrow=ni, ncol=n.vals)
      
      #run projection models
      for(i in 1:n.vals){ 
        growth.ex <- lambda.calc(parms=parms.scen.t(), sd.parms=sd.parms.scen.t(), eta.a=eta.a.vec[i], n.ex=1, 
                                 n.iter=ni)
        lambda.mat.ex[,i] <- growth.ex$lambda
        abans.ex[,i] <- growth.ex$aban.counts.ex
        incProgress(1/n.vals) #progress bar for simulations 
      }
      #summarize data  
      ref.line <- lambda.calc(parms=parms.scen.t(), sd.parms=sd.parms.scen.t(),   
                              n.ex=0, n.iter=ni) #unexclosed reference
      aban.ex.means <- colMeans(abans.ex)
      prob.decline <- prob.decline.ref <- rep(NA, n.vals)
      prob.decline.ref <- mean(ref.line$lambda[]<1)
      for (i in 1:n.vals){
        prob.decline[i] <- length(which(lambda.mat.ex[,i]<1))/ni 
      }
      prob.curv <- smooth.spline(aban.ex.means, prob.decline)
      diffs <- abs(prob.curv$y-prob.decline.ref)
      aban.tolerance <- round(prob.curv$x[min(which(diffs[]==min(diffs)))],0)
      #package data for output  
      list(prob.decline=prob.decline, prob.decline.ref=prob.decline.ref, aban.ex.means=aban.ex.means, prob.curv=prob.curv,
           aban.tolerance=aban.tolerance)
    }) #scen.thresholds withProgress
  }) #scen.thresholds 
  #########################################################################################################
  
  scen.thresh.summary <- reactiveValues(prob.decline=NA, prob.decline.ref=NA, 
                                        aban.ex.means=NA, prob.curv=NA,aban.tolerance=NA)  
  
  observeEvent(input$Scen.threshold.data, {
    if(input$Scen.threshold.data) {
      scen.thresh.summary$prob.decline <- scen.thresholds()$prob.decline; scen.thresh.summary$prob.decline.ref <- scen.thresholds()$prob.decline.ref;
      scen.thresh.summary$aban.ex.means <- scen.thresholds()$aban.ex.means; scen.thresh.summary$prob.curv <- scen.thresholds()$prob.curv;
      scen.thresh.summary$aban.tolerance <- scen.thresholds()$aban.tolerance
    }
  })  
    
###############################################################################################################
  
#PACKAGE OUTPUT#### 
  #example data download####
  #note; this feature does not work in Rstudio window, must open in browser first
  output$example <- downloadHandler(
    filename = "exampleData.csv",
    content = function(file) {
      write.csv(exampleDat, file, row.names=F)
    }
  )

  #download cleaned data####  
  output$downloadClean <- downloadHandler(
    filename = "cleanedData.csv",
    content = function(file){
      write.csv(nestdata$nestdata, file, row.names=F)
    }
  )
  
  #input summary####
  input.summary<- reactive({
    if(is.null(nestdata$nestdata)) {nest.count <- 0; hatches <- 0; preds <- 0; abans <- 0; floods <- 0;
                        uknfail <- 0; uknfate <- 0; otherfail <- 0}
    else {nest.count<-length(unique(as.factor(nestdata$nestdata$Nest.ID))) 
    
    hatches <- sum(aggregate(nestdata$nestdata$Status, list(nestdata$nestdata$Nest.ID),max)$x==6)
    preds <- sum(nestdata$nestdata$Status==2, na.rm=T)
    abans <- sum(nestdata$nestdata$Status==4, na.rm=T)
    floods <- sum(nestdata$nestdata$Status==3, na.rm=T)
    uknfail <- sum(nestdata$nestdata$Status==5, na.rm=T)
    uknfate <- sum(nestdata$nestdata$Status==7, na.rm=T)
    otherfail <- sum(nestdata$nestdata$Status==8, na.rm=T)}
    input.summary <- as.data.frame(cbind(nest.count,  hatches, preds, abans, floods, uknfail, uknfate, otherfail))
    input.summary
  })

  output$nests<-renderText({paste("You have entered data for ", input.summary()$nest.count, "nests")})
  output$hatches<-renderText({paste("Number of hatches = ", input.summary()$hatches)})  
  output$preds<-renderText({paste("Number of depredations = ", input.summary()$preds)})
  output$abans<-renderText({paste("Number of abandonments = ", input.summary()$abans)})
  output$floods<-renderText({paste("Number of tidal/weather-caused failures = ", input.summary()$floods)})
  output$uknfail<-renderText({paste("Number of failures due to unknown cause = ", input.summary()$uknfail)})
  output$uknfate<-renderText({paste("Number of unknown fates = ", input.summary()$uknfate)})
  output$otherfail<-renderText({paste("Number of failures from other causes =", input.summary()$otherfail)})
  
#editable data table####
  output$hot <- renderRHandsontable({
    DF = nestdata$nestdata
    if (!is.null(DF)) 
      rhandsontable(DF, rowheaders = NULL, height = 900, width = 1000) %>%
     hot_col(col="Nest.ID", type="dropdown", source=unique(nestdata$nestdata$Nest.ID), strict=F, allowInvalid=T) %>%
      hot_col(col="Exclosed",type="dropdown", source=c("Y","N"),allowInvalid=F) %>%
      hot_col(col="Status", type="dropdown", source=c("1","2","3","4","5","6","7","8"), allowInvalid=F) %>%
      hot_col(col='Interval', readOnly=T) %>%
      hot_cols(columnSorting=T) %>%
      hot_context_menu(allowRowEdit = T, allowColEdit=F)  #allows user to add or delete rows
    
  })
  

  
#Bayesian nest fate analysis summary output####  
  output$hatch.ex<-renderText({paste(round(survival()$Fate.est$Hatch.Ex, digits=2)*100,"(",
                                     round(survival()$Fate.SD$Hatch.Ex, 2)*100, ") %")})
  output$hatch.un<-renderText({paste(round(survival()$Fate.est$Hatch.Un,2)*100, "(",
                                    round(survival()$Fate.SD$Hatch.Un,2)*100,") %")})
  output$aban.ex<-renderText({paste(round(survival()$Fate.est$Aban.Ex,2)*100, "(",
                                    round(survival()$Fate.SD$Aban.Ex, 2)*100, ") %")})
  output$aban.un<-renderText({paste(round(survival()$Fate.est$Aban.Un,2)*100,"(",
                                    round(survival()$Fate.SD$Aban.Un, 2)*100, ") %")})
  output$pred.ex<-renderText({paste(round(survival()$Fate.est$Pred.Ex,2)*100,"(",
                                    round(survival()$Fate.SD$Pred.Ex, 2)*100, ") %")})
  output$pred.un<-renderText({paste(round(survival()$Fate.est$Pred.Un,2)*100,"(",
                                    round(survival()$Fate.SD$Pred.Un, 2)*100, ") %" )})
  output$flood.ex<-renderText({paste(round(survival()$Fate.est$Flood,2)*100,"(",
                                     round(survival()$Fate.SD$Flood, 2)*100, ") %")}) 
  output$flood.un<-renderText({paste(round(survival()$Fate.est$Flood,2)*100,"(",
                                     round(survival()$Fate.SD$Flood, 2)*100, ") %")}) #can't show same output in two places
  
  Fate.summary <- reactiveValues(
    mat = matrix(rep(NA,8), nrow=4, ncol = 2, byrow=T, dimnames = list(c("Hatch", "Predation", "Abandonment", "Tide/Weather"),
                                                          c("Exclosed", "Not Exclosed")))
  )
  

  observeEvent(input$analyze.go, {
    if(!is.null(input$analyze.go)) {
      Fate.summary$mat[1,1] <- paste(round(survival()$Fate.est$Hatch.Ex, digits=2)*100,"(",
                                 round(survival()$Fate.SD$Hatch.Ex, 2)*100, ") %")
      Fate.summary$mat[1,2] <- paste(round(survival()$Fate.est$Hatch.Un,2)*100, "(",
                                 round(survival()$Fate.SD$Hatch.Un,2)*100,") %")
      Fate.summary$mat[2,1] <- paste(round(survival()$Fate.est$Pred.Ex,2)*100,"(",
                                  round(survival()$Fate.SD$Pred.Ex, 2)*100, ") %")
      Fate.summary$mat[2,2] <- paste(round(survival()$Fate.est$Pred.Un,2)*100,"(",
                                 round(survival()$Fate.SD$Pred.Un, 2)*100, ") %" )
      Fate.summary$mat[3,1] <- paste(round(survival()$Fate.est$Aban.Ex,2)*100, "(",
                                 round(survival()$Fate.SD$Aban.Ex, 2)*100, ") %")
      Fate.summary$mat[3,2] <- paste(round(survival()$Fate.est$Aban.Un,2)*100,"(",
                                 round(survival()$Fate.SD$Aban.Un, 2)*100, ") %")
      Fate.summary$mat[4,1] <- paste(round(survival()$Fate.est$Flood,2)*100,"(",
                                                      round(survival()$Fate.SD$Flood, 2)*100, ") %")
      Fate.summary$mat[4,2] <- paste(round(survival()$Fate.est$Flood,2)*100,"(",
                                     round(survival()$Fate.SD$Flood, 2)*100, ") %")
    }
  })
  
  #3D heat plot####
  site.pt <- c(0,0)
  output$ThreeDplot <- renderPlot({

    graph <- persp(x,y,z, expand=0.9,mar=c(10,1,0,2),
                   zlab="", xlab="",ylab="", border=NA,
                   phi=20,theta=30, ticktype="detailed", col=color[facetcol])
    x.label <- c("Abandonment Probability")
    x.label2 <- c("with Exclosures")
    y.label <- c("Predation Probability")
    y.label2 <- c("without Exclosures")
    x.label.pos <- trans3d(0.1, -0.1, -1.1, pmat=graph)
    x.label2.pos <- trans3d(0.1,-0.18,-1.15,pmat=graph)
    y.label.pos <- trans3d(0.6, 0.1, -1.05, pmat=graph)
    y.label2.pos <- trans3d(0.64, 0.1, -1.1, pmat=graph)
    text(x.label.pos$x,x.label.pos$y,labels=x.label, srt=-25, cex=1.2, adj=c(0,NA), xpd=T)
    text(x.label2.pos$x,x.label2.pos$y,labels=x.label2,srt=-25,cex=1.2,adj=c(0,NA),xpd=T)
    text(y.label.pos$x, y.label.pos$y, labels=y.label, srt=60, cex=1.2, adj=c(0,NA), xpd=T)
    text(y.label2.pos$x, y.label2.pos$y, labels=y.label2, srt=60, cex=1.2, adj=c(0,NA), xpd=T)
    z.label <- ("Gain in Growth Rate")
    z.label2 <- ("with Exclosures")
    z.label.pos <- trans3d(-0.2, 0, 0.2, pmat=graph)
    z.label2.pos <- trans3d(-0.25,-0.03,0.2, pmat=graph)
    text(z.label.pos$x, z.label.pos$y, labels=z.label, srt=-82, cex=1.2, adj=c(0,NA), xpd=T)
    text(z.label2.pos$x, z.label2.pos$y, labels=z.label2, srt=-82, cex=1.2, adj=c(0,NA), xpd=T)
    no.effect.line <- contourLines(x, y, z, nlevels=1, level=0)
    no.effect <- trans3d(no.effect.line[[1]]$x, no.effect.line[[1]]$y, rep(0, length(no.effect.line[[1]]$x)), pmat=graph)
    lines(no.effect)
    
    if(input$analyze.go) { #add site-specific estimates
    site.pt <- c(survival()$Fate.est$Aban.Ex, survival()$Fate.est$Pred.Un)
    xval<-which.min(abs(x-site.pt[1]))
    yval<-which.min(abs(y-site.pt[2]))
    zval<-z[xval,yval]
    add.pts<-trans3d(x[xval],y[yval],zval,pmat=graph)
    points(add.pts, pch=19,col="black")
    site.pos <- trans3d(site.pt[1]+0.05, site.pt[2]+0.05, zval, pmat=graph)
    text(site.pos,labels=c("Your Site"), adj=c(0,NA), cex=2)
    #variance ellipse
    var.points<-ellipse(site.pt,shape=matrix(c(survival()$Fate.SD$Aban.Ex^2,0,0,survival()$Fate.SD$Pred.Un^2),nrow=2,byrow=T), 
                        radius=1)
    var.z <- predict(lambda.loess2, newdata = var.points) #get z values for ellipse
    var.points.trans <- trans3d(var.points[,1], var.points[,2], z=(var.z+0.01), pmat=graph)
    lines(var.points.trans)
    }
  })
  
  #lambda plot####
  output$lambdaplot<-renderPlot({ 
    plot(c(1,3), c(mean(trajectory()$lambda.ex), mean(trajectory()$lambda.un)), xlim=c(0,4), 
         ylim=c(
           min(quantile(trajectory()$lambda.ex,0.025),quantile(trajectory()$lambda.un, 0.025))-0.05, 
           max(quantile(trajectory()$lambda.ex,0.975),quantile(trajectory()$lambda.un, 0.975)+0.05)
         ),
         ylab="Population Growth Rate", xlab="", pch=c(15,0), lwd=2, cex=2, xaxt="n", cex.lab=1.5)
    axis(1, at=c(1,3), labels=c("Exclosures", "No Exclosures"), cex.axis=1.5)
    lines(c(1,1), c(quantile(trajectory()$lambda.ex,0.975), quantile(trajectory()$lambda.ex, 0.025)), lwd=2)
    lines(c(3,3), c(quantile(trajectory()$lambda.un,0.975), quantile(trajectory()$lambda.un, 0.025)), lwd=2)
    lines(c(1+0.05,1-0.05), c(quantile(trajectory()$lambda.ex,0.975), quantile(trajectory()$lambda.ex, 0.975)), lwd=2)
    lines(c(1+0.05,1-0.05), c(quantile(trajectory()$lambda.ex,0.025), quantile(trajectory()$lambda.ex, 0.025)), lwd=2)
    lines(c(3+0.05,3-0.05), c(quantile(trajectory()$lambda.un,0.975), quantile(trajectory()$lambda.un, 0.975)), lwd=2)
    lines(c(3+0.05,3-0.05), c(quantile(trajectory()$lambda.un,0.025), quantile(trajectory()$lambda.un, 0.025)), lwd=2)
    abline(h=1)
  })
  
  #lambda probabilities####
  output$declineEx <- renderText({paste(round(trajectory()$decline.ex*100,0),"%")})
  output$declineUn <- renderText({ paste(round(trajectory()$decline.un*100,0), "%")})
  output$declineSteepEx <- renderText({paste(round(trajectory()$decline.steep.ex*100,0), "%")})
  output$declineSteepUn <- renderText({paste(round(trajectory()$decline.steep.un*100,0), "%")})
  output$growthEx <- renderText({paste(round(trajectory()$growth.ex*100,0), "%")})
  output$growthUn <- renderText({paste(round(trajectory()$growth.un*100,0), "%")})
  output$growthSteepEx <- renderText({paste(round(trajectory()$growth.steep.ex*100,0), "%")})
  output$growthSteepUn <- renderText({paste(round(trajectory()$growth.steep.un*100,0), "%")})
  
  
  
  traj.probs <- reactiveValues(
    mat = matrix(rep(NA,8), nrow=4, ncol = 2, byrow=T, dimnames = list(c("Rapid Growth", "Growth", 
                                                                         "Decline", "Rapid Decline"),
                                                                       c("With Exclosures", "Without Exclosures")))
  )
  
  observeEvent(input$button, {
    if(input$button) {
      traj.probs$mat[1,1] <- paste(round(trajectory()$growth.steep.ex*100,0), "%")
      traj.probs$mat[1,2] <- paste(round(trajectory()$growth.steep.un*100,0), "%")
      traj.probs$mat[2,1] <- paste(round(trajectory()$growth.ex*100,0), "%")
      traj.probs$mat[2,2] <- paste(round(trajectory()$growth.un*100,0), "%")
      traj.probs$mat[3,1] <- paste(round(trajectory()$decline.ex*100,0),"%")
      traj.probs$mat[3,2] <- paste(round(trajectory()$decline.un*100,0), "%")
      traj.probs$mat[4,1] <- paste(round(trajectory()$decline.steep.ex*100,0), "%")
      traj.probs$mat[4,2] <- paste(round(trajectory()$decline.steep.un*100,0), "%")
    }
  })
  
  
  #thresholds plot####
  
  output$thresh.abans <- renderPlot({
    par(mar=c(5,5,4,2)+0.1)
    plot(x=NULL, y=NULL, xlab="Number of Abandonments",
         ylab="Probability of Population Decline", ylim=c(0,1), xlim=c(0, max(thresholds()$aban.ex.means)), 
         xaxt = "n")
    axis(1, at=seq(0,max(thresholds()$aban.ex.means),by=1))
    prob.curv <- smooth.spline(thresholds()$aban.ex.means, thresholds()$prob.decline)
    lines(thresholds()$prob.curv, lwd=2)
    if (length(thresholds()$prob.curv$y)>=n.vals) {y.coord <- thresholds()$prob.curv$y[1:200]} else 
    {y.coord <- c(thresholds()$prob.curv$y[1:(n.vals-length(thresholds()$prob.curv$y))], thresholds()$prob.curv)}
    polygon(c(thresholds()$aban.ex.means[1], thresholds()$aban.ex.means, thresholds()$aban.ex.means[n.vals]), 
            c(0,y.coord,0), density=10, angle=45)
    abline(h=thresholds()$prob.decline.ref, lty=2, lwd=2)
    par(xpd=T)
    legend("topright", c("Exclosed", "Unexclosed Reference"), xpd=T, lty=c(1,2), lwd=c(2,2), adj=c(0,NA)) 
    
  })
  
  output$reassess <- renderText({paste("Reassess or pull exclosures after ", thresholds()$aban.tolerance, 
                                       "observed nest abandonments")})
  
  #scenario modeling output####
  #resettable slider values####
  output$resettableScenarioValues <- renderUI({
    times <- input$reset_input
    div(id=letters[(times %% length(letters)) + 1],
        sliderInput("predRisk", "Predation Probability without Exclosures", min = 0, max = 99, value = 36), #average values for defaults
        sliderInput("abanRisk", "Abandonment Probability with Exclosures", min = 0, max = 99, value = 6),
        sliderInput("mortality", "Mortality Probability Given Abandonment", min = 0, max = 100, value = 70),
        sliderInput("fledge","Chick Survival Probability (Hatch to Fledge)", min=1, max=99, value=40))
  })
  
   
  #plot
  output$scenLambdaPlot<-renderPlot({ 
    plot(c(1,3), c(mean(scenario()$lambda.ex), mean(scenario()$lambda.un)), xlim=c(0,4), 
         ylim=c(
           min(quantile(scenario()$lambda.ex,0.025),quantile(scenario()$lambda.un, 0.025))-0.05, 
           max(quantile(scenario()$lambda.ex,0.975),quantile(scenario()$lambda.un, 0.975)+0.05)
         ),
         ylab="Population Growth Rate", xlab="", pch=c(15,0), lwd=2, cex=2, xaxt="n", cex.lab=1.5)
    axis(1, at=c(1,3), labels=c("Exclosures", "No Exclosures"), cex.axis=1.5)
    lines(c(1,1), c(quantile(scenario()$lambda.ex,0.975), quantile(scenario()$lambda.ex, 0.025)), lwd=2)
    lines(c(3,3), c(quantile(scenario()$lambda.un,0.975), quantile(scenario()$lambda.un, 0.025)), lwd=2)
    lines(c(1+0.05,1-0.05), c(quantile(scenario()$lambda.ex,0.975), quantile(scenario()$lambda.ex, 0.975)), lwd=2)
    lines(c(1+0.05,1-0.05), c(quantile(scenario()$lambda.ex,0.025), quantile(scenario()$lambda.ex, 0.025)), lwd=2)
    lines(c(3+0.05,3-0.05), c(quantile(scenario()$lambda.un,0.975), quantile(scenario()$lambda.un, 0.975)), lwd=2)
    lines(c(3+0.05,3-0.05), c(quantile(scenario()$lambda.un,0.025), quantile(scenario()$lambda.un, 0.025)), lwd=2)
    abline(h=1)
  })
  
  scen.summary <- reactiveValues(lambda.ex=NA, lambda.un=NA)

  observeEvent(input$scenario, {
    if(input$scenario) {
      scen.summary$lambda.un <- scenario()$lambda.un; scen.summary$lambda.ex <- scenario()$lambda.ex
    }
  })
  
  #growth probabilities
  output$SCENdeclineEx <- renderText({paste(round(scenario()$decline.ex*100,0),"%")})
  output$SCENdeclineUn <- renderText({ paste(round(scenario()$decline.un*100,0), "%")})
  output$SCENdeclineSteepEx <- renderText({paste(round(scenario()$decline.steep.ex*100,0), "%")})
  output$SCENdeclineSteepUn <- renderText({paste(round(scenario()$decline.steep.un*100,0), "%")})
  output$SCENgrowthEx <- renderText({paste(round(scenario()$growth.ex*100,0), "%")})
  output$SCENgrowthUn <- renderText({paste(round(scenario()$growth.un*100,0), "%")})
  output$SCENgrowthSteepEx <- renderText({paste(round(scenario()$growth.steep.ex*100,0), "%")})
  output$SCENgrowthSteepUn <- renderText({paste(round(scenario()$growth.steep.un*100,0), "%")})
  
  scen.traj.probs <- reactiveValues(
    mat = matrix(rep(NA,8), nrow=4, ncol = 2, byrow=T, dimnames = list(c("Rapid Growth", "Growth", 
                                                 "Decline", "Rapid Decline"), c("With Exclosures", "Without Exclosures")))
  )
  
  observeEvent(input$scenario, {
    if(input$scenario) {
      scen.traj.probs$mat[1,1] <- paste(round(scenario()$growth.steep.ex*100,0), "%")
      scen.traj.probs$mat[1,2] <- paste(round(scenario()$growth.steep.un*100,0), "%")
      scen.traj.probs$mat[2,1] <- paste(round(scenario()$growth.ex*100,0), "%")
      scen.traj.probs$mat[2,2] <- paste(round(scenario()$growth.un*100,0), "%")
      scen.traj.probs$mat[3,1] <- paste(round(scenario()$decline.ex*100,0),"%")
      scen.traj.probs$mat[3,2] <- paste(round(scenario()$decline.un*100,0), "%")
      scen.traj.probs$mat[4,1] <- paste(round(scenario()$decline.steep.ex*100,0), "%")
      scen.traj.probs$mat[4,2] <- paste(round(scenario()$decline.steep.un*100,0), "%")
    }
  })
  
  
  #scenario thresholds plot####
  
  output$scen.thresh.abans <- renderPlot({
    par(mar=c(5,5,4,2)+0.1)
    plot(x=NULL, y=NULL, xlab="Number of Abandonments",
         ylab="Probability of Population Decline", ylim=c(0,1), xlim=c(0, max(scen.thresholds()$aban.ex.means)), 
         xaxt = "n")
    axis(1, at=seq(0,max(scen.thresholds()$aban.ex.means),by=1))
    prob.curv <- smooth.spline(scen.thresholds()$aban.ex.means, scen.thresholds()$prob.decline)
    lines(scen.thresholds()$prob.curv, lwd=2)
    if (length(scen.thresholds()$prob.curv$y)>=n.vals) {y.coord <- scen.thresholds()$prob.curv$y[1:200]} else 
    {y.coord <- c(scen.thresholds()$prob.curv$y[1:(n.vals-length(scen.thresholds()$prob.curv$y))], scen.thresholds()$prob.curv)}
    polygon(c(scen.thresholds()$aban.ex.means[1], scen.thresholds()$aban.ex.means, scen.thresholds()$aban.ex.means[n.vals]), 
            c(0,y.coord,0), density=10, angle=45)
    abline(h=scen.thresholds()$prob.decline.ref, lty=2, lwd=2)
    par(xpd=T)
    legend("topright", c("Exclosed", "Unexclosed Reference"), xpd=T, lty=c(1,2), lwd=c(2,2), adj=c(0,NA)) 
    
  })
  
  output$SCENreassess <- renderText({paste("Reassess or pull exclosures after ", scen.thresholds()$aban.tolerance, 
                                       "observed nest abandonments")})
  
  
  
#REPORT####
  output$Report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "Report.Rmd")
      file.copy("Report.Rmd", tempReport, overwrite = TRUE)
      
    #  if(input$NestSum) {sum.include = 1} else {sum.include = 0}
        summary<-as.data.frame(input.summary())
     # if(input$NestFate) {fate.include=1} else {fate.include=0}
        
      if(!is.na(lambda.summary$lambda.ex)) {plot.include=1} else {plot.include=0}
        
      if(!is.na(thresh.summary$prob.curv)) {thresh.include=1} else {thresh.include=0}
 
      params <- list(summary=summary, lambda.ex=lambda.summary$lambda.ex, lambda.un=lambda.summary$lambda.un, 
                     plot.include=plot.include, thresh.include=thresh.include, thresholds=thresh.summary,
                      Fate.summary=Fate.summary$mat,  traj.probs=traj.probs$mat, N0=input$Pairs)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  ) #report
  
  output$ReportWord <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "Report.docx",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "ReportWord.Rmd")
      file.copy("ReportWord.Rmd", tempReport, overwrite = TRUE)
      
      #  if(input$NestSum) {sum.include = 1} else {sum.include = 0}
      summary<-as.data.frame(input.summary())
      # if(input$NestFate) {fate.include=1} else {fate.include=0}
      
      if(!is.na(lambda.summary$lambda.ex)) {plot.include=1} else {plot.include=0}
      
      if(!is.na(thresh.summary$prob.curv)) {thresh.include=1} else {thresh.include=0}
      
      params <- list(summary=summary, lambda.ex=lambda.summary$lambda.ex, lambda.un=lambda.summary$lambda.un, 
                     plot.include=plot.include, thresh.include=thresh.include, thresholds=thresh.summary,
                     Fate.summary=Fate.summary$mat,  traj.probs=traj.probs$mat, N0=input$Pairs)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  ) #report word
  
  output$ScenReport <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "ScenReport.html",
    content = function(file) {
      tempReport <- file.path(tempdir(), "ScenReport.Rmd")
      file.copy("ScenReport.Rmd", tempReport, overwrite = TRUE)

      if(!is.na(scen.thresh.summary$prob.curv)) {scen.thresh.include=1} else {scen.thresh.include=0}
      
      if(!is.na(scen.summary$lambda.ex)) {scen.include=1} else {scen.include=0}
      
      params <- list(aban.choose=input$abanRisk, fledge.choose=input$fledge,
                     pred.choose=input$predRisk, mort.choose = input$mortality, scen.include=scen.include,
                     scenario=scen.summary, scen.traj.probs=scen.traj.probs$mat, N0=input$ScenPairs,
                     thresh.include=scen.thresh.include, thresholds=scen.thresh.summary)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  ) #scenreport
  
  output$ScenReportWord <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "ScenReport.docx",
    content = function(file) {
      tempReport <- file.path(tempdir(), "ScenReportWord.Rmd")
      file.copy("ScenReportWord.Rmd", tempReport, overwrite = TRUE)
      
      if(!is.na(scen.thresh.summary$prob.curv)) {scen.thresh.include=1} else {scen.thresh.include=0}
      
      if(!is.na(scen.summary$lambda.ex)) {scen.include=1} else {scen.include=0}
      
      params <- list(aban.choose=input$abanRisk, fledge.choose=input$fledge,
                     pred.choose=input$predRisk, mort.choose = input$mortality, scen.include=scen.include,
                     scenario=scen.summary, scen.traj.probs=scen.traj.probs$mat, N0=input$ScenPairs,
                     thresh.include=scen.thresh.include, thresholds=scen.thresh.summary)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  ) #scenario report word
  
  
 }
)#shinyServer
