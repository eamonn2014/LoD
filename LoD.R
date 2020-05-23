#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R code 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls()) 
set.seed(333) # reproducible


list.of.packages <- c("directlabels", "ggplot2" , "xtable", "doBy", "VCA", "reshape", "nlme", "vcd","car",
                      "MASS","R2wd","tables","gtools", "rtf", "binom", "coin", 
                      "lmec", "coxme", "lme4",  
                      "arm", "rms", "plyr",
                      "directlabels","shiny","shinyWidgets","shinythemes","DT","shinyalert",
                      "Hmisc","reshape","rms","ggplot2","tidyverse","digest","DT")

lapply(list.of.packages, require, character.only = TRUE)


#options(mc.cores = parallel::detectCores())
#rstan_options(auto_write = TRUE)
options(max.print=1000000)    
fig.width7 <- 600
fig.height7 <- 570
fig.width1 <- 1340

## convenience functions
p0 <- function(x) {formatC(x, format="f", digits=1)}
p1 <- function(x) {formatC(x, format="f", digits=1)}
p2 <- function(x) {formatC(x, format="f", digits=2)}
p3 <- function(x) {formatC(x, format="f", digits=3)}
p5 <- function(x) {formatC(x, format="f", digits=5)}
logit <- function(p) log(1/(1/p-1))
expit <- function(x) 1/(1/exp(x) + 1)
inv_logit <- function(logit) exp(logit) / (1 + exp(logit))
is.even <- function(x){ x %% 2 == 0 } # function to id. odd maybe useful
options(width=200)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  the data
# LoD_data <- read.delim("~/LoD/LoD_data.txt")
# write.table(LoD_data, file = 'LoD_data.txt', sep = ' ', row.names = FALSE)

file <- "https://raw.githubusercontent.com/eamonn2014/LoD/master/LoD_data.txt"

d99 <- read.csv(file, sep="")

assa <- unique(d99$assay)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

assay <- "Assay02"
MODEL <- "logit"
MODEL2 <- "Buckley James"
LoB <- 39
hit <- .95
jitt <- "yes"
jitt1 <- .3
knots <- 5
agger <- "yes"

#--------------------------------------------------------------------------------------------------------------------

dd <- d99[d99$assay %in% assay,]

#create a hit variable to identify detected and non-detected (censored) observations
dd$count<-ifelse( dd$CT<LoB,1,0)      
dd$count[is.na(dd$count)] <- 0


LoD.count <- as.data.frame(ftable(xtabs(count ~ pool + run + dil, data =dd)))
LoD.countn <- as.data.frame(ftable(xtabs( ~ pool + run + dil, data =  dd)))
names(LoD.countn )[names(LoD.countn) == "Freq"]<-c("N")
LoD.count <- merge(LoD.count, LoD.countn, by.x = c("pool", "run", "dil"))
LoD.count$dil <- as.numeric(LoD.count$dil)

### Model detection count against dilution using logit model (probit gives slightly lower LoD for some assays)
fp.l <- tryCatch(glm(Freq/N ~ dil, weights=N, data=LoD.count, family=binomial(link=MODEL)), error=function(e) NULL)


if (agger %in% "yes") {
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  we can show we get the sam results even if aggregate, as our model does not include covariates
  LoD.count<- LoD.count[,c("dil","Freq","N")]
  detach(package:plyr)    
  library(dplyr)
  LoD.count <- LoD.count %>%
    group_by(dil) %>%
    summarize(Freq=sum(Freq), N=sum(N))
  
  fp.l <- tryCatch(glm(Freq/N ~ dil, weights=N, data=LoD.count, family=binomial(link=MODEL)), error=function(e) NULL)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
}

# use this info for the logit probit model plot
x =seq(1, 8, 0.001)
data <- tryCatch(data.frame(x , pred.fp.l = predict(fp.l, data.frame(dil=x), type="resp")), error=function(e) NULL)


# use model coefficients to estimate dilution at hit rate of 95% 

if (MODEL %in% "logit") {
  
  # use the logit model to obtain the dilution for required hit rate (probability)
  intercept <-  fp.l$coefficients[1][[1]]
  beta <-       fp.l$coefficients[2][[1]]
  d.fp.l <- (qlogis(hit) - intercept) / beta
  
} else if  ( MODEL %in% 'probit') {
  
  #  use the prbit model to obtain the dilution for required hit rate (probability)
  intercept <-  fp.l$coefficients[1][[1]]
  beta <-       fp.l$coefficients[2][[1]]
  d.fp.l <- (qnorm(hit) - intercept) / beta
  
}

## get LoD Dilution Point to predict Cq value for
dPred <- data.frame(dil  = d.fp.l)

#-------------------------------------------------------------------------------------------------------------
#create a hit variable to identify detected and non-detected (censored) observations
ddist <<- datadist(dd)
options(datadist='ddist')#

#--------------------------------------------------------------------------------------------------------------
if (MODEL2 %in% "Buckley James") {
  
  tag <- "BJ"
  
  fbj <<- bj(Surv(CT, count) ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE, link="identity",
             control=list(iter.max =250))
  
} else if  ( MODEL2 %in% 'Ordinary Least Squares') {
  
  tag <- "Ols"  
  
  fbj <<- ols(CT ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE)
}

pj<-as.data.frame(cbind(dPred, predicted = predict(fbj, type="lp", newdata=dPred, se.fit=T)))

pj$lower<-pj[2][[1]] + c(-1) * 1.96*pj[3][[1]]
pj$upper<-pj[2][[1]] + c(1) * 1.96*pj[3][[1]]
names(pj) <-c("dil", "LoD Cq", "SE", "LoD Lower 95% CI", "LoD Upper 95% CI")
print(pj)
# for some reason rms::Predict is not working so I have to use predict see lot 2
Pre <- rms::Predict(fbj)

#--------------------------------------------------------------------------------------------------------------


### Plot count vs dilution point 
plot1 <- plot(x=LoD.count$dil, y=LoD.count$Freq/LoD.count$N, pch="o"
              ,xlab=paste0("Dilution Series (LoD dil. series estimate ",dPred,", denoted by red vertical dashed line)"), ylab="Hit Rate"
              , main=paste("Plot of Hit Rates against Dilution Series for Assay \n", assay, 
                           " with Line of Estimated Dilution from the Fitted ",MODEL," Model\n using limit of blank of ", LoB,"Cq", sep= ""), cex.main=1.0, bty='n')
lines(data$x, data$pred.fp.l, type="l", lwd=1, col="blue") 
abline(h=hit, lwd=1, col="gray") 
legend('bottomleft','groups', c(paste0(MODEL, " fit"), paste0(hit*100,"% Hit Rate")),lty = c(1),
       col=c('blue', 'gray'), ncol=1, bty ="n", cex=0.8)
text(1.75, 0.80, paste("LoD Dil.Series =", dPred, sep=" "), cex = 0.95)
abline(v=dPred, lwd=1, col="red", lty=2) 

#--------------------------------------------------------------------------------------------------------------

pj <- Pre

pl <-  ggplot(Pre , anova=NULL, pval=FALSE, ylim.=c(20,45)) 

plot2 <- pl + geom_point(data=dd, aes(x = dil, y = CT )) +
  xlab(paste0("Dilution Series (LoD dil. series estimate ",dPred,", denoted by vertical dashed line)")) + ylab("Cq value") +
  geom_hline(yintercept=LoB, colour="gray", linetype="solid") +
  geom_hline(yintercept=pj[1,2], colour="blue", linetype="dashed") +
  ggtitle(paste("Assay ",assay," Plot of Cq vs. Dilution Series with the Fitted BJ Model \nLoD Dilution Point ="
                , p2(d.fp.l), "; LoD Cq Value = ", p2(pj[1,2]), ", 95%CI (", p2(pj[1,4]), ",", p2(pj[1,5])     ,")", sep=" "))



#__________________________________________________________________________________________
# it was difficult to get the ribbon added to the data.
# The approach I used before works and is implemented below.

ddx <- dd[,c("CT","dil")]#

if (jitt %in% "yes") {
  
  nn <- length(ddx$CT >=40)
  ddx$CT <- ifelse(ddx$CT >=40, ddx$CT+runif(nn,0,7),ddx$CT)  
  
} 
#-------------------------------------------------------------------------------------


ddx <- merge(ddx, Pre)
zz <- ddx

WW <-jitt1
#-------------------------------------------------------------------------------------

plot3 <- ggplot() +
  
  geom_blank(data=zz, aes(x=dil, y=CT)) +
  
  geom_point(data=zz ,
             aes(x=dil,y=CT) ,
             colour="black",  size=.3,   
             position = position_jitter(w = WW , h = 0)) +
  
  geom_line(data=zz, aes(x= dil,y=yhat), colour="blue") + 
  scale_x_continuous(breaks=1:9, labels=1:9, limits=c(.7,8.2)) +
  scale_y_continuous(breaks=c(20,25,30,35,40,45,50), labels=c(20,25,30,35,40,45,50), limits=c(20,50)) +
  
  EnvStats::stat_n_text( data=zz,
                         aes(x=dil,y=CT) ,
                         size = 5, y.pos =  20, y.expand.factor=.0, 
                         angle = 0, hjust = 0.5, family = "mono", fontface = "bold")  + 
  
  ggtitle(paste("N=",
                length(!is.na(zz$dil)),  
                ":",assay,"Cq vs. Dilution Series and fitted BJ Model \nLoD Dilution Point ="
                , p2(d.fp.l), "; LoD Cq = ", p2(pj[1,2]), ", 95%CI (", p2(pj[1,4]), ",", p2(pj[1,5])     ,")", sep=" ")) +
  
  geom_ribbon(data=zz, 
              aes(x= dil,
                  ymin=lower, ymax=upper, alpha=0.2,fill="purple") ) +
  
  labs(caption = paste0("- Cq values do not exceed machine upper limit of 40\n- LoD dilution series estimate ",dPred,", denoted by vertical dashed line")) +
  
  geom_hline(yintercept=LoB, colour="lightgray", linetype="dashed") +
  geom_hline(yintercept=pj[1,2], colour="blue", linetype="dashed") +
  geom_vline(xintercept=d.fp.l, colour="blue", linetype="dashed") 


plot3 <- plot3 +labs(y=paste0(assay, " Cq"),  x=paste("Dilution Series")) +    
  
  theme(legend.position="none") +
  theme(#panel.background=element_blank(),
    # axis.text.y=element_blank(),
    # axis.ticks.y=element_blank(),
    # https://stackoverflow.com/questions/46482846/ggplot2-x-axis-extreme-right-tick-label-clipped-after-insetting-legend
    # stop axis being clipped
    plot.title=element_text(size = 18), plot.margin = unit(c(5.5,12,5.5,5.5), "pt"),
    legend.text=element_text(size=14),
    legend.title=element_text(size=14),
    legend.position="none",
    axis.text.x  = element_text(size=15),
    axis.text.y  = element_text(size=15),
    axis.line.x = element_line(color="black"),
    axis.line.y = element_line(color="black"),
    plot.caption=element_text(hjust = 0, size = 11),
    strip.text.x = element_text(size = 16, colour = "black", angle = 0),
    axis.title.y = element_text(size = rel(1.5), angle = 90),
    axis.title.x = element_text(size = rel(1.5), angle = 0),
    strip.background = element_rect(colour = "black", fill = "#ececf0"),
    panel.background = element_rect(fill = '#ececf0', colour = '#ececf0'),
    plot.background = element_rect(fill = '#ececf0', colour = '#ececf0'),#
  ) 


print(plot3)

foo<-  d99

namez <- c("Assay", "Sample Pool", "Dilution", "Run", "Replicate", "Cq" )


names(foo) <- namez
rownames(foo) <- NULL


datatable(foo,   
          
          rownames = TRUE,
          #           
          options = list(
            searching = TRUE,
            pageLength = 15,
            paging=TRUE,
            lengthMenu = FALSE ,
            lengthChange = FALSE,
            autoWidth = FALSE,
            #  colReorder = TRUE,
            #deferRender = TRUE,
            # scrollY = 200,
            scroller = T
          ) ) #%>%

# formatRound(
#    columns= namez,   
#    digits=c(0,1,1,1,1,1,1,2,2)  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# run all assays in one go using a loop, basically repeating above code inside a loop

assay <- "Assay02"
MODEL <- "logit"
MODEL2 <- "Buckley James"
LoB <- 39
hit <- .95
knots <- 5

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

assa <- unique(d99$assay)

J<-unique(assa)
k<-length(J)

dnam = list( Assay =J, component=c("dilution", "LoD Cq",  "LoD Lower 95% CI", "LoD Upper 95% CI") )
pow1 <- array( NA, dim=sapply(dnam, length), dimnames=dnam )


for ( i in 1:k) { 
  
  dd<-d99[(d99$assay %in% J[i]),] 
  dd$count<-ifelse( dd$CT<LoB,1,0)      
  dd$count[is.na(dd$count)] <- 0
  
  LoD.count <- as.data.frame(ftable(xtabs(count ~ pool + run + dil, data =dd)))
  LoD.countn <- as.data.frame(ftable(xtabs( ~ pool + run + dil, data =  dd)))
  names(LoD.countn )[names(LoD.countn) == "Freq"]<-c("N")
  LoD.count <- merge(LoD.count, LoD.countn, by.x = c("pool", "run", "dil"))
  LoD.count$dil <- as.numeric(LoD.count$dil)
  
  dd$assay <-factor(dd$assay)
  
  # set verified and established lob values
  fp.l <- tryCatch(glm(Freq/N ~ dil, weights=N, data=LoD.count, family=binomial(link=MODEL)), 
                   error=function(e) NULL)
  
  if (MODEL %in% "logit") {
    
    # use the logit model to obtain the dilution for required hit rate (probability)
    intercept <-  fp.l$coefficients[1][[1]]
    beta <-       fp.l$coefficients[2][[1]]
    d.fp.l <- (qlogis(hit) - intercept) / beta
    
  } else if  ( MODEL %in% 'probit') {
    
    #  use the probit model to obtain the dilution for required hit rate (probability)
    intercept <-  fp.l$coefficients[1][[1]]
    beta <-       fp.l$coefficients[2][[1]]
    d.fp.l <- (qnorm(hit) - intercept) / beta
    
  }
  
  # get LoD Dilution Point to predict Cq value for
  dPred <- data.frame(dil  = d.fp.l)
  # Establishing LoD Cq values 
  # Investigating appropriate Models for Cq against dilution; reflecting Non detected observations as censored
  ddist <<- datadist(dd)
  options(datadist='ddist')
  
  fbj <- pj <- NULL
  
  ## <<- assignemnt is not necessary?
  
  if (MODEL2 %in% "Buckley James") {
    
    # invisible(capture.output(   # stops convergence error appearing in Rshiny app
    fbj <- (bj(Surv(CT, count) ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE, link="identity",
               control=list(iter.max =250)))   
    #))
    
  } else if  ( MODEL2 %in% 'Ordinary Least Squares') {
    
    fbj <- (ols(CT ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE))
    
  }
  
  # 
  #### PREDICT LoD Cq Values
  pj <- as.data.frame(cbind(dPred, predicted = predict(fbj, type="lp", newdata=dPred, se.fit=T)))
  pj$lower<-pj[2][[1]] + c(-1) * 1.96*pj[3][[1]]
  pj$upper<-pj[2][[1]] + c( 1) * 1.96*pj[3][[1]]
  names(pj) <-c("dil", "LoD Cq", "SE", "LoD Cq Lower 95% CI", "LoD Cq Upper 95% CI")
  pow1[i,] <-  c(p2(d.fp.l), p2(pj[1,2]), p2(pj[1,4]), p2(pj[1,5]))
  
}

#-----------------------------------------------------------------------------------------
pow1 <- plyr::adply(pow1, 1, c)  # convert to numeric

print(pow1)

warnings()