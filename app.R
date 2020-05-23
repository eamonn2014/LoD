#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rshiny ideas from on https://gallery.shinyapps.io/multi_regression/  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())  
set.seed(333) # reproducible


list.of.packages <- c("directlabels", "ggplot2" , "xtable", "doBy", "VCA", "reshape", "nlme", "vcd","car",
                      "MASS","R2wd","tables","gtools", "rtf", "binom", "coin", 
                      "lmec", "coxme", "lme4",  
                      "arm", "rms", "plyr",
                      "directlabels","shiny","shinyWidgets","shinythemes","DT","shinyalert",
                      "Hmisc","reshape","rms","ggplot2","tidyverse","digest")

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
file <- "https://raw.githubusercontent.com/eamonn2014/LoD/master/LoD_data.txt"

d99 <- read.delim(file)

assa <- unique(d99$assay)
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ui <- fluidPage(theme = shinytheme("journal"), #https://www.rdocumentation.org/packages/shinythemes/versions/1.1.2
                # paper
                useShinyalert(),  # Set up shinyalert
                setBackgroundColor(
                  color = c( "#2171B5", "#F7FBFF"), 
                  gradient = "linear",
                  direction = "bottom"
                ),
                
                h2("Establishing Limits of Detection (LoD) for quantitative polymerase chain reaction (qPCR) assays"), 
                
                h4("
                
The objective of this app is to show an analysis approach that can be used to establish the limit of detection for qPCR assays designed to quantify RNA expression 
from a biological sample when the total amount of RNA
input is constant but the percentage of mutant RNA content is varied [1,2]. An LoD experiment was simulated for nine 
qPCR mutation assays. The outcome of this analysis will identify the concentration of mutant at which 95% of observations are detected, 
however, the Cq value corresponding to the identified RNA mutant
content will be predicted from a fitted model to determine a Cq LoD threshold for each assay."),
                
                h4("
Methods: A study was designed in accordance with CLSI EP17-A2, with the exception we assume a single lot 
of reagents [3]. RNA is derived from three independent pools of cell 
lines that are 'lowly positive' determined by a preliminary experiment to identify the samples for each of the nine assay mutations. For each assay, each biomarker positive pool sample is
then assayed 7 times on each of 3 testing days (run) to generate a total of 21 measurements and a total of 42 Cq values (2 Cq values per measurment) per dilution value. 
The mutant RNA content is subject to a 2-fold serial dilution encompassing 8 dilution points. It can be seen the more dilute, the higher the Cq 
                   value, note the qPCR machine has a technical upper limit of 40 Cq."),
                
                h4("
Data Analysis: Each assay LoD is established using a two step approach; binary regression modelling (logit/probit as in EP17-A2) followed by linear regression modelling. 
The Cq value was used to determine the assay
detection calls based on the assay specific verified or established limit of blank (LoB) thresholds. The hit rate (percentage of detection calls)
was calculated from the total replicates tested at each dilution concentration, and modelled against the dilution using the logit/probit model. 
It should be noted that the logit and probit models are essentially the same. The LoD will be established as the Cq value at which the lowest 
concentration of analyte can be routinely detected; more concretely defined as Cq value corresponding to the dilution concentration at which 95% 
of the observed Cq values for samples tested are below the LoB. 
The corresponding Cq value can be predicted from the fit of the Buckley-James (BJ) censored data 
regression model that takes censoring of unobserved or undetected observations into consideration rather than imputing them as 40 Cq. 
There is an option to use ordinary least squares instead of the BJ model as a comparison.              
"), 
                
              #  h4("test  "), 
                
                
                sidebarLayout(
                  
                  sidebarPanel( width=3 ,
                               # h4("Options for modelling and presentation:"),
                                tags$style(type="text/css", ".span8 .well { background-color: #00FFFF; }"),
                                
                                
                                actionButton(inputId='ab1', label="R Shiny ",   icon = icon("th"),   
                                             onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/LoD/master/app.R', '_blank')"), 
                                actionButton(inputId='ab1', label="R code",   icon = icon("th"),   
                                             onclick ="window.open('https://raw.githubusercontent.com/eamonn2014/LoD/master/LoD.R', '_blank')"),  
                                #actionButton("resample", "Simulate a new sample"),
                                br(),  
                                tags$style(".well {background-color:#b6aebd ;}"), 
                                
                             
                                div(
                                  
                                  tags$head(
                                    tags$style(HTML('#ab1{background-color:orange}'))
                                  ),
                                  
                                  tags$head(
                                    tags$style(HTML('#resample{background-color:orange}'))
                                  ),
                                  
                                  #  textInput('Assay', 
                                  #           div(h5(tags$span(style="color:blue", "Assay"))), "BAIA"),
                                  
                                  selectInput("Assay",
                                              div(h5(tags$span(style="color:blue", "Select assay"))),
                                              choices=assa, selected = assa[2]),
                                  
                                  selectInput("MODEL",
                                              div(h5(tags$span(style="color:blue", "Select modelling approach step 1"))),
                                              choices=c("logit","probit"), selected = "logit"),
                                  
                             
                                  
                              #    tags$hr(),
                                  textInput('LoB', 
                                            div(h5(tags$span(style="color:blue", "Limit of Blank (Cq) used in step 1"))), "40"),
                                  
                                 # tags$hr(),
                                  textInput('Hitrate', 
                                            div(h5(tags$span(style="color:blue", "Hit rate used in step 1"))), "0.95"),
                              textInput('agger', 
                                        div(h5(tags$span(style="color:blue", "Enter 'yes' explicitly to aggregate over pools and runs only affects presentation : figures 1 & 3 in step 1"))), "ye"),
                              
                              selectInput("MODEL2",
                                          div(h5(tags$span(style="color:blue", "Select Cq prediction modelling approach step 2"))),
                                          choices=c("Buckley James","Ordinary Least Squares"), selected = "Buckley James"),
                               #   tags$hr(),
                                  textInput('knots', 
                                            div(h5(tags$span(style="color:blue", "Number of restricted cubic spline knots in Cq prediction model"))), "5"),
                                  
                               #   tags$hr(),
                                  textInput('jitt', 
                                            div(h5(tags$span(style="color:blue", "Enter 'yes' explicitly to add vertical jitter to 40 Cq data points in Cq prediction model plot (only affects presentation)"))), "ye"),
                                  
                                  textInput('jitt1', 
                                            div(h5(tags$span(style="color:blue", "Enter magnitude of jitter for Cq prediction model plot (only affects presentation)"))), "0.3"),
                                  
                               #   tags$hr(),
                         
                                  
                                   div(h5("References:")),  
                              tags$a(href = "https://en.wikipedia.org/wiki/Detection_limit",  tags$span(style="color:blue", "[1] LoD wiki"),),   
                              div(p(" ")),
                                   tags$a(href = "https://en.wikipedia.org/wiki/Real-time_polymerase_chain_reaction", tags$span(style="color:blue", "[2] qPCR wiki"),),   
                                   div(p(" ")),
                                 
                                  tags$a(href = "https://community.clsi.org/media/1430/ep17a2_sample.pdf", tags$span(style="color:blue", "[3] EP17-A2 Guidance"),),
                                   div(p(" ")),
                               
                                )
                                
                                
                  ),
                  
                  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~tab panels
                  mainPanel(width=9,
                            
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            navbarPage(       
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
                              tags$style(HTML("
                            .navbar-default .navbar-brand {color: orange;}
                            .navbar-default .navbar-brand:hover {color: blue;}
                            .navbar { background-color: #b6aebd;}
                            .navbar-default .navbar-nav > li > a {color:black;}
                            .navbar-default .navbar-nav > .active > a,
                            .navbar-default .navbar-nav > .active > a:focus,
                            .navbar-default .navbar-nav > .active > a:hover {color: pink;background-color: purple;}
                            .navbar-default .navbar-nav > li > a:hover {color: black;background-color:yellow;text-decoration:underline;}
                            .navbar-default .navbar-nav > li > a[data-value='t1'] {color: red;background-color: pink;}
                            .navbar-default .navbar-nav > li > a[data-value='t2'] {color: blue;background-color: lightblue;}
                            .navbar-default .navbar-nav > li > a[data-value='t3'] {color: green;background-color: lightgreen;}
                   ")),
                              
                              
                              tabPanel("1 Establishing LoD", value=7, 
                              
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                
                                                div(plotOutput("plot1",  width=fig.width7, height=fig.height7)),
                                                h4(paste("Figure 1. Step 1, Estimating the dilution that satisfies the hit rate")), 
                                         ) ,
                                         
                                         
                                         fluidRow(
                                           column(width = 5, offset = 0, style='padding:1px;',
                                                  
                                                  div(plotOutput("plot3x",  width=fig.width7, height=fig.height7)) ,
                                                  h4(paste("Figure 2. Step 2, Fitted Prediction model, predicting the LoD Cq at the dilution that satisfied the required hit rate.\n Horizontal jitter added to improve visualisation.")),
                                                  
                                                     
                                           ))),
                                      
                                       
                              ) ,
                              
                              tabPanel("2 Buckley James (BJ) or Ols Model", value=3, 
                                       
                                         div( verbatimTextOutput("datx") ),
                                        
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4(paste("Table 1. Model regression table")), 
                                                
                                         )),
                                       
                                       
                              ),
                              
                              tabPanel("3 BJ or Ols Model Predictions", value=3, 
                                       
                                       div( verbatimTextOutput("dat") ),
                                      
                                      
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4(paste("Table 2. Model regression predictions")),
                                                
                                         )),
                                       
                                       
                              ),
                              
                              tabPanel("4 LoD", value=7, 
                                       
                                       fluidRow(
                                         column(width = 6, offset = 0, style='padding:1px;',
                                                
                                                div( verbatimTextOutput("daty") ),
                                                h4(paste("Table 3. Limit of detection estimate")),
                                         ) ,
                                         
                                         
                                         
                                        # h4("Table 2 xxxxxxx"),
                                         fluidRow(
                                           column(width = 6, offset = 0, style='padding:1px;',
                                               
                                           ))),
                                       
                              ) ,
                 
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              tabPanel("5 Probit/Logit Plot", value=3, 
                                       
                                       div(plotOutput("plot4", width=fig.width1, height=fig.height7)),  
                                       
                                       
                                       fluidRow(
                                         column(width = 12, offset = 0, style='padding:1px;',
                                                h4(paste("Figure 3. Finding the dilution that satisfies the hit rate (estimated from either logistic or probit regression)
                                                          ")),
                                                
                                         )),
                                       
                                       
                              ),
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              tabPanel("6 BJ or Ols Model plot", value=3, 
                                       
                                       div(plotOutput("plot3", width=fig.width1, height=fig.height7)),  
                                
                                       
                                       fluidRow(
                                         column(width = 7, offset = 0, style='padding:1px;',
                                                h4(paste("Figure 4. Predicting the LoD Cq at the dilution that satisfied the required hit rate; horizontal jitter added to data points to aid visualisation")), 
                                         )),
                                       
                                       
                              ),
                              
                         
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              tabPanel("7 Conclusion", value=3, 
                                        
                                        h4(paste("We have shown how to establish LoD values for the assays. 
                                        It is advisable to perform summary plots and investigate missing data before the main analysis.
                                        Briefly, for each dilution, score 0 or 1 if Cq is above or below LoB. 
                                        Fit a logit/probit model to the aformentioned hit rate 
                                        information and calculate using the estimated intercept and regression coefficients the
                                        dilution that corresponds to the 'hit rate', that is probability of 
                                        detection, typically 95%. Finally, armed with the dilution estimate, fit a flexible
                                        model to data (one that accounts for censoring preferably and non linearity) and predict 
                                         the Cq and its standard error. 
  Notice the hit rate is 100% for most of the dilution points 
in the case of Assay01 and Assay08. The model parameters are unstable and result 
in occurrence of fitted probabilities being numerically 0 or 1. The near 100% hit rate for most of the dilution points for Assay01 and Assay08 
implies that the experimental data for the two assays do not meet the desirable CLSI EP17 A2 
hit rate criteria. Therefore the LoD estimation would benefit from further dilutions for these two assays. In our example the modelling choices do not substantially affect the final estimates.")), 
                                         

                                        fluidRow(
                                          column(width = 7, offset = 0, style='padding:1px;',
                                              
                                                 
                                          )),
                               ),
                              
                              tabPanel("8 Data", value=3, 
                                       DT::dataTableOutput("table1"),
                                       
                              ) # ,
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                              
                              # tabPanel("9 BJ Model plot alt.", value=3, 
                              #          
                              #          div(plotOutput("plot2", width=fig.width1, height=fig.height7)),  
                              #          
                              #          
                              #          fluidRow(
                              #            column(width = 7, offset = 0, style='padding:1px;',
                              #                   h4(paste("Figure 5. Finding the dilution that satisfies the hit rate")), 
                              #            )),
                              #          
                              #)
                              
                              
                              #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   END NEW   
                            )
                            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                  )
                ) 
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end tab panels 
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- shinyServer(function(input, output   ) {
  
  shinyalert("Welcome! \nEver heard of the Buckley-James model?",
             "Establishing Assay Limit of Detection", 
             type = "info")
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # analyse
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Loddp <- reactive({
 
    assay <- (unlist(strsplit(input$Assay,",")))
    
    MODEL <- (unlist(strsplit(input$MODEL,",")))
    MODEL2 <- (unlist(strsplit(input$MODEL2,",")))
    
    LoB <- as.numeric(unlist(strsplit(input$LoB,",")))
    
    hit <- as.numeric(unlist(strsplit(input$Hitrate,",")))
    
    jitt <- (unlist(strsplit(input$jitt,",")))
    jitt1 <- as.numeric(unlist(strsplit(input$jitt1,",")))
    
    knots <- as.numeric(unlist(strsplit(input$knots,",")))
    
    agger <- (unlist(strsplit(input$agger,",")))
   
 
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
      
      LoD.count <- LoD.count %>%
      group_by(dil) %>%
      summarize(Freq=sum(Freq), N=sum(N))
    
      fp.l <- tryCatch(glm(Freq/N ~ dil, weights=N, data=LoD.count, family=binomial(link=MODEL)), error=function(e) NULL)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    }
    
    # Predicting LoD (dilution point)....I use model coefficients instead see below
    
     # use this info for the logit/probit model plot
     x =seq(1, 8, 0.001)
     data <- tryCatch(data.frame(x , pred.fp.l = predict(fp.l, data.frame(dil=x), type="resp")), error=function(e) NULL)
    
    # find nearest dilution at which prob=0.95 ######ESTIMATE 
    # it<-which.min(abs(data$pred.fp.l - hit))
    # d.fp.l<-data[it,]$x 
    # d.fp.l
    
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
    ## <<- assignemnt is necessary!!

    if (MODEL2 %in% "Buckley James") {
  
    tag <- "BJ"
      
    fbj <<- bj(Surv(CT, count) ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE, link="identity",
              control=list(iter.max =250))
    
    } else if  ( MODEL2 %in% 'Ordinary Least Squares') {

    tag <- "Ols"  
      
    fbj <<- ols(CT ~ rcs(dil,knots) , data=dd, x=TRUE, y=TRUE)
    }
    #--------------------------------------------------------------------------------------------------------------
  
    pj<-as.data.frame(cbind(dPred, predicted = predict(fbj, type="lp", newdata=dPred, se.fit=T)))
    
    pj$lower<-pj[2][[1]] + c(-1) * 1.96*pj[3][[1]]
    pj$upper<-pj[2][[1]] + c(1) * 1.96*pj[3][[1]]
    names(pj) <-c("dil", "LoD Cq", "SE", "LoD Lower 95% CI", "LoD Upper 95% CI")
    # for some reason rms::Predict is not working so I have to use predict see lot 2
    Pre <- rms::Predict(fbj)
    #--------------------------------------------------------------------------------------------------------------

    ### Plot count vs dilution point 
   plot1 <- plot(x=LoD.count$dil, y=LoD.count$Freq/LoD.count$N, pch="o"
         ,xlab=paste0("Dilution Series (LoD dil. series estimate ",p2(d.fp.l),", denoted by red vertical dashed line)"), ylab="Hit Rate"
         , main=paste("Plot of Hit Rates against Dilution Series for \n", assay, 
                      " with Line of Estimated Dilution from the Fitted ",MODEL," Model\n using limit of blank of ", LoB,"Cq", sep= ""), cex.main=1.0, bty='n')
    lines(data$x, data$pred.fp.l, type="l", lwd=1, col="blue")    # the prediction line
    abline(h=hit, lwd=1, col="gray") 
    legend('bottomleft','groups', c(paste0(MODEL, " fit"), paste0(hit*100,"% Hit Rate")),lty = c(1),
           col=c('blue', 'gray'), ncol=1, bty ="n", cex=0.8)
    text(1.75, 0.80, paste("LoD Dil.Series =", p2(d.fp.l), sep=" "), cex = 0.95)
    abline(v=(d.fp.l), lwd=1, col="red", lty=2) 
  
    #--------------------------------------------------------------------------------------------------------------
    #pj <- Pre
    
    pl <-  ggplot(Pre , anova=NULL, pval=FALSE, ylim.=c(20,45)) 
    
    plot2 <- pl + geom_point(data=dd, aes(x = dil, y = CT )) +
      xlab(paste0("Dilution Series (LoD dil. series estimate ",dPred,", denoted by vertical dashed line)")) + ylab("Cq value") +
      geom_hline(yintercept=LoB, colour="gray", linetype="solid") +
      geom_hline(yintercept=pj[1,2], colour="blue", linetype="dashed") +
      ggtitle(paste("Assay ",assay," Plot of Cq vs. Dilution Series with the Fitted",tag,"Model \nLoD Dilution Point ="
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
                    ":",assay,", Cq vs. Dilution Series and fitted",tag,"Model \nLoD Dilution Point ="
                    , p2(d.fp.l), "; LoD Cq = ", p2(pj[1,2]), ", 95%CI (", p2(pj[1,4]), ",", p2(pj[1,5])     ,")", sep=" ")) +
    
      geom_ribbon(data=zz, 
                  aes(x= dil,
                      ymin=lower, ymax=upper, alpha=0.2,fill="purple") ) +
      
      labs(caption = paste0("- Cq values do not exceed machine upper limit of 40\n- LoD dilution series estimate ",p2(d.fp.l),", denoted by vertical dashed line")) +
 
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
    
    #-----------------------------------------------------------------------------------------
    return(list(plot1=plot1, plot2=plot2, fbj=fbj,  pj=pj, dd=dd, plot3=plot3 , Pre=Pre ))  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  })

  #________________________________________________________________________________________________
  #________________________________________________________________________________________________
 
  #https://stackoverflow.com/questions/44205137/plotting-the-same-output-in-two-tabpanels-in-shiny
  output$plot1 <- output$plot4  <- renderPlot({         
   
    x<- Loddp()
    (x$plot1)
    
  })
  
  output$plot2<- renderPlot({         
    
    x<- Loddp()
    plot(x$plot2)
    
  })
  
  output$plot3x <- output$plot3 <- renderPlot({         
    
    x<- Loddp()
    plot(x$plot3)

  })
  

  #________________________________________________________________________________________________
  
  output$dat <- renderPrint({
    
    d <- Loddp()$Pre
    
    return(print(d, digits=4))
  })
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  output$datx <- renderPrint({

    d <- Loddp()$fbj

    return(print(d, digits=4))
  })
  
  output$daty <- renderPrint({
    
    d <- Loddp()$pj
    
    return(print(d, digits=4))
  })
   
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  output$table1 <- DT::renderDataTable({
    
    foo<-  d99
    
    namez <- c("Assay", "Biological Sample", "Sample Pool", "Dilution", "Run", "Replicate", "Cq" )
    
    
    names(foo) <- namez
    rownames(foo) <- NULL
    library(DT)
    
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
  })
})

# Run the application 
shinyApp(ui = ui, server = server)