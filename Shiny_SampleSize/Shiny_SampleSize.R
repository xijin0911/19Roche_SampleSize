# **********************************
# Required sample size based on two methods by Riley Richard and van Smeden
# **********************************

library(shiny)
library(ggplot2)
library(markdown)
library(speff2trial)
library(shinyjs)
library(ggpubr)
library(rmarkdown)
#library(scatterplot3d)
#library(rgl)
#library("plot3D")
library(plotly)
#library(tidyr)


server <- function(input, output,session)
{
  output$currentTime <- renderText({
    invalidateLater(1000, session)
    paste( Sys.time())
  })
 
  
  ### functions (prediction error)
  MSPE <- function(ef,P,N){
    exp(-0.59+-1.06*log(N)+0.36*log(ef)+0.94*log(P))
  }
  MAPE <- function(ef,P,N){
    exp(-0.48+-0.53*log(N)+0.31*log(ef)+0.47*log(P))
  }
  ### functions (prediction error)
##########################################
######Riley###############################
##########################################
  
# Outcome of three sample sizes corresponding to three criterion and the smallest one among those three
# either by option 1 or option 2
  output$Size1 <- renderText({
    
      n1 <- input$p/((input$S_VH-1)*log(1-input$csR2_adj/input$S_VH))
      return(ceiling(n1))
  })
    
  output$Size2 <- renderText({
       max_csR2 <- 1- (input$Ef^input$Ef*(1-input$Ef)^(1-input$Ef))^2
        S_VH <- input$csR2_adj/(input$csR2_adj+input$d*max_csR2)
        n2 <- input$p/((S_VH-1)*log(1-input$csR2_adj/S_VH))
        return(ceiling(n2))
    })
    
  output$Size3 <- renderText({
        n3 <- ((1.96/input$error_mar)^2)*(input$Ef)*(1-(input$Ef))
        return(ceiling(n3))
    })
    
  output$n <- renderText({
        n1 <- input$p/((input$S_VH-1)*log(1-input$csR2_adj/input$S_VH))
        max_csR2 <- 1- (input$Ef^input$Ef*(1-input$Ef)^(1-input$Ef))^2
        S_VH <- input$csR2_adj/(input$csR2_adj+input$d*max_csR2)
        n2<- input$p/((S_VH-1)*log(1-input$csR2_adj/S_VH))
        n3<- (1.96/input$error_mar)^2*(input$Ef)*(1-(input$Ef))
        n <- max(ceiling(n1),
                 ceiling(n2),
                 ceiling(n3))
        return(n)
    })
  
  

      csR <-  eventReactive(input$goButton,{
        switch(input$Reported_information,
               "LR"= 1-exp(-input$corresponds_stat/input$n), 
               "Nagelkerke" = input$corresponds_stat*(1-exp((input$n*input$Ef*log(input$Ef)+(input$n-input$n*input$Ef)*log(1-input$Ef))*2/input$n)),
              "McFadden" = 1-exp(- (-2*((input$n*input$Ef*log(input$Ef)+(input$n-input$n*input$Ef)*log(1-input$Ef)) - (input$n*input$Ef*log(input$Ef)+(input$n-input$n*input$Ef)*log(1-input$Ef))*(1- input$corresponds_stat))) /input$n),
               "O'Quigley" = 1-exp(log((1-input$corresponds_stat)*log(input$n*input$Ef))/input$n)
        )
      })
    output$adjR2 <- renderPrint({
      csR()
    }) 
    
    csR_nagel <-  eventReactive(input$goButton2,{
      switch(input$NagelR2,
             "0.15"= 0.15*(1-((input$Ef^input$Ef)*(1-input$Ef)^(1-input$Ef))^2), 
             "0.5" = 0.5*(1-((input$Ef^input$Ef)*(1-input$Ef)^(1-input$Ef))^2 )
      )
    })
    output$adjR2_nagel <- renderPrint({
      csR_nagel()
    }) 
    
   #plot for three criteria by Riley method
    output$plot_Riley <- renderPlot({
        R <- seq(from=0.05,to=0.5,length=500)     #adjusted cox snell R2 as x values
        y <- ceiling(seq(from=10,length=500))
        d <- data.frame(x=R, y=y)
        max_csR2 <- 1- (input$Ef^input$Ef*(1-input$Ef)^(1-input$Ef))^2
        S <- R/(R+input$d*max_csR2)  # S is calculated by criterion 2
        
        a <- ggplot(d, aes(x, y)) + #geom_line(color="blue", size=2) +
            xlab("Adjusted Cox-Snell R-squared") + ylab("Sample Size") +
            ggtitle("Sensitivity of three criteria calculated size on Adjusted Cox-Snell R-squared") +
            theme(legend.background = element_rect(fill="white", 
                                                 size=0.5, linetype="solid"))+
            theme(text=element_text(size=20)) +
             scale_y_log10() +
          geom_line(aes(x,y=input$p/((input$S_VH-1)*log(1-x/input$S_VH)),colour = "1. Relative drop"))+
          geom_line(aes(x,y= input$p/((S-1)*log(1-x/S)),colour = "2. Absolute difference"))+
          geom_line(aes(x,y=(1.96/input$error_mar)^2*(input$Ef)*(1-(input$Ef)), 
                          colour = "3. Precise estimation"))+
           scale_colour_manual("Three criteria",values=c('1. Relative drop'="red",
                                                        '2. Absolute difference'="black",
                                                        '3. Precise estimation'="blue"))+
          geom_vline(xintercept = 0.15/max_csR2,linetype = "dashed")
        print(a)
    })
    
    library(plotly)
    output$threedplot <- renderPlotly({
      S_VH <- seq(0.9,0.95,length=20)   #z-axis
      csR2_adj <- seq(0.001,0.5,length=50)  #x-axis
      d <- c(0.1,0.075,0.05) #color
      
      n <- rep(0,length(S_VH)*length(csR2_adj)*length(d))
      n1 <- n
      n2 <- n
      n3<- n
      S <- n # shrinkage factor calculated by criterion 2
      
      dat <- data.frame(cbind(
        csR2_adj=rep(rep(csR2_adj,length(S_VH)),length(d)),
        S_VH=rep(rep(S_VH,each=length(csR2_adj)),length(d)),
        logn=n,
        d=rep(d,each=length(S_VH)*length(csR2_adj))))
      
      for (i in 1:nrow(dat)){
        n1[i] <- input$p/((dat$S_VH[i]-1)*log(1-dat$csR2_adj[i]/dat$S_VH[i]))
        max_csR2 <- 1- ( input$Ef^ input$Ef*(1- input$Ef)^(1- input$Ef))^2
        S[i] <-   dat$csR2_adj[i]/(dat$csR2_adj[i]+dat$d[i]* dat$csR2_adj[i])
        n2[i]<-  input$p/(( S[i]-1)*log(1-  dat$csR2_adj[i]/ S[i]))
        n3[i]<- (1.96/ input$error_mar)^2*( input$Ef)*(1-( input$Ef))
        n[i] <- max(ceiling(n1[i]),
                    ceiling(n2[i]),
                    ceiling(n3[i]))
      }
      dat$logn <- log(n)
      dat$d <- as.factor(dat$d)
      plot_ly(x = dat$csR2_adj, y = dat$S_VH, z = dat$logn, type = "scatter3d",color = dat$d)
      
#      scatter3d(dat[,1:3], pch = 16, groups=dat$d,surface=FALSE,fit = "smooth",surface = FALSE)
 #     rglwidget()
    })
    
    # Dynamic report download
    output$Report_Riley <- downloadHandler(
      filename = function() {
        paste('Report-Three-Criteria', sep = '.', switch(
          input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
        ))
      },  
      content = function(file) {
        src <- normalizePath('Report_Riley.Rmd')
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        file.copy(src, 'Report_Riley.Rmd', overwrite = TRUE)
        params <- list(Ef = input$Ef,
                       p =  input$p, 
                       csR2_adj = input$csR2_adj,
                       S_VH = input$S_VH,
                       d = input$d,
                       error_mar = input$error_mar
        )
        
        out <- render('Report_Riley.Rmd', switch(
          input$format,
          PDF = pdf_document(), HTML = html_document(), Word = word_document()
        ),params = params,
        envir = new.env(parent = globalenv()))
        file.rename(out, file)
      }
    )
    
    output$Report_Smeden <- downloadHandler(
      filename = function() {
      paste('Report-Metamodel', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
      },  
      content = function(file) {
        src <- normalizePath('Report_Smeden.Rmd')
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        file.copy(src, 'Report_Smeden.Rmd', overwrite = TRUE)
      params <- list(Ef_van = input$Ef_van,
                     p_van =  input$p_van,
                    Modeling = input$modelling_strategies,
                     MAPE = input$MAPE,
                       MSPE = input$MSPE,
                     Brier = input$Brier
    
        )
        
        out <- render('Report_Smeden.Rmd', switch(
          input$format,
          PDF = pdf_document(), HTML = html_document(), Word = word_document()
        ),params = params,
        envir = new.env(parent = globalenv()))
        file.rename(out, file)
      }
    )
                                                                                                 
  
    ##########################################
    ######van Smeden##########################
    ##########################################
    
    output$plot_vanSmeden_ML <- renderPlot({
      metrics <- seq(from=0.01,to=0.5,length=100)
      size_smeden <- ceiling(seq(from=300,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
         xlab("MSPE") + 
        xlim(0.01,0.5) +
        ylim(0,250)+
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.59 - log(x) + 0.36*log((input$Ef_van)) + 0.94*log(input$p_van))/1.06),colour="blue"))+
        ylab("Sample Size") +
        theme(text=element_text(size=20))
      b <-  ggplot(dat, aes(x, y))+
        xlab(" MAPE") + 
        xlim(0.01,0.5) +
        ylim(0,300)+
        geom_line(aes(x,y=exp((-0.48 - log(x) + 0.31*log((input$Ef_van)) + 0.48*log(input$p_van))/0.53),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")+
        theme(text=element_text(size=20))
      
      c <-  ggplot(dat, aes(x, y))+
        xlab("Brier") + 
        xlim(0.125,0.155) +
        ylim(0,300)+
        geom_line(aes(x,y=exp((-0.91 - log(x) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")+
        theme(text=element_text(size=20))
        ggarrange(a,b,c,
                 # labels = c("MSPE","MAPE","Brier"),
        ncol = 3, nrow = 1)
    })
    
    
    output$plot_vanSmeden_Firth <- renderPlot({
      metrics <- seq(from=0.01,to=0.5,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
         xlab(" MSPE") + 
        xlim(0.01,0.5) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.86 - log(x) + 0.33*log((input$Ef_van)) + 0.93*log(input$p_van))/1.03),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab("MAPE") + 
        xlim(0.01,0.5) +
        scale_y_log10(limits = c(2,500)) + #,
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.61 - log(x) + 0.29*log((input$Ef_van)) + 0.47*log(input$p_van))/0.51),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab("Brier") + 
        xlim(0.125,0.15) +
        scale_y_log10(limits = c(2,500)) + #,
        geom_line(aes(x,y=exp((-0.95 - log(x) + 0.62*log((input$Ef_van)) + 0.03*log(input$p_van))/0.03),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    output$plot_vanSmeden_HS <- renderPlot({
      metrics <- seq(from=0.01,to=0.8,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab("MSPE") + 
        xlim(0.01,0.8)+
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.75 - log(x) + 0.44*log((input$Ef_van)) + 0.74*log(input$p_van))/0.49),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab(" MAPE") + 
        xlim(0.55,0.75) +
        scale_y_log10(limits = c(2,500)) + #,
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.56 - log(x) + 0.33*log((input$Ef_van)) + 0.39*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab(" Brier") + 
        xlim(0.125,0.15) +
        scale_y_log10(limits = c(2,500)) + #,
        geom_line(aes(x,y=exp((-0.93 - log(x) + 0.62*log((input$Ef_van)) + 0.02*log(input$p_van))/0.03),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    output$plot_vanSmeden_Lasso <- renderPlot({
      metrics <- seq(from=0.01,to=0.7,length=100)
      size_smeden <- ceiling(seq(from=3000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab(" MSPE") + 
        xlim(0.01,1.2) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y= exp((-0.86 - log(x) + 0.46*log((input$Ef_van)) + 0.68*log(input$p_van))/0.48),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab(" MAPE") + 
        xlim(0.5,0.75) +
        scale_y_log10(limits = c(2,500)) + #,
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.59 - log(x) + 0.34*log((input$Ef_van)) + 0.35*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab("Brier") + 
        xlim(0.12,0.145) +
        scale_y_log10(limits = c(2,500)) + #,
        geom_line(aes(x,y=exp((-0.96 - log(x) + 0.62*log((input$Ef_van)) + 0.02*log(input$p_van))/0.03),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
              #  labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    
    output$plot_vanSmeden_Ridge <- renderPlot({
      metrics <- seq(from=0.01,to=0.5,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab("MSPE") + 
        xlim(0.01,0.5) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y= exp((-0.93 - log(x) + 0.50*log((input$Ef_van)) + 0.49*log(input$p_van))/0.45),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab("MAPE") + 
        xlim(0.4,0.5) +
        scale_y_log10(limits = c(2,500)) + #,
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.61 - log(x) + 0.36*log((input$Ef_van)) + 0.26*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab("Brier") + 
        xlim(0.12,0.145) +
        scale_y_log10(limits = c(2,500)) + #,
        geom_line(aes(x,y=exp((-0.98 - log(x) + 0.62*log((input$Ef_van)) + 0.01*log(input$p_van))/0.02),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
      
    })
    
    
    output$plot_vanSmeden_ML_p <- renderPlot({
      metrics <- seq(from=0.01,to=1.0,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab(" MSPE") + 
        xlim(0.01,1.2) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.57 - log(x) + 0.40*log((input$Ef_van)) + 0.96*log(input$p_van))/0.52),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab(" MAPE") + 
        xlim(0.8,1) +
        scale_y_log10(limits = c(2,500)) + #,
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.45 - log(x) + 0.31*log((input$Ef_van)) + 0.49*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab(" Brier") + 
        xlim(0.8,1) +
        scale_y_log10(limits = c(2,500)) + #,
        geom_line(aes(x,y=exp((-0.45 - log(x) + 0.31*log((input$Ef_van)) + 0.49*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    
    output$plot_vanSmeden_ML_aic <- renderPlot({
      metrics <- seq(from=0.01,to=1.0,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab(" MSPE") + 
        xlim(0.01,0.5) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.59 - log(x) + 0.38*log((input$Ef_van)) + 0.95*log(input$p_van))/0.53),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab("MAPE") + 
        xlim(0.8,1) +
        scale_y_log10(limits = c(2,500)) + 
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.48 - log(x) + 0.31*log((input$Ef_van)) + 0.49*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab("Brier") + 
        xlim(0.13,0.14) +
        scale_y_log10(limits = c(2,500)) + 
        geom_line(aes(x,y=exp((-0.91 - log(x) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
  
    output$plot_vanSmeden_Firth_p <- renderPlot({
      metrics <- seq(from=0.01,to=1.0,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab(" MSPE") + 
        xlim(0.01,1.0) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.66 - log(x) + 0.39*log((input$Ef_van)) + 0.95*log(input$p_van))/0.52),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab(" MAPE") + 
        xlim(0.8,1.0) +
        scale_y_log10(limits = c(2,500)) + 
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.49 - log(x) + 0.30*log((input$Ef_van)) + 0.50*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab(" Brier") + 
        xlim(0.125,0.155) +
        scale_y_log10(limits = c(2,500)) + 
        geom_line(aes(x,y=exp((-0.90 - log(x) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    output$plot_vanSmeden_Firth_aic <- renderPlot({
      metrics <- seq(from=0.01,to=1.0,length=100)
      size_smeden <- ceiling(seq(from=4000,to=10 ,length=100))
      dat <- data.frame(x=metrics, y=size_smeden)
      a <-  ggplot(dat, aes(x, y))+
        xlab("MSPE") + 
        xlim(0.01,1.0) +
        scale_y_log10(limits = c(2,500)) +
        theme(legend.position = "none")+
        geom_line(aes(x,y= exp((-0.74 - log(x) + 0.37*log((input$Ef_van)) + 0.95*log(input$p_van))/0.52),colour="MSPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      
      b <-  ggplot(dat, aes(x, y))+
        xlab("MAPE") + 
        xlim(0.75,1.0) +
        scale_y_log10(limits = c(2,500)) + 
        theme(legend.position = "none")+
        geom_line(aes(x,y=exp((-0.55 - log(x) + 0.30*log((input$Ef_van)) + 0.49*log(input$p_van))/0.04),colour="MAPE"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      c <-  ggplot(dat, aes(x, y))+
        xlab(" Brier") + 
        xlim(0.205,0.25) +
        scale_y_log10(limits = c(2,500)) + 
        geom_line(aes(x,y=exp((-0.93 - log(x) + 0.31*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04),colour="Brier"))+
        ylab("Sample Size") +
        theme(legend.position = "none")
      ggarrange(a,b,c,
               # labels = c("MSPE","MAPE","Brier"),
                ncol = 3, nrow = 1)
    })
    
    
    #output of meta-model depends on modelling strategies by van Smeden
    strategy <- reactive({
      switch(input$modelling_strategies,
             "ML: Maximum likelihood (full model)" = 
               rbind("Size by MSPE" = ceiling (exp((-0.59 - log(input$MSPE) + 0.36*log((input$Ef_van)) + 0.94*log(input$p_van))/1.06)),
                          "Size by MAPE" = ceiling(exp((-0.48 - log(input$MAPE) + 0.31*log((input$Ef_van)) + 0.48*log(input$p_van))/0.53)),
                          "Size by Brier"  = ceiling (exp((-0.91 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04))
             ),
             "Firth: Firth's penalized likelihood (full model)" = 
               rbind("Size by MSPE" = ceiling(exp((-0.86 - log(input$MSPE) + 0.33*log((input$Ef_van)) + 0.93*log(input$p_van))/1.03)),
                             "Size by MAPE" =ceiling (exp((-0.61 - log(input$MAPE) + 0.29*log((input$Ef_van)) + 0.47*log(input$p_van))/0.51)),
                             "Size by Brier"  =  ceiling(exp((-0.95 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.03*log(input$p_van))/0.03))
             ),
             "HS: Heuistic shrinkage" = 
               rbind("Size by MSPE" =ceiling( exp((-0.75 - log(input$MSPE) + 0.44*log((input$Ef_van)) + 0.74*log(input$p_van))/0.97)),
                          "Size by MAPE" = ceiling(exp((-0.56 - log(input$MAPE) + 0.33*log((input$Ef_van)) + 0.39*log(input$p_van))/0.49)),
                          "Size by Brier"  =  ceiling(exp((-0.93 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.02*log(input$p_van))/0.03))
             ),
             "Lasso: Lasso penalized likelihood" = 
               rbind("Size by MSPE" = ceiling(exp((-0.86 - log(input$MSPE) + 0.46*log((input$Ef_van)) + 0.68*log(input$p_van))/0.93)),
                             "Size by MAPE" = ceiling(exp((-0.59 - log(input$MAPE) + 0.34*log((input$Ef_van)) + 0.35*log(input$p_van))/0.48)),
                             "Size by Brier"  = ceiling (exp((-0.96 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.02*log(input$p_van))/0.03))
             ),
             "Ridge: Ridge penalized likelihood" = 
               rbind("Size by MSPE" = ceiling(exp((-0.93 - log(input$MSPE) + 0.50*log((input$Ef_van)) + 0.49*log(input$p_van))/0.88)),
                             "Size by MAPE" =ceiling (exp((-0.61 - log(input$MAPE) + 0.36*log((input$Ef_van)) + 0.26*log(input$p_van))/0.45)),
                             "Size by Brier"  =ceiling ( exp((-0.98 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.01*log(input$p_van))/0.02))
             ),
             "ML_p: Maximum likelihood (backward 1)" = 
               rbind("Size by MSPE" = ceiling(exp((-0.57 - log(input$MSPE) + 0.40*log((input$Ef_van)) + 0.96*log(input$p_van))/1.03)),
                            "Size by MAPE" =ceiling(exp((-0.45 - log(input$MAPE) + 0.31*log((input$Ef_van)) + 0.49*log(input$p_van))/0.52)),
                            "Size by Brier"  = ceiling (exp((-0.89 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04))
             ),
             "ML_aic: Maximum likelihood (backward 2)" = 
               rbind("Size by MSPE" = ceiling(exp((-0.59 - log(input$MSPE) + 0.38*log((input$Ef_van)) + 0.95*log(input$p_van))/1.05)),
                              "Size by MAPE" =ceiling (exp((-0.48 - log(input$MAPE) + 0.31*log((input$Ef_van)) + 0.49*log(input$p_van))/0.53)),
                              "Size by Brier"  = ceiling (exp((-0.91 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04))
             ),
             "Firth_p: Firth's penalized likelihood (backward 1)" = 
               rbind("Size by MSPE" = ceiling(exp((-0.66 - log(input$MSPE) + 0.39*log((input$Ef_van)) + 0.95*log(input$p_van))/1.01)),
                               "Size by MAPE" =ceiling (exp((-0.49 - log(input$MAPE) + 0.30*log((input$Ef_van)) + 0.50*log(input$p_van))/0.52)),
                               "Size by Brier"  = ceiling( exp((-0.90 - log(input$Brier) + 0.62*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04))
             ),
             "Firth_aic: Firth's penalized likelihood (backward 2)" = 
               rbind("Size by MSPE" =ceiling (exp((-0.74 - log(input$MSPE) + 0.37*log((input$Ef_van)) + 0.95*log(input$p_van))/1.03)),
                                 "Size by MAPE" =ceiling (exp((-0.55 - log(input$MAPE) + 0.30*log((input$Ef_van)) + 0.49*log(input$p_van))/0.52)),
                                 "Size by Brier"  = ceiling (exp((-0.93 - log(input$Brier) + 0.31*log((input$Ef_van)) + 0.04*log(input$p_van))/0.04))
             )
      )
    })
    output$Smeden <- renderTable( digits = NULL,{
      table <- strategy()
      table <- cbind(c("MSPE","MAPE","Brier"),table)
      colnames(table) <- c("Predictive error metrics","Required sample size")
      table
    })
    
    
    
    #explanation of three criteria explanation
    observeEvent(input$btn1, {
      show("element1")
    })
    observeEvent(input$btn2, {
      show("element2")
    })
    observeEvent(input$btn3, {
      show("element3")
    })
  #  showModal(modalDialog(
  #    HTML('<img src="https://upload.wikimedia.org/wikipedia/commons/1/17/Roche_Logo.svg">'),
      #HTML('<img src="https://en.wikipedia.org/wiki/Hoffmann-La_Roche#/media/File:Hoffmann-La_Roche_logo.svg">'),
  ##     easyClose = TRUE,
  #    footer = NULL
  #  ))
    
    ##########################################
    ###### Guidence ##########################
    ##########################################    
    
    
    
    output$metrics_value <- renderTable({
      max_csR2 <- 1- (input$Ef_compare^input$Ef_compare*(1-input$Ef_compare)^(1-input$Ef_compare))^2
      S <- input$csR2_compare/(input$csR2_compare+input$d_compare*max_csR2) 
      s1 <- input$p_compare/((input$S_VH_compare-1)*log(1-input$csR2_compare/input$S_VH_compare))
      s2 <- input$p_compare/((S-1)*log(1-input$csR2_compare/S))
      s3 <- (1.96/input$error_mar_compare)^2*(input$Ef_compare)*(1-(input$Ef_compare))
      Riley_size  <- apply(cbind(s1,s2,s3),1,max)
      mspe_Riley_size <- MSPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size)
      mape_Riley_size <- MAPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size)
      cbind(MSPE=mspe_Riley_size,MAPE=mape_Riley_size)
    },digits = 6)
    
    
    output$plot_Compare <- renderPlot({
      max_csR2 <- 1- (input$Ef_compare^input$Ef_compare*(1-input$Ef_compare)^(1-input$Ef_compare))^2
      S <- input$csR2_compare/(input$csR2_compare+input$d_compare*max_csR2) 
      s1 <- input$p_compare/((input$S_VH_compare-1)*log(1-input$csR2_compare/input$S_VH_compare))
      s2 <- input$p_compare/((S-1)*log(1-input$csR2_compare/S))
      s3 <- (1.96/input$error_mar_compare)^2*(input$Ef_compare)*(1-(input$Ef_compare))
      Riley_size  <- apply(cbind(s1,s2,s3),1,max)
      mspe_Riley_size_show <- round(MSPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size),6)
      mape_Riley_size_show <- round(MAPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size),6)
      
      
      R <- seq(from=0,to=0.3,length=500) 
      
      ###computation of Sample size by Riley three-criteria method based on inputs
      max_csR2 <- 1- (input$Ef_compare^input$Ef_compare*(1-input$Ef_compare)^(1-input$Ef_compare))^2
      S <- R/(R+input$d_compare*max_csR2) 
      s1 <- input$p_compare/((input$S_VH_compare-1)*log(1-R/input$S_VH_compare))
      s2 <- input$p_compare/((S-1)*log(1-R/S))
      s3 <- (1.96/input$error_mar_compare)^2*(input$Ef_compare)*(1-(input$Ef_compare))
      Riley_size  <- apply(cbind(s1[1:length(R)],s2[1:length(R)],rep(s3,length(R))),1,max)
      
      mspe_Riley_size <- MSPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size)
      mape_Riley_size <- MAPE(ef=input$Ef_compare,P=input$p_compare,N=Riley_size)
      
      y_metrics <- seq(0,1,length=500)
      metrics_plot <- ggplot(data.frame(x=R, y=y_metrics),aes(x, y))+
        xlab("Anticipated adjusted Cox-Snell R-squared") + ylab("Prediction error metrics values") +
        ggtitle("Prediction Error for consistency based on Adjusted Cox-Snell R-squared") +
        theme(legend.background = element_rect(fill="white", 
                                               size=0.5, linetype="solid"))+
        theme(text=element_text(size=20)) +
        geom_line(aes(x,y=mspe_Riley_size,colour = "MSPE"))+
        geom_line(aes(x,y=mape_Riley_size,colour = "MAPE"))+
        scale_colour_manual("Prediction error metrics",values=c('MSPE'="red",
                                                                'MAPE'="black"))+
        geom_vline(xintercept = input$csR2_compare,linetype = "dashed")+
        annotate("text", x = input$csR2_compare-0.05, y = mspe_Riley_size_show, label = mspe_Riley_size_show)+
        annotate("text", x = input$csR2_compare-0.05, y = mape_Riley_size_show, label = mape_Riley_size_show)
      
      print(metrics_plot)
      
      
    })
    
    
}


ui <- fluidPage(
    headerPanel("Sample Size computation for a prediction model (Logistic Regression)"),
    withMathJax(),
    tabsetPanel(
      ##########################################
      ######Riley###############################
      ##########################################
        tabPanel("Riley",
                 sidebarLayout(
                     sidebarPanel(
                         fluidRow(
                             
                                    wellPanel(
                                      h3("Data characteristics"),
                             sliderInput(inputId="p", "Number of Predictors",
                                         value=3,min=0,max=30,
                                         animate = TRUE),
                             sliderInput(inputId="Ef", "Events fraction",
                                         min = 0, max = 0.5,
                                         value = 0.25,animate = TRUE)
                             )
                             ,
                          
                                  wellPanel(
                                    h3("Specification of desired 
                                       values"),
                                    sliderInput("S_VH", "Desired Shrinkage Factor:",
                                                min = 0.9, max = 1,
                                                value = 0.9),
                                    sliderInput("d", "Absolute Difference",
                                                min = 0.01, max = 0.1,
                                                value = 0.05,animate = TRUE),
                                    sliderInput("error_mar", "Margin of Error",
                                                min = 0.05, max = 0.1,
                                                value = 0.05,animate = TRUE),
                                    sliderInput("csR2_adj", "Adjusted Cox-Snell R-squared (prior or calculated below):",
                                                min = 0, max = 0.5,
                                                value = 0.30,animate = TRUE))
                                    ,
                         
                                 wellPanel(h3("Specification of adjusted Cox-Snell R-squared"),
                                           p("Calculated adjusted Cox-Snell R-squared value is", em("signal to noise ratio"), "reflect data characteristics"),
                                           wellPanel(p(code("With Reported Information")),
                                                     selectInput(inputId = "Reported_information",
                                                                 label = "Reported Information",
                                                                 choices = c(
                                                                   "LR","Nagelkerke","McFadden","McFadden","O'Quigley")),
                                                     numericInput("corresponds_stat", "Corresponding Statistical Value:",10,step = 0.001),
                                                     p("You would also the total number of patients if you want to use those reported information"),
                                                     numericInput(inputId="n", "Number of Patients", 100),
                                                     br(),
                                                     actionButton("goButton", "Go!"),
                                                     verbatimTextOutput("adjR2"),
                                                     div(strong("Note:"), "Reported information could be from related studies", style = "color:blue")                                  
                                           ),
                                           wellPanel(p(code("With Recommended Value"),"(vertical dashed line)"),
                                                     div(strong("Note:"), "With no information avaliable, there is indirect recommendation by $$R^2_{Nagelkerke}=0.15$$ ", style = "color:blue"),
                                                     div("Recommended value in some special cases (predictors include 'direct' measurements) is 0.5", style = "color:blue"),
                                                     selectInput(inputId = "NagelR2",
                                                                 label = "Recommendation",
                                                                 choices = c(0.15,0.5)),
                                                     br(),
                                                     actionButton("goButton2", "Go!"),
                                                     verbatimTextOutput("adjR2_nagel")
                                           )
                                           ),
                                          
                                   p("Calculated adjusted Cox-Snell R-squared value by reported information if there is no prior information about it."),
                                    p("Click the button to and then put the output value in the slider for anticipated Cox-Snell R-squared.")
                                  )
                         ),
                     mainPanel(
                       wellPanel(
                         # Set up shinyjs
                         useShinyjs(),
                         h3(helpText('Sample sizes calculated based on three criteria')),
                         div(style="display:inline-block",p("Sample size with criterion 1 (Relative drop)")),
                       actionButton("btn1", "Know more!"),
                         hidden(
                           div(style="color:orange",id = "element1", " reduce the distance between adjusted and apparent model performance 
                               by desired shrinkage factor, $$S_{VH}=1+\\frac{p}{nln(1-R^2_{CS\\_{app}})}.$$ ")
                         ),
                         
                         textOutput(outputId="Size1"),   #checkboxInput("chk", label = "NULL: ", value = T),
                         tags$br(),
                         div(style="display:inline-block",p("Sample size with criterion 2 (Absolute difference)")),
                      actionButton("btn2","Know more!"),
                      hidden(
                        div(style="color:orange",id = "element2", "reduce the distance between adjusted and apparent model performance by absolute difference of R-squared values,
                             $$R^2_{Nagelkerkes\\_app}-R^2_{Nagelkerkes\\_adj} \\le \\delta.$$")
                      ),   
                      textOutput(outputId="Size2"),
                         tags$br(),
                         div(style="display:inline-block",p("Sample size with criterion 3 (Precise estimation)")),
                      actionButton("btn3", "Know more!"),
                      hidden(
                        div(style="color:orange",id = "element3", "ensure estimation of outcome proportion in a limit of margin of error,
                             $$1.96\\sqrt{\\frac{\\hat\\phi(1-\\hat\\phi)}{n}}\\le\\delta.$$ ")
                      ),   
                      textOutput(outputId="Size3"),
                         tags$br(),
                         div("Smallest sample size requried by all three criteria", style = "color:blue"),
                         textOutput(outputId="n"),
                      div(strong(strong("Note:"), ""), "If the calculated sample size is not achievable, please do data reduction at first to reduce the number of predictors.", style = "color:blue")
                       ),
                       wellPanel(plotOutput("plot_Riley"),
                                 p(strong(strong("Note:"), ""),"The vertical dashed line indicated the recommended Cox-Snell R-squared value when there is no avaiable information."),
                                 p(HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'),HTML('&nbsp;'), "This recommended value is given Nagelkerke R-squred (0.15) in general case."),
                                 h3("Sensitivity of final sample size on desired values"),
                                 plotlyOutput("threedplot",  width = 800, height = 600),
                                 p("x: Cox-Snell R-squared"),
                                 p("y: Shrinkage factor"),
                                 p("z: sample size on log scale"),
                                 p("Three colors correspond to three values for absolute difference")
                                 ),
                      wellPanel(
                        p("This report is based on your input of your data set"),
                        radioButtons('format', 'Report', c('PDF', 'HTML', 'Word'),
                                     inline = TRUE),
                        downloadButton('Report_Riley') ,
                        p("If you want to know more about this method, there is a document about it in detail"),
                        a("click on me",target="_blank",href="Document_Riley.pdf")
                      )
                     )
                 )
                 
        ),
      ##########################################
      ######van Smeden##########################
      ##########################################
        tabPanel("vanSmeden",
                 sidebarLayout(
                   sidebarPanel(
                     wellPanel(h3("Data characteristics"),
                     sliderInput(inputId="p_van", "Number of Predictors",
                                 min = 4, max = 12,step = 1,
                                 value = 8,animate = TRUE),
                     sliderInput(inputId="Ef_van", "Events Fraction",
                                 min = 1/16, max = 1/2,step = 0.001,
                                 value = 0.2,animate = TRUE)),
                     wellPanel(h3("Specification of desired values"),
                               sliderInput("MSPE", "MSPE: Mean squared prediction error, y denotes the true probability,
                                          $$\\frac{1}{n}\\sum(y_{true}−\\hat{p})^2$$",
                                           min = 0, max = 0.5,step = 0.0001,
                                           value = 0.02,
                                           animate = TRUE),
                               sliderInput("MAPE", "MAPE: Mean absolute prediction error, y denotes the true probability, 
                                           $$\\frac{1}{n}\\lvert(y_{true}-\\hat{p})\\lvert$$",
                                           min = 0, max = 0.5,step = 0.0001,
                                           value = 0.35,
                                           animate = TRUE),
                               sliderInput("Brier", "Brier: Brier score, y denotes the binary observations, 
                                              $$\\frac{1}{n}\\sum(y_{obs}-\\hat{p})^2$$",
                                           min = 0, max = 0.5,step = 0.0001,
                                           value = 0.125,
                                           animate = TRUE),
                               div(style="color:red", " Brier score is not recommended by van Smeden even though it could be explained by data characteristics as other two metrics ")
                     ),
                    
                     wellPanel(
                     selectInput(inputId = "modelling_strategies",
                                 label = "Modelling Strategies",
                                 choices = c("ML: Maximum likelihood (full model)",
                                             "Firth: Firth's penalized likelihood (full model)",
                                             "HS: Heuistic shrinkage",
                                             "Lasso: Lasso penalized likelihood",
                                             "Ridge: Ridge penalized likelihood",
                                             "ML_p: Maximum likelihood (bwd 1)",
                                             "ML_AIC: Maximum likelihood (bwd 2)",
                                             "Firth_p: Firth's penalized likelihood (bwd 1)",
                                             "Firth_AIC: Firth's penalized likelihood (bwd 2)")),
                     h6("Two backward elimination predictor selection:"),
                     p("bwd1: conventional p=.05 as stopping rule "),
                     p("bwd2: conventional p=.157 (AIC) as stopping rule "))
                     
                
                   ),
                   mainPanel(
                     fluidRow(
                       column(10,
                              wellPanel(
                              h3(helpText('Sample sizes calculated based on three predictive error metrics seperately')),
                              tableOutput(outputId="Smeden"),
                              div(strong("Note:"), " please take the one (or ones) that satisfy your requirement for predictive error",style = "color:blue")),
                              tags$br(),
                              wellPanel(
                                p("Abbreviations corresponding to 9 modeling strategies on the left side"),
                                tabsetPanel(type = "tabs",
                                            tabPanel("ML", plotOutput("plot_vanSmeden_ML")),
                                            tabPanel("Firth", plotOutput("plot_vanSmeden_Firth")),
                                            tabPanel("Lasso", plotOutput("plot_vanSmeden_Lasso")),
                                            tabPanel("Ridge", plotOutput("plot_vanSmeden_Ridge")),
                                            tabPanel("HS", plotOutput("plot_vanSmeden_HS")),
                                            tabPanel("ML_p", plotOutput("plot_vanSmeden_ML_p")),
                                            tabPanel("ML_aic", plotOutput("plot_vanSmeden_ML_aic")),
                                            tabPanel("Firth_p", plotOutput("plot_vanSmeden_Firth_p")),
                                            tabPanel("Firth_aic", plotOutput("plot_vanSmeden_Firth_aic"))
                                ),
                                div(style="color:blue", strong("Note:"), " Please keep in mind that these three metrics are not comparable!
                        Computation of sample size is totally seperated.")
                              ),
                    wellPanel(
                      p("This report is based on your input of your data set"),
                      radioButtons('format', 'Report', c('PDF', 'HTML', 'Word'),
                                   inline = TRUE),
                      downloadButton('Report_Smeden') 
                    )
                        )
                       )
                     )
                   )
                 ),
      ##########################################
      ######Guidence ##########################
      ##########################################
      tabPanel("Guidance ",
               sidebarLayout(
                 sidebarPanel(
                   wellPanel(h3("Data characteristics"),
                     sliderInput(inputId="p_compare", "Number of Predictors",
                               min = 4, max = 12,step = 1,
                               value = 8,animate = TRUE),
                   sliderInput(inputId="Ef_compare", "Events Fraction",
                               min = 1/16, max = 1/2,step = 0.001,
                               value = 0.2,animate = TRUE)),
                   
                   wellPanel(h3("Specification of desired values for Riley method"),
                   sliderInput("S_VH_compare", "Desired Shrinkage Factor:",
                               min = 0.9, max = 1,
                               value = 0.9),
                   sliderInput("d_compare", "Absolute Difference",
                               min = 0.01, max = 0.1,
                               value = 0.05,animate = TRUE),
                   sliderInput("error_mar_compare", "Margin of Error",
                               min = 0.05, max = 0.1,
                               value = 0.05,animate = TRUE)),
                   wellPanel(h3("Model performance measure by Riley"),
                   sliderInput("csR2_compare", "Anticipated adjusted Cox-Snell R-squared",
                               min = 0.1, max = 0.3,
                               value = 0.15,animate = TRUE),
                   div(strong("Note:"), " Anticipated adjusted Cox-Snell R-squared would decide the prediction error metrics values with consistent sample size", style = "color:blue")
                   )
                 ),
                 mainPanel(
                   fluidRow(
                     column(10,
                            wellPanel(h3("Consistent model performance measures by van Smeden"),
                              p("Prediction errors wtih respect to consistent sample size based on Cox-Snell R-squared value"),
                              tableOutput(outputId="metrics_value"),
                              plotOutput(outputId="plot_Compare"),
                              div(strong("Note:"), " The vertical dashed line is the input value of Cox-Snell R-squared and corresponds to values of 
                                  predictions error metrics which could give us the same result of calculated sample size ",style = "color:blue"))
                     )
                   )
                 )
               )
      )
    ),

h4("Developed by Xijin Chen"),
h3(textOutput("currentTime"))
)

shinyApp(ui = ui, server = server)
