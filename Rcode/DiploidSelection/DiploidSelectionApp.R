shinyUI(fluidPage(
  titlePanel("driftR: Population Genetic Simulations in R"),
  sidebarLayout(
    sidebarPanel(
      numericInput("p","Starting allele frequency",value=0.5,min=0,max=1),
      sliderInput("Uab","Mutation Rate",value=0,min=0,max=0.1),
      sliderInput("Waa","Fitness of genotype AA",value=1,min=0,max=1),
      sliderInput("Wab","Fitness of genotype AB",value=1,min=0,max=1),
      sliderInput("Wbb","Fitness of genotype BB",value=1,min=0,max=1),
      sliderInput("m","Migration Rate",0,min=0,max=0.35),
      numericInput("n","Population Size",100,min=1,max=1e5),
      numericInput("nPop","Number of Populations",2,min=1,max=100),
      numericInput("gen","Number of Generations",100,min=1,max=5000),
      checkboxInput("infinitePop","Infinite Population (no drift)",value = F),
      checkboxGroupInput(inputId="plotStats",label="plot:",choices=c("p","He","Hs","Ht","Fst"),inline=T,selected="p"),
      checkboxInput("legend","Legend",value = F),
      actionButton("go","go",width="100%"),
      div(helpText("driftR simulates allele and genotype frequencies for a single biallelic locus in biological populations. 
                   Core functions adapted from the Java program popG (http://evolution.gs.washington.edu/popgen/popg.html). 
                   Full code available on github: https://github.com/cjbattey/driftR"),style="font-size:75%")
      ),
    
    mainPanel(
      plotOutput("plot"),
      textOutput("nLost"),
      tableOutput("endStateTable"),
      div(style="border-top:1px solid black;"),
      helpText("Click the button below to run 100, 100-generation simulations of 2 populations using the current 
               parameters."),
      actionButton("runSim","Run Replicate Simulations"),
      tableOutput("meanTable"),
      tableOutput("varTable"),
      div(tableOutput("sumTable"), style = "font-size: 75%; width: 75%;")
      )
      )
  
  ))

shinyServer(function(input,output,session){
  
  source("./dev.R")
  
  sim.data <- eventReactive(input$go,ignoreNULL = F,{
    validate(
      need(input$gen<5001,"Please select < 5000 generations."),
      need(input$nPop<101,"Please select < 100 populations"),
      need(input$n<1000001,"Please select n < 1,000,000"),
      need(input$plotStats!="","Select a variable to plot.")
    )
    runPopSim(gen=input$gen,p=input$p,Waa=input$Waa,Wab=input$Wab,Wbb=input$Wbb,n=input$n,
              nPop=input$nPop,m=input$m,stats=input$plotStats,infinitePop=input$infinitePop,Uab=input$Uab,Uba=input$Uab)
  })
  
  plot.data <- eventReactive(sim.data(),{
    meltPlotData(allele.freq.df = sim.data(),gen=input$gen,nPop=input$nPop,stats=input$plotStats)
  })
  
  output$plot <- renderPlot({
    plotSingleRun(plot.data(),nPop=input$nPop,gen=input$gen,legend=input$legend)
  })
  
  nLost.text <- eventReactive(sim.data(),{
    p <- sim.data()[input$gen+1,1:input$nPop]
    nFixed <- length(p[p==1])
    nLost <- length(p[p==0])
    paste0("Fixed: ",nFixed,".  Lost: ",nLost,".")
  })
  
  output$nLost <- renderText({
    nLost.text()
  })
  
  output$endStateTable <- renderTable({
    endState <- sim.data()[input$gen,c("Fis","Hs","Ht","Fst")]
  })
  
  sumTable <- eventReactive(input$runSim,{
    validate(
      need(input$n<=100000,"Please select n <= 100,000")
    )
    sumTable <- data.frame(matrix(ncol=14))
    withProgress(message="simulating populations...",value=0,{
      for(i in 1:100){
        df <- runPopSim2(gen=100,p=input$p,Waa=input$Waa,Wab=input$Wab,Wbb=input$Wbb,n=input$n,
                         nPop=2,m=input$m,infinitePop=input$infinitePop,Uab=input$Uab,Uba=input$Uab)
        names(sumTable) <- names(df)
        sumTable[i,] <- df[nrow(df),]
        incProgress(1/100)
      }
    })
    sumTable
  })
  
  meanTable <- reactive({
    tbl <- colMeans(sumTable(),na.rm=T)
    tbl <- tbl[c("Fis","Hs","Ht","Fst")]
    t(tbl)
  })
  
  varTable <- reactive({
    tbl <- apply(sumTable(),2,function(e) var(e,na.rm=T))
    tbl <- tbl[c("Fis","Hs","Ht","Fst")]
    t(tbl)
  })
  
  
  
  output$meanTable <- renderTable(meanTable(),colnames = T,digits=4,caption = "Mean state at final generation:",
                                  caption.placement = getOption("xtable.caption.placement", "top"),
                                  caption.width = getOption("xtable.caption.width", NULL))
  
  output$varTable <- renderTable(varTable(),colnames = T,digits=4,caption = "Variance:",
                                 caption.placement = getOption("xtable.caption.placement", "top"),
                                 caption.width = getOption("xtable.caption.width", NULL))
  
  output$sumTable <- renderTable(sumTable(),caption = "Final Generation States:",
                                 caption.placement = getOption("xtable.caption.placement", "top"),
                                 caption.width = getOption("xtable.caption.width", NULL))
})

#development script for popgen simulations (+/- R translation of PopG, adding more visualization and summary stats)
library(plyr);library(reshape);library(ggplot2);library(magrittr)

#binomial draw for new genotype frequencies
binomialDraw <- function(n,p){ 
  draws <- runif(n)
  bnl <- length(draws[draws<p])
  return(bnl)
}

#main simulation function
runPopSim <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst"),Uab=0,Uba=0,infinitePop=F){
  withProgress(message="simulating populations...",value=0,{
    allele.freq <- data.frame(matrix(ncol=3*nPop))
    allele.freq[1,(1:nPop)] <- rep(p,nPop) #starting allele freqs
    for(i in 1:gen){ 
      mean.p <- as.numeric(rowMeans(allele.freq[i,(1:nPop)]))
      for(j in 1:nPop){
        p <- allele.freq[i,j]
        p <- (1 - Uab) * p + Uba * (1 - p) #mutation
        p <- p*(1-m)+m*mean.p # migration
        q <- 1-p
        if(p>0 && p<1){ #if alleles are not fixed
          w <- p*p*Waa+2*p*q*Wab+q*q*Wbb #population average fitness
          freq.aa <- (p*p*Waa)/w #post-selection genotype frequencies (weighted by relative fitness)
          freq.ab <- (2*p*q*Wab)/w
          if(infinitePop==F){ 
            Naa <- binomialDraw(n,freq.aa) #binomial draw for new genotype counts (drift)
            if(freq.aa<1){ 
              Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
            }
            else {
              Nab <- 0
            }
            p <- ((2*Naa)+Nab)/(2*n)
            q <- 1-p
            allele.freq[(i+1),j] <- p #new p after drift in columns 1:nPop
            allele.freq[(i+1),(j+nPop)] <- Nab/n #Ho in columns (nPop+1):(nPop*2)
            allele.freq[(i+1),(j+2*nPop)] <- 2*p*q #He in columns (nPop*2+1):nPop*3
          } 
          else { #no drift (infinite population) conditions
            p <- freq.aa+(freq.ab/2)
            q <- 1-p
            allele.freq[(i+1),j] <- p
            allele.freq[(i+1),(j+nPop)] <- freq.ab 
            allele.freq[(i+1),(j+2*nPop)] <- 2*p*q
          }
        } else { #if alleles are fixed
          if(p<=0){
            p <- 0
            
          } else {
            p <- 1
            
          }
          allele.freq[(i+1),j] <- p
          allele.freq[(i+1),(j+nPop)] <- 0
          allele.freq[(i+1),(j+2*nPop)] <- 0
        }
      } #end populations loop
      incProgress(1/gen)
    } #end generations loop
    #summary stats
    names <- c()
    for(i in 1:nPop){names[i]<-paste0("p",i)}
    for(i in (nPop+1):(2*nPop)){names[i]<-paste0("Ho",i-nPop)}
    for(i in (nPop+nPop+1):(3*nPop)){names[i]<-paste0("He",i-2*nPop)}
    colnames(allele.freq) <- names
    allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
    allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
    allele.freq$Fis <- abs(1-(allele.freq$meanHo/allele.freq$meanHe))
    allele.freq$mean.p <- rowMeans(allele.freq[1:nPop]) 
    allele.freq$Hs <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
    allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
    allele.freq$Fst <- (allele.freq$Ht-allele.freq$Hs)/allele.freq$Ht
    allele.freq$Fst[allele.freq$Fst<0] <- 0
    allele.freq$gen <- 0:gen
    return(allele.freq)
  })
}

#format for plotting
meltPlotData <- function(allele.freq.df=allele.freq.df,gen=100,nPop=2,stats=c("p","Fst")){
  df <- melt(allele.freq.df,id.vars = "gen")
  df$dataType <- c(rep("p",(nPop*(gen+1))),rep("Ho",nPop*(gen+1)),rep("He",nPop*(gen+1)),rep("meanHo",(gen+1)),rep("meanHe",(gen+1)),
                   rep("Fis",(gen+1)),rep("mean.p",(gen+1)),rep("Hs",(gen+1)),rep("Ht",(gen+1)),rep("Fst",(gen+1)))
  df <- subset(df,dataType %in% stats)
  return(df)
}

#plotting function
plotSingleRun <- function(df,nPop,gen,legend){
  if(legend==T){
    print(ggplot(df,aes(x=gen,y=value,col=variable))+ylim(0,1)+facet_wrap(~dataType)+geom_line()+xlab("Generations")+ylab(""))
  } else {
    print(ggplot(df,aes(x=gen,y=value,col=variable))+theme(legend.position="none")+ylim(0,1)+facet_wrap(~dataType)+geom_line()+xlab("Generations")+ylab(""))
  }
  
}

#simulation function w/o progress bar for debug+replicate runs
runPopSim2 <- function(gen=100,p=0.5,Waa=1,Wab=1,Wbb=1,n=100,nPop=2,m=0,stats=c("p","Fst"),Uab=0,Uba=0,infinitePop=F){
  allele.freq <- data.frame(matrix(ncol=3*nPop))
  allele.freq[1,(1:nPop)] <- rep(p,nPop) #starting allele freqs
  for(i in 1:gen){ 
    mean.p <- as.numeric(rowMeans(allele.freq[i,(1:nPop)]))
    for(j in 1:nPop){
      p <- allele.freq[i,j]
      p <- (1 - Uab) * p + Uba * (1 - p) #mutation
      p <- p*(1-m)+m*mean.p # migration
      q <- 1-p
      if(p>0 && p<1){ #if alleles are not fixed
        w <- p*p*Waa+2*p*q*Wab+q*q*Wbb #population average fitness
        freq.aa <- (p*p*Waa)/w #post-selection genotype frequencies (weighted by relative fitness)
        freq.ab <- (2*p*q*Wab)/w
        if(infinitePop==F){ 
          Naa <- binomialDraw(n,freq.aa) #binomial draw for new genotype counts (drift)
          if(freq.aa<1){ 
            Nab <- binomialDraw((n-Naa),(freq.ab/(1-freq.aa)))
          }
          else {
            Nab <- 0
          }
          p <- ((2*Naa)+Nab)/(2*n)
          q <- 1-p
          allele.freq[(i+1),j] <- p #new p after drift in columns 1:nPop
          allele.freq[(i+1),(j+nPop)] <- Nab/n #Ho in columns (nPop+1):(nPop*2)
          allele.freq[(i+1),(j+2*nPop)] <- 2*p*q #He in columns (nPop*2+1):nPop*3
        } 
        else { #no drift (infinite population) conditions
          p <- freq.aa+(freq.ab/2)
          q <- 1-p
          allele.freq[(i+1),j] <- p
          allele.freq[(i+1),(j+nPop)] <- freq.ab 
          allele.freq[(i+1),(j+2*nPop)] <- 2*p*q
        }
      } else { #if alleles are fixed
        if(p<=0){
          p <- 0
          
        } else {
          p <- 1
          
        }
        allele.freq[(i+1),j] <- p
        allele.freq[(i+1),(j+nPop)] <- 0
        allele.freq[(i+1),(j+2*nPop)] <- 0
      }
    } #end populations loop
  } #end generations loop
  #summary stats
  names <- c()
  for(i in 1:nPop){names[i]<-paste0("p",i)}
  for(i in (nPop+1):(2*nPop)){names[i]<-paste0("Ho",i-nPop)}
  for(i in (nPop+nPop+1):(3*nPop)){names[i]<-paste0("He",i-2*nPop)}
  colnames(allele.freq) <- names
  allele.freq$meanHo <- rowMeans(allele.freq[(nPop+1):(nPop*2)])
  allele.freq$meanHe <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Fis <- 1-(allele.freq$meanHo/allele.freq$meanHe)
  allele.freq$mean.p <- rowMeans(allele.freq[1:nPop]) 
  allele.freq$Hs <- rowMeans(allele.freq[(nPop*2+1):(nPop*3)])
  allele.freq$Ht <- 2*allele.freq$mean.p*(1 - allele.freq$mean.p)
  allele.freq$Fst <- (allele.freq$Ht-allele.freq$Hs)/allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst<0] <- 0
  allele.freq$gen <- 0:gen
  return(allele.freq)
}

# Run the application 
shinyApp(ui = ui, server = server)

