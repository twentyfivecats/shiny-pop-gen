if (! "shiny" %in% rownames(installed.packages())) {
	install.packages("shiny")
}

library(shiny)

ui <- fluidPage(

  # App title ----
  titlePanel("Genetic Differentiation"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the population size ----
      sliderInput(inputId = "popSize",
                  label = "Number of individuals in each pair of populations:",
                  min = 1,
                  max = 500,
                  value = 100,
                  step =10), 
                  
      # Input: Slider for the initial allele frequency ----
      sliderInput(inputId = "p",
                  label = "Initial allele frequency:",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.1), 
 
       # Input: Slider for the number of generations to simulate ----
      sliderInput(inputId = "gens",
                  label = "Number of generations:",
                  min = 1,
                  max = 500,
                  value = 1),
      sliderInput(inputId = "reps",
                  label = "Number of loci to simulate:",
                  min = 2,
                  max = 100,
                  value = 5)
    ),

    # Main panel for displaying outputs ----
    mainPanel(
		#fluidRow(column(8, "Allele frequencies through time",plotOutput(outputId = "driftPlot") ), column(4,  plotOutput(outputId = "final1"))			 ),
     
      #fluidRow (column(8, "Genetic differentiation through time", plotOutput(outputId = "differentiation")), column(4,plotOutput(outputId = "final2")))
      
      fluidRow("Allele frequencies through time",plotOutput(outputId = "driftPlot")),
      fluidRow("Genetic differentiation through time", plotOutput(outputId = "differentiation"))
    )
  )
)



server <- function(input, output) {

 drift = function (i,p,popSize,gens)
	{
	ps = c(p)
	for (i in 1:gens)
		{
		ps=c(ps, rbinom(n=1, size=(2*popSize), prob = ps[i])/(2*popSize)) 
		}
	return(ps)	
	}


drift_outcomes<- function(reps, p, popSize, gens)
{
	outcomes  =matrix(data=NA, nrow = 2*(input$reps), ncol = input$gens+1)
	outcomes = apply(outcomes,1,drift, p, popSize, gens)
	return(outcomes)	
}


outcomes <- reactive(drift_outcomes(input$reps, input$p, input$popSize, input$gens))

output$driftPlot <- renderPlot ({

	plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "Generations", ylab = "Allele Frequency")
	abline(h=input$p, lty=3)
	subpop1 =outcomes()[,1:input$reps]
	subpop2 = outcomes()[,(1+input$reps):(2*input$reps)]
	apply(subpop1, 2,points, x=0:input$gens, pch = 19, col = "powderblue")
	apply(subpop1, 2,lines, x=0:input$gens, col = "powderblue")
	apply(subpop2, 2,points, x=0:input$gens, pch = 19, col = "tomato3")
	apply(subpop2, 2,lines, x=0:input$gens, col = "tomato3")


	})


output$final1 <-
	renderPlot({

		subpop1 =outcomes()[,1:input$reps]

		par(lwd = 3)
		thisCol=  "powderblue"#colors()[250]  
		bs = ifelse(length(unique(outcomes()[input$gens+1,]))<5, "custom", "Sturges") 
		if (bs == "custom") {bs = seq(0,1,0.1)}
		hist(subpop1[input$gens+1,], xlim = c(0,1),  border = "gray", col = thisCol, main ="", ylab="Number of populations with this outcome", xlab="Allele Frequency", breaks = bs)
		mline = paste("Allele Frequency for first subpopulation in a pair at generation", input$gens)
		mtext(side = 3, text = mline, line = 3, cex = 1.15, adj = -0.1)
		})
		
output$final2 <-
	renderPlot({
		subpop2 = outcomes()[,(1+input$reps):(2*input$reps)]
		par(lwd = 3)
		thisCol=  "tomato3"#colors()[250]  
		bs = ifelse(length(unique(outcomes()[input$gens+1,]))<5, "custom", "Sturges") 
		if (bs == "custom") {bs = seq(0,1,0.1)}
		hist(subpop2[input$gens+1,], xlim = c(0,1),  border = "gray", col = thisCol, main ="", ylab="Number of populations with this outcome", xlab="Allele Frequency", breaks = bs)
		mline = paste("Allele Frequency for second subpopulation in a pair  at generation", input$gens)
		mtext(side = 3, text = mline, line = 3, cex = 1.15, adj = -0.1)
		})
		
		
output$differentiation <-
	renderPlot({
		het = function(af){return(2.0*af*(1-af))}
		fst =function(Hs, Ht){return(ifelse(Ht==0, 0, (Ht-Hs)/Ht))}		
		subpop1 =outcomes()[,1:input$reps]
		subpop2 = outcomes()[,(1+input$reps):(2*input$reps)]
		afT = (subpop1+ subpop2)/2
		Ht = het(afT)
		Hs1 = het(subpop1)
		Hs2 = het(subpop2)
		Hs = (Hs1 + Hs2)/2	
		obsFst = fst(Hs,Ht)			
		if (!is.null(dim(obsFst))) {
	 		meanD = apply(obsFst,1,mean)
			sD = apply(obsFst,1,sd)
		} else {
			meanD = obsFst
			sD = 0
			
		}
 		
		mline = "Mean genetic differentiation across populations"
		plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "", ylab = "")
		abline(h=2*(1-input$p)*input$p, lty=3)
		lines(x=0:input$gens, meanD+sD, lwd= 2, lty=2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD-sD, lwd= 2, lty=2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD, lwd=5, col ="steelblue")
		mtext(side=1, line = 3, adj=0, "Mean Fst (dark line) ± 1 standard deviation (light lines)")
		mtext(side = 1, text = "Generations", line =2, cex =1.2)
		mtext(side = 2, text = "Fst", line =2, cex =1.2)
		})

}


shinyApp(ui = ui, server = server)

