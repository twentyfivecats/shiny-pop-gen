if (! "shiny" %in% rownames(installed.packages())) {
	install.packages("shiny")
}

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Genetic Drift"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Slider for the population size ----
      sliderInput(inputId = "popSize",
                  label = "Number of individuals:",
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
                  label = "Number of populations to simulate:",
                  min = 1,
                  max = 100,
                  value = 1)
    ),

    # Main panel for displaying outputs ----
    mainPanel(
		fluidRow(column(8, "Allele frequencies through time",plotOutput(outputId = "driftPlot") ),
			 column(4,  plotOutput(outputId = "final"))
			 ),
     
      fluidRow (column(8, "Genetic diversity through time", plotOutput(outputId = "diversity")), column(4,""))
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
	outcomes  =matrix(data=NA, nrow = input$reps, ncol = input$gens+1)
	outcomes = apply(outcomes,1,drift, p, popSize, gens)
	return(outcomes)	
}


outcomes <- reactive(drift_outcomes(input$reps, input$p, input$popSize, input$gens))

output$driftPlot <- renderPlot ({

	plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "Generations", ylab = "Allele Frequency")
	abline(h=input$p, lty=3)
	apply(outcomes(),2,points, x=0:input$gens, pch = 19, col = "powderblue")
	apply(outcomes(),2,lines, x=0:input$gens, col = "powderblue")
	})


output$final <-
	renderPlot({
		par(lwd = 3)
		thisCol= colors()[250]  
		bs = ifelse(length(unique(outcomes()[input$gens+1,]))<5, "custom", "Sturges") 
		if (bs == "custom") {bs = seq(0,1,0.1)}
		hist(outcomes()[input$gens+1,], xlim = c(0,1),  border = "gray", col = thisCol, main ="", ylab="Number of populations with this outcome", xlab="Allele Frequency", breaks = bs)
		mline = paste("Allele Frequency at generation", input$gens)
		mtext(side = 3, text = mline, line = 3, cex = 1.15, adj = -0.1)
		})
		
		
output$diversity <-
	renderPlot({
		het = function(af){return(2.0*af*(1-af))}
		hetO = het(outcomes())
		meanD = apply(hetO,1,mean)
		sD = apply(hetO,1,sd)

		mline = "Mean genetic diversity across populations"
		plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "", ylab = "")
		abline(h=2*(1-input$p)*input$p, lty=3)
		lines(x=0:input$gens, meanD+sD, lwd= 2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD-sD, lwd= 2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD, lwd=5, col ="steelblue")
		mtext(side=1, line = 3, adj=0, "Mean diversity (dark line) ± 1 standard deviation (light lines)")
		mtext(side = 1, text = "Generations", line =2, cex =1.2)
		mtext(side = 2, text = "Genetic Diversity (heterozygosity, 2pq)", line =2, cex =1.2)
		})

}


shinyApp(ui = ui, server = server)

