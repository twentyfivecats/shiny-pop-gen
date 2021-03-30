if (! "shiny" %in% rownames(installed.packages())) {
	install.packages("shiny")
}

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Directional Selection and Genetic Drift"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
                
      # Input: Slider for the initial allele frequency ----
      sliderInput(inputId = "p",
                  label = "Initial frequency of the A allele:",
                  min = 0,
                  max = 1,
                  value = 0.1,
                  step = 0.01), 
 
      sliderInput(inputId = "s",
                  label = "Selection coefficient acting on A:",
                  min = -0.5,
                  max = 0.5,
                  value = 0.1,
                  step = 0.01), 
      sliderInput(inputId = "h",
                  label = "Dominance:",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.1), 
       # Input: Slider for the number of generations to simulate ----
      sliderInput(inputId = "popSize",
                  label = "Number of individuals:",
                  min = 1,
                  max = 500,
                  value = 1,
                  step =10), 
      sliderInput(inputId = "gens",
                  label = "Number of generations:",
                  min = 1,
                  max = 500,
                  value = 5),
                  
      sliderInput(inputId = "reps",
                  label = "Number of populations to simulate:",
                  min = 1,
                  max = 100,
                  value = 10),
   
  		fluidRow("Fitness of Genotypes"), tableOutput(outputId = "fitnessScheme")
  		),

    

    # Main panel for displaying outputs ----
    mainPanel(

     	fluidRow("", plotOutput(outputId = "frequencyPlot")),
      # Output: Histogram ----
   		fluidRow("", plotOutput(outputId = "final")),
		#tableOutput("selMat")
    )
  )

)


server <- function(input, output) {



newP<-function(p, h, s)
{
	num <- (p**2 * (1+s)) + p*(1-p)*(1+h*s)
	denom = (p**2 * (1+s)) + 2*p*(1-p)*(1+h*s) +(1-p)**2
	return(num/denom)
}

selection_drift = function (i,p,h, s,popSize,gens)
	{
	ps = c(p)
	for (i in 1:gens)
		{
		effectiveP = newP(ps[i],h,s)
		nextP =rbinom(n=1, size=(2*popSize), prob = effectiveP)/(2*popSize)
		ps=c(ps,nextP) 
		}
	return(ps)	
	}



selection_drift_outcomes<- function()
{
	outcomes  =matrix(data=NA, ncol = input$reps, nrow = input$gens+1)
	outcomes = apply(outcomes,2, selection_drift, input$p, input$h,input$s, input$popSize, input$gens)
	return(outcomes)	
}

outcomes <- reactive(selection_drift_outcomes())

output$selMat = renderTable({outcomes()})
output$frequencyPlot <- renderPlot ({

	plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "Generations", ylab = "Allele Frequency")
	abline(h=input$p, lty=3)
	mtext (text="Trajectory of selected allele over time", font = 2, cex = 1.5,side=3, adj=0.1, line=2, col="darkslategray")
	sdO = selection_drift_outcomes()
	apply(sdO,2,points, x=0:input$gens, pch = 19, col = "deepskyblue4")
	apply(sdO,2,lines, x=0:input$gens, col = "deepskyblue4")
	})

output$final <-
	renderPlot({
		par(lwd = 3) 
		thisCol= "cornsilk"  
		bs = ifelse(length(unique(outcomes()[input$gens+1,]))<5, "custom", "Sturges") 
		if (bs == "custom") {bs = seq(0,1,0.1)}
		hist(outcomes()[input$gens+1,], xlim = c(0,1), main ="", border = "cornsilk3", col = thisCol, ylab="Number of populations with this outcome", xlab="Allele Frequency", breaks = bs) 
		mtext(col="darkslategray", cex = 1.5,side=3, adj=0.1, font = 2,line=2,text=paste("Allele frequency at generation", input$gens, " across ", input$reps, " populations"))
		})
	
output$fitnessScheme <-renderTable({data.frame(wAA=1+input$s, wAa = 1+input$h*input$s, waa = 1)}, align="c")

}


shinyApp(ui = ui, server = server)

