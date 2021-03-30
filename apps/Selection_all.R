if (! "shiny" %in% rownames(installed.packages())) {
	install.packages("shiny")
}

library(shiny)

ui <- fluidPage(

  # App title ----
  titlePanel("General selection model with or without genetic drift"),

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
       sliderInput(inputId = "gens",
                  label = "Number of generations:",
                  min = 1,
                  max = 500,
                  value = 5),

     radioButtons(inputId="modeldrift", label="Should the model be deterministic (with no genetic drift) or stochastic (with genetic drift)?", choices=c('deterministic', 'stochastic')), 
	 conditionalPanel("input.modeldrift == 'stochastic'", 
      	sliderInput(inputId = "popSize",
                  label = "Number of individuals:",
                  min = 1,
                  max = 500,
                  value = 1,
                  step =10),
                  
      sliderInput(inputId = "reps",
                  label = "Number of populations to simulate:",
                  min = 1,
                  max = 100,
                  value = 1)
                   
			),

	    sliderInput(inputId = "wAA",
                  label = "Fitness of AA homozygote, wAA:",
                  min = 0,
                  max = 2,
                  value = 1,
                  step = 0.01), 
	    sliderInput(inputId = "wAa",
                  label = "Fitness of Aa heterozygote, wAa:",
                  min = 0,
                  max = 2,
                  value = 1,
                  step = 0.01), 
	 sliderInput(inputId = "waa",
                  label = "Fitness of aa homozygote, waa:",
                  min = 0,
                  max = 2,
                  value = 1,
                  step = 0.01),

      plotOutput(outputId = "fitnessGraph")     	
     ),
	
    # Main panel for displaying outputs ----
    mainPanel(
    
   		fluidRow(column(8, "",plotOutput(outputId = "frequencyPlot")),
		 column(4,  plotOutput(outputId = "final"))),

		#tableOutput("selMat")
      fluidRow (column(8,"", plotOutput(outputId = "diversity")))

		
    )
  )
  )
  
  
  

server <- function(input, output) {

eqP <- reactive(
	{
		s = abs(wAA() -1)
		t = abs(waa() -1)
		return(t/(s+t))		
	}
)

wAA <-reactive({
	ifelse(input$model == "under", input$wAA_under, input$wAA_over)
})

waa <-reactive({
	ifelse(input$model == "under", input$waa_under, input$waa_over)
})


newP<-function(p)
{
	num <- (p**2 * input$wAA) + p*(1-p)* input$wAa 
	denom = (p**2 * input$wAA) + 2*p*(1-p)*input$wAa + ((1-p)**2)*input$waa
	return(num/denom)
}

selection_drift = function (i,p,h, s,popSize,gens)
	{
	ps = c(p)
	for (i in 1:gens)
		{
		effectiveP = newP(ps[i])
		
		if (input$modeldrift == "stochastic") 
			{	
			ps=c(ps,rbinom(n=1, size=(2*popSize), prob = effectiveP)/(2*popSize))
				} else {
			nextP = ps=c(ps,effectiveP)
				}
		}
	return(ps)	
	}



outcomes<- reactive(
{
	outcomes  =matrix(data=NA, ncol = input$reps, nrow = input$gens+1)
	outcomes = apply(outcomes,2, selection_drift, input$p, input$h,input$s, input$popSize, input$gens)
	return(outcomes)	
})



output$spacer <- renderText({rep('-',40)})


output$frequencyPlot <- renderPlot ({
	plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n",  cex.lab = 1.5, cex.axis = 1.2, xlab = "Generations", ylab = "Allele Frequency")
	abline(h=eqP(), lty=3)
	mtext (text="Trajectory of selected allele over time", font = 2, cex = 1.5,side=3, adj=0.1, line=2, col="darkslategray")
	sdO = outcomes()
	apply(sdO,2,points, x=0:input$gens, pch = 19, col = "cadetblue4")
	apply(sdO,2,lines, x=0:input$gens, col = "deepskyblue4")
	})

output$final <-
	renderPlot({
		par(lwd = 3) 
		thisCol= "cadetblue4"  
		bs = ifelse(length(unique(outcomes()[input$gens+1,]))<5, "custom", "Sturges") 
		if (bs == "custom") {bs = seq(0,1,0.1)}
		hist(outcomes()[input$gens+1,], xlim = c(0,1), main ="", border = "cadetblue4", cex.axis = 1.2, cex.lab = 1.5, col = thisCol, ylab="Number of populations with this outcome", xlab="Allele Frequency", breaks = bs) 
		mtext(col="darkslategray", cex = 1.5,side=3, adj=0.1, font = 2,line=-1,text=paste("Allele frequency at generation", input$gens, "\n across ", input$reps, " populations"))
		})

output$fitnessGraph <-renderPlot({
	barplot(names = c("wAA", "wAa", "waa"), height=c(input$wAA, input$wAa, input$waa), col = "cadetblue3", border ="gray", ylim=c(0,2), cex.axis = 1.2, cex.lab = 1.6, font.lab =2)
	abline(h=1, lty=2, lwd =2)
})

output$diversity <-
	renderPlot({
		het = function(af){return(2.0*af*(1-af))}
		hetO = het(outcomes())
		meanD = apply(hetO,1,mean)
		sD = apply(hetO,1,sd)

		plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "", ylab = "", cex.axis = 1.2)
		abline(h=2*(1-input$p)*input$p, lty=3)
		lines(x=0:input$gens, meanD+sD, lwd= 2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD-sD, lwd= 2, col = adjustcolor("steelblue", 0.5))
		lines(x=0:input$gens, meanD, lwd=5, col ="steelblue")
		mtext(col="darkslategray", cex = 1.5,side=3, adj=0.1, font = 2,line=2,text= "Mean genetic diversity across populations")
		mtext(side=1, line = 3, adj=0, "Mean diversity (dark line) ± 1 standard deviation (light lines)", cex=1.2)
		mtext(side = 1, text = "Generations", line =2, cex =1.5)
		mtext(side = 2, text = "Genetic Diversity (2pq)", line =3, cex =1.5)
		})
	
}


shinyApp(ui = ui, server = server)

