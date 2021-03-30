if (! "shiny" %in% rownames(installed.packages())) {
	install.packages("shiny")
}

library(shiny)

# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Directional Selection"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(
                
      # Input: Slider for the initial allele frequency ----
      sliderInput(inputId = "p",
                  label = "Initial frequency of the A allele:",
                  min = 0,
                  max = 1,
                  value = 0.01,
                  step = 0.01), 
 
      sliderInput(inputId = "s",
                  label = "Selection coefficient acting on A:",
                  min = -0.5,
                  max = 0.5,
                  value = 0.0,
                  step = 0.01), 
      sliderInput(inputId = "h",
                  label = "Dominance:",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step = 0.1), 
       # Input: Slider for the number of generations to simulate ----
      sliderInput(inputId = "gens",
                  label = "Number of generations:",
                  min = 1,
                  max = 500,
                  value = 50),
    

    ),

    # Main panel for displaying outputs ----
    mainPanel(

    fluidRow("", plotOutput(outputId = "frequencyPlot")),
    fluidRow(htmlOutput("fitnessSchemeLabel")),
	fluidRow (column(4,fluidRow(""), fluidRow(""), tableOutput(outputId = "fitnessScheme")),column(8, plotOutput(outputId = "fitnessGraph")))
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


newP<-function(p, h, s)
{
	num <- (p**2 * (1+s)) + p*(1-p)*(1+h*s)
	denom = (p**2 * (1+s)) + 2*p*(1-p)*(1+h*s) +(1-p)**2
	return(num/denom)
}



sel_trajectory<- function(p, h, s, gens)
{
	pnow= p
	traj = c(pnow)

	for (i in 1:gens)
		{
		pnow= newP(pnow,h, s)
		traj = c(traj,pnow)			
		}
	return(traj)	
}


trajectory <- reactive({sel_trajectory(input$p, input$h, input$s, input$gens)})

output$frequencyPlot <- renderPlot ({

	plot(0,0, xlim=c(0,input$gens), ylim = c(0,1), type = "n", xlab = "", ylab = "")
	mtext(text = "Generations", line = 2, side = 1, cex = 1.2)
	mtext(text = "Allele Frequency", side = 2, line = 3, cex = 1.2)
	abline(h=input$p, lty=3)
	abline(h=0, col="gray",lty=3)
	abline(h=1, col="gray", lty=3)
	points(x=0:input$gens,y=trajectory(), pch = 19, col = "deepskyblue4")
	lines(x=0:input$gens,y=trajectory(), col = "deepskyblue4")
	mtext(col="darkslategray", cex = 1.5,side=3, adj=0, font = 2,line=2,text="Trajectory of selected allele over time")

	})


#output$fitnessSchemeLabel <-renderText({paste("style=color:blue;font-size:18px; <b>Fitness scheme table and plot</b>")})
output$fitnessSchemeLabel <-renderText({"<b>Fitness scheme table and plot</b>"})

output$fitnessScheme <-renderTable({data.frame(row.names = c("fitness", "initial.frequency") , AA=c(1+input$s, input$p**2), Aa = c(1+input$h*input$s, 2*(input$p)*(1-input$p)), aa = c(1, (1-input$p)**2))}, rownames=TRUE)
output$fitnessGraph <-renderPlot({
	barplot(names = c("wAA", "wAa", "waa"), height=c(1+input$s,1+input$h*input$s, 1), col = "deepskyblue4", border ="gray", cex.axis = 1.2, cex.lab = 1.6, font.lab =2)
	abline(h=1, lty=2, lwd =2)
	mtext(col="darkslategray", cex = 1.5,side=3, adj=-5, font = 2,line=2,text="Fitness scheme")

})
}


shinyApp(ui = ui, server = server)



