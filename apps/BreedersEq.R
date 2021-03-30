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

      sliderInput(inputId = "N",
                  label = "How many individuals in the breeding population?",
                  min = 10,
                  max = 1000,
                  value = 100,
                  step =0.1), 

      sliderInput(inputId = "h2",
                  label = "Heritability, h2",
                  min = 0,
                  max = 1,
                  value = 0.5,
                  step =0.1), 
                  
       # Input: Slider for the number of generations to simulate ----
      sliderInput(inputId = "thresh",
                  label = "Selection threshold:",
                  min = -3,
                  max = 3,
                  value = 0.5, 
                  step = 0.1),
     radioButtons(inputId="direction", label="Breed from parents with trait values above or below the threshhold value:", choices=c('above', 'below'))
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      plotOutput(outputId = "parents"),
      plotOutput(outputId = "selectedParents"),
      plotOutput(outputId = "offspring"),

      # Output: Histogram ----
		
    )
  )
)



server <- function(input, output) {


parents <- reactive(rnorm(n=input$N))
	
	
selected <-reactive(	{
	if (input$direction == "above")
		{
		parents()[which(parents()>input$thresh)]			
		} else if (input$direction == "below")
		{
		parents()[which(parents()<input$thresh)]			
		}
	})
kids <-reactive({
	rnorm(n=input$N,mean=mean(selected())*input$h2 )
	})

	
selDiff <- reactive({mean(selected()) - mean(parents())})
selRes <- reactive({mean(kids()) - mean(parents())})
	
output$parents <- renderPlot ({
	x=hist(parents(), col="cadetblue3", border="gray", main = "",xlab = "", ylab = "Frequency", breaks =20, cex.axis=1.5, cex.lab=1.5, xlim =c (-3, 3))
	mtext(line=2, "Distribution of trait values for 100 parents", cex = 1.8, font= 2,col ="darkgray", adj=-.1)
	abline(v=mean(parents()), lwd=5, lty=5, col="dodgerblue4")
	mtext(side = 3, adj=0, paste("Trait Value for all parents: mean = ", format(mean(parents()),digits=2)), line = 0.5, cex =1.3, font=2,  col ="darkgray")
	if (input$direction == "above") {	rect(xleft = input$thresh, xright = max(parents()) , ybottom = 0, ytop = max(x$counts), lty = 1, lwd = 2, border = "darkolivegreen4" ) } else {
		rect(xright = input$thresh, xleft = min(parents()) , ybottom = 0, ytop = max(x$counts), lty = 2, lwd = 2, border = "darkolivegreen4") }
	})


output$selectedParents <- renderPlot ({
	hist(selected(), col="cadetblue3", border="gray", xlab = "", ylab = "Frequency", main = "", xlim=c(-3,3), cex.axis=1.5, cex.lab=1.5)
	if (length(selected())<=1) 
		{mtext("Selection is too strong. There are no parents to breed from; try changing the selection threshhold.", cex = 1.5)} else
		{
		mm=paste("Distribution of trait values from the ", format(length(selected())/length(parents())*100, digits=2), " % of parents selected for breeding", sep="")
		mtext(line=2, mm, font= 2,cex=1.8, col ="darkgray", adj=-.1)
		mtext(side = 3, adj=0, paste("Trait Value for all parents: mean = ", format(mean(selected()),digits=2)), line = 0.5, cex =1.3, font=2,  col ="darkgray")
		abline(v= mean(selected()), lty=5, col="darkolivegreen4", lwd=5)

		}
	})

output$offspring <- renderPlot ({
		hist(kids(), col="cadetblue3", border="gray", xlab = "", ylab = "Frequency", main = "", xlim=c(-3,3), cex.axis = 1.5, cex.lab=1.5)
		mtext(line=2, "Distribution of trait values for offspring of selected parents", font= 2,cex=1.8, col ="darkgray", adj=-.1)
		mtext(side = 3, adj=0, paste("Trait Value for offspring: mean = ", format(mean(kids()),digits=2)), line = 0.5, cex =1.3, font=2,  col ="darkgray")
		abline(v=mean(kids()), lty=1, lwd=5, col = "firebrick4")
		abline(v=mean(parents()), lwd=5, lty=5, col="dodgerblue4")
		abline(v= mean(selected()), lty=5,lwd=5, col="darkolivegreen4")
		co=max(hist(kids(), plot=FALSE)$counts)
		arrows(lwd=3, y0=co/3,y1=co/3,x0=mean(parents()), x1=mean(selected()), code=3, angle = 90, col="dimgrey")
		arrows(lwd=3, y0=co/2,y1=co/2,x0=mean(parents()), x1=mean(kids()), code=3, angle = 90, col="dimgrey")
		text(paste("Expected response = h2 * S = ",  format(input$h2, digits=2), "*", format(selDiff(), digits=2)," = ", format(input$h2*selDiff(), digits=2)), pos = 4, x= min(kids()), y = co*0.8, cex = 1.2)
		text(paste("Actual response = ", format(selRes(),digits=2)), pos = 4,  x= min(kids()), y = co*0.7, cex = 1.2)
		
		text(paste("S = ", format(selDiff(), digits=2)), x = mean(parents())+(0.5* selDiff()), y= co/3, pos = 3, cex = 1.2, font=2, adj = 0.5)  
		text(paste("R = ", format(selRes(),digits=2)), x=mean(parents())+(0.5* selRes()), y=co/2 , pos = 3, cex = 1.2, font=2) 
		})

}


shinyApp(ui = ui, server = server)

