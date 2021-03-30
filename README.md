# shiny-pop-gen

This is a suite of little shiny apps and worksheets I wrote to teach population genetics to undergraduates at the University of Liverpool. The code is no doubt rough, and mostly uncommented. Suggestions welcome; corrections rec'd with gratitude.

The easiest way to run them is to copy and paste the whole app into R.  Slightly more elegantly, load the app as source, and then run the last line in the console "shinyApp(ui = ui, server = server)".  In either case, the message "Listening on http://(some numbers)" should appear, and the app should open a page in your default web browser, and show you some sliders and plots. To quit the app, press 'escape' in the R console.

You can find the apps in the apps folder, and the worksheets in the worksheets folder.  I did not make an answer key for the worksheets.

Selection.R --illustrates deterministic directional selection models

GeneticDriftDiv.R-- illustrates genetic drift in finite populations, showing allele frequency changes and declines in genetic diversity

Selection_and_Drift.R-- illustrates directional selection in finite populations

Balancing_Selection_and_Drift.R-- illustrates models of under- and over-dominant selection in finite populations

Selection_all.R-- a catch-all app; can be used to show any of the above models

BreedersEq.R-- the odd one out; illustrates the Breeders' Equation in quantitative genetics

Worksheets of the same name give instructions and a small number of questions to work through.





