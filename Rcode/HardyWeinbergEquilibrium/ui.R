library(shiny)

shinyUI(
  fluidPage(
    h2("Hardy-Weinberg Equilibrium Calculator",align="left"),
    sidebarPanel(
      p('Observed genotype frequencies'),
      sliderInput('gaa', 'AA',25, min = 1, max = 100, step = 1),
      sliderInput('gab', 'AB',50, min = 1, max = 100, step = 1),
      sliderInput('gbb', 'BB',25, min = 1, max = 100, step = 1),
      submitButton('calculate'),
     br(),
     p('')
      ),
  # Summary of input, results, graph and decision
    mainPanel(
      tableOutput('values'),
      verbatimTextOutput("hwp"),
      verbatimTextOutput("decision"),
      plotOutput('GenoFreqplot')

    )
  )
)