library("shiny")
library("gridExtra")

ou.simulate = function(alpha = 0, mu0 = 0, mu = 0, sigma, n_steps = 100) {

  time_series = rep(mu0, n_steps)
  sd          = sqrt(sigma)

  for(i in 2:n_steps) {
    time_series[i] = { alpha * {mu - time_series[i-1]} } + rnorm(1, time_series[ i-1 ], sd)
  }

  return(time_series)
}


ui = fluidPage(

  titlePanel("OU process"),

  # UI Controls
  sidebarLayout(

    sidebarPanel(
      sliderInput(inputId = "alpha_m",
                  label   = "alpha -- 'selection'",
                  min     = 0,
                  max     = 0.07,
                  value   = 0,
                  step    = 0.0005,
                  ticks   = FALSE),

      sliderInput(inputId = "sigma_m",
                  label   = "sigma -- rate of evolution",
                  min     = 0,
                  max     = 1,
                  value   = 0.001,
                  step    = 0.001,
                  ticks   = FALSE),

      sliderInput(inputId = "mu_m",
                  label   = "Optimum trait",
                  min     = -40,
                  max     = 40,
                  value   = 0,
                  step    = 0.5,
                  ticks   = FALSE),

      sliderInput(inputId = "mu0_m",
                  label   = "Ancestral trait",
                  min     = -40,
                  max     = 40,
                  value   = 0,
                  step    = 0.5,
                  ticks   = FALSE),

      sliderInput(inputId = 'replicates_m',
                  label   = 'Number of Replicates',
                  value   = 3,
                  min     = 1,
                  max     = 100),

      sliderInput(inputId = 'nsteps_m',
                  label   = 'Time steps',
                  value   = 50,
                  min     = 10,
                  max     = 200)
    ),

    # Render
    mainPanel(
      plotOutput("mainplot")
    ),
    fluid = TRUE
  )
)


server = function(input, output){

  output$mainplot = renderPlot(
    {
      y = matrix(nrow = input$nsteps_m, ncol = input$replicates_m)

      for(i in 1:ncol(y) ){
        y[,i] = ou.simulate(alpha = input$alpha_m, mu0 = input$mu0_m, mu = input$mu_m, sigma = input$sigma_m, n_steps =  input$nsteps_m)
      }

      plot(x = nrow(y), y = input$mu_m, type = "n", pch = 16, col = "blue", cex = 1,
           xlim = c(1, nrow(y)), ylim = range( c(y, input$mu_m ) ), xlab = "Time", ylab = "trait")

      points(x = nrow(y), y = input$mu_m, col = "red", las = 2, cex = 3, pch = 16)


      for(i in 1:ncol(y) ){
        lines(y[, i], col = "grey40", lwd = 1.5)
      }

    }, width = 600 ,height = 600
  )

}


shinyApp(ui = ui, server = server)

