library("shiny")
library("phytools")
library("nlme")


ui = fluidPage(

  titlePanel("Phylogenetic GLS"),

  # UI Controls

  sidebarLayout(
    sidebarPanel(
      actionButton(inputId = "refresh", label = "New data", ),
      shiny::selectInput(inputId = "tree", label = NULL,
                         choices = c("Tree 1", "Tree 2" ),
                         selected = "Tree 1"), width = 2
      ),
      mainPanel(
        plotOutput("mainplot")
        ),
    fluid = TRUE
  )
)

server = function(input, output){

  output$mainplot = renderPlot(
    {

      set.seed(2)
      n     = 100
      tree1 = pbtree(n = n, scale = 1)
      tree1 = ladderize(tree1)

      #set.seed(1222)
      set.seed(30)  # 30 42 14
      tree2 = rcoal(n)
      tree2 = rescaleSimmap(tree2, 1)
      tree2 = ladderize(tree2)

      if(input$tree == "Tree 1"){
        tree = tree1
      } else {
        tree = tree2
      }

      input$refresh


      ## simulate evolution along each edge
      set.seed(NULL)


      traits = data.frame(trait1 = fastBM(tree = tree, a = 1, sig2 = 0.25),
                          trait2 = fastBM(tree = tree, a = 1, sig2 = 0.25))

      fitlm   = lm(trait1 ~ trait2, data = traits)
      fitpgls = gls(trait1 ~ trait2, correlation = corBrownian(phy = tree),
                    data = traits, method = "ML")

      ################################################################################
      # Plot
      ################################################################################

      par(mfrow = c(1, 2))

      plot.phylo(tree, show.tip.label = FALSE, edge.width = 2)
      plot(trait1 ~ trait2, data = traits, pch = 16, col = "black")
      abline(a = coef(fitlm)[1], b = coef(fitlm)[2], col = "blue", lwd = 2)
      abline(a = coef(fitpgls)[1], b = coef(fitpgls)[2], col = "orange", lwd = 2)
      legend("topright", legend = c("ls", "pgls"), fill = c("blue", "orange"), bty = "n", cex = 1.2)



    }, width = 1000 ,height = 600
  )

}


shinyApp(ui = ui, server = server)

