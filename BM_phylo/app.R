library("shiny")
library("phytools")


ui = fluidPage(

  titlePanel("BM on Tree"),

  # UI Controls

  fluidRow(
    column(1,
           actionButton(inputId = "refresh", label = "Evolve!")
    ),
    column(11,
           plotOutput("mainplot")
    )
  )
)


server = function(input, output){

  output$mainplot = renderPlot(
    {
      input$refresh

      set.seed(12220)

      t    = 100  # total time
      n    = 5    # total taxa
      b    = (log(n) - log(2)) / t
      sig2 = 1

      tree = pbtree(b = b, n = n, t = t, type = "discrete")
      tree = ape::ladderize(tree)

      ## reorder the edges of the tree for pre-order traversal
      cw   = reorder(tree)

      ################################################################################
      # Sim data
      ################################################################################

      ## simulate evolution along each edge
      set.seed(NULL)

      X = lapply(cw$edge.length, function(x) c(0, cumsum(rnorm(n = x, sd = sqrt(sig2)))))

      ## now simulate on the tree
      ll = cw$edge.length + 1

      for(i in 1:nrow(cw$edge)) {
        pp=which(cw$edge[, 2] == cw$edge[i, 1])

        if (length(pp) > 0){
          X[[i]] = X[[i]] + X[[pp]][ll[pp]]
        } else {
          X[[i]] = X[[i]] + X[[1]][1]
        }
      }

      ## get the starting and ending points of each edge for plotting
      H <- nodeHeights(cw)

      ################################################################################
      # Plot
      ################################################################################

      par(mfrow = c(1,2))

      col = c("goldenrod1", "deeppink", "red", "firebrick", "orange", "dark green", "deepskyblue", "blue")
      lwd = 2

      plot(cw, edge.width = lwd, show.tip.label = FALSE, edge.color = col)
      axis(1)
      mtext("time", 1, line = 3)

      ## plot the simulation
      plot(H[1, 1], X[[1]][1], ylim = range(X), xlim = range(H),
           xlab = "time", ylab = "phenotype", cex = 0.6, pch = 16, col = "black")


      for (i in 1:length(X)){
        lines(H[i, 1]:H[i, 2], X[[i]], lwd = lwd, col = col[i])
        points(H[i, 1], X[[i]][1], cex = 0.6, pch = 16, col = "black")
      }


      ## Add tip labels
      yy = sapply(1:length(cw$tip.label),
                  function(x, y){
                    which(x == y)
                  },
                  y = cw$edge[ , 2])

      yy = sapply(yy,
                  function(x, y){
                    y[[x]][length(y[[x]])]
                  },
                  y = X)

      #text(x = max(H), y = yy, tree$tip.label, col = "red", font = 3)
      points(x = rep(max(H), length(yy)), y = yy, cex = 1, pch = 16, col = "black")

    }, width = 1000 ,height = 600
  )

}


shinyApp(ui = ui, server = server)

