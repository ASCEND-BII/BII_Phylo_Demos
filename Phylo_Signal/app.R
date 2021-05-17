library("shiny")
library("phytools")
library("geiger")

# Define UI
ui <- fluidPage(

    # Application title
    titlePanel("Phylogenetic Signal"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId  = "n_sp",
                        label    = "Number of species:",
                        choices  = c(10, 25, 50, 100),
                        multiple = FALSE),

            sliderInput(inputId = "lambda",
                        label   = "Lambda",
                        min     = 0,
                        max     = 1,
                        value   = 1,
                        step    = 0.1),
            checkboxInput(inputId = "phylosignal", value = FALSE,
                          label = "Compute Phylogenetic Signal?")
        ),
        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("plot")
        )
    )
)


server = function(input, output) {

    v_phylosig = function(tree, x, method){
        a = apply(x, MARGIN = 2, FUN = phytools::phylosig,
                  tree = tree, method = method, test = FALSE)
        unlist(sapply(a, `[`, 1))
    }

    seed   = 1234
    set.seed(seed)

    tree   = phytools::pbtree(n = 100)
    tree   = ladderize(tree)

    output$plot = renderPlot({

        t_drop       = sample(seq.int( 100 - as.numeric(input$n_sp) ))
        p_tree       = drop.tip(tree, t_drop)
        l_fun        = geiger::rescale(p_tree, "lambda")


        if(input$lambda > 0){
            l_tree = l_fun(input$lambda)
        } else {
            l_tree = l_fun(input$lambda + 1e-24)
        }

        set.seed(NULL)
        p_trait      = fastBM(l_tree, a = 0, sig2 = 1, nsim = 21)
        p_null_trait = fastBM(l_fun(0), a = 0, sig2 = 1, nsim = 20)


        if(input$phylosignal){

            lambda_est      = v_phylosig(tree   = p_tree,
                                         x      = p_trait,
                                         method = "lambda")
            lambda_null_est = v_phylosig(tree   = p_tree,
                                         x      = p_null_trait,
                                         method = "lambda")

            l = data.frame(type = rep(c("Estimate","Sim Lambda", "Star Tree"),
                                      c(1, 20, 20)),
                           lambda = c(lambda_est, lambda_null_est))

            k_est      = v_phylosig(tree   = p_tree,
                                    x      = p_trait,
                                    method = "K")
            k_null_est = v_phylosig(tree   = p_tree,
                                    x      = p_null_trait,
                                    method = "K")

            k = data.frame(type = rep(c("Estimate","Sim K", "Star Tree"),
                                      c(1, 20, 20)),
                           K = c(k_est, k_null_est))
        }


        if(input$phylosignal){
            m = matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
        } else {
            m = matrix(c(1, 2), nrow = 2, ncol = 1)
        }
        layout(m)

        plotTree.wBars(tree   = p_tree,
                       x      = p_trait[ , 1],
                       col    = sapply(p_trait[ , 1],
                                       function(x) if(x>=0) "blue" else "orange"),
                       mar    = c(3, 3, 5, 3))
        title("Tree and Traits")

        plotTree(l_tree, ftype = "off", mar    = c(3, 3, 5, 3))
        title("Tree scaled by lambda")

        if(input$phylosignal){
            par(mar = c(5, 4, 10, 3) )

            boxplot(lambda ~ type, data = l, main = "Pagel's Lambda", las = 2,
                    xlab = NULL, ylab = "lambda")
            abline(h = input$lambda, lty = 2, col = "red")
            boxplot(K ~ type, data = k, main = "Blomberg's K", las = 2,
                    xlab = NULL, ylab = "K")
        }


    }, width = 500, height = 800)
}

# Run the application
shinyApp(ui = ui, server = server)
