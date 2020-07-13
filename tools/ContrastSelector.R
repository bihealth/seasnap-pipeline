#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

read_cov <- function(my_config) {
    my_cov <- NULL
    if (!is.null(my_config$experiment$covariate_file$salmon)) {
        my_cov <- read.table(my_config$experiment$covariate_file$salmon, sep="\t", header=1, stringsAsFactors=FALSE)
    } else {
        if (!is.null(my_config$experiment$covariate_file$star)) {
            my_cov <- read.table(my_config$experiment$covariate_file$star, sep="\t", header=1, stringsAsFactors=FALSE)
        }
    }
    if (is.null(my_cov)) return(NULL)
    columns <- names(my_config$experiment$columns)
    my_cov <- my_cov[,columns]
    for (i in names(my_cov)) {
        my_cov[,i] <- factor(as.character(my_cov[,i]), levels=my_config$experiment$columns[[i]])
    }
    my_cov$group <- interaction(my_cov, drop=TRUE, sep=":")
    my_cov
}

build_dummy_dds <- function(X, my_cov, N=1000) {
    dummy <- DESeq2::DESeqDataSetFromMatrix(
        countData=matrix(sample.int(1000, size=N*nrow(my_cov), replace=TRUE),
                         nrow=N, ncol=nrow(my_cov), dimnames=list(sprintf("feature_%d", 1:N), my_cov$sample_id)),
        colData=my_cov,
        design=X
    )
    DESeq2::DESeq(dummy, betaPrior=FALSE, test="Wald", quiet=TRUE, fitType="local")
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    # Application title
    titlePanel("Model & constrast selection"),
   
    sidebarLayout(
        sidebarPanel(
            fileInput(
                "config",
                "Choose config file:",
                accept=c("text/yaml", ".yaml", ".yml", ".csv", ".tsv", ".txt")
            ),
            tags$hr(),
            textInput(
                "design",
                "Experimental design",
                value="~ group"
            ),
            actionButton(
                "compute",
                "Compute contrasts")
        ),
        mainPanel(
            tableOutput("cov"),
            tableOutput("contrasts")
            # verbatimTextOutput("dump", placeholder=TRUE)
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   config <- reactive({
       fn <- input$config
       if (is.null(fn)) return(NULL)
       my_config <- yaml::read_yaml(fn$datapath)
       my_config
   })

   Cov <- reactive({
       my_config <- config()
       if (is.null(my_config)) return(NULL)
       read_cov(my_config)
   })
   
   Design <- reactive({
       my_cov <- Cov()
       if (is.null(my_cov)) return(NULL)
       my_design <- input$design
       model.matrix(as.formula(my_design), data=my_cov)
   })
   
   Contrasts <- eventReactive(input$compute, {
       my_cov <- Cov()
       if (is.null(my_cov)) return(NULL)
       X <- Design()
       if (is.null(X)) return(NULL)

       my_cov$group <- factor(my_cov$group)

       dummy <- build_dummy_dds(X, my_cov)
          
       result_names <- DESeq2::resultsNames(dummy)
         
       # In Z: mapping between group conditions and test model
       Z <- model.matrix(~ -1 + group, data=my_cov)
       Z <- solve(t(Z) %*% Z) %*% t(Z) %*% X

       n <- (nrow(Z)*(nrow(Z)-1))/2
       tmp <- data.frame(Numerator=rep(as.character(NA), n), Denominator=rep(as.character(NA), n),
                         Coefficient=rep(as.character(NA), n),
                         Contrast=rep(as.character(NA), n),
                         stringsAsFactors=FALSE)
       n <- 0
       for (j in 2:nrow(Z)) for (i in 1:(j-1)) {
           n <- n+1
           tmp$Numerator[n]   <- levels(my_cov$group)[j]
           tmp$Denominator[n] <- levels(my_cov$group)[i]
           tmp$Coefficient[n] <- ""
           d <- Z[j,] - Z[i,]
           if (sum(abs(d)<1.0e-7) == length(d)-1) tmp$Coefficient[n] <- result_names[which(abs(d)>1.0e-7)]
           if (all(abs(d)<1.0e-7)) d <- "Same model value"
           tmp$Contrast[n]    <- paste("[ ", paste(d, collapse=", "), " ]", sep="")
       }

       tmp
   })
   
   output$dump <- renderText(({
       my_config <- config()
       if (is.null(my_config)) return(NULL)
       yaml::as.yaml(my_config)
   }))
   
   output$cov <- renderTable({
       my_cov <- Cov()
       if (is.null(my_cov)) return(NULL)
       explanatory_variables <- colnames(my_cov)[-ncol(my_cov)]
       df <- data.frame("Explanatory variable"="All", Value="", Reference=FALSE, "Number of experiment"=nrow(my_cov),
                        stringsAsFactors=FALSE, check.names=FALSE)
       for (i in c(explanatory_variables, "group")) {
           tmp <- table(my_cov[,i])
           tmp <- data.frame("Explanatory variable"=i, Value=names(tmp), Reference=(1:length(tmp))==1, "Number of experiment"=as.integer(tmp),
                             stringsAsFactors=FALSE, check.names=FALSE)
           df <- rbind(df, tmp)
       }
       df
   })
   
   output$contrasts <- renderTable({
       my_table <- Contrasts()
       if (is.null(my_table)) return(NULL)
       my_table
   })
}

# Run the application 
app <- shinyApp(ui = ui, server = server)

