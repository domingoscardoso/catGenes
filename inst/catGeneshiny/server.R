shinyServer <- function(input, output, session) {

    rv <- reactiveValues(outgroup = NULL)

    observeEvent(input$addoutgroup, {
        rv$outgroup <- c(rv$outgroup, input$outgroups)
        updateTextInput(session, "outgroups", value="")
    })

    output$text <- renderUI({
        #<br> is a break; <b> is bold; </b> ends bold
        str1 <- paste("<b>Total outgroups selected:</b>", length(rv$outgroup))
        str2 <- paste("<b>Outgroups selected:<br></b>", paste(rv$outgroup, collapse = "<br>"))
        HTML(paste(str1, str2, sep = "<br/>"))
    })

    catGenes_result <- eventReactive(input$run_search, {

        req(input$datafiles)

        genedatasets = list()
        for(i in 1:length(input$datafiles[, 1])){
            genedatasets[[i]] <- read.nexus.data(file = input$datafiles[[i, 'datapath']])
        }

        # Getting the name of the files stored within the input datafiles
        data_temp <- input$datafiles
        names(genedatasets) <- gsub("[.].*", "", data_temp$name)


        missdat <- input$missdata
        shortax <- input$shortaxlabel
        maxispp <- input$maximizespp

        if(input$outbutton == 1) {
            outgrps <- NULL
        } else {
            outgrps <- rv$outgroup
        }

        if(input$multiaccessions == 1) {
            catfullGenes(genedatasets,
                         shortaxlabel = shortax,
                         missdata = missdat,
                         outgroup = outgrps)
        } else {
            catmultGenes(genedatasets,
                         maxspp = maxispp,
                         shortaxlabel = shortax,
                         missdata = missdat,
                         outgroup = outgrps)
        }
    })

    # Enabling the donwload button only after the concatenation in done
    observe({
        if (input$run_search > 0) {
            Sys.sleep(1)
            # enable the download button
            shinyjs::enable("down_results")
            # change the html of the download button
            shinyjs::html("down_results",
                          sprintf("<i class='fa fa-download'></i>
                              Download (file size: %s)",
                                  round(runif(1, 1, 10000))
                          )
            )
            shinyjs::enable("down_partition")
            shinyjs::html("down_partition",
                          sprintf("<i class='fa fa-download'></i>
                              Download RAxML partition file (file size: %s)",
                                  round(runif(1, 1, 10000))
                          )
            )
        }
    })

    # This part is not working
    # Progress bar of the concatenation
    # output$bar <- renderPlot({
    #     withProgress(message = 'Concatenating DNA alignments...', {
    #         catGenes_result()
    #     })
    # })

    nexus <- reactive({
        nexus <- input$nexus
    })
    interl <- reactive({
        interl <- input$interleave
    })
    charset <- reactive({
        charset <- input$charset
    })
    phylip <- reactive({
        phylip <- input$phylip
    })

    output$down_results <- downloadHandler(
        filename = function() {
            if(nexus() == 2 & phylip() == 1){
                paste0('concatenated_dataset', '.nex')
            }else{
                paste0('concatenated_dataset', '.phy')
            }
        },
        content = function(file) {
            if(nexus() == 2 & phylip() == 1){
                writeNexus(catGenes_result(), file,
                           bayesblock = charset(),
                           interleave = interl())
            }else{
                writePhylip(catGenes_result(), file,
                            catalignments = TRUE,
                            partitionfile = FALSE)
            }
        }
    )

    output$down_partition <- downloadHandler(
        filename = function() {
            if(phylip() == 2){
                paste0('RAxML_partition_file', '.txt')
            }
        },
        content = function(file) {
            if(phylip() == 2){
                writePhylip(catGenes_result(), file,
                            catalignments = FALSE,
                            partitionfile = TRUE)
            }
        }
    )
    shinyjs::disable("down_results")
    shinyjs::disable("down_partition")
}
