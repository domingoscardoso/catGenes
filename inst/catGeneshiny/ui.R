shinyUI <- shinydashboard::dashboardPage(skin = "red",
                                         shinydashboard::dashboardHeader(title = "catGenes"),
                                         shinydashboard::dashboardSidebar(
                                           shinydashboard::sidebarMenu(
                                             shinydashboard::menuItem('Search Options', tabName = 'search',
                                                                      icon = icon('info-circle')),
                                             shinydashboard::menuItem('Result', tabName = 'results',
                                                                      icon = icon('table')),
                                             mainPanel(img(src = "DBOS2_logo.png", height = 120, width = 170))
                                           )
                                         ),
                                         shinydashboard::dashboardBody(
                                           shinydashboard::tabItems(
                                             shinydashboard::tabItem(tabName = 'search',
                                                                     h2('Enter concatenation and writing parameters'),
                                                                     #h3('Input concatenation and saving parameters'),
                                                                     fluidRow(
                                                                       shinydashboard::box(title = 'Select DNA alignments and concatenation parameters',
                                                                                           status = 'primary', width = 4,
                                                                                           fileInput('datafiles', 'Choose DNA alignments files',
                                                                                                     multiple = TRUE,
                                                                                                     accept = c('.nex', '.txt')),
                                                                                           p('You must select more than one DNA alignment'),
                                                                                           p('DNA alignments files must be named simply as "ITS.nex", "rbcL.nex", "COX1.nex", etc'),
                                                                                           hr(),
                                                                                           radioButtons('multiaccessions', 'DNA alignments include any species with multiple accessions?',
                                                                                                        choices = c('No' = 1, 'Yes' = 2),
                                                                                                        selected = 1),
                                                                                           uiOutput('multiaccessions'),
                                                                                           conditionalPanel(condition = 'input.multiaccessions == 2',
                                                                                                            radioButtons('maximizespp', 'Maximize species coverage?',
                                                                                                                         choices = c('No' = FALSE, 'Yes' = TRUE),
                                                                                                                         selected = TRUE)),
                                                                                           uiOutput('maximizespp'),
                                                                                           radioButtons('missdata', 'Concatenation should keep taxa with any missing gene?',
                                                                                                        choices = c('No' = FALSE, 'Yes' = TRUE),
                                                                                                        selected = FALSE),
                                                                                           uiOutput('missdata'),
                                                                                           radioButtons('shortaxlabel', 'Remove collector or GenBank numbers from concatenated dataset?',
                                                                                                        choices = c('No' = FALSE, 'Yes' = TRUE),
                                                                                                        selected = TRUE),
                                                                                           uiOutput('shortaxlabel'),
                                                                                           radioButtons('outbutton', 'Define outgroup taxa?',
                                                                                                        choices = c('No' = 1, 'Yes' = 2),
                                                                                                        selected = 1),
                                                                                           conditionalPanel(condition = 'input.outbutton == 2',
                                                                                                            textInput('outgroups', h5(strong("Enter each outgroup taxa separately")), ""),
                                                                                                            actionButton("addoutgroup","Add"),
                                                                                                            mainPanel(htmlOutput("text")))
                                                                       ),
                                                                       shinydashboard::box(title = 'Choose file format to write the concatenated dataset',
                                                                                           status = 'primary', width = 3,
                                                                                           hr(),
                                                                                           radioButtons('nexus', 'Save nexus format?',
                                                                                                        choices = c('No' = 1, 'Yes' = 2),
                                                                                                        selected = 2),
                                                                                           uiOutput('nexus'),
                                                                                           conditionalPanel(condition = 'input.nexus == 2',
                                                                                                            radioButtons('interleave', 'Interleaved matrix?',
                                                                                                                         choices = c('No' = FALSE, 'Yes' = TRUE),
                                                                                                                         selected = TRUE),
                                                                                                            uiOutput('interleave'),
                                                                                                            radioButtons('charset', 'Write charset block?',
                                                                                                                         choices = c('No' = FALSE, 'Yes' = TRUE),
                                                                                                                         selected = TRUE),
                                                                                                            uiOutput('charset')),
                                                                                           radioButtons('phylip', 'Save phylip format?',
                                                                                                        choices = c('No' = 1, 'Yes' = 2),
                                                                                                        selected = 1),
                                                                                           uiOutput('phylip'))
                                                                     ),
                                                                     fluidRow(
                                                                       shinyjs::useShinyjs(),
                                                                       shinydashboard::box(title = 'Run Concatenation', status = 'warning',
                                                                                           width = 12,
                                                                                           actionButton('run_search', 'Run catGenes on uploaded datasets')
                                                                       )
                                                                     )
                                             ),
                                             shinydashboard::tabItem(tabName = 'results',
                                                                     h2('Download result'),
                                                                     fluidRow(
                                                                       shinydashboard::box(title = 'Download Concatenated Dataset',
                                                                                           status = 'info',
                                                                                           width = 5,
                                                                                           downloadButton('down_results', 'Download'),
                                                                                           h2(''),
                                                                                           conditionalPanel(condition = 'input.phylip == 2',
                                                                                                            downloadButton('down_partition', 'Download RAxML partition file'))
                                                                       )
                                                                     )

                                             )
                                           )
                                         )#,
                                         #plotOutput('bar')
)
