library(shiny)
library(shinyjs)
library(tidyverse)
library(rlang)
library(DT)
library(shinycssloaders)
library(profvis)


result <- read_tsv("result.tsv") %>% mutate_if(is.integer, funs(as.double(.)))
chembl <- read_tsv("chembl.tsv", trim_ws = F)
drugb <- read_tsv("drugb.tsv", na = "NA")
ipro <- read_tsv("ipro.tsv", col_types = "cccn", na = "NA")
prot_types <- read_tsv("prot_types.tsv")
screenS <- read_tsv("screenS.tsv")
shdf <- read_tsv("shdf.tsv")
numprot <- shdf %>% select(MS3_id) %>% count() %>% pull()
wght <- rep(1, ncol(shdf)+2)
names(wght) <- c(colnames(shdf), "one_cmpd", "more_cmpds")
more_cmpds <- 5
ptdd <- c("All identifiers",paste(prot_types$ipro_id, prot_types$ipro_descr, sep = " "))
equalWeight <- 5

#set preselected cut-offs for blast coverage and similarity to determine homology
simil_cutoff <- 80
cov_cutoff <- 50

cutOff <- function(df, cut=0)
{
if_else(df >= cut, 1, 0)
}

# Define UI for Shiny app ----

ui <- tagList(useShinyjs(),
              
              div(id = "navbar", navbarPage(title = "S. haematobium drug target ranking",
                                            
                                            selected = "Transcription",
                                            
                                            tabPanel("Transcription",
                                                     {
                                                       fluidRow(
                                                         column(5,
                                                                
                                                                sliderInput(inputId = "trAdS", "Transcribed in adults?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "trAdRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "require", inline = T)
                                                         )
                                                       )
                                                       
                                                     }
                                            ),
                                            
                                            tabPanel("Orthology", 
                                                     {
                                                       list(
                                                       fluidRow(
                                                         column(12,
                                                                
                                                                numericInput(inputId = "sim_cut", label = "Similarity cut-off (%)", value = simil_cutoff, min = 0, max = 100),
                                                                
                                                                numericInput(inputId = "cov_cut", label = "Coverage cut-off (%)", value = cov_cutoff, min = 0, max = 100)
                                                         )
                                                       ), 
                                                       
                                                       fluidRow(
                                                         column(3,
                                                                sliderInput(inputId = "HomCelS", "Ortholog in C. elegans?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "HomCelRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "HomDmelS", "Ortholog in D. melanogaster?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "HomDmelRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "HomMusS", "Ortholog in M. musculus?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "HomMusRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "exclude", inline = T)
                                                                
                                                         ),
                                                         
                                                         column(3,
                                                                
                                                                sliderInput(inputId = "HomSjS", "Ortholog in S. japonicum?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 7, step = 1),
                                                                
                                                                radioButtons(inputId = "HomSjRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "HomSmS", "Ortholog in S. mansoni?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 8, step = 1),
                                                                
                                                                radioButtons(inputId = "HomSmRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T)
                                                         ),
                                                         
                                                         column(3,   
                                                                
                                                                sliderInput(inputId = "HomCsinS", "Ortholog in C. sinensis?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 6, step = 1),
                                                                
                                                                radioButtons(inputId = "HomCsinRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "HomOvivS", label = "Ortholog in O. viverrini?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 6, step = 1),
                                                                
                                                                radioButtons(inputId = "HomOvivRB", label = "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "HomFhepS", label = "Ortholog in F. hepatica?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 8, step = 1),
                                                                
                                                                radioButtons(inputId = "HomFhepRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T)
                                                         ),
                                                         
                                                         column(3,
                                                                
                                                                sliderInput(inputId = "HomSPS", "Ortholog in SwissProt?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = 7, step = 1),
                                                                
                                                                radioButtons(inputId = "HomSPRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                sliderInput(inputId = "quant_slider", "No ortholog or similarity <75th percentile of H. sapiens?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "QuantHsRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "require", inline = T)
                                                         ) 
                                                       )
                                                      )
                                                     }
                                            ),
                                            
                                            tabPanel("Essentiality", 
                                                     {  fluidRow(
                                                       column(6,
                                                              
                                                              sliderInput(inputId = "letCelS", "Lethal phenotype C. elegans?",
                                                                          min = 1, max = 10, width = "100%",
                                                                          value = 9, step = 1),
                                                              
                                                              radioButtons(inputId = "letCelRB", "",
                                                                           c("WEIGHTED" = "weighted",
                                                                             "IGNORE" = "ignore",
                                                                             "REQUIRE" = "require",
                                                                             "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                              
                                                              
                                                              sliderInput(inputId = "letDmelS", "Lethal phenotype D. melanogaster?",
                                                                          min = 1, max = 10, width = "100%",
                                                                          value = 9, step = 1),
                                                              
                                                              radioButtons(inputId = "letDmelRB", "",
                                                                           c("WEIGHTED" = "weighted",
                                                                             "IGNORE" = "ignore",
                                                                             "REQUIRE" = "require",
                                                                             "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                              
                                                              
                                                              sliderInput(inputId = "letMusS", "Lethal phenotype M. musculus?",
                                                                          min = 1, max = 10, width = "100%",
                                                                          value = equalWeight, step = 1),
                                                              
                                                              radioButtons(inputId = "letMusRB", "",
                                                                           c("WEIGHTED" = "weighted",
                                                                             "IGNORE" = "ignore",
                                                                             "REQUIRE" = "require",
                                                                             "EXCLUDE" = "exclude"), selected = "exclude", inline = T
                                                              )
                                                              
                                                       ),
                                                       
                                                       column(6,
                                                              
                                                              sliderInput(inputId = "chokeS", "Unique KEGG term ('choke-point')?",
                                                                          min = 1, max = 10, width = "100%",
                                                                          value = 8, step = 1),
                                                              
                                                              radioButtons(inputId = "chokeRB", "",
                                                                           c("WEIGHTED" = "weighted",
                                                                             "IGNORE" = "ignore",
                                                                             "REQUIRE" = "require",
                                                                             "EXCLUDE" = "exclude"),selected = "weighted", inline = T),
                                                              
                                                              
                                                              sliderInput(inputId = "uniqIproS", "Unique InterPro identifier?",
                                                                          min = 1, max = 10, width = "100%",
                                                                          value = 7, step = 1),
                                                              
                                                              radioButtons(inputId = "uniqIproRB", "",
                                                                           c("WEIGHTED" = "weighted",
                                                                             "IGNORE" = "ignore",
                                                                             "REQUIRE" = "require",
                                                                             "EXCLUDE" = "exclude"), selected = "weighted", inline = T)
                                                              
                                                       )
                                                     )
                                                     }
                                            ),
                                            
                                            tabPanel("Annotation", 
                                                     {
                                                       fluidRow(
                                                         column(12,
                                                                selectInput("prot_type", label = "Select InterPro protein identifier(s)", selected = ptdd[1], choices = ptdd, multiple = T, width = '100%'),
                                                                radioButtons(inputId = "groupcheck", "", choices = c("ANY OF SELECTED" = "requireAny", "ALL OF SELECTED" = "requireAll", "NONE OF SELECTED" = "exclude"), selected = "requireAny", inline = T)
                                                         )
                                                       )
                                                     }
                                            ),
                                            
                                            tabPanel("Drugs", 
                                                     {
                                                       fluidRow(
                                                         column(6,
                                                                
                                                                h2("ChEMBL"),
                                                                
                                                                sliderInput(inputId = "ChCmpdOneS", "At least one associated compound in ChEMBL?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "ChCmpdOneRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "require", inline = T),
                                                                
                                                                sliderInput(inputId = "ChCmpdMoreS", paste(more_cmpds, " or more associated compounds in ChEMBL?"),
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "ChCmpdMoreRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                
                                                                sliderInput(inputId = "chPhase", "Drug phase",
                                                                            min = 0, max = 4, width = "100%",
                                                                            value = c(2,4), step = 1),
                                                                
                                                                
                                                                sliderInput(inputId = "Ch_ro5_viol", "Max. no. of rule-of-5 violations",
                                                                            min = 0, max = 3, width = "100%",
                                                                            value = 0, step = 1),
                                                                
                                                                h5("Other requirements"),
                                                                
                                                                checkboxInput(inputId = "Ch_ro3_pass", label = "Require rule-of-3 pass"),
                                                                
                                                                checkboxInput(inputId = "InNatProd", label = "Require compounds to be natural products"),
                                                                
                                                                checkboxInput(inputId = "InScreenedSchisto", label = "Require drug screening data for schistosomes", value = TRUE)
                                                                
                                                         ),
                                                         
                                                         
                                                         column(6,
                                                                
                                                                h2("DrugBank"),
                                                                
                                                                sliderInput(inputId = "DBCmpdOneS", "At least one associated compound in DrugBank?",
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "DBCmpdOneRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "require", inline = T),
                                                                
                                                                sliderInput(inputId = "DBCmpdMoreS", paste(more_cmpds, " or more associated compounds in DrugBank?"),
                                                                            min = 1, max = 10, width = "100%",
                                                                            value = equalWeight, step = 1),
                                                                
                                                                radioButtons(inputId = "DBCmpdMoreRB", "",
                                                                             c("WEIGHTED" = "weighted",
                                                                               "IGNORE" = "ignore",
                                                                               "REQUIRE" = "require",
                                                                               "EXCLUDE" = "exclude"), selected = "weighted", inline = T),
                                                                
                                                                checkboxGroupInput(inputId = "DBstatus", 
                                                                                   choiceValues = c("approved",
                                                                                   "experimental",
                                                                                   "illicit",
                                                                                   "investigational",
                                                                                   "nutraceutical",
                                                                                   "withdrawn"),
                                                                                   choiceNames = c("Approved",
                                                                                                   "Experimental",
                                                                                                   "Illicit",
                                                                                                   "Investigational",
                                                                                                   "Nutraceutical",
                                                                                                   "Withdrawn"),
                                                                                   selected = "approved", label = "Include drugs with status/properties"
                                                                                   
                                                                                   ),
                                                                
                                                                h5("Other requirements"),
                                                                
                                                                checkboxInput(inputId = "InMDDR", label = "Require MDDR-likeness"),
                                                                
                                                                checkboxInput(inputId = "InRo5", label = "Require rule-of-5 pass", value = TRUE)

                                                              
                                                         )
                                                       )
                                                     }
                                            ),
                                            
                                            tabPanel("Summary and Analysis", 
                                                     {
                                                       list(
                                                         fluidRow(
                                                           column(4,
                                                                  actionButton("reset_input", "Reset weights"),
                                                                  actionButton("calc", "Calculate ranking"),
                                                                  tableOutput(outputId = "outTable")
                                                           ),
                                                           column(8,
                                                                  textOutput(outputId = "numRes"),
                                                                  plotOutput(outputId = "scoredistr") %>% withSpinner(),
                                                                  uiOutput("downloads"),
                                                                  uiOutput(outputId = "DBversions")
                                                                  #,
                                                                  #dataTableOutput(outputId = "drugSelTable")
                                                           )
                                                         ),
                                                         
                                                         fluidRow(
                                                           column(12,
                                                                  dataTableOutput(outputId = "protSelTable") %>% withSpinner()
                                                           )
                                                         )
                                                       )
                                                     }
                                            )
                                            ,
                                            
                                            navbarMenu("Drug Results",
                                                       tabPanel("ChEMBL",
                                                                
                                                                fluidRow(
                                                                  column(12,
                                                                         uiOutput("CHout"),
                                                                         uiOutput("CHdl"),
                                                                         dataTableOutput(outputId = "drugSelTableCH") %>% withSpinner()
                                                                  )
                                                                )
                                                                
                                                       ),
                                                       
                                                       tabPanel("DrugBank",
                                                                
                                                                fluidRow(
                                                                  column(12,
                                                                         uiOutput("DBout"),
                                                                         uiOutput("DBdl"),
                                                                         dataTableOutput(outputId = "drugSelTableDB") %>%  withSpinner()
                                                                  )
                                                                )
                                                                
                                                       )
                                            )
  
              )
              )
)

# Define server side ----
server <- function(input, output, session) {
  
#Gather data for selected proteins  
  gatherProtSelect <- eventReactive(eventExpr = input$calc,
                                    {
                                      if("All identifiers" %in% input$prot_type)
                                      {
                                        selectProt <- c("",ptdd[-1])
                                      }
                                      else
                                      {
                                        selectProt <- as.vector(input$prot_type)
                                      }
                                      
                                      #Get just the IDs (minus the actual description) of the selected protein identifiers
                                      selectProt <- str_sub(selectProt, 1,9)
                                      
                                      if(input$groupcheck == "requireAny")
                                      {
                                        
                                        id_tbl <- ipro %>%
                                          filter(ipro_id %in% selectProt) %>%
                                          select(MS3_id) %>%
                                          distinct()
                                      }
                                      else
                                      {
                                        
                                        if(input$groupcheck == "requireAll")
                                        {
                                          #this statement returns MS3_id
                                          #find all MS3_id that have a row for all of the values in selectProt
                                          id_tbl <- ipro %>% 
                                            group_by(MS3_id) %>%
                                            filter(., all(selectProt %in% ipro_id)) %>%
                                            select(MS3_id) %>%
                                            distinct() %>% ungroup()
                                        }
                                        else
                                        {
                                          
                                          if(input$groupcheck == "exclude")
                                          {
                                            id_tbl <- ipro %>%
                                              filter(ipro_id %in% selectProt) %>%
                                              select(MS3_id) %>%
                                              distinct() %>%
                                              anti_join(ipro, ., by = "MS3_id") %>%
                                              select(MS3_id) %>%
                                              distinct()
                                          }
                                        }
                                      }
                                      
                                      #subset the result based on the selection made by the user
                                      subsResult <- inner_join(id_tbl, result, by="MS3_id") %>% 
                                        mutate_at(vars(ends_with("_simil")), funs(replace_na(., 0))) %>% 
                                        mutate_at(vars(ends_with("_cov")), funs(replace_na(., 0))) %>% 
                                        #if the cov or simil value in the result table is bigger than the selected cut-off, set it to 1
                                        # else set it to 0
                                        mutate_at(vars(ends_with("_simil")), funs(cutOff(., cut = input$sim_cut))) %>% 
                                        mutate_at(vars(ends_with("_cov")), funs(cutOff(., cut = input$cov_cut/100))) %>% 
                                        gather(type, value, union(ends_with("_simil"), ends_with("_cov"))) %>%
                                        separate(type, into = c("Species", "simCov")) %>%
                                        arrange(MS3_id, Species, simCov) %>%
                                        group_by(MS3_id, Species) %>%
                                        mutate(value = max(value)) %>%
                                        ungroup %>%
                                        filter(simCov == "cov") %>%
                                        unite(Species, simCov, col = "column1", sep = "_") %>%
                                        spread(column1, value) %>% 
                                        select(MS3_id,
                                               male_tpm,
                                               cel_leth,
                                               dmel_leth,
                                               mus_leth,
                                               chokepoint,
                                               uniqIpro,
                                               cel_cov,
                                               dmel_cov,
                                               mus_cov,
                                               sjap_cov,
                                               sman_cov,
                                               closi_cov,
                                               opvi_cov,
                                               fahe_cov,
                                               sp_cov,
                                               hsa_lower_75_qtl,
                                               ChCmpdOne,
                                               ChCmpdMore,
                                               DBCmpdOne,
                                               DBCmpdMore
                                              )
                                        
                                      tibOut <- tibble(
                                        
                                        name = c(
                                          "male_tpm",
                                          "cel_leth",
                                          "dmel_leth",
                                          "mus_leth",
                                          "chokepoint", 
                                          "uniqIpro",
                                          "cel_cov",
                                          "dmel_cov",
                                          "mus_cov",
                                          "sjap_cov",
                                          "sman_cov",
                                          "closi_cov",
                                          "opvi_cov",
                                          "fahe_cov",
                                          "sp_cov",
                                          "hsa_lower_75_qtl",
                                          "ChCmpdOne",
                                          "ChCmpdMore",
                                          "DBCmpdOne",
                                          "DBCmpdMore"
                                        ),
                                        
                                        
                                        weighted = c(
                                          input$trAdS, 
                                          input$letCelS,
                                          input$letDmelS,
                                          input$letMusS,
                                          input$chokeS,
                                          input$uniqIproS,
                                          input$HomCelS,
                                          input$HomDmelS,
                                          input$HomMusS,
                                          input$HomSjS,
                                          input$HomSmS,
                                          input$HomCsinS,
                                          input$HomOvivS,
                                          input$HomFhepS,
                                          input$HomSPS,
                                          input$quant_slider,
                                          input$ChCmpdOneS,
                                          input$ChCmpdMoreS,
                                          input$DBCmpdOneS,
                                          input$DBCmpdMoreS
                                        ),
                                        
                                        
                                        
                                        mode = as.character(c(
                                          input$trAdRB,
                                          input$letCelRB,
                                          input$letDmelRB,
                                          input$letMusRB,
                                          input$chokeRB,
                                          input$uniqIproRB,
                                          input$HomCelRB,
                                          input$HomDmelRB,
                                          input$HomMusRB,
                                          input$HomSjRB,
                                          input$HomSmRB,
                                          input$HomCsinRB,
                                          input$HomOvivRB,
                                          input$HomFhepRB,
                                          input$HomSPRB,
                                          input$QuantHsRB,
                                          input$ChCmpdOneRB,
                                          input$ChCmpdMoreRB,
                                          input$DBCmpdOneRB,
                                          input$DBCmpdMoreRB
                                        ))
                                        
                                        
                                      )
                                      
                                      
                                      newWeights <- filter(tibOut ,mode == "weighted") %>% select(., name, weighted) %>% spread(.,name,weighted)
                                      
                                      if(!is_empty(newWeights)){
                                        for (nm in names(newWeights)) {
                                          subsResult[nm] <- subsResult[nm]*as.double(newWeights[1,nm])
                                        }
                                      }
                                      
                                      
                                      required <- select(tibOut, name, mode) %>% filter(.,mode == "require") %>% pull(., var="name")
                                      if(!is_empty(required)){
                                        subsResult <- filter_at(subsResult, vars(required), all_vars(. > 0))  
                                        subsResult <- select(subsResult, -(one_of(required)))
                                      }
                                      
                                      
                                      excluded <- select(tibOut, name, mode) %>% filter(.,mode == "exclude") %>% pull(., var="name")
                                      if(!is_empty(excluded)){
                                        subsResult <- filter_at(subsResult, vars(excluded), all_vars(. == 0))
                                        subsResult <- select(subsResult, -(one_of(excluded)))
                                      }
                                      
                                      
                                      ignored <- select(tibOut, name, mode) %>% filter(.,mode == "ignore") %>% pull(., var="name")
                                      if(!is_empty(ignored)){
                                        subsResult <- select(subsResult, -(one_of(ignored)))
                                      }
                                      
                                      
                                      subsResult <- mutate(subsResult, score = rowSums(subsResult[,-1])) %>% arrange(.,desc(score))
                                      
                                      plotout <- ggplot(subsResult, aes(score)) + geom_bar(fill = "seagreen4") + labs(x = "Score", y = "Frequency") + theme_minimal() + theme(text = element_text(size = 16))
                                      
                                      NonDistRes <- as.character(subsResult %>% select(MS3_id) %>% count())
                                      
                                      list(table = subsResult, plot = plotout, text = paste(NonDistRes, sep = ""))
                                    }                                  
  )
  
  
  # Define output server side ----
  
  # Create table of selected parameters as an output on the UI
  output$outTable <- renderTable(
    makeParamsTable()
  )
   
  # Create table of selected parameters
  makeParamsTable <- reactive( 
    {
      
      outtbl <- tibble(
        
        Criterion = c(
          #"male_tpm",
          "Transcribed in adults?",
          #"cel_leth",
          "Lethal phenotype in C. elegans?",
          #"dmel_leth",
          "Lethal phenotype in D. melanogaster?",
          #"mus_leth",
          "Lethal phenotype in M. musculus?",
          #"chokepoint", 
          "KEGG 'choke-point'?",
          #"uniqIpro",
          "Unique InterPro identifier?",
          #"cel_cov",
          "C. elegans ortholog?",
          #"dmel_cov",
          "D. melanogaster ortholog?",
          #"mus_cov",
          "M. musculus ortholog?",
          #"sjap_cov",
          "S. japonicum ortholog?",
          #"sman_cov",
          "S. mansoni ortholog?",
          #"closi_cov",
          "C. sinensis ortholog?",
          #"opvi_cov",
          "O. viverrini ortholog?",
          #"fahe_cov",
          "F. hepatica ortholog?",
          #"sp_cov",
          "SwissProt ortholog?",
          #"hsa_lower_75_qtl",
          "Similarity to H. sapiens ortholog <75th percentile?",
          #"ChCmpdOne",
          "At least one associated compound in ChEMBL",
          #"ChCmpdMore",
          paste(more_cmpds, " or more associated compounds in ChEMBL?"),
          #"DBCmpdOne",
          "At least one associated compound in DrugBank",
          #"DBCmpdMore"
          paste(more_cmpds, " or more associated compounds in DrugBank?")
        ),
        
        Weight = c(
          input$trAdS, 
          input$letCelS,
          input$letDmelS,
          input$letMusS,
          input$chokeS,
          input$uniqIproS,
          input$HomCelS,
          input$HomDmelS,
          input$HomMusS,
          input$HomSjS,
          input$HomSmS,
          input$HomCsinS,
          input$HomOvivS,
          input$HomFhepS,
          input$HomSPS,
          input$quant_slider,
          input$ChCmpdOneS,
          input$ChCmpdMoreS,
          input$DBCmpdOneS,
          input$DBCmpdMoreS
        ),
        
        "Ranking mode" = factor(c(
          input$trAdRB,
          input$letCelRB,
          input$letDmelRB,
          input$letMusRB,
          input$chokeRB,
          input$uniqIproRB,
          input$HomCelRB,
          input$HomDmelRB,
          input$HomMusRB,
          input$HomSjRB,
          input$HomSmRB,
          input$HomCsinRB,
          input$HomOvivRB,
          input$HomFhepRB,
          input$HomSPRB,
          input$QuantHsRB,
          input$ChCmpdOneRB,
          input$ChCmpdMoreRB,
          input$DBCmpdOneRB,
          input$DBCmpdMoreRB
        ))
        
        
      )
      
      arrange(outtbl, UQ(sym("Ranking mode")))
      
      
    }
  )
  
  # Summary sentence of how many proteins were found that match the user criteria
  output$numRes <- renderText(
    {
      out <- gatherProtSelect()
      paste("Of ", numprot, " proteins, ", unlist(out$text), " match your protein filtering criteria, with the following score distribution:")
    }
  )


# Summary sentence of how many proteins were found that match the user criteria
  output$DBversions <- renderUI(
    {
      out <- gatherProtSelect()
      HTML(paste(
      "",
      "Results are calculated using the following database versions/releases:",
      "WormBase: WS262",
      "FlyBase: FB2017_06",
      "Ensembl: 91.38",
      "WormBase Parasite: WBPS9",
      "Swiss-Prot: 05/2017",
      "KEGG: 05/2017",
      "InterProScan: 5.15.54",
      "ChEMBL: 23",
      "DrugBank: 5-0-11",
      sep="<br/>"))
    }
  )

  
  # Make plot of the score distribution
  output$scoredistr <- renderPlot(
    {
      out <- gatherProtSelect()
      out$plot
    }
  )
  
  
  # Clean up and create a table for viewing on UI
  createTable <- reactive(
    {
      out <- gatherProtSelect()
      
      out$table <- out$table %>%
        select(MS3_id) %>%
        left_join(., ipro) %>%
        select(MS3_id, ipro_descr) %>%
        distinct() %>% 
        group_by(MS3_id) %>% 
        summarise(anno = paste0(ipro_descr, collapse = "; ")) %>% 
        ungroup() %>%
        select(MS3_id, anno) %>% 
        right_join(., out$table, by = "MS3_id") %>% 
        left_join(., select(shdf, MS3_id, gene_name_species, uniprot_id), by = "MS3_id") %>% 
        mutate(gene_name_species = paste(gene_name_species, uniprot_id, sep = " ")) %>% 
        mutate_at(., vars(gene_name_species), funs(if_else(. == "NA NA", "", .))) %>% 
        select(-uniprot_id) %>% 
        select(
          MS3_id,
          score,
          anno,
          gene_name_species,
          everything()
          )
      
    }
  )

# Create two different drug tables on click
  DrugCHTable <- eventReactive(eventExpr = input$show_drug_chembl,
                               {
                                 tab <- createTable() %>% slice(1:isolate(input$num_dl_chembl)) %>%
                                   select(MS3_id) %>%
                                   left_join(., chembl) %>% 
                                   select(
                                     MS3_id,
                                     acc_no,
                                     chembl_target_id,
                                     ro5_viol,
                                     ro3_pass,
                                     target_name,
                                     cmpd_name,
                                     cmpd_syn,
                                     chembl_cmpd_id,
                                     smiles,
                                     phase,
                                     nat_prod,
                                     first_approv,
                                     pat_no,
                                     pat_exp
                                   )
                                 
                                 
                                tab <-  tab %>% select(-c(cmpd_name, cmpd_syn, pat_no, pat_exp)) %>% distinct() %>% 
                                   left_join(.,
                                              {
                                                tab %>% select(chembl_cmpd_id, cmpd_name, cmpd_syn, pat_no, pat_exp) %>% 
                                                  distinct() %>% 
                                                  mutate(cmpd_name = str_to_title(cmpd_name)) %>% 
                                                  mutate(cmpd_syn = str_to_title(cmpd_syn)) %>% 
                                                  group_by(chembl_cmpd_id) %>% 
                                                  mutate(cmpd_name = paste0(unique(cmpd_name), collapse = "; ")) %>% distinct() %>% 
                                                  mutate(cmpd_syn = paste0(unique(cmpd_syn), collapse = "; ")) %>% distinct() %>% 
                                                  mutate(pat_no = str_replace(pat_no, "^None$" , NA_character_)) %>% 
                                                  mutate(pat_exp = str_replace(pat_exp, "^None$" , NA_character_)) %>% 
                                                  mutate(pat_no = paste(pat_no[!is.na(pat_no)], " (exp. ",  pat_exp[!is.na(pat_exp)], ")", sep =  "", collapse = "; ")) %>% 
                                                  ungroup() %>% select(-pat_exp) %>% distinct() %>% 
                                                  mutate(pat_no = str_replace(pat_no, "[:blank:]\\(exp\\.[:blank:]\\)" , "")) %>% 
                                                  mutate(cmpd_name = str_replace(cmpd_name, "^NA$" , "")) %>% 
                                                  mutate(cmpd_syn = str_replace(cmpd_syn, "^NA$" , ""))
                                              }
                                    , by = "chembl_cmpd_id"
                                    ) %>% 
                                   select(
                                     MS3_id,
                                     acc_no,
                                     chembl_target_id,
                                     ro5_viol,
                                     ro3_pass,
                                     target_name,
                                     cmpd_name,
                                     cmpd_syn,
                                     chembl_cmpd_id,
                                     smiles,
                                     phase,
                                     nat_prod,
                                     first_approv,
                                     pat_no
                                   )
                                
                                 tab <- tab %>% filter(ro5_viol <= input$Ch_ro5_viol | is.na(ro5_viol)) %>% 
                                 filter(phase %in% range(input$chPhase) | is.na(phase))
                                 
                                 if(input$Ch_ro3_pass){
                                   tab <- tab %>% filter(ro3_pass == "Y")
                                 }
                                 
                                 if(input$InNatProd){
                                   tab <- tab %>% filter(nat_prod == "yes")
                                 }
                                 
                                 
                                 if(input$InScreenedSchisto){
                                   tab <- tab %>% filter(chembl_cmpd_id %in% screenS$chembl_cmpd_id)
                                 }
                                 
                                 tab <- tab %>% left_join(screenS, by = "chembl_cmpd_id")
                                 
                                 tab
                               }
  )


DrugDBTable <- eventReactive(eventExpr = input$show_drug_db,
                             {
                               
                              tab <- createTable() %>%
                                slice(1:isolate(input$num_dl_db)) %>% 
                                select(MS3_id) %>%
                                left_join(., drugb)
                              
                              tab <- tab %>% filter(status %in% input$DBstatus | is.na(status))
                              
                               if(input$InMDDR){
                                 tab <- tab %>% filter(mddr.like == 1)
                               }
                               
                               if(input$InRo5){
                                 tab <- tab %>% filter(ro5_pass == 1)
                               }

                              tab <- tab %>% 
                                mutate_at(.,
                                          vars(ro5_pass, mddr.like),
                                          funs(
                                            if_else(. == 1, "yes", "no")
                                          )
                                )
                              
                              tab
                              }
)

# Create link for drug IDs to add to data table
createLinkCH <- function(val) {
  sprintf('<a href="https://www.ebi.ac.uk/chembl/compound/inspect/%s" target="_blank" class="btn btn-primary">%s</a>',val, val)
}

createLinkDB <- function(val) {
  sprintf('<a href="https://www.drugbank.ca/drugs/%s" target="_blank" class="btn btn-primary">%s</a>',val, val)
}

# Create datatables for drug output
output$drugSelTableCH <- renderDataTable(
                                            {
                                              
                                              tableCH <- DrugCHTable()
                                              tableCH$chembl_cmpd_id <- createLinkCH(tableCH$chembl_cmpd_id)
                                              tableCH$chembl_cmpd_id[tableCH$chembl_cmpd_id == '<a href=\"https://www.ebi.ac.uk/chembl/compound/inspect/NA\" target=\"_blank\" class=\"btn btn-primary\">NA</a>'] <- ""
                                              columns <- tibble(name = colnames(tableCH))
                                              
                                              lookup <- as.tibble(matrix(
                                                
                                                c(
                                                  "MS3_id", "Identifier",
                                                  "acc_no", "Target accession number",
                                                  "chembl_target_id", "Target ID",
                                                  #"target_type", "Target type",
                                                  "ro5_viol", "Ro5 violations",
                                                  "ro3_pass", "Ro3 pass",
                                                  "target_name", "Target name",
                                                  "cmpd_name", "Compound name(s)",
                                                  "cmpd_syn", "Compound synonym(s)",
                                                  "chembl_cmpd_id", "Compound ID",
                                                  "smiles", "SMILES",
                                                  "phase", "Drug phase",
                                                  "nat_prod", "Natural product?",
                                                  "first_approv", "Approval date",
                                                  "pat_no", "Patent number(s) (expiry date)",
                                                  "link", "Screening data for Schistosoma"
                                                ),
                                                ncol = 2, byrow = T
                                              )
                                              
                                              )
                                              
                                              colnames(lookup) <- c("name", "cleanname")
                                              
                                              cleannames <- columns %>% left_join(., lookup, by = "name") %>% select(cleanname) %>% pull()
                                              
                                              
                                              datatable(
                                                data = tableCH
                                                , colnames = cleannames
                                                , options = list(scrollX = T, lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')))
                                                , escape = F
                                              )
                                            }
                                        )

output$drugSelTableDB <- renderDataTable(
                                         {
                                           
                                           tableDB <- DrugDBTable()
                                           columns <- tibble(name = colnames(tableDB))
                                           
                                           tableDB$cmpd_id <- createLinkDB(tableDB$cmpd_id)
                                           tableDB$cmpd_id[tableDB$cmpd_id == '<a href=\"https://www.drugbank.ca/drugs/NA\" target=\"_blank\" class=\"btn btn-primary\">NA</a>'] <- ""
                                           
                                           lookup <- as.tibble(matrix(
                                             
                                             c(
                                               "MS3_id", "Identifier",
                                               "target_id", "Target ID",
                                               "cmpd_id", "Compound ID",
                                               "cmpd_name", "Compound name",
                                               "status", "Approval status",
                                               "smiles", "SMILES",
                                               "ro5_pass", "Ro5 pass",
                                               "mddr.like", "MDDR-like"
                                             ),
                                             ncol = 2, byrow = T
                                           )
                                           
                                           )
                                           
                                           colnames(lookup) <- c("name", "cleanname")
                                           
                                           cleannames <- columns %>% left_join(., lookup, by = "name") %>% select(cleanname) %>% pull()
                                           
                                           
                                           dt <- datatable(
                                             data = tableDB
                                             , colnames = cleannames
                                             , options = list(scrollX = T, lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')))
                                             , escape = F
                                           )
                                           
}
                                         
                                         )

  
 # Create datatable of selected proteins
  output$protSelTable <- renderDataTable(
    {
      columns <- tibble(name = colnames(createTable()))
      
      lookup <- as.tibble(matrix(
        
        c("MS3_id" , "Identifier",
          "score" , "Score",
          "male_tpm" , "Adult transcription",
          "chokepoint"	, "Choke-point",
          "mus_leth"	, "Lethal phenotype mouse",
          "dmel_leth"	, "Lethal phenotype D. melanogaster",
          "cel_leth" , "Lethal phenotype C. elegans",
          "hsa_lower_75_qtl"	, "Similarity to human ortholog <75th percentile",
          "sp_cov"	, "SwissProt Ortholog",
          "sman_cov"	, "S. mansoni Ortholog",
          "sjap_cov"	, "S. japonicum Ortholog",
          "opvi_cov"	, "O. viverrini Ortholog",
          "closi_cov"	, "C. sinensis Ortholog",
          "fahe_cov"	, "F. hepatica Ortholog",
          "cel_cov"	, "C. elegans Ortholog",
          "dmel_cov"	, "D. melanogaster Ortholog",
          "mus_cov"	, "Mouse Ortholog",
          "ChCmpdOne"	, ">= 1 compound in ChEMBL",
          "ChCmpdMore"	, ">= 5 compound in ChEMBL",
          "DBCmpdOne"	, ">= 1 compound in DrugBank",
          "DBCmpdMore"	, ">= 5 compound in DrugBank",
          "uniqIpro" , "Unique InterPro identifier",
          "anno", "Annotation",
          "gene_name_species", "SwissProt Ortholog"),
        ncol = 2, byrow = T
      )
      
      )
      
      colnames(lookup) <- c("name", "cleanname")
      
      cleannames <- columns %>% left_join(., lookup, by = "name") %>% select(cleanname) %>% pull()
     
      
      datatable(
        data = createTable()
        , colnames = cleannames
        , options = list(scrollX = T, lengthMenu = list(c(10, 25, 50, 100, -1), list('10', '25', '50', '100', 'All')))
      ) 
      
    }
  )
  
  # Create UI (buttons) for protein download
  output$downloads <- renderUI(
    {
      out <- gatherProtSelect()
      top10 <- ceiling((count(out$table) %>% pull())*0.1)
      list(
        numericInput(inputId = "num_dl", label = "Download data for top n proteins:", value = top10),
        downloadButton(outputId = "dl_prot", label = "Download protein data"), 
        downloadButton(outputId = "dl_rank", label = "Download ranking parameters")
      )
    }
  )
    
  # Create UI (buttons) for drug downloads
  output$CHout <- renderUI(
    {
      out <- gatherProtSelect()
      top10 <- ceiling((count(out$table) %>% pull())*0.1)
      list(
        numericInput(inputId = "num_dl_chembl", label = "Show ChEMBL data for top n proteins:", value = top10),
        actionButton(inputId = "show_drug_chembl", label =  "Show ChEMBL data")
      )
    }
  )
  
  output$CHdl <- renderUI(
    {
      DrugCHTable()
      downloadButton(outputId = "dl_drug_chembl", label =  "Download ChEMBL data") 
    }
  )
  
  output$DBout <- renderUI(
    {
      out <- gatherProtSelect()
      top10 <- ceiling((count(out$table) %>% pull())*0.1)
      list(
        numericInput(inputId = "num_dl_db", label = "Show DrugBank data for top n proteins:", value = top10),
        actionButton(inputId = "show_drug_db", label =  "Show DrugBank data")
      )
    }
  )
  
  output$DBdl <- renderUI(
    {
      DrugDBTable()
      downloadButton(outputId = "dl_drug_db", label =  "Download DrugBank data")    
    }
  )

# Download ranking parameters  
output$dl_rank <- downloadHandler(
  filename = function() {
    paste('ranking-parameters-', Sys.Date(), '.tsv', sep='')
    },
  content = function(file) {
    write_tsv(makeParamsTable(), file)
  }
)

# Download protein table
output$dl_prot <- downloadHandler(
  filename = function() {
    paste('protein-data-', Sys.Date(), '.tsv', sep='')
  },
  content = function(file) {
    write_tsv((createTable() %>% slice(1:input$num_dl)), file)
  }
)
  
# Download drug tables
output$dl_drug_chembl <- downloadHandler(
  filename = function() {
    paste('drug-data-chembl-top-', input$num_dl_chembl ,"-", Sys.Date(), '.tsv', sep='')
  },
  content = function(file) {
    write_tsv(
      (
        DrugCHTable()
      ), file
    )
  }
)

output$dl_drug_db <- downloadHandler(
  filename = function() {
    paste('drug-data-drugbank-top-', input$num_dl_db ,"-", Sys.Date(), '.tsv', sep='')
  },
  content = function(file) {
    write_tsv(
      (
        DrugDBTable()
      ), file
    )
  }
)

 # Reset input values
  observeEvent(input$reset_input, {
    reset(id = "navbar")
  })
  
  #Grey out sliders that are not used for weighting
  observeEvent(input$trAdRB, {toggleState("trAdS", input$trAdRB == "weighted")})
  observeEvent(input$letCelRB, {toggleState("letCelS", input$letCelRB == "weighted")})
  observeEvent(input$letDmelRB, {toggleState("letDmelS", input$letDmelRB == "weighted")})
  observeEvent(input$letMusRB, {toggleState("letMusS", input$letMusRB == "weighted")})
  observeEvent(input$chokeRB, {toggleState("chokeS", input$chokeRB == "weighted")})
  observeEvent(input$uniqIproRB, {toggleState("uniqIproS", input$uniqIproRB == "weighted")})
  observeEvent(input$HomCelRB, {toggleState("HomCelS", input$HomCelRB == "weighted")})
  observeEvent(input$HomDmelRB, {toggleState("HomDmelS", input$HomDmelRB == "weighted")})
  observeEvent(input$HomMusRB, {toggleState("HomMusS", input$HomMusRB == "weighted")})
  observeEvent(input$HomSjRB, {toggleState("HomSjS", input$HomSjRB == "weighted")})
  observeEvent(input$HomSmRB, {toggleState("HomSmS", input$HomSmRB == "weighted")})
  observeEvent(input$HomCsinRB, {toggleState("HomCsinS", input$HomCsinRB == "weighted")})
  observeEvent(input$HomOvivRB, {toggleState("HomOvivS", input$HomOvivRB == "weighted")})
  observeEvent(input$HomFhepRB, {toggleState("HomFhepS", input$HomFhepRB == "weighted")})
  observeEvent(input$HomSPRB, {toggleState("HomSPS", input$HomSPRB == "weighted")})
  observeEvent(input$QuantHsRB, {toggleState("quant_slider", input$QuantHsRB == "weighted")})
  observeEvent(input$ChCmpdOneRB, {toggleState("ChCmpdOneS", input$ChCmpdOneRB == "weighted")})
  observeEvent(input$ChCmpdMoreRB, {toggleState("ChCmpdMoreS", input$ChCmpdMoreRB == "weighted")})
  observeEvent(input$DBCmpdOneRB, {toggleState("DBCmpdOneS", input$DBCmpdOneRB == "weighted")})
  observeEvent(input$DBCmpdMoreRB, {toggleState("DBCmpdMoreS", input$DBCmpdMoreRB == "weighted")})
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)