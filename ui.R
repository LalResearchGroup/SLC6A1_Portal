################################################
################# SLC6A1 PORTAL ##################
################################################
################## LIBRARIES ################

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(tippy)
library(r3dmol)
library(DT)
library(readr)
library(tidyverse)
library(vembedr)
library(RColorBrewer)
library(shinyhelper)
library(plyr)
library(rsconnect)

# CSS/STYLE/INFO #
landing_panel <- "color: #333333;
      height: 200px;
      width: 260px"

spinner_color <- "#2196f3"

sub_style <- "color:gray; font-style: italic; font-size: 14px;"

#Layout colors 

schema_color <- "#85C1E9"

schema_color_strong <- "#2E86C1"

schema_color_light <- "gainsboro"

##### General Variable #####

contact_us <- "mailto:lald@ccf.org; stefana4@ccf.org"
contact_katrine <- "mailto:kamaa@filadelfia.dk"


##### Landing Page Variable #####
#General
landing_portal_title =  "SLC6A1 Portal"

landing_bannername = "SLC6A1-Portal_Banner_inks_v3.png"

#Tabs 

landing_tab1 = HTML(paste0(em("SLC6A1"), ", its function and associated disorders"))
landing_tab2 = "Foundations, family groups, links to resources, and more"
landing_tab3 = "Comprehensive information on variant interpretation"
landing_tab4 = "Filter and select a subset of variants for research"
landing_tab5 = "Information on where to register with a SLC6A1-related disorder"

#####Basic Information Variable #####
#General 
basic_title = HTML(paste0(em("SLC6A1"),"-related disorders"))

##Genes 
#Gene 1
gene1 = "SLC6A1"
gene1_info_header = div(HTML(paste0("Information about ",em("SLC6A1"))),style = "font-size:17px;")

basic_text_title_1 <- HTML(paste0(em("SLC6A1"),"-related disorders"))

basic_mechanism_png <- "SLC6A1_mechanism_white-bg_v2.png"

#Visualisations
basic_visual_title = "Clinical Information of our patients"

basic_phenotype_title1 = "Phenotype overview"
basic_phenotype_title2 = "Seizure onset"
basic_phenotype_title3 = "Epilepsy syndromes" 
basic_phenotype_title4 = "Development" 

basic_plot_title_fac1 = "Epilepsy"
basic_plot_title_fac2 = "Autism"
basic_plot_title_fac3 = "Epilepsy syndromes"

basic_plot_title_num1 = ""
basic_plot_title_num2 = ""

basic_abbreviations1 = ""

basic_abbreviations2 = "PTV: Protein truncating variants"

basic_abbreviations3 = "EMAS = Epilepsy with myoclonic-atonic seizures; GGE = Genetic generalized epilepsy;
NAFE = Non-acquired focal epilepsy; EOAE = Early onset absence epilepsy; DEE = Developmental epileptic encephalopathy; CAE = Childhood absence epilepsy; TLE = Temporal lobe epilepsy; LGS = Lennox-Gastaut syndrome"

basic_abbreviations4 = "ID = Intellectual disability; DD = Developmental delay"

#tippy
epi_syndrome_tippy <- h5(HTML(paste0("<ul><li>CAE = Childhood Absence Epilepsy</li>")), 
                         HTML(paste0("<li>DEE = Developmental Epileptic Encephalopathy</li>")), 
                         HTML(paste0("<li>EMAS = Epilepsy with Myoclonic-Atonic Seizures</li>")),
                         HTML(paste0("<li>GGE = Genetic Generalized Epilepsy</li>")), 
                         HTML(paste0("<li>LGS = Lennox-Gastaut Syndrome</li>")), 
                         HTML(paste0("<li>NAFE = Non-Acquired Focal Epilepsy</li>")), 
                         # HTML(paste0("<li>EOAE = Early onset absence epilepsy</li>")), 
                         HTML(paste0("<li>TLE = Temporal Lobe Epilepsy</li>")), 
                         # HTML(paste0("<li>No SZ = No seizures</li>")), 
                         # HTML(paste0("<li>UG = Unclassified general</li>")), 
                         align = "left")

dd_id_tippy <- h5(HTML(paste0("<ul><li>DD = Developmental delay</li>")),
                  HTML(paste0("<li>ID = Intellectual disability</li>")), align = "left")

Gene1_basic_text <- p(HTML(paste0("",em("SLC6A1"),"-related disorders are a group of genetic disorders characterized by disease-causing variants in the ",em("SLC6A1")," gene. Clinical manifestations of ",em("SLC6A1")," disorders include mild to severe intellectual disability, behavioral disturbances such as autism and ADHD and seizures (mean onset 3.7 years), often absences, myoclonic and atonic seizures. One or more of these features can be present in each individual affected by an ",em("SLC6A1"),"-related disorder.")))



#####Families Variable #####  
family_title <- HTML(paste0("What are ",em("SLC6A1"),"-related disorders"))

family_video_url <- "https://www.youtube.com/watch?v=iNmf1drJzo4"

family_video_sub_description <- "Subtitles available in English"

family_video_description_title <- HTML(paste0("What are ",em("SLC6A1"),"-related disorders"))

family_video_description_content <- p("If your family has received a diagnosis of ",em("SLC6A1"),"-related disorders, you’ve probably never heard of it before. That’s okay. Your doctors probably haven’t either! But the Educational Science community is here to help.", 
                                      br(),
                                      br(),
                                      "This 4-minute “What are ",em("SLC6A1"),"-related disorders?” explainer video helps families learn about these disorders in plain language. You’ll learn:",
                                      br(),
                                      HTML(paste0("<ul><li> Scientific terms (such as genes, proteins, and variants/mutations) that will help you understand ",em("SLC6A1"),"-related disorders 
                                      </li><li> What causes ",em("SLC6A1"),"-related disorders 
                                      </li><li> How ",em("SLC6A1"),"-related disorders are diagnosed 
                                      </li><li> The symptoms of ",em("SLC6A1"),"-related disorders 
                                      </li><li> Where to find support")))



##### Variant Analysis variable #####

master.df <- read_delim("data/master_table.txt", delim = "\t") 

var_possible_genes_title <- "Select gene"

var_possible_genes <- "SLC6A1"

var_possible_phenotype <- c("Epilepsy + ID", "Epilepsy without ID",
                            "Autism", "Other")

var_patient_info_title <- p("Individuals with the same variant in the SLC6A1 Portal",
                            tippy(icon("question-circle"), tooltip = h5(HTML(paste0("Variant collection:")),
                                                                        HTML(paste0("<ul><li>Goodspeed et al.,2020</li>")), 
                                                                        HTML(paste0("<li>Internal variant database (unpublished)</li>")), align = "left"),
                                                           animation = "scale", theme = "light"))
var_patient_info_abb <- ""


var_paralog_info_abb <- ""


variants_compareplots_abb <- "DEE: Developmental and Epileptic Encephalopathies; BNFS: Benign Familial Neonatal Seizures; EOEE: Early-Onset Epileptic Encephalopathies; ASD: Autism Spectrum Disorder "

#####Research Variable ####
research_pheno_abb <- "DD: Developmental delay, ID: Intellectual disability"
research_pheno_abb2 <- "EMAS = Epilepsy with myoclonic-atonic seizures; GGE = Genetic generalized epilepsy; UG = Unclassified general;
NAFE = Non-acquired focal epilepsy; EOAE = Early onset absence epilepsy; DEE = Developmental epileptic encephalopathy; CAE = Childhood absence epilepsy; TLE = Temporal lobe epilepsy; LGS = Lennox-Gastaut syndrome; No SZ = No seizures; NA = Not available"

research_geno_transcripts <- p("The following transcript was used:",em("SLC6A1"),": ENST00000287766",
                               br(),
                               "H: Helix")

#Datasets required for research tab 
Patient_data.df <- read_delim("data/Patient_variants_SLC6A1_v8.txt", delim = "\t") %>%
  dplyr::rename(Sz_onset = "Age at seizure onset (months)",
         Epilepsy = "Epilepsy") %>% 
  dplyr::rename(Autism = "Autistic traits",
         Epilepsy_syndrome = "Epilepsy Syndrome Classification",
         ID_after_sz_onset = "Cognitive Level AFTER Seizure Onset")

Patient_data_missense.df <- read_delim("data/Patient_variants_SLC6A1_v8.txt", delim = "\t") %>% 
  dplyr::rename(Sz_onset = "Age at seizure onset (months)") %>% 
  filter(Vartype == "Missense ")

##### About Variable #####
about_terms_of_use <- p("All data presented here are publicly available for the benefit of the wider biomedical community. 
                               You can freely explore the data, and we encourage the use and publication of results generated 
                               from these data. However, we encourage you to contact us before embarking on analyses to 
                               check if your proposed analysis overlaps with work currently underway by our team. Further, 
                               we request that you actively acknowledge and give attribution to the SLC6A1 Portal project, and 
                               link back to the relevant page, wherever possible. All users of our data agree to not attempt 
                               to reidentify participants. Our data set has been subjected to extensive quality control, 
                               but may be imperfect so errors may remain. If you spot any results that seem impossible, 
                               or have suggestions for SLC6A1 Portal improvements: ")

about_data_generation <- "Data has been curated in a community effort."

## Tutorials variables
tut_basic_info <- "https://youtu.be/q9Usxkj_wBs"

tut_family <- "https://youtu.be/FdVgO6bWaWc"

tut_variant_analysis <- "https://youtu.be/gCPjOzyH5u8"

tut_research <- "https://youtu.be/i2w8Z2qH9kQ"



################# UI #################
shinyUI(
    ##### LANDING PAGE #####
    ui <- 
      navbarPage(#fluid = FALSE, 
        windowTitle = "SLC6A1 Portal",  
        id = "TabDisplay",
        theme = "mytheme2.css",
        title = p(icon("dna"), landing_portal_title), 
        header = tagList(useShinydashboard()),
        
        tabPanel(
          tags$style(HTML(paste0("
                      .navbar-default .navbar-brand {
                      background-color: ",schema_color_light,";
                      color: black
                    }"))),
          tags$style(HTML(paste0("
                      .navbar-default {
                      background-color: ",schema_color_light,";
                    }"))),
          
          title = "Welcome",
          value = "welcomeTab",
          tabName = "welcomeTab",
          
          div(style = "font-size:100%", 
              fluidRow(
                column(12, 
                   style = "background-color: white; color: #676767",
                   align = "center",
                   br(), br(),
                   img(src = landing_bannername, width = "100%")   #insert your banner, saved in the www-folder 
            )
          ),
          
          div(style = 'background-color: #F8FCFE',
              fluidRow(
                
                # Basic Information Tab
                div(width = "100%", 
                    column(2, offset = 2, align = "center",
                           br(),
                           panel(width = 12, 
                                 status = "success",
                                 heading = "",
                                 h2(tags$i(class = "fa fa-dna", style = "color:#676767")),
                                 br(),
                                 div(landing_tab1,  #c
                                     br(), br(),
                                     actionBttn(
                                       inputId = "infoBtn",
                                       label = "Basic Information",
                                       color = "success",
                                       block = TRUE,
                                       size = "md",
                                       style = "stretch")),
                                 style =  "background-color: #f3faf4;",   #set color of the tab
                           ), 
                    ),
                    
                    #Family tab / Educational resources tab
                    column(2,align = "center",
                           br(),
                           div(panel(
                             status = "warning",
                             heading = "",
                             h2(tags$i(class = "fa fa-child", style = "color:#676767")),
                             br(),
                             div(landing_tab2,
                                 br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "familyBtn",
                                   label = "Educational resources",
                                   color = "warning",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #fff8ef;",   #set color of the tab
                           ))
                    ),
                    
                    # Variant Analysis Tab
                    column(2, align = "center",
                           br(),
                           div(panel(
                             status = "info",
                             heading = "",
                             h2(tags$i(class = "fa fa-code-branch", style = "color:#676767")),
                             br(),
                             div(landing_tab3,br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "variantBtn",
                                   label = "Variant Analysis",
                                   color = "royal",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #f9f1fa;",    #set color of the tab
                           ))
                    ),
                    
                    # Research Tab
                    column(2, align = "center",
                           br(),
                           div(panel(
                             status = "danger",
                             heading = "",
                             h2(tags$i(class = "fa fa-microscope", style = "color:#676767")),
                             br(),
                             div(landing_tab4,br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "researchBtn",
                                   label = "Research",
                                   color = "danger",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #fdf0f1;",
                             
                           ))
                    ),
                )),
              fluidRow(
                column(4, offset = 2, align = "center",
                       p("Visit our other Portals:", style = "font-size:14px;", align = "center", style = "font-size:14px;"),
                       p(shiny::a("SCN-Portal", href="http://scn-portal.broadinstitute.org/", target = '_blank'), "& ",
                         shiny::a("GRIN-Portal", href="http://grin-portal.broadinstitute.org/", target = '_blank'),
                         style = "font-size:12px;")
                ),
                column(4, align = "center",
                       
                       p("You want to join the project or provide feedback?", align = "center",style = "font-size:14px;"),
                       p("Please ", shiny::a("contact us!", href=contact_us), style = "font-size:12px;", align = "center"),
                       br()
                ),
                column(12, align = "center",offset = 4,
                      fluidRow(
                        div(width = "100%",valueBoxOutput("disclaimer")))
               )
            ))
          )), # end tab Panel 
        ##### BASIC INFORMATION #####
        tabPanel(title = "Basic Information", value = "infoTab",

                 fluidRow(
                   column(10, offset = 1,
                     panel(heading = basic_title,    
                           status = "success",
                       fluidRow(
                         column(4,
                           box(width = 12,
                             title = p("History of ",em(gene1), " research"),   
                             timelineBlock(
                                 reversed = FALSE,
                                 width = 12,

                                 #Add as many blocks as desired for new publications
                                 timelineLabel(2015, color = "teal"),
                                 timelineItem(
                                   title = div(strong("Discovery of",em("SLC6A1"),"as cause for epilepsy with myoclonic-atonic seizures")),
                                   icon = icon("dna"),
                                   color = "olive",
                                   #ADD link to the publication
                                   time = shiny::a("Carvill et al.", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4570550/", target = '_blank'),
                                   border = FALSE,
                                 ),
                                 timelineLabel(2018, color = "teal"),
                                 timelineItem(
                                     title = div("Delineation of the phenotypic spectrum"),
                                     icon = icon("user"),
                                     color = "aqua",
                                     #ADD link to the publication
                                     time = shiny::a("Johannesen et al.", href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5912688/", target = '_blank'),
                                     border = FALSE,
                                 ),
                                 timelineLabel(2020, color = "teal"),
                                 timelineItem(
                                   title = HTML(paste0("Current knowledge of ",em("SLC6A1"),"-related neurodevelopmental disorders")),
                                   icon = icon("user"),
                                   color = "aqua",
                                   #ADD link to the publication
                                   time = shiny::a("Goodspeed et al.", href="https://academic.oup.com/braincomms/article/2/2/fcaa170/5922604?login=true", target = '_blank'),
                                   border = FALSE,
                                 ),
                                 timelineLabel(2021, color = "teal"),
                                 timelineItem(
                                   title = div("Common molecular mechanisms of",em("SLC6A1")),
                                   icon = icon("user"),
                                   color = "aqua",
                                   #ADD link to the publication
                                   time = shiny::a("Mermer et al.", href="https://pubmed.ncbi.nlm.nih.gov/34028503/", target = '_blank'),
                                   border = FALSE,
                                 ),
                             )
                           )),
                           column(4, 
                            box(title=basic_text_title_1, width = 12, 
                                  Gene1_basic_text)),
                           column(4,
                             box(title =  HTML(paste0(em("SLC6A1") ," mechanism")), width = 12,
                                 img(src = basic_mechanism_png, width = "100%"))),

                           column(8,
                             box(title = p("Summary of clincal information in our registry",
                                           tippy(icon("question-circle"),
                                 tooltip = h5(HTML(paste0("Variant collection:")),
                                              HTML(paste0("<ul><li>Goodspeed et al.,2020</li>")), 
                                              HTML(paste0("<li>Internal variant database (unpublished)</li>")), align = "left"),
                                 animation = "scale", 
                                 theme = "light")),
                                 width = 12,

                              tabsetPanel(
                              #Phenotype overview
                                tabPanel(title = basic_phenotype_title1,
                                         br(),
                                         div(p("N =" ,Patient_data.df %>% filter(Gene == gene1,!is.na(Epilepsy)) %>% nrow(), " patients"),style = "font-size:15px;color:black;"),
                                         br(),
                                 fluidRow(align = "center",
                                    div(title = paste0("Phenotypes available for " ,Patient_data.df %>% filter(Gene == gene1) %>% nrow(), " patients"), width = 12),
                                     column(6, 
                                            box(title=div(basic_plot_title_fac1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac1_gene1"))),
                                     column(6, 
                                            box(title=div(basic_plot_title_fac2, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac2_gene1"))),
                                     column(12, 
                                            align="center", 
                                            p(basic_abbreviations1, style=sub_style)))),
                                #Seizure onset 
                                 tabPanel(title = basic_phenotype_title2,  
                                          br(),
                                          div(p("N =" ,Patient_data.df %>% filter(Gene == gene1,!is.na(Sz_onset)) %>% nrow(), " patients"),style = "font-size:15px;color:black;"),
                                          br(),
                                  fluidRow(
                                           align = "center",
                                    column(12, 
                                           align = "center",
                                           plotOutput(outputId = "basic_legend_plot1", height = 20, width = 300)),
                                    column(6, 
                                           box(title=div(basic_plot_title_num1, style = "font-size: 10px"), width=12, plotlyOutput(outputId = "Phenotype_num1_gene1"))),
    
                                    column(12, align="center", p(basic_abbreviations2, style=sub_style)))),
                                #Phenotype tab 3
                                 tabPanel(title = basic_phenotype_title3,   
                                          br(),
                                          div(p("N =" ,Patient_data.df %>% filter(Gene == gene1, !is.na(Epilepsy_syndrome)) %>% nrow(), " patients"),style = "font-size:15px;color:black;"),
                                          br(),
                                  fluidRow(align = "center",
                                          div(title = paste0("Phenotypes available for " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Epilepsy_syndrome)) %>% nrow(), " patients"), 
                                              width = 12),
                                          column(12, 
                                                 box(title=div(style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac3_Missense", height = 500))),
                                          column(12, align="center", 
                                                 p(basic_abbreviations3, style=sub_style)),
                                          column(12, 
                                                 box(title=div(style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac3_PTV", height = 500))),
                                          column(12, align="center", 
                                                 p(basic_abbreviations3, style=sub_style)))),
                                #Phenotype tab 4
                                tabPanel(title = basic_phenotype_title4,    #c
                                         br(),
                                         div(p("N =" ,Patient_data.df %>% filter(Gene == gene1,!is.na(ID_after_sz_onset)) %>% nrow(), " patients"),style = "font-size:15px;color:black;"),
                                         br(),
                                 fluidRow(align = "center",
                                          div(title = paste0("Phenotypes available for " ,Patient_data.df %>% filter(Gene == gene1,!is.na(ID_after_sz_onset)) %>% nrow(), " patients"), 
                                          width = 12),
                                  column(12, 
                                         box(title=div(style = "font-size: 15px"), 
                                             width=12, plotlyOutput(outputId = "Phenotype_fac4_gene1"))),
                                  column(12, align="center", 
                                         p(basic_abbreviations4, style=sub_style))))
                          ))))
                     )))),#end basic info tab panel

        ##### FAMILIES ##### Educational resources

        tabPanel(title = "Educational resources", value = "familyTab",
                 tags$style(HTML(paste0("
                      .box.box-solid.box-danger>.box-header {
                      background-color: ",schema_color_light,";
                      border: black
                      }"))),
                 tags$style(HTML(paste0("
                      .box.box-solid.box-danger {
                      border: 1px solid black
                      }"))),
         fluidRow(
           column(10, offset = 1,
            panel(heading = "Educational resources and patient advocacy groups", status = "warning",
              tabsetPanel(
                tabPanel(title = "Videos",
                  fluidRow(
                    column(5, offset = 1,
                      br(), 
                      br(),
                      box(title= HTML(paste0("What are ",em("SLC6A1"),"-related disorders?")),
                          width = 12,
                          embed_url(family_video_url) %>% use_rounded() %>% use_align("center") %>% use_bs_responsive(),
                          footer = p(div(align="center", family_video_sub_description, style=sub_style))
                    )),
                    column(5, offset=0, br(), align = "justify",  br(), br(), br(),
                     box(title= div(family_video_description_title, style = "color:black;"),
                     width = 12,
                     status = "danger", solidHeader = TRUE,
                     family_video_description_content
                    ))
                 )),

              #Logos linked to family organizations
               tabPanel(title = "Family and other organizations",
                panel(status="default",heading = "Family organizations", width = 12,
                 fluidRow(
                   column(12, align = "center", 
                          offset = 0,
                          div(width = 12, shiny::a(
                            img(src = "SLC6A1ConnectCircle.jpg", width = "20%"),
                            href = "https://slc6a1connect.org",
                            target = '_blank'))
                   )
                  )
                 ),
                 panel(status="default",heading = "Related organizations", width = 12,
                   fluidRow(
                     column(3, align = "center", offset = 1,
                            div(width = 12, shiny::a(
                               img(src = "aes.png", width = "70%"),
                               href = "https://www.aesnet.org/",
                               target = '_blank'))
                     ),
                     column(3, align = "center",
                            div(width = 12, shiny::a(
                               img(src = "nord.png", width = "70%"),
                               href = "https://rarediseases.org/",
                               target = '_blank'))
                     ),
                     column(3, align = "center",
                            div(width = 12, shiny::a(
                               img(src = "epi_foundation.png", width = "70%"),
                               href = "https://www.epilepsy.com/",
                               target = '_blank'))
                     ),
                     column(3, align = "center", offset = 1,
                            div(width = 12, shiny::a(
                              img(src = "combined.png", width = "70%"),
                              href = "https://combinedbrain.org/",
                              target = '_blank'))
                     ),
                     column(3, align = "center",
                            div(width = 12, shiny::a(
                               img(src = "global_genes.png", width = "70%"),
                               href = "https://globalgenes.org/",
                               target = '_blank'))
                     ),
                     column(3, align = "center",
                            div(width = 12, shiny::a(
                              img(src = "ren.png", width = "70%"),
                              href = "https://www.rareepilepsynetwork.org/",
                              target = '_blank'))
                     ),
                     column(3, align = "center",
                            div(width = 12, shiny::a(
                               img(src = "epi_council.png", width = "70%"),
                               href = "https://www.epilepsyleadershipcouncil.org/",
                               target = '_blank'))
                     )
                   ))
         ),
         tabPanel(title = "Tutorials", value = "tutorialTab", #### Tutorials
                  fluidRow(
                    column(10, offset = 1,
                           panel(heading = "Tutorials", status = "primary",
                                 div(align = "center",
                                     tabsetPanel(
                                       tabPanel("Basic Information",
                                                column(8, offset = 2,
                                                       box(title= "",
                                                           align = "right",
                                                           width = 12,
                                                           embed_url(tut_basic_info) %>%
                                                             use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                       ))
                                       ),
                                       tabPanel("Educational resources",
                                                column(8, offset = 2,
                                                       box(title= "",
                                                           align = "right",
                                                           width = 12,
                                                           embed_url(tut_family) %>%
                                                             use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                       ))
                                       ),
                                       tabPanel("Variant Analysis",
                                                column(8, offset = 2,
                                                       box(title= "",
                                                           align = "right",
                                                           width = 12,
                                                           embed_url(tut_variant_analysis) %>%
                                                             use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                       ))
                                       ),
                                       tabPanel("Research",
                                                column(8, offset = 2,
                                                       box(title= "",
                                                           align = "right",
                                                           width = 12,
                                                           embed_url(tut_research) %>%
                                                             use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                       )))
                                     )))))
         )
         ))))

        ), # end families tab / Educational resources tab

        # ##### VARIANT ANALYSIS #####

        tabPanel(title = "Variant Analysis", value = "variantTab",
                 
         fluidRow(
           column(10, offset = 1,
                  panel(heading = "Analyse your variants", status = "info",
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default", heading = "Enter variant",
                   column(3,
                     pickerInput(
                         inputId = 'var_gene',
                         label =  h5(strong(var_possible_genes_title)),
                         width = "100%",
                         choices = var_possible_genes,
                         options = list(`style` = "default")
                     )
                   ),
                   column(3,
                     numericInputIcon(
                         inputId = "search_cDNA_pos",
                         label = h5(strong("cDNA Position")),
                         min = 1,
                         max = 10000,
                         value = 1024,
                         width = "100%"
                     ),
                     pickerInput(
                         inputId = "search_Allele",
                         label = "Ref",
                         choices = c("G", "A", "C", "T"),
                         selected = "G"
                     ),
                     pickerInput(
                         inputId = "search_cDNA_alt",
                         label = "Alt",
                         choices = c("G", "A", "C", "T", "null"),
                         selected = "A"
                     ),
                     actionButton(inputId = "search_var_c", label = "Search")),
                   column(3,
                     numericInputIcon(
                         inputId = "search_AA_pos",
                         label = h5(strong("Amino Acid Position")),
                         min = 1,
                         max = 10000,
                         value = 342,
                         width = "100%"
                     ),
                     pickerInput(
                         inputId = "search_AA_ref",
                         label = "Ref",
                         choices = sort(unique(master.df$AA_ref)),
                         selected = "Val"
                     ),
                     pickerInput(
                         inputId = "search_AA_alt",
                         label = "Alt",
                         choices = sort(unique(master.df$AA_alt)),
                         selected = "Met"
                     ),
                     actionButton(inputId = "search_var_p", label = "Search")),
                   column(2,
                     pickerInput(
                       inputId = "get_var_type",
                       label = h5(strong("Variant Consequence")),
                       choices = "Missense",
                       selected = "Missense"))
           )))),
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default", 
                       heading = "Variant Information",
                   fluidRow(
                     div(width = "100%",valueBoxOutput("geneBox1")),
                     div(width = "100%",valueBoxOutput("geneBox2")),
                     div(width = "100%",valueBoxOutput("geneBox3"))
                   ),
                   
                   # panel(heading = "Clinical Significance (Preliminary Use)", ##add report
                   #       fluidRow(
                   #         column(2,
                   #                radioButtons(
                   #                  inputId = "var_origin",
                   #                  label = h5(strong("Origin")), 
                   #                  choices = c("de novo (confirmed)", "de novo (assumed)", "Inherited", "Unknown/Untested"),
                   #                  selected = "Unknown/Untested"
                   #                )),
                   #         
                   #         column(2, 
                   #                radioButtons(
                   #                  inputId = "var_phenotype",
                   #                  label = h5(strong("Phenotype")),
                   #                  choices = c("Neurodevelopmental Disorder (NDD) + Cortical Malformation", "NDD + Epilepsy or Cerebral Visual Impairment", 
                   #                              "Unspecified NDD", "Other"),
                   #                  selected = "Other"
                   #                )),
                   #         br(),
                   #         column(8,
                   #                fluidRow(
                   #                  uiOutput(outputId = "ACMGUI_crit"),
                   #                  box(title="Explanation", collapsible = TRUE, collapsed = TRUE, width = 12,
                   #                      
                   #                      div(style="background-color:#f7f7f7",
                   #                          h4(style=" margin:5px; padding:10px","Very strong evidence"),  
                   #                          div(style=" margin:5px; padding:10px",strong("PVS1"),textOutput("PVS1_text")),
                   #                          
                   #                          br(),
                   #                          
                   #                          h4(style=" margin:5px; padding:10px","Strong evidence"),  
                   #                          div(style="display: flex;",
                   #                              div(style="width:33%; margin:5px; padding:10px",strong("PS2/PM6"),textOutput("PS2PM6_text")),
                   #                              div(style="width:33%; margin:5px; padding:10px",strong("PS3"), textOutput("PS3_text")),
                   #                              div(style="width:33%; margin:5px; padding:10px",strong("PS4"), textOutput("PS4_text"))),
                   #                          
                   #                          br(),
                   #                          
                   #                          h4(style=" margin:5px; padding:10px","Moderate evidence"), 
                   #                          div(style="display: flex;",
                   #                              div(style="width:50%; margin:5px; padding:10px",strong("PM1"), textOutput("PM1_text")),
                   #                              div(style="width:50%; margin:5px; padding:10px",strong("PM2"), textOutput("PM2_text"))),
                   #                          
                   #                          br(),
                   #                          
                   #                          h4(style=" margin:5px; padding:10px","Supporting evidence"), 
                   #                          div(style="display: flex;",
                   #                              div(style="width:50%; margin:5px; padding:10px;",
                   #                                  strong("PP2"), textOutput("PP2_text")),
                   #                              
                   #                              div(style="width:50%; margin:5px; padding:10px",strong("PP3"), textOutput("PP3_text"))))
                   #                      
                   #                  )))),
                   #       fluidRow(column(4,
                   #                       downloadButton("report", label = "Generate Report"))), br(), 
                   #       fluidRow(box(width=12,status = "danger", strong("Preliminary/incomplete use."), "Criteria for GRIN Registry variants 
                   #                             were automatically 
                   #                             annotated according to standard ", shiny::a(href="https://www.acmg.net/docs/Standards_Guidelines_for_the_Interpretation_of_Sequence_Variants.pdf", 
                   #                                                                         "ACMG Variant Interpretation Guidelines") , "(GRIN1, GRIN2A, GRIN2D) or the ClinGen GRIN Expert Panel 
                   #                                                               Specifications to the ACMG/AMP Variant Interpretation Guidelines (GRIN2B, unpublished). They are yet incomplete and not to be 
                   #                                                               used in a medical context.", style = sub_style))
                   # )
                   )))
             ),
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default", 
                       heading = "Patient information",
                       box(title = var_patient_info_title, width = 12,
                       DT::dataTableOutput(outputId = "patientTable"),
                       br(), p(var_patient_info_abb, style=sub_style, align = "center")))#,
             ))),
           fluidRow(
             column(12,style='padding:30px;',
              fluidRow(
                panel(status="default", 
                      heading = "Custom variant analysis",
                  tabsetPanel(
                    tabPanel("Comparison of phenotypes",
                             
                             br(),
                             h4("Explore what phenotypes have been associated with patients in our cohort that have a similar variant."),
                             br(),
                             radioGroupButtons(inputId = "compareButtons",
                                               label = "Variants with the same:",
                                               choices = c("Amino Acid Position", "Protein region", "Variant Type"),
                                               justified = TRUE,
                                               status = "default",
                                               checkIcon = list(yes = icon("ok",lib = "glyphicon"))
                             ),
                    br(),
                    div(width = "100%", plotlyOutput("comparePlot"))
                  ),
                  tabPanel("Comparison of variant location and molecular activity",
                           column(6,
                                  addSpinner(color = spinner_color,
                                             r3dmolOutput(
                                               outputId = "Var_analysis_compare_var",
                                               width = "100%",
                                               height = "400px"
                                             ))),
                           column(6,
                                  addSpinner(color = spinner_color,
                                             r3dmolOutput(
                                               outputId = "Var_analyis_hotzone",
                                               width = "100%",
                                               height = "400px"
                                             ))),
                           column(12,
                                  div(width = "100%", plotlyOutput("compare_act")))
                           
                  ),
                  
                  tabPanel("In silico scores",
                           br(),
                           #column(5,plotOutput(outputId = "paraz_legend", height = 30, width = 450)),
                           #column(12,div(width = "100%", plotlyOutput("paraz_legend"))),
                           column(4,
                                  "Paralog conservation score",
                                  align = "center",
                                  div(width = "100%", plotlyOutput("Var_analyis_paraz"))
                           ),
                           column(4,
                                  "Missense constraint score",
                                  align = "center",
                                  div(width = "100%", plotlyOutput("Var_analyis_mtr"))
                           ),
                           column(4,
                                  "Pathogenic variant enrichment",
                                  align = "center",
                                  div(width = "100%", plotlyOutput("Var_analyis_per"))
                           ),
                           column(12,plotOutput(outputId = "paraz_legend", height = 30, width = 600))
                           
                  )
  ))
             )))
                  )))
         ), # end variant analysis tab

         ##### RESEARCH #####

        tabPanel(title = "Research", value = "researchTab",
         fluidRow(
           column(10, offset = 1,
                  panel(heading = "Analyse your variants", status = "danger",
           fluidRow(
            column(12,style='padding:15px;',
              panel(heading="Filter Registry", status="default",
                fluidRow(
                  column(12,
                    selectizeGroupUI(
                      id = "research-filters",
                      btn_label = "Reset filters",
                      params = list(
                         varcons = list(
                             inputId = "Vartype",
                             title = p(strong("Variant Type")),
                             placeholder="all",
                             choices = c("Missense", "Other","PTV","Synonymous")
                         ),
                         aachange = list(
                             inputId = "AA_alt",
                             title = p(strong("Amino Acid Change")),
                             placeholder="all",
                             choices = unique(Patient_data.df$AA_alt)
                         ),
                         domain = list(
                           inputId = "Domain",
                           title = p(strong("Protein region")),
                           placeholder="all",
                           choices = unique(Patient_data.df$Domain)
                         ),
                         epilepsy = list(
                             inputId = "Epilepsy_syndrome",
                             title = p(strong("Epilepsy Syndromes"),tippy(icon("question-circle"),
                                                                          tooltip = epi_syndrome_tippy,
                                                                          animation = "scale", theme = "light")),
                             placeholder="all",
                             choices = unique(Patient_data.df$Epilepsy_syndrome)
                         ),
                         autism = list(
                           inputId = "Autism",
                           title = p(strong("Autism")),
                           placeholder="all",
                           choices = unique(Patient_data.df$Autism)
                         ),
                         id = list(
                           inputId = "ID_after_sz_onset",
                           title = p(strong("DD/ID"),tippy(icon("question-circle"),
                                                          tooltip = dd_id_tippy,
                                                          animation = "scale", theme = "light")),
                           placeholder="all",
                           choices = unique(Patient_data.df$ID_after_sz_onset)
                         ),
                         Hotzone_2D= list(
                           inputId = "Hotzone_2D",
                           title = p(strong("PER"),
                                     tippy(icon("question-circle"),
                                           tooltip = h6(HTML(paste0("PER: Pathogenic enriched regions<br>",
                                                                    "Reference: Pérez-Palma et al., 2019, PMID:31871067" )),
                                                        align = "left"),
                                           animation = "scale",
                                           theme = "light")),
                           placeholder="all",
                           choices = c(" No", "Mild","Moderate","Severe/Profound","Unkno"),
                           multiple = TRUE
                         ),
                         PER3D= list(
                           inputId = "PER3D",
                           title = p(strong("PER 3D"),
                                     tippy(icon("question-circle"),
                                           tooltip = h6(HTML(paste0("PER3D: Pathogenic enriched regions in 3D<br>",
                                                                    "Reference: Iqbal et al., to be published" )),
                                                        align = "left"),
                                           animation = "scale",
                                           theme = "light")),
                           placeholder="all",
                           choices = c(" No", "Mild","Moderate","Severe/Profound","Unkno"),
                           multiple = TRUE
                         ))
                    )),
           )))),
           fluidRow(
             column(12,style='padding:15px;',
               panel(status = "default",heading = "Custom variant exploration",
                     div(textOutput(outputId = "filtered_n"),
                         br(),
                         "  "),
                 tabsetPanel(
                   tabPanel(
                     "Genotype Interface",
                     fluidRow(
                       column(
                         12,
                         br(),
                         p(
                           strong("Selected variants are displayed in 1D (lolliplot) and 3D (protein structure).")
                         ),
                         fluidRow(
                           column(7,
                             materialSwitch(
                               inputId = "gnomad_m",
                               label = "Reference population missense variants (gnomAD)",
                               status = "primary",
                               right = T,
                               inline = T
                             ),
                             materialSwitch(
                               inputId = "PER",
                               label = "Pathogenic enriched regions (PER)",
                               status = "primary",
                               right = T,
                               inline = T
                             )
                         )),
                         fluidRow(
                           column(5,plotOutput(outputId = "Genotype_legend_plot", height = 30, width = 450) ),
                         ),
                         fluidRow(
                           column(7,
                             addSpinner(plotlyOutput(outputId = "Genotype_overview_plot"), color =
                                          spinner_color),
                             br(),
                             p(research_geno_transcripts , align="center", style=sub_style),
                             br()
                          ),
                           column(5,
                             addSpinner(color = spinner_color,
                               r3dmolOutput(
                                 outputId = "threeDmolGene_all_new",
                                 width = "100%",
                                 height = "400px"
                               )),
                             div("UniProt:",
                                 align="center",
                                 style=sub_style,
                                 shiny::a("P30531", href="https://www.uniprot.org/uniprot/P30531", target="_blank"),
                                 br(),
                                 div("Predicted structure", align="center", style=sub_style,
                                 shiny::a("Alpha fold", href="https://www.alphafold.ebi.ac.uk/entry/P30531", target="_blank")))),
                         column(6,align = "justify",plotlyOutput("domain_enrichment_pat_pop"))
             )))),
             tabPanel("Phenotype Interface",br(),
               fluidRow(
                 column(12,align="justify", plotlyOutput("research_phenotype1")
               )),
               fluidRow(
                 column(12,align="justify", plotlyOutput("research_phenotype6")
                 ),
                 column(12, align="center", p(research_pheno_abb2, style=sub_style))),
               fluidRow(
                 column(5, offset = 1,plotlyOutput("research_phenotype2"), align="center",
                        ),
                 column(5, plotlyOutput("research_phenotype3"))
               ),
               fluidRow(
                 column(5, offset = 1,plotlyOutput("research_phenotype4"), align="center",
                 ),
                 column(5,plotlyOutput("research_phenotype5"))),
               fluidRow(
                 column(12,
                        p(research_pheno_abb, style=sub_style, align = "center")))
               ),
             tabPanel("Functional Interface",
                      br(),
                      fluidRow(
                        column(6,
                               materialSwitch(
                                 inputId = "pat_only",
                                 label = "Show only functional data for patient variants",
                                 status = "primary",
                                 right = T,
                                 inline = T
                               ),
                               align="justify",
                               plotlyOutput("research_functional1")),
                        column(6,
                               br(),
                               br(),
                               plotlyOutput("research_functional2")),
                        column(4,
                               addSpinner(color = spinner_color,
                                          r3dmolOutput(
                                            outputId = "research_var_map1",
                                            width = "100%",
                                            height = "400px"
                                          ))),
                        column(4,
                               addSpinner(color = spinner_color,
                                          r3dmolOutput(
                                            outputId = "research_var_map2",
                                            width = "100%",
                                            height = "400px"
                                          ))),
                        column(4,
                               addSpinner(color = spinner_color,
                                          r3dmolOutput(
                                            outputId = "research_var_map3",
                                            width = "100%",
                                            height = "400px"
                                          ))),
                        column(6,
                               addSpinner(color = spinner_color,
                                          r3dmolOutput(
                                            outputId = "research_hotzone",
                                            width = "100%",
                                            height = "400px"
                                          ))),
                        column(6,
                               addSpinner(color = spinner_color,
                                          r3dmolOutput(
                                            outputId = "research_dist_struc",
                                            width = "100%",
                                            height = "400px"
                                          ))),
                        
                        column(6,
                               plotlyOutput("research_functional3"),
                               plotlyOutput("research_functional5")),
                        column(6,
                               plotlyOutput("research_functional4"),
                               plotlyOutput("research_functional6")),
                        column(12,align = "center",
                               div(h3("Other measured parameter"))),
                        column(6,
                               plotlyOutput("research_functional7")),
                        column(6,
                               plotlyOutput("research_functional8")),
                        
                      )
               # fluidRow(
               #   column(4,align="justify",plotlyOutput("research_functional1")),
               #   column(4,align="justify",plotlyOutput("research_functional2")),
               #   column(4,align="justify",plotlyOutput("research_functional3")),
               #   column(4,align="justify",plotlyOutput("research_functional4")),
               #   column(4,align="justify",plotlyOutput("research_functional5"))
               # )
               )
           ))))
                  )))
      ), #end research tab
##### ABOUT #####
      tabPanel(title = "About", 
               value = "aboutTab",
       fluidRow(
         column(10, offset = 1,
                panel(heading = "About", status = "primary",
            tabsetPanel(
              tabPanel("General Information",
                panel(heading ="Portal version", status = "default",
                      p(strong("This is the alpha version of the  SLC6A1 Portal")),
                      p("The SLC6A1 Portal is a coalition of investigators seeking to aggregate and harmonize data generated to study
                            SLC6A1-related disorders, and to make summary data interactively accessible for the wider scientific community,
                            while providing educational resources for everyone. The goals of this project are: "),
                      br(),
                  fluidRow(
                    column(6,
                           p(HTML(paste0("<ul><li> Providing information on ", em("SLC6A1"),"-related disorders</li>"))),
                           p(HTML(paste0("<li> Supporting research on ", em("SLC6A1"),"-related disorders</li>"))),
                           p(HTML(paste0("<li> Facilitating recruitment of individuals to the global ", em("SLC6A1")," registries</li>")))),
                    column(6,
                           p(HTML("<ul><li> Providing support in variant interpretation and classification</li>")),
                           p(HTML(paste0("<li> Visualizing data from the global ", em("SLC6A1")," registries</li>"))),
                           p(HTML("<li> Linking researchers, clinicians and families</li>")))),
                  br(),
                  footer = div("The SLC6A1 Portal is an ongoing project of the scientific community in collaboration with SLC6A1Connect. Interested collaborators are invited to reach out to join the project.")),
                panel(heading = "Teams and People", status = "default",
                      p("The current version of the SLC6A1 Portal has been developed by an international team of researchers and clinicians:"),
                  fluidRow(
                    column(10, 
                      panel(heading = "Team Leaders", 
                        p(strong("Dennis Lal"), "(Cleveland, US): General concept, web development, bioinformatics, video production"),
                        p(strong("Rikke Møller"), "(Dianalund, DK): Clinical & genetic data"),
                        p(strong("Kimberly Goodspeed"), "(Dallas, US): Clinical & genetic data"),
                        p(strong("Katty (Jing-Qiong) Kang"), "(Vanderbilt, US): Molecular data"),
                    )),
                    column(12, p("")),
                    column(2, 
                     panel(style="height: 230px;",heading = "Clinical & Genetic Data",
                           div(style="height: 100%;",
                           p(strong("Katrine Johannesen ")))
                    )),
                    
                    column(2, 
                      panel(style="height: 230px;",heading = "Molecular Data",
                            div(style="height: 100%",p(strong("Felicia Mermer")))
                    )),
                    column(2, 
                      panel(style="height: 230px;",heading = "Web Development",
                            div(style="height: 100%",p(strong("Arthur Stefanski")),
                            p(strong("Tobias Brünger")),
                            p("Eduardo Perez-Palma"),
                            p("Marie Macnee"),
                            p("Chiara Klöckner"))
                    )),
                    column(2, 
                      panel(style="height: 230px;",heading = "Bioinformatics",
                            div(style="height: 100%",
                                p(strong("Tobias Brünger")),
                                p("Eduardo Perez-Palma"),
                                p("Marie Mcnee"),
                                p("Patrick May"),
                                p("Chiara Klöckner"),
                                p("Johannes Lemke"))
                    )),
                    column(2, 
                      panel(style="height: 230px;", heading = "Video",
                            div(style="height: 100%",p(strong("Arthur Stefanski")),
                            p("Amber Freed"))
                    ))
                )),
                panel(heading = "Imprint", 
                      status = "default", 
                      p("We object to any commercial use and disclosure of data."),
                      p(strong("Copyright and use:"), "The authors grants you the right of use to make a private copy for personal purposes.
                        However, you are not entitled to change and/or pass on the materials or publish them yourself.
                        Upon request, you will receive free information about the personal data stored about you.
                        To do so, please contact the administrator."),
                      p(strong("No liability:"), "The contents of this web project have been carefully checked and created to
                        the best of our knowledge. But for the information presented here is no claim to completeness,
                        timeliness, quality and accuracy. No responsibility can be accepted for any damage caused by reliance
                        on or use of the contents of this website."))
                     ),
                tabPanel("Terms and Data Information",
                  panel(heading = "Terms of Use", status = "default",
                        p(about_terms_of_use,
                        shiny::a(href=contact_us, "Contact us"),
                        "that we can improve.")),
                  panel(heading = "Data Generation", status = "default",
                        about_data_generation)
                  )
 ))))))) # end ui
