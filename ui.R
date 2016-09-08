# Copyright 2016 Nicholas J. Seewald and Peng Liao
# This file is part of MRT-SS Calculator.
# 
# MRT-SS Calculator is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MRT-SS Calculator is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with MRT-SS Calculator.  If not, see <http://www.gnu.org/licenses/>.
# 
# Contact the authors at nseewald@umich.edu or pengliao@umich.edu


library(shiny)
library(shinyBS)
source("server.R")

shinyUI(fluidPage(
  titlePanel(HTML("<strong>MRT-SS Calculator</strong>: A Sample Size Calculator for Micro-Randomized Trials"), 
             windowTitle = "MRT-SS Calculator"),    
  
  ####Sample Size Calculator Simple version####
  
  ### Print error validation in "red" ###
  tags$head(
    tags$style(HTML("
      .shiny-output-error-validation {
        color: red;
      }
    ")),
    includeScript("www/google-analytics.js")
  ),
  
  ### CSS HEADER ###
  ### Apply style attributes across the application
  tags$head( 
    tags$style(type = 'text/css',"input[type='number'] {width:90px}"),   #set width of numericInputs 
    tags$link(rel = "stylesheet", href = "//fonts.googleapis.com/css?family=Roboto|Roboto+Condensed"), 
    tags$style("body {font-family: 'Roboto', sans-serif;}  
                                         h1 {font-family: 'Roboto Condensed', sans-serif;}  
                                         h2 {font-family: 'Roboto Condensed', sans-serif;} 
                                         h3 {font-family: 'Roboto Condensed', sans-serif;}  
                                         h4 {font-family: 'Roboto Condensed', sans-serif;}  
                                         h5 {font-family: 'Roboto Condensed', sans-serif;}  
                                         h6 {font-family: 'Roboto Condensed', sans-serif;}")    #apply font styles 
  ),
  
  ### Introduction of MRT-SS Calculator on left-hand side panel ###
  sidebarPanel(
    includeHTML("www/sidebar.html")
  ),
  
  ### Main Panel on the right-hand side ###
  mainPanel(
    tags$hr(),
    
    ### Study Setup ###
    
    h3("Study Setup"),
    fluidRow(
      column(3,
             numericInput("days",label = "Duration of the Study (Days)",value = 42)
      ),
      column(4,
             numericInput("occ_per_day",label = "Number of Decision Time Points per Day",value = 5)
      ),
      column(5,
             textOutput("setting_warning")  ### output warnings when you type in wrong format for
      )                                     ### "Duration of the study" and "Number of decision time
    ),                                      ### per day"
    
    tags$hr(),
    
    #### time-varying randomization probability ####
    h3("Randomization Probability"), 
    tabsetPanel(id = "ranPro",
                
                ### Input of constant randomization probability ###
                tabPanel("Constant",
                         br(),
                         fluidRow(
                           column(6,
                                  numericInput("P_intervene",label = "Constant Randomization Probability",value = 0.4 )
                           ),
                           column(6,
                                  textOutput("setting_warning_ranPro") ### Output warnings when you type in wrong format for 
                           )                                           ### "Value of Constant Randomization Probability"
                         )
                ),
                ### Input of time-varying randomization probability ###
                tabPanel("Time-varying",
                         br(),
                         fluidRow(helpText(textOutput("timevar_prob_text"))),
                         hr(), 
                         fluidRow(
                           column(4,
                                  # textOutput("timevar_prob_text"),
                                  radioButtons("numbers", "Number of Inputs:",
                                               c("One per Day" = "re_days",
                                                 "One per Decision Time" = "re_dec")), ### two choices for uploading files ###
                                  ### uploading file respect to decision times ###
                                  conditionalPanel(condition="input.numbers =='re_dec'",
                                                   fileInput('file1', 'Choose a .csv file containing time-varying randomization probabilities (one per decision time) to upload',
                                                             accept = c('.csv')
                                                   )
                                  ),
                                  ###uploading file respect to dyas ###
                                  conditionalPanel(condition="input.numbers =='re_days'",
                                                   fileInput('file2', 'Choose a .csv file containing time-varying randomization probabilities (one per day) to upload',
                                                             accept = c('.csv')
                                                   )
                                  )
                           ),
                           column(4,
                                  ### downloading the sample file for "respecting to days" and
                                  ### "respecting to number of days"
                                  p("If you would like to use a template, you can download one here:"),
                                  downloadButton("timevar_prob_template", "Download Template"),
                                  br(),
                                  textOutput("download_template_caption")
                           ),
                           column(4,
                                  ### showing the first 5 rows of the uploaded file ###
                                  conditionalPanel(condition="input.numbers =='re_dec'",
                                                   p('Showing the first 5 rows of the uploaded file. '),
                                                   tableOutput('P_inter_table_dec')
                                  ),
                                  conditionalPanel(condition="input.numbers =='re_days'",
                                                   p('Showing the first 5 rows of the uploaded file. '),
                                                   tableOutput('P_inter_table_days')
                                  )
                           )
                         )  
                )
    ),
    
    tags$hr(),
    
    #### Expected Availability ####
    h3("Expected Availability"),
    fluidRow(
      column(5,
             ### Three patterns of expected availability to choose from 
             ### quadratic, constant and linear
             selectizeInput("tau_choices", label = "Select one of the following patterns for the expected availability", 
                            choices=list("Quadratic" = "quadratic","Constant"="constant",
                                         "Linear"="linear"),
                            options = list(
                              placeholder = "Please select a pattern",
                              onInitialize = I('function() { this.setValue(0); }')
                            )
             ),
             ### Inputs for constant pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices =='constant' ",
                              sliderInput("tau_mean_a",label="Average of Expected Availability",min = 0, max = 1,value = 0.5)),
             
             ### Inputs for linear pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices == 'linear' ",
                              sliderInput("tau_mean_c",label="Average of Expected Availability",min = 0, max = 1,value = 0.5),
                              sliderInput("tau_initial_c", label = "Initial Value of Expected Availability",
                                          min = 0, max = 1,value = 0.2)),
             
             ### Inputs for quadratic pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices == 'quadratic' ",
                              sliderInput("tau_mean_e",label="Average of Expected Availability",min = 0, max = 1,value = 0.5),
                              sliderInput("quadratic_initial", label = "Initial value of Expected Availability",min =0, max = 1, value = 0.7),
                              numericInput("quadratic_max", label = "Changing Point of Availability", value = 42)),
             
             ### Comments on constant pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices == 'constant'",
                              p(em("Notes: ")," A simplistic constant availability pattern.")
                              
             ),
             ### Comments on linear pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices == 'linear'",
                              p(em("Notes: "), "A linearly increasing pattern of expected availability might be used  if participants
                                will find the intervention useful and thus more likely to turn the intervention on"),
                              p("A linearly decreasing pattern of expected availability might be used if participants learn more about the intervetion
                                and get bored through the course of the study and thus getting less likely to turn on the invervention.")
             ),
             ### Comments on quadratic pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices == 'quadratic'",
                              p(em("Notes: "),"A quadratic pattern of availability. Here the changing point of availability refers to day of either maximal of minimal availability,
                                depending on the input values of initial and average availability")
             ),
             
             ### invisible checkboxInput that will not show on the screen
             ### This invisible checkboxinput is used to check whether a specific 
             ### actionbutton is pressed.
             
             conditionalPanel(condition="input.tau_choices == 'invisible'", ### this condition is nver true, which make it invisible ###
                              
                              ### invisible_a_size: to check whether I pressed the action button for calculating the sample size
                              # for quadratic trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a_size",label = "Invisible_a_size",value = FALSE ),
                              
                              ### invisible_c_size: to check whether I pressed the action button for calculating the sample size
                              # for quadratic trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c_size",label = "Invisible_c_size",value = FALSE ),
                              
                              ### invisible_e_size: to check whether I pressed the action button for calculating the sample size
                              # for quadratic trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e_size",label = "Invisible_e_size",value = FALSE ),
                              
                              ### invisible_a1_size: to check whether I pressed the action button for calculating the sample size
                              # for constant trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a1_size",label = "Invisible_a1_size",value = FALSE ),
                              
                              ### invisible_c1_size: to check whether I pressed the action button for calculating the sample size
                              # for constant trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c1_size",label = "Invisible_c1_size",value = FALSE ),
                              
                              ### invisible_e1_size: to check whether I pressed the action button for calculating the sample size
                              # for constant trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e1_size",label = "Invisible_e1_size",value = FALSE ),
                              
                              ### invisible_a2_size: to check whether I pressed the action button for calculating the sample size
                              # for linear trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a2_size",label = "Invisible_a2_size",value = FALSE ),
                              
                              ### invisible_c2_size: to check whether I pressed the action button for calculating the sample size
                              # for linear trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c2_size",label = "Invisible_c2_size",value = FALSE ),
                              
                              ### invisible_e2_size: to check whether I pressed the action button for calculating the sample size
                              # for linear trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e2_size",label = "Invisible_e2_size",value = FALSE ),
                              
                              ### invisible_a_power: to check whether I pressed the action button for calculating the power
                              # for quadratic trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a_power",label = "Invisible_a_power",value = FALSE ),
                              
                              ### invisible_c_power: to check whether I pressed the action button for calculating the power
                              # for quadratic trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c_power",label = "Invisible_c_power",value = FALSE ),
                              
                              ### invisible_e_power: to check whether I pressed the action button for calculating the power
                              # for quadratic trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e_power",label = "Invisible_e_power",value = FALSE ),
                              
                              ### invisible_a1_power: to check whether I pressed the action button for calculating the power
                              # for constant trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a1_power",label = "Invisible_a1_power",value = FALSE ),
                              
                              ### invisible_c1_power: to check whether I pressed the action button for calculating the power
                              # for constant trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c1_power",label = "Invisible_c1_power",value = FALSE ),
                              
                              ### invisible_e1_power: to check whether I pressed the action button for calculating the power
                              # for constant trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e1_power",label = "Invisible_e1_power",value = FALSE ),
                              
                              ### invisible_a2_power: to check whether I pressed the action button for calculating the power
                              # for linear trend of proximal effect and constant pattern for expected availability
                              checkboxInput("invisible_a2_power",label = "Invisible_a2_power",value = FALSE ),
                              
                              ### invisible_c2_power: to check whether I pressed the action button for calculating the power
                              # for linear trend of proximal effect and linear pattern for expected availability
                              checkboxInput("invisible_c2_power",label = "Invisible_c2_power",value = FALSE ),
                              
                              ### invisible_e2_power: to check whether I pressed the action button for calculating the power
                              # for linear trend of proximal effect and quadratic pattern for expected availability
                              checkboxInput("invisible_e2_power",label = "Invisible_e2_power",value = FALSE )
             )
             
      ),
      column(1),
      column(6,
             
             ### Output the graph of constant pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices=='constant'",
                              
                              plotOutput("cst_graph")
             ),
             
             ### Output the graph of linear pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices=='linear'", 
                              
                              plotOutput("linear_graph")
             ),
             
             ### Output the graph of quadratic pattern of expected availability ###
             conditionalPanel(condition="input.tau_choices=='quadratic'", 
                              
                              plotOutput("quadratic_graph")
             )  
      )
    ),
    tags$hr(),
    
    #### Specifying trend for Proximal Treatment Effect ####
    h3("Proximal Treatment Effect"),
    
    fluidRow(
      column(5,
             ### Three trends of proximal treatment effect to choose from
             ### quadratic, constant, linear
             selectizeInput("beta_choices", label = "Select one of the following trends for the proximal treatment effect", 
                            choices=list("Quadratic" = "beta_quadratic","Constant"="beta_constant",
                                         "Linear"="beta_linear"),
                            options = list(
                              placeholder = "Please select a trend",
                              onInitialize = I('function() { this.setValue(0); }')
                            )
             ),
             ### Inputs for quadratic trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices =='beta_quadratic'",
                              sliderInput("beta_quadratic__mean",label="Average of Proximal Treatment Effect",min = 0, max = 0.2,value = 0.1),
                              numericInput("beta_quadratic__max", label = "Day of Maximal Proximal Treatment Effect", value = 28),
                              numericInput("beta_quadratic__initial", label = "Initial value of Proximal Treatment Effect",value = 0),
                              p(em("Notes"),": The quadratic form of a proximal treatment effect might be used if you expect that 
                                                        initially participants will enthusiastically engage in the apps and thus the 
                                                         proximal effect will get higher. Then, as the study goes on, some participants are 
                                                         likely to disengage or begin to ignore the activity suggestions and hence a downward trend.")
             ),
             ### Inputs for constant trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices =='beta_constant'",
                              sliderInput("beta_constant_mean", label = "Average of Proximal Treatment Effect",min = 0, max = 0.2, value = 0.1),
                              p(em("Notes"),": The proximal treatment effect stays constant over the study.")
             ),
             ### Inputs for linear trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices == 'beta_linear' ",
                              sliderInput("beta_linear_mean",label="Average of Proximal Treatment Effect",min = 0, max = 0.2,value = 0.1),
                              sliderInput("beta_linear_initial", label = "Initial Value of Proximal Treatment Effect",
                                          min = 0, max = 0.2,value = 0),
                              p(em("Notes"),": The linearly increasing form of a proximal treatment effect might be used if participants
                                                         will get more enthusiastically engage in the apps and thus the proximal effect will increase as the 
                                                         study goes."),
                              p("The linearly decreasing form of a proximal treatment effect might be used if participants
                                                        are likely to disengage the activity suggestionss and thus the proximal effect will decrease as the 
                                                         study goes.")
             )
      ),
      column(1),
      column(6,
             ### output graph for quadratic trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices =='beta_quadratic'",
                              plotOutput("beta_graph_quadratic")),
             ### output graph for constant trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices =='beta_constant'",
                              plotOutput("beta_graph_constant")),
             ### output graph for linear trend of proximal treatment effect ###
             conditionalPanel(condition="input.beta_choices =='beta_linear'",
                              plotOutput("beta_graph_linear"))
      )
      
    ),
    
    tags$hr(),
    
    
    ### Specifying whether you are interested in calculating the sample size or the power ###
    fluidRow(
      column(4,
             
             ### Choices to choose from "sample size" or "power" ###
             radioButtons("radio_choices", label = "Are you interested in finding sample size or power?", 
                          choices=list("Sample Size"="choice_sample_size","Power"="choice_power"),
                          selected = "choice_sample_size")
      ),
      column(4,
             
             ### type in the desired power if you want to calculate the sample size ###
             conditionalPanel(condition="input.radio_choices=='choice_sample_size'",
                              numericInput("power", label = HTML("Desired Power"), value = 0.8)
             ),
             
             ### type in the sample size if you want to calculate the power attained ###
             conditionalPanel(condition="input.radio_choices=='choice_power'",
                              numericInput("sample_size", label = HTML("Number of Participants"), value = 40)
             ),
             
             ### type in significance level for both cases ###
             numericInput("sigLev", label = HTML("Significance Level"), value = 0.05)
             
      ),
      column(4,
             
             ### Output warnings if you type in wrong format for "desired power" ###
             conditionalPanel(condition="input.radio_choices=='choice_sample_size'",
                              textOutput("choice_sample_size_warning")
             ),
             
             ### Output warnings if you type in wrong format for "Number of participants" ###
             conditionalPanel(condition="input.radio_choices=='choice_power'",
                              textOutput("choice_power_warning")
             ),
             
             ### Output warnings if you type in wrong format for "Significance level" ###
             textOutput("significance_warning")
      )
    ),
    
    tags$hr(),
    
    ### Output warnings if you didn't choose any patterns for expected availability ###
    conditionalPanel(condition="input.tau_choices !='constant'
                                               && input.tau_choices != 'quadratic'
                                               && input.tau_choices != 'linear'",
                     uiOutput("result_warning")
    ),
    ### Output warnings if you didn't choose any patterns for proximal treatment effect ###
    conditionalPanel(condition="input.beta_choices !='beta_constant'
                                               &&input.beta_choices != 'beta_linear'
                                               && input.beta_choices != 'beta_quadratic'",
                     uiOutput("result_warning_beta")
    ),
    
    
    ##### choice to calculate sample size(action buttons) #####
    
    ### action button to calculate sample size for constant expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                                               && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_a_size","Get Result")
    ),
    
    ### action button to calculate sample size for linear expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                                               && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_c_size","Get Result")
    ),
    
    ### action button to calculate sample size for quadratic expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                                               && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_e_size","Get Result")
    ),
    
    ### action button to calculate sample size for constant expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                                               && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_a1_size","Get Result")
    ),
    
    ### action button to calculate sample size for linear expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                                               && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_c1_size","Get Result")
    ),
    
    ### action button to calculate sample size for quadratic expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                                               && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_e1_size","Get Result")
    ),
    
    ### action button to calculate sample size for constant expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_a2_size","Get Result")
    ),
    
    ### action button to calculate sample size for linear expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_c2_size","Get Result")
    ),
    
    ### action button to calculate sample size for quadratic expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                     actionButton("getResult_e2_size","Get Result")
    ),
    
    ##### choice to calculate power(action buttons) #####
    
    ### action button to calculate power for constant expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_a_power","Get Result")
    ),
    
    ### action button to calculate power for linear expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_c_power","Get Result")
    ),
    
    ### action button to calculate power for quadratic expected availability and 
    ### quadratic proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_e_power","Get Result")
    ),
    
    ### action button to calculate power for constant expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_a1_power","Get Result")
    ),
    
    ### action button to calculate power for linear expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_c1_power","Get Result")
    ),
    
    ### action button to calculate power for quadratic expected availability and 
    ### constant proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                                               && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_e1_power","Get Result")
    ),
    
    ### action button to calculate power for constant expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='constant' 
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_a2_power","Get Result")
    ),
    
    ### action button to calculate power for linear expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices =='linear'
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_c2_power","Get Result")
    ),
    
    ### action button to calculate power for quadratic expected availability and 
    ### linear proximal treatment effect
    conditionalPanel(condition="input.tau_choices=='quadratic'
                                               && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                     actionButton("getResult_e2_power","Get Result")
    ),
    
    br(),
    conditionalPanel(condition="(input.tau_choices =='constant'
                                || input.tau_choices == 'quadratic'
                                || input.tau_choices == 'linear')
                                && (input.beta_choices =='beta_constant'
                                ||input.beta_choices == 'beta_linear'
                                || input.beta_choices == 'beta_quadratic')",
                     tabsetPanel(
                       
                       ### Current result ###
                       tabPanel("Current Result",
                                
                                ################choice to calculate sample size(results) ###############
                                
                                ###current result of sample size for constant expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_a")
                                ),
                                
                                ###current result of sample size for linear expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_c")
                                ),
                                
                                ###current result of sample size for quaratic expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_e")
                                ),
                                
                                ###current result of sample size for constant expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_a1")
                                ),
                                
                                ###current result of sample size for linear expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_c1")
                                ),
                                
                                ###current result of sample size for quaratic expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_e1")
                                ),
                                
                                ###current result of sample size for constant expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_a2")
                                ),
                                
                                ###current result of sample size for linear expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_c2")
                                ),
                                
                                ###current result of sample size for quadratic expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_sample_size' ",
                                                 uiOutput("result_e2")
                                ),
                                
                                ############# choice to calculate the power(results) ####################
                                
                                
                                ###current result of power for constant expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_a_power")
                                ),
                                
                                ###current result of power for linear expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_c_power")
                                ),
                                
                                ###current result of power for quadratic expected availability and
                                ### quadratic proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_quadratic'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_e_power")
                                ),
                                
                                ###current result of power for constant expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_a1_power")
                                ),
                                
                                ###current result of power for linear expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_c1_power")
                                ),
                                
                                ###current result of power for quadratic expected availability and
                                ### constant proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_constant'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_e1_power")
                                ),
                                
                                ###current result of power for constant expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices =='constant' 
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_a2_power")
                                ),
                                
                                ###current result of power for linear expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices =='linear'
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_c2_power")
                                ),
                                
                                ###current result of power for quadratic expected availability and
                                ### linear proximal effect
                                conditionalPanel(condition="input.tau_choices=='quadratic'
                     && input.beta_choices == 'beta_linear'
                     && input.radio_choices == 'choice_power' ",
                                                 uiOutput("result_e2_power")
                                )
                       ),
                       
                       ################result of history table ######################
                       tabPanel("History",         
                                
                                ### download the result table of sample size ###
                                conditionalPanel(condition="input.radio_choices=='choice_sample_size'
                      && (input.tau_choices == 'constant'
                     || input.tau_choices =='linear'
                     || input.tau_choices == 'quadratic' )
                     &&(input.beta_choices == 'beta_constant'
                        || input.beta_choices =='beta_linear'
                        || input.beta_choices == 'beta_quadratic')",
                                                 fluidRow(
                                                   column(2,
                                                          radioButtons("filetype_size", "File type:",
                                                                       choices = c("csv", "tsv"))
                                                   ),
                                                   column(3,
                                                          br(),
                                                          downloadButton('downloadData_size', 'Download')
                                                   )
                                                 ),
                                                 tableOutput("result_table_size") ### Output the table result for sample size ###
                                ),
                                
                                ### download the result table of power ###
                                conditionalPanel(condition="input.radio_choices=='choice_power'
                     && (input.tau_choices == 'constant'
                     || input.tau_choices =='linear'
                     || input.tau_choices == 'quadratic' )
                     &&(input.beta_choices == 'beta_constant'
                        || input.beta_choices =='beta_linear'
                        || input.beta_choices == 'beta_quadratic')",
                                                 fluidRow(
                                                   column(2,
                                                          radioButtons("filetype_power", "File type:",
                                                                       choices = c("csv", "tsv"))
                                                   ),
                                                   column(3,
                                                          br(),
                                                          downloadButton('downloadData_power', 'Download')
                                                   )
                                                 ),
                                                 tableOutput("result_table_power") ### Output the table result for power ###
                                )
                       )
                     )
    ),
    br(),
    br(),
    br(),
    br(),
    br()
    
  ),
  
  
  #### A TOOLTIPS: to give you tips for typeing in numeric inputs. ####
  
  
  bsTooltip(id = "days", title = "integer greater than 0",placement="right", trigger = "focus"),
  bsTooltip(id = "occ_per_day", title = "integer greater than 0",placement="right", trigger = "focus"),
  bsTooltip(id = "P_intervene", title = "Input must range from 0-1" ,placement="right", trigger = "focus"),
  bsTooltip(id = "sigLev", title = "Input must range from 0-1" ,placement="right", trigger = "focus"),
  bsTooltip(id = "power", title = "Input must range from 0-1" ,placement="right", trigger = "focus"),
  bsTooltip(id = "beta_quadratic__initial", title = "Input must be greater than or equal to 0, and less than or equal to Average Standardized Effect",placement="right", trigger = "focus"),
  bsTooltip(id = "beta_quadratic__max", title = "Input must be integer greater than 0",placement="right", trigger = "focus"),
  bsTooltip(id = "sample_size", title = "Input must be integer greater than 1",placement="right", trigger = "focus"),
  bsTooltip(id = "quadratic_max", title = "Input must be integer greater than 0",placement="right", trigger = "focus"),
  
  collapsable = TRUE,
  footer = HTML("<p style = 'font-size:12px'> Please direct correspondence to <a href='mailto:pengliao@umich.edu'>pengliao@umich.edu</a></p>")
  
))
