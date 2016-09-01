# Copyright 2016 Nicholas J. Seewald and Peng Liao
# This file is part of MRT-SS Calculator.

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

# Contact the authors at nseewald@umich.edu or pengliao@umich.edu


library(shiny)
source("SampleSizeCalculator.R")

### Non-reactive function declarations ###
true.beta <- function(q1 = 0, q2 = 0.3, q3 = 28, days){
  
  ### output (coefficients of) standardized treatment effect (quadratic
  ### form) given q1, q2, q3 ###
  ### only used when q3 (in days) is larger than half of the study
  
  beta = matrix(0, 3)
  beta[1] = q1
  a = sum(c(1:(days-1)))
  b = sum(c(1:(days-1))^2)
  mat = t(matrix(c(1, 2*(q3-1),a,b), 2, 2))
  beta[2:3] = solve(mat) %*% matrix(c(0, q2 * days))
  
  return(beta)
}

### Begin Shiny Server ###
shinyServer(function(input,output,session){
  
  ### Render helper text for time-varying randomization probability ###
  output$timevar_prob_text <-
    renderText({paste("To specify time-varying randomization probabilities, upload a .csv file containing
        (index, probability) pairs. Depending on how frequently you want to vary the
        randomization probability, you should provide either",
        input$days,
        "(one per day)
        or",
        input$days * input$occ_per_day,
        "(one per decision time) pairs. "
      )
    })
  
  downloadParams <- reactive({
    filetag <- ifelse(input$numbers == "re_dec",
                      "DecTimes.csv",
                      "Days.csv")
    numrows <- ifelse(input$numbers == "re_dec",
                      input$days * input$occ_per_day,
                      input$days)
    header.text <- ifelse(input$numbers == "re_dec",
                          "decision.time",
                          "day")
    return(list("filetag" = filetag, "numrows" = numrows, "header" = header.text))
  })
  
  ### Download handlers for reactively-created randomization probability templates ###
  output$timevar_prob_template <- downloadHandler(
    filename = function() {paste0("MRT-SS-Randomization-Probability-", downloadParams()$filetag)},
    content = function(file) {
      write.csv(
        as.data.frame(matrix(c(1:downloadParams()$numrows, rep(NA, downloadParams()$numrows)),
                             nrow = downloadParams()$numrows,
                             dimnames = list(NULL, c(downloadParams()$header, "randomization.probability")))),
        file = file, na = "", row.names = F
      )
    },
    contentType = "text/csv"
  )
  
  output$download_template_caption <- renderText({
    paste("The template will contain one row per", 
          ifelse(input$numbers == "re_dec",
                 "decision point.", "day."), 
          "Just fill in your desired randomization probabilities, and upload the file.")
  })
  
  ### Reading the file for time-varying randomization probability###
  #### Reading the file of decision times ###
  P_inter_dec <- reactive({    
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE, sep = ",",
             col.names = c("Dec.Times", "Randomization.Probability"))
  })
  
  #### Output the first five rows of the table reading from the file with respect to decision times 
  ### and output warnings if the format of the file is not correct
  output$P_inter_table_dec <- renderTable({       
    data <- as.vector(P_inter_dec()$Randomization.Probability)
    validate(
      need(!is.null(input$file1), "Please upload a file."),
      need(is.null(input$file1) || length(data) == input$days * input$occ_per_day, 
           "Error: The number of rows in the uploaded file does not match the number of decision times provided in Study Setup. Either upload a corrected file, or adjust your Study Setup."),
      need(is.null(input$file1) || max(data) <= 1, 
           "Error: The provided randomization probability is greater than 1 for one or more decision times. Please upload a corrected file."),
      need(is.null(input$file1) || min(data) >= 0, 
           "Error: The provided randomization probability is less than 0 for one or more decision times. Please upload a corrected file.")
    )
    head(P_inter_dec(), n = 5)
  })
  
  ### Reading the file with respect to days ###
  P_inter_days <- reactive({     
    inFile <<- input$file2
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE, sep = ",",
             col.names = c("Days", "Randomization.Probability"))
  })
  
  #### Output the first five rows of the table reading from the file with respect to days
  ### and output warnings if the format of the file is not correct
  output$P_inter_table_days <- renderTable({
    #### Output the first five rows of the table reading from the file for days
    data <- as.vector(P_inter_days()$Randomization.Probability)
    validate(
      need(!is.null(input$file2), "Please upload a file."),
      need(is.null(input$file2) || length(data) == input$days , 
           "Error: The number of rows in the uploaded file does not match the number of days provided in Study Setup. Either upload a corrected file, or adjust your Study Setup."),
      need(is.null(input$file2) || max(data) <= 1, 
           "Error: The provided randomization probability is greater than 1 for one or more days. Please upload a corrected file."),
      need(is.null(input$file2) || min(data) >= 0, 
           "Error: The provided randomization probability is less than 0 for one or more days. Please upload a corrected file.")
    )
    head(P_inter_days(), n = 5)
  })

                 #### Data generating for proximal treatment effect ####
  
  ### Quadratic class of proximal treatment effect ###
  beta_quadratic__input <- reactive({               
    N = input$days*input$occ_per_day
    H = input$beta_quadratic__max * input$occ_per_day
    M = input$beta_quadratic__mean
    I = input$beta_quadratic__initial
    
    beta1 <- true.beta(q1 = 0, q2 = M-I, q3 = input$beta_quadratic__max , days = input$days)
    b <- beta1[2]; c <- beta1[3]
    sequence = c(0:(input$days-1))
    input <- I + b*sequence + c * sequence^2
  })
  
  ### Constant class of proximal treatment effect ###
  beta_constant_input <- reactive({                    
    replicate(input$days, input$beta_constant_mean)
  })
  
  ### Linear class of proximal treatment effect ###
  beta_linear_input <- reactive({                  
    initial = input$beta_linear_initial
    mean = input$beta_linear_mean
    days = input$days
    
    range = (initial - mean) * 2
    
    if(days == 1) {
      input = mean
    } else if (range == 0) {
      replicate(days, mean)
    } else {
      num = range / (days - 1)
      input = seq(from = mean + range / 2,
                  to = mean - range / 2,
                  by = -num)
    }
  })
  
  ### plot of the graphs for the proximal treatment effect ###
  
  ### plot the graph of constant trend of proximal treatment effect. ###
  output$beta_graph_constant <- renderPlot({               
    validate(
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0")
    )
    up = 0.2
    if(max(beta_constant_input()) > 0.2)
    {
      up = max(beta_constant_input())
    }
    plot(beta_constant_input(),xlab = "Days", ylab = "Proximal Effect", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_constant_input()), y = rep(0, length(beta_constant_input())), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=input$beta_constant_mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")
  })
  
  ### plot the graph of linear trend of proximal treatment effect. ###
  output$beta_graph_linear <- renderPlot({
    validate(
      need(min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0")
    )
    
    up = 0.2
    if(max(beta_linear_input()) > 0.2)
    {
      up = max(beta_linear_input())
    }
    plot(beta_linear_input(),xlab = "Days", ylab = "Proximal Effect", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_linear_input()), y = rep(0, length(beta_linear_input())), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=input$beta_linear_mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")
  })
  
  ### plot the graph of quadratic trend of proximal treatment effect. ###
  output$beta_graph_quadratic <-renderPlot({
    
    validate(
      need(input$beta_quadratic__mean > 0,"Error: Please specify the average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify the average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0 ,"Warning: Some values of proximal treatment effect are less than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify a Standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify Standardized Initial Effect less than or equal to Average Standardized Effect"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0")
      
    )
    
    up = 0.2
    if(max(beta_quadratic__input()) > 0.2)
    {
      up = max(beta_quadratic__input())
    }
    
    plot(beta_quadratic__input(), xlab = "Days", ylab = "Effects", ylim = c(0,up), type = "o",
         pch = 16, cex = 0.8, col = 2)
    points(x = 1:length(beta_quadratic__input()), y = rep(0, length(beta_quadratic__input())), type = 'o', pch = 16, col = 4, cex = 0.8)
    abline(h=input$beta_quadratic__mean, lty = 2, col = 1)
    legend("topleft", cex = 1, legend=c('Null Hypothesis','Alternate Hypothesis','Average Effect'), col = c(4,2,1),lty = c(1,1,2), pch=c(16,16,NA),bty = "n")
    
  })
  
  
           #### Expected Availability ####
  
  ### Constant class of expected availability ###
  constant_input <- reactive({                      
    replicate(input$days, input$tau_mean_a)
  })
  
  ### Linear class of expected availability ###
  linear_input <- reactive({                      
    initial = input$tau_initial_c
    mean = input$tau_mean_c
    
    range = (initial - mean) * 2
    days = input$days
    if(days == 1)
    {
      input <- mean
    }
    else if(range == 0)
    {
      
      input <- replicate(days,mean)
    }
    else
    {
      num <- range/(days - 1)
      input <- seq(from = mean+range/2, to = mean-range/2, by = -num)
    }
    
    return(input)
  })
  
  
  ### data generation for quadratic class of expected availability ###
  quadratic_input <- reactive({           
    
    N = input$days*input$occ_per_day
    H = input$quadratic_max * input$occ_per_day
    M = input$tau_mean_e
    I = 2*M - input$quadratic_initial
    days = input$days
    
    
    true.beta <- function(q1 = 0, q2 = 0.3, q3 = 28, days){
      
      ### output (coefficients of) standardized treatment effect (quadratic
      ### form) given q1, q2, q3 ###
      ### only used when q3 (in days) is larger than half of the study
      
      beta = matrix(0, 3)
      beta[1] = q1
      a = sum(c(1:(days-1)))
      b = sum(c(1:(days-1))^2)
      mat = t(matrix(c(1, 2*(q3-1),a,b), 2, 2))
      beta[2:3] = solve(mat) %*% matrix(c(0, q2 * days))
      
      return(beta)
    }
    
    beta1 <- true.beta(q1 = 0, q2 = M-I, q3 = input$quadratic_max , days = input$days)
    b <- beta1[2]; c <- beta1[3]
    sequence = c(0:(input$days-1))
    input <- I + b*sequence + c * sequence^2
    input <- 2*M-input
    
    return(input)
  })
          ### Plot the graph for expected availability ###
  
  ### plot the graph of constant pattern of expected availability ###
  output$cst_graph <- renderPlot({      
    validate(
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0")
      
    )
    plot(constant_input(),xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=input$tau_mean_a, lty = 2)
    legend("topleft", legend=c('Availability','Average Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")
  })
  
  ### plot the graph of linear pattern of expected availability ###
  output$linear_graph <- renderPlot({
    validate(
      need(min(linear_input()) >= 0, "Warning: Some values of the availability are less than 0."),
      need(max(linear_input()) <= 1, "Warning: Some values of the availability are bigger than 1."),
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0")
    )
    plot(linear_input(),xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=input$tau_mean_c, lty = 2)
    legend("topleft", legend=c('Availability','Average of Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")
    
  })
  
  ### plot the graph of quadratic pattern of expected availability ###
  output$quadratic_graph <- renderPlot({
    validate(
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0"),
      need(min(quadratic_input()) >= 0,"Warning: Some values of Availability are less than 0."),
      need(max(quadratic_input()) <= 1 ,"Warning: Some values of Availability are bigger than 1")
    )
    plot(quadratic_input(),xlab = "Days", ylab = "Availability", ylim = c(0,1), type = "o",
         pch = 16, cex = 0.8, col = 4)
    abline(h=input$tau_mean_e, lty = 2)
    legend("topleft", legend=c('Availability','Average of Availability'), col = c(4,1),lty = c(1,2), pch=c(16,NA),bty = "n")
  })
  
  ### generating the warnings for the numeric input: "Duration of the Study (Days)"
  ### and "Number of Decision Time Points per day"
  output$setting_warning <-renderText({                         
    validate(
      need(input$days == round(input$days),"Error: Please enter integer values for the number of days"),
      need(input$days > 0 ,"Error: Please specify the number of days greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error:Please enter integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0")
    )
  })
  
  ### generating the warnings for numeric input: "Average of Randomization Probability". ###
  output$setting_warning_ranPro <-renderText({ 
    validate(
      need(input$P_intervene > 0,"Error: Please specify the randomization probability greater than 0"), 
      need(input$P_intervene < 1,"Error: Please specify the randomization probability less than 1")
    )
  })
  
  ### generating the warnings for numeric input: "Desired Power" ###
  output$choice_sample_size_warning <-renderText({      
    validate(
      need(input$power >= 0,"Error: Please specify the power greater than or equal to 0"), 
      need(input$power <= 1,"Error: Please specify the power less than or equal to 1")
    )
  })
  
  ### generating the warnings for numeric input: "Number of Participants" ###
  output$choice_power_warning <-renderText({      
    
    validate(
      need(input$sample_size > 0,"Error: Please specify Number of Participants greater than 0"), 
      need(input$sample_size == round(input$sample_size),"Error: Please enter integer value for Number of Participants")
    )
  })
  
  ### generating warnings for numeric input: "Significance Level" ###
  output$significance_warning <- renderText({     
    validate(
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"), 
      need(input$sigLev <= 1,"Error: Please specify the significance level less than or equal to 1")
    )
  })
  
  ### generating warnings if you don't choose a pattern of expected availability. ###
  output$result_warning <- renderUI({
    validate(
      need(input$tau_choices == "constant"
           || input$tau_choices =="linear"
           || input$tau_choices == "quadratic",
           "Please select a pattern for expected availability")
    )
  })
  
  ### generating warnings if you don't choose a trend of proximal treatment effect. ###
  output$result_warning_beta <- renderUI({
    
    validate(
      need(input$beta_choices == "beta_constant"
           || input$beta_choices =="beta_linear"
           || input$beta_choices == "beta_quadratic",
           "Please select a pattern for proximal treatment effect")
    )
  })
                       ### Output the current results of sample size ### 
  
  ### Generate the current result of sample size for Quadratic proximal treatment effect and 
  ### constant expected availability
  resultA <- eventReactive(input$getResult_a_size, {  ### Generate this current result of sample size if the corresponding 
                                                      ### action button is pressed
    days = input$days
    occ = input$occ_per_day
    Total = days*occ;
    
    
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    
    ### proximal treatment effect is consant on each day ###
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    ### expected availability is constant on each day ###
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        N <- SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of sample size for quadratic proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,
           "Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0,"Warning: Some values of proximal treatment effect are less than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_size()
    resultA()
  })
  
  ### Generate the current result of sample size for quadratic proximal treatment effect and 
  ### linear expected availability
  resultC <- eventReactive(input$getResult_c_size, { ### Generate this current result of sample size if the corresponding 
                                                     ### action button is pressed
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    
    ### We assume that the proximal treatment effect is consant on each day ###
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    ### We assume that the expected availability is constant on each day ###
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        N <- SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
  })
  
  ### Output the current result of sample size for quadratic proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0,"Warning: Some values of proximal treatment effect are less than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need( min(linear_input()) >= 0, "Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1, "Warning: Some values of Availability are bigger than 1")
    )
    data_get_size()
    resultC()
    
  })
  
  ### Generate the current result of sample size for quadratic proximal treatment effect and 
  ### quadratic expected availability
  resultE<- eventReactive(input$getResult_e_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    })
  
  ### Output the current result of sample size for quadratic proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0 ,"Warning: Some values of proximal treatment effect are less than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are greater than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1,"Warning: Some values of Availability are greater than 1")
      
    )
    data_get_size()
    resultE()
    
  })
  
  ### Generate the current result of sample size for constant proximal treatment effect and 
  ### constant expected availability
  resultA1<- eventReactive(input$getResult_a1_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of sample size for constant proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a1 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_size()
    resultA1()
    
  })
  
  ### Generate the current result of sample size for constant proximal treatment effect and 
  ### linear expected availability
  resultC1<- eventReactive(input$getResult_c1_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of sample size for constant proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c1 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(linear_input()) >= 0,"Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_size()
    resultC1()
    
  })
  
  ### Generate the current result of sample size for constant proximal treatment effect and 
  ### quadratic expected availability
  resultE1<- eventReactive(input$getResult_e1_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
    
  })
  
  ### Output the current result of sample size for constant proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e1 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1  ,"Warning: Some values of Availability are greater than 1")
      
    )
    data_get_size()
    resultE1()
    
  })
  
  ### Generate the current result of sample size for linear proximal treatment effect and 
  ### constant expected availability
  resultA2<- eventReactive(input$getResult_a2_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_a, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
  })
  
  ### Output the current result of sample size for linear proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a2 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      
      need(min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1 , "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_size()
    resultA2()
    
  })
  
  ### Generate the current result of sample size for linear proximal treatment effect and 
  ### linear expected availability
  resultC2<- eventReactive(input$getResult_c2_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_c, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of sample size for linear proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c2 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      need(min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(linear_input()) >= 0,"Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_size()
    resultC2()
    
  })
  
  ### Generate the current result of sample size for linear proximal treatment effect and 
  ### quadratic expected availability
  resultE2<- eventReactive(input$getResult_e2_size, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    
    if(input$ranPro == "Constant"){    
      
      ### If the randomization probability is constant ###
      
      delta <- input$P_intervene;
      
      if(input$P_intervene > 0 && input$P_intervene < 1){ 
        
        
        
        N <- SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
        
        if(N > 10){
          HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
        }else{
          ### if the calculated sample size is less than 10, we won't output the exact sample size ###
          HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
          
        }
      }else{
        ### Output the warnings if the randomization probability is not correctly specified ###
      }
    }else{   
      
      ### If the randomization probability is time-varying ###
      
      if(input$numbers == "re_days"){  
        
        ### if the input file is respect to days ###
        
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }
            else
            {      ### if the calculated sample size is less than 10, we won't output the exact sample size ###
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{   
        
        ### if the input file is respect to decision times ###
        
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            
            N <- SampleSize (b_input, input_e, delta, alpha0=input$sigLev, beta0=input$power, setup=list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
            
            if(N > 10){
              HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            }else{
              HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of sample size for linear proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e2 <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$power >= 0 ,"Error: Please specify power greater than or equal to 0"),
      need(input$power <= 1 ,"Error: Please specify the power less than or equal to 1"),
      need(min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1  ,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_size()
    resultE2()
    
  })
  
  
             ### Output the current results of power ### 
  
  ### Generate the current result of power for quadratic proximal treatment effect and 
  ### constant expected availability
  resultA_power <- eventReactive(input$getResult_a_power, { ### Generate this current result of power if the corresponding 
                                                            ### action button is pressed
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    
    size = input$sample_size
    
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size, "when the significance level is", input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for quadratic proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a_power <- renderUI({
    
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0,"Warning: Some values of proximal treatment effect are les than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_power()
    resultA_power()
    
  })
  
  ### Generate the current result of power for quadratic proximal treatment effect and 
  ### linear expected availability
  resultC_power<- eventReactive(input$getResult_c_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
  })
  
  ### Output the current result of power for quadratic proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0,"Warning: Some values of proximal treatment effect are les than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(linear_input()) >= 0,"Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_power()
    resultC_power()
    
  })
  
  ### Generate the current result of power for quadratic proximal treatment effect and 
  ### quadratic expected availability
  resultE_power<- eventReactive(input$getResult_e_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
  })
  
  ### Output the current result of power for quadratic proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_quadratic__mean > 0,"Error: Please specify an average standardized effect greater than 0"),
      need(input$beta_quadratic__mean < 1,"Error: Please specify an average standardized effect less than 1"),
      need(min(beta_quadratic__input()) >= 0,"Warning: Some values of proximal treatment effect are less than 0"),
      need(max(beta_quadratic__input()) <= 1,"Warning: Some values of proximal treatment effect are bigger than 1"),
      need(input$beta_quadratic__initial >= 0, "Error: Please specify the standardized Initial Effect greater than or equal to 0"),
      need(input$beta_quadratic__initial <= input$beta_quadratic__mean, "Error: Please specify the standardized Initial Effect less than or equal to average standardized effect"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1 ,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_power()
    resultE_power()
    
  })
  
  ### Generate the current result of power for constant proximal treatment effect and 
  ### constant expected availability
  resultA1_power<- eventReactive(input$getResult_a1_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    size = input$sample_size
    
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for constant proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a1_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_power()
    resultA1_power()
    
  })
  
  ### Generate the current result of power for constant proximal treatment effect and 
  ### linear expected availability
  resultC1_power<- eventReactive(input$getResult_c1_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
  })
  
  ### Output the current result of power for constant proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c1_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(linear_input()) >= 0,"Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_power()
    resultC1_power()
    
  })
  
  ### Generate the current result of power for constant proximal treatment effect and 
  ### quadratic expected availability
  resultE1_power<- eventReactive(input$getResult_e1_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for constant proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e1_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(input$beta_constant_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need( min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1 ,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_power()
    resultE1_power()
    
  })
  
  ### Generate the current result of power for linear proximal treatment effect and 
  ### constant expected availability
  resultA2_power<- eventReactive(input$getResult_a2_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_a = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
    }
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for linear proximal treatment effect and 
  ### constant expected availability 
  ### Output the validation errors
  output$result_a2_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need( min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_a > 0, "Error: Please specify the mean of availability greatert than 0")
      
    )
    data_get_power()
    resultA2_power()
    
  })
  
  ### Generate the current result of power for linear proximal treatment effect and 
  ### linear expected availability
  resultC2_power<- eventReactive(input$getResult_c2_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_c = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
    }
    
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for linear proximal treatment effect and 
  ### linear expected availability 
  ### Output the validation errors
  output$result_c2_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need(min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_c > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(linear_input()) >= 0,"Warning: Some values of Availability are less than 0"),
      need(max(linear_input()) <= 1,"Warning: Some values of Availability are bigger than 1")
    )
    data_get_power()
    resultC2_power()
    
  })
  
  ### Generate the current result of power for linear proximal treatment effect and 
  ### quadratic expected availability
  resultE2_power<- eventReactive(input$getResult_e2_power, {
    days = input$days
    occ = input$occ_per_day
    Total = days*occ
    input_e = vector('numeric', length = Total)
    
    b_input = vector('numeric',Total)
    for(k in 1:days)
    {
      b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
    }
    
    for(k in 1:input$days)
    {
      input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
    }
    
    
    size = input$sample_size
    if(input$ranPro == "Constant"){
      
      delta <- input$P_intervene;
      
      power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
      
      if(power >= 0.5){
        
        HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."));
        
      }
      else
      { 
        ### if the calculated power is less than 50%, output an warning ###
        HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
        
      }
    }else{
      if(input$numbers == "re_days"){
        if(!is.null(input$file2)){
          
          delta <- as.vector(P_inter_days()$Randomization.Probability)
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }else{
        if(!is.null(input$file1)){
          
          delta <- as.vector(P_inter_dec()$Randomization.Probability);
          
          if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
            
            power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
            
            if(power >= 0.5){
              HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
            else
            {
              HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            }
          }else{
            HTML(paste("<h5 style = 'color:red';> Wrong format with the uploaded file."))
          }
        }else{
          HTML(paste("<h5 style = 'color:red';> No uploaded file."))
        }
      }
    }
    
  })
  
  ### Output the current result of power for linear proximal treatment effect and 
  ### quadratic expected availability 
  ### Output the validation errors
  output$result_e2_power <- renderUI({
    validate(
      need(input$days == round(input$days),"Error: Please enter an integer value for duration of the study"),
      need(input$days > 0 ,"Error: Please specify the duration of the study greater than 0"),
      need(input$occ_per_day == round(input$occ_per_day),"Error: Please enter an integer for the number of occasions per day"),
      need(input$occ_per_day > 0 ,"Error: Please specify the number of occasions per day greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene > 0,"Error: Please specify randomization probability greater than 0"),
      need(!input$ranPro == "Constant" || input$P_intervene < 1,"Error: Please specify randomization probability less than 1"),
      need(input$sigLev >= 0,"Error: Please specify the significance level greater than or equal to 0"),
      need(input$sigLev <= 1,"Error: Please specify the significance Level less than or equal to 1"),
      need(input$sample_size > 0 ,"Error: Please specify the number of participants greater than 0"),
      need(input$sample_size == round(input$sample_size) ,"Error: Please enter an integer value for the number of participants"),
      
      need( min(beta_linear_input()) >= 0, "Warning: Some values of the proximal treatment effect are less than 0."),
      need(max(beta_linear_input()) <= 1, "Warning: Some values of the proximal treatment effect are bigger than 1."),
      need(input$beta_linear_mean > 0, "Error: Please specify the mean of proximal treatment effect greatert than 0"),
      need(input$tau_mean_e > 0, "Error: Please specify the mean of availability greatert than 0"),
      need(min(quadratic_input()) >= 0 ,"Warning: Some values of Availability are less than 0"),
      need(max(quadratic_input()) <= 1 ,"Warning: Some values of Availability are bigger than 1")
      
    )
    data_get_power()
    resultE2_power()
    
  })
  
  
  
  ###################calculate the sample size (history table output) ##################
  
  
  ### Whenever the action button is pressed, the corresponding checkBox is set to be true,
  ### which says that we can add the element to the history table.  
  ### These invisible checkBox is for calculating the sample size
  observeEvent(input$getResult_a_size,{
    updateCheckboxInput(session, "invisible_a_size", value = TRUE)
  })
  observeEvent(input$getResult_c_size,{
    updateCheckboxInput(session, "invisible_c_size", value = TRUE)
  })
  observeEvent(input$getResult_e_size,{
    updateCheckboxInput(session, "invisible_e_size", value = TRUE)
  })
  observeEvent(input$getResult_a1_size,{
    updateCheckboxInput(session, "invisible_a1_size", value = TRUE)
  })
  observeEvent(input$getResult_c1_size,{
    updateCheckboxInput(session, "invisible_c1_size", value = TRUE)
  })
  observeEvent(input$getResult_e1_size,{
    updateCheckboxInput(session, "invisible_e1_size", value = TRUE)
  })
  observeEvent(input$getResult_a2_size,{
    updateCheckboxInput(session, "invisible_a2_size", value = TRUE)
  })
  observeEvent(input$getResult_c2_size,{
    updateCheckboxInput(session, "invisible_c2_size", value = TRUE)
  })
  observeEvent(input$getResult_e2_size,{
    updateCheckboxInput(session, "invisible_e2_size", value = TRUE)
  })

  ### the global vectors which contain the information added to the history table. ###
  data_effect_size = vector("character")
  data_avail_size=vector("character")
  data_size = vector("integer")
  data_power = vector("numeric")
  data_sig_size = vector('numeric')
  data_rsize = vector("integer")
  data_rpower = vector('numeric')
  ave_effect_size = vector('numeric')
  ave_avail_size = vector('numeric')
  data_ranPro_size = vector("character")
  
  
  ### data generation for the history table of sample size ###
  data_get_size <- reactive({
    
      if(input$invisible_a_size){ ### if the statement is true, we can add the result for calculating
                                  ### sample size for quadratic proximal effect and constant expected availability
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1 && min(beta_quadratic__input()) >= 0 
         && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0 && input$beta_quadratic__initial <= input$beta_quadratic__mean
         && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {   
          ### If te randomization probability is constant ###
          delta <- input$P_intervene
          N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            
            ### Once the result is added to the history table
            ### the value of checkBox is set to be false to prevent
            ### duplicate adding. The result will be added again until
            ### the corresponding action button is pressed to set the value 
            ### of checkBox to be true. 
            updateCheckboxInput(session, "invisible_a_size", value = FALSE)                                                   
          }                                                                 
          else
          { ### the case where the calculated sample size is less than 10 ###
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_a_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Quadratic")
          data_avail_size <<- c(data_avail_size, "Constant")
          ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){   ### If the randomiztion probability is time-varying. ###
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){  ### If the uploaded file is respect to days ###
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a_size", value = FALSE)
                }
                else          
                {   ### If the calculated sample size is less than or equal to 10 ###
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                ### If the uploaded file with respect to days
                ### is not corrected formatted, the result won't be added to the history table.
                updateCheckboxInput(session, "invisible_a_size", value = FALSE)
              }
            }else{
              ### If you don't upload the file with respect to days, 
              ### the result won't be added to the history table.
              updateCheckboxInput(session, "invisible_a_size", value = FALSE)
            }
          }else{   
            ### If the uploaded file is respect to decision times ###
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a_size", value = FALSE)
                }else{
                  ### If the calculated sample size is less than or equal to 10 ###
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{ ### If the uploaded file with respect to decision times 
                     ### is not corrected formatted, the result won't be added to the history table.
                updateCheckboxInput(session, "invisible_a_size", value = FALSE)  
              }
            }else{  ### If you don't upload the file with respect to decision times, 
                    ### the result won't be added to the history table.
              updateCheckboxInput(session, "invisible_a_size", value = FALSE)
            }
          }
        }else{
          updateCheckboxInput(session, "invisible_a_size", value = FALSE)
        }
      }else{ 
        
        ### If some numeric input is not correct or the graph of proximal treatment effect is not correct 
        ### or the graph of expected availability is not correct, 
        ### the result won't be added to the history table.
        updateCheckboxInput(session, "invisible_a_size", value = FALSE)
      }
    }
    
    if(input$invisible_c_size){  ### if the statement is true, we can add the result for calculating
                                 ### sample size for quadratic proximal effect and linear expected availability
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1 && min(beta_quadratic__input()) >= 0 
         && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0 && input$beta_quadratic__initial <= input$beta_quadratic__mean
         && input$tau_mean_c > 0 && min(linear_input()) >= 0 && max(linear_input()) <= 1){
        
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene
          N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Quadratic")
          data_avail_size <<- c(data_avail_size, "Linear")
          ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
          data_sig_size <<- c(data_sig_size, input$sigLev)
          
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
                
              }else{
                updateCheckboxInput(session, "invisible_c_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
                
              }else{
                updateCheckboxInput(session, "invisible_c_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c_size", value = FALSE)
      }
    }
    
    if(input$invisible_e_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for quadratic proximal effect and quadratic expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1 && min(beta_quadratic__input()) >= 0 
         && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0 && input$beta_quadratic__initial <= input$beta_quadratic__mean
         && input$tau_mean_e > 0 && max(quadratic_input()) <= 1 && min(quadratic_input()) >= 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene
          N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Quadratic")
          data_avail_size <<- c(data_avail_size, "Quadratic")
          ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=3, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Quadratic")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_quadratic__mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e_size", value = FALSE)
      }
      
      
    }
    
    if(input$invisible_a1_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for constant proximal effect and constant expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_constant_mean > 0 && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene
          N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Constant")
          data_avail_size <<- c(data_avail_size, "Constant")
          ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_a1_size", value = FALSE)
      }
      
    }
    
    if(input$invisible_c1_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for constant proximal effect and linear expected availability
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_constant_mean > 0 && input$tau_mean_c > 0 && min(linear_input()) >= 0 && max(linear_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene
          N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Constant")
          data_avail_size <<- c(data_avail_size, "Linear")
          ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000);
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c1_size", value = FALSE)
      }
    }
    
    
    if(input$invisible_e1_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for constant proximal effect and quadratic expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && input$beta_constant_mean > 0 && input$tau_mean_e > 0 && min(quadratic_input()) >= 0 && max(quadratic_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene; 
          N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
          
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Constant")
          data_avail_size <<- c(data_avail_size, "Quadratic")
          ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=1, q=3, Nmax=1000);
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Constant")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_constant_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e1_size", value = FALSE)
      }
      
    }
    
    if(input$invisible_a2_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for linear proximal effect and constant expected availability
      
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1 && input$beta_linear_mean > 0 && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene; 
          N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
          
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Linear")
          data_avail_size <<- c(data_avail_size, "Constant")
          ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                N = SampleSize (b_input, input_a, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Constant")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_a)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_a2_size", value = FALSE)
      }
      
    }
    if(input$invisible_c2_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for linear proximal effect and linear expected availability
      
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1 && input$beta_linear_mean > 0
         && input$tau_mean_c > 0 && min(linear_input()) >= 0 && max(linear_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene; 
          N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
          
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Linear")
          data_avail_size <<- c(data_avail_size, "Linear")
          ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                
                N = SampleSize (b_input, input_c, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Linear")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_c)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c2_size", value = FALSE)
      }
      
    }
    if(input$invisible_e2_size)
    {### if the statement is true, we can add the result for calculating
      ### sample size for linear proximal effect and quadratic expected availability
      
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$power >= 0 && input$power <= 1
         && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1 && input$beta_linear_mean > 0 && input$tau_mean_e > 0
         && min(quadratic_input()) >= 0 && max(quadratic_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene; 
          N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
          
          if(N > 10){
            HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
            data_size <<- c(data_size,N)
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
            
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
            data_size <<- c(data_size, "less than or equal to 10")  
            data_rpower <<- c(data_rpower, input$power)
            data_ranPro_size <<- c(data_ranPro_size, input$P_intervene)
            updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
          }
          data_effect_size <<- c(data_effect_size, "Linear")
          data_avail_size <<- c(data_avail_size, "Quadratic")
          ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
          ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
          data_sig_size <<- c(data_sig_size, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
                  
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                
                N = SampleSize (b_input, input_e, delta = delta, alpha0=input$sigLev, beta0=input$power, setup = list(days = days, occ.per.day = occ), p=2, q=3, Nmax=1000)
                
                if(N > 10){
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is ", N, "to attain", input$power*100,"% power when the significance level is",input$sigLev,".")) 
                  data_size <<- c(data_size,N)
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The required sample size is less than or equal to 10 to attain", input$power*100,"% power when the significance level is",input$sigLev,". Please refer to the result section in the left column for suggestions.")) 
                  data_size <<- c(data_size, "less than or equal to 10")  
                  data_rpower <<- c(data_rpower, input$power)
                  data_ranPro_size <<- c(data_ranPro_size, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
                }
                data_effect_size <<- c(data_effect_size, "Linear")
                data_avail_size <<- c(data_avail_size, "Quadratic")
                ave_effect_size <<- c(ave_effect_size, input$beta_linear_mean)
                ave_avail_size <<- c(ave_avail_size, input$tau_mean_e)
                data_sig_size <<- c(data_sig_size, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e2_size", value = FALSE)
      }
    }
    
    ### Generate the data frame of sample size which contains all the vectors of information we have added through the above procedure ### 
    data = data.frame("Sample Size" = data_size, "Treatment Effect" = data_effect_size, "Average Effect" = ave_effect_size, "Availability" = data_avail_size, "Average Availability" = ave_avail_size,"Type I Error" = data_sig_size,"Power" = data_rpower,
                      "Randomization Probability" = data_ranPro_size)
    
    
  })
  
  
  ###################calculate the power (history table output) ##################
  
  
  ### Whenever the action button is pressed, the corresponding checkBox is set to be true,
  ### which says that we can add the element to the history table.  
  ### These invisible checkBox is for calculating the power
  observeEvent(input$getResult_a_power,{
    updateCheckboxInput(session, "invisible_a_power", value = TRUE)
  })
  observeEvent(input$getResult_c_power,{
    updateCheckboxInput(session, "invisible_c_power", value = TRUE)
  })
  observeEvent(input$getResult_e_power,{
    updateCheckboxInput(session, "invisible_e_power", value = TRUE)
  })
  observeEvent(input$getResult_a1_power,{
    updateCheckboxInput(session, "invisible_a1_power", value = TRUE)
  })
  observeEvent(input$getResult_c1_power,{
    updateCheckboxInput(session, "invisible_c1_power", value = TRUE)
  })
  observeEvent(input$getResult_e1_power,{
    updateCheckboxInput(session, "invisible_e1_power", value = TRUE)
  })
  observeEvent(input$getResult_a2_power,{
    updateCheckboxInput(session, "invisible_a2_power", value = TRUE)
  })
  observeEvent(input$getResult_c2_power,{
    updateCheckboxInput(session, "invisible_c2_power", value = TRUE)
  })
  observeEvent(input$getResult_e2_power,{
    updateCheckboxInput(session, "invisible_e2_power", value = TRUE)
  })
  
  ### the global vectors which contain the information added to the history table. ###
  data_effect_power = vector("character")
  data_avail_power = vector("character")
  ave_effect_power = vector("numeric")
  ave_avail_power = vector("numeric")
  data_sig_power = vector("numeric")
  data_ranPro_power = vector("character")
  
  ### data generation for the history table of power ###
  data_get_power <- reactive({
    
    if(input$invisible_a_power){ ### if the statement is true, we can add the result for calculating
      ### the power for quadratic proximal effect and constant expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      
      size = input$sample_size
      
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1
         && min(beta_quadratic__input()) >= 0 && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0
         && input$beta_quadratic__initial <= input$beta_quadratic__mean && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        { ### If te randomization probability is constant ###
          
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            
            ### Once the result is added to the history table
            ### the value of checkBox is set to be false to prevent
            ### duplicate adding. The result will be added again until
            ### the corresponding action button is pressed again to set the value 
            ### of checkBox to be true. 
            updateCheckboxInput(session, "invisible_a_power", value = FALSE)
          }
          else
          {### If the calculated power is less than 50% ###
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_a_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Quadratic")
          data_avail_power <<- c(data_avail_power, "Constant")
          ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){ ### If the randomiztion probability is time-varying. ###
          if(input$numbers == "re_days"){  ### If the uploaded file is respect to days ###
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a_power", value = FALSE)
                }
                else
                {### If the calculated power is less than 50% ###
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                ### If the uploaded file with respect to days
                ### is not corrected formatted, the result won't be added to the history table.
                updateCheckboxInput(session, "invisible_a_power", value = FALSE)
              }
            }else{
              ### If you don't upload the file with respect to days, 
              ### the result won't be added to the history table.
              updateCheckboxInput(session, "invisible_a_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){  ### If the uploaded file is respect to decision times ###
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a_power", value = FALSE)
                }
                else
                {### If the calculated power is less than 50% ###
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{### If the uploaded file with respect to decision times 
                ### is not corrected formatted, the result won't be added to the history table.
                updateCheckboxInput(session, "invisible_a_power", value = FALSE)
              }
            }else{
              ### If you don't upload the file with respect to decision times, 
              ### the result won't be added to the history table.
              updateCheckboxInput(session, "invisible_a_power", value = FALSE)
            }
          }
        }else{
          updateCheckboxInput(session, "invisible_a_power", value = FALSE)
        }
      }else{
        ### If some numeric input is not correct or the graph of proximal treatment effect is not correct 
        ### or the graph of expected availability is not correct, 
        ###the result won't be added to the history table.
        updateCheckboxInput(session, "invisible_a_power", value = FALSE)
      }
    }
    
    
    
    if(input$invisible_c_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for quadratic proximal effect and linear expected availability
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1
         && min(beta_quadratic__input()) >= 0 && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0
         && input$beta_quadratic__initial <= input$beta_quadratic__mean && input$tau_mean_c > 0 && min(linear_input()) >= 0
         && max(linear_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Quadratic")
          data_avail_power <<- c(data_avail_power, "Linear")
          ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c_power", value = FALSE)
      }
    }
    
    if(input$invisible_e_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for quadratic proximal effect and quadratic expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_quadratic__input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_quadratic__mean > 0 && input$beta_quadratic__mean < 1
         && min(beta_quadratic__input()) >= 0 && max(beta_quadratic__input()) <= 1 && input$beta_quadratic__initial >= 0
         && input$beta_quadratic__initial <= input$beta_quadratic__mean && input$tau_mean_e > 0 && min(quadratic_input()) >= 0
         && max(quadratic_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Quadratic")
          data_avail_power <<- c(data_avail_power, "Quadratic")
          ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=3, q=3);
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Quadratic")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_quadratic__mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e_power", value = FALSE)
      }
    }
    
    if(input$invisible_a1_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for constant proximal effect and constant expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_constant_mean > 0 && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Constant")
          data_avail_power <<- c(data_avail_power, "Constant")
          ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                
                
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_a1_power", value = FALSE)
      }
      
      
    }
    if(input$invisible_c1_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for constant proximal effect and linear expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_constant_mean > 0 && input$tau_mean_c > 0 
         && min(linear_input()) >= 0 && max(linear_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Constant")
          data_avail_power <<- c(data_avail_power, "Linear")
          ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c1_power", value = FALSE)
      }
    }
    
    if(input$invisible_e1_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for constant proximal effect and quadratic expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_constant_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && input$beta_constant_mean > 0 && input$tau_mean_e > 0
         && min(quadratic_input()) >= 0 && max(quadratic_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Constant")
          data_avail_power <<- c(data_avail_power, "Quadratic")
          ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=1, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Constant")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_constant_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e1_power", value = FALSE)
      }
    }
    if(input$invisible_a2_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for linear proximal effect and constant expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_a = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_a[(occ*k-occ+1):(occ*k)] = replicate(occ, constant_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1
         && input$beta_linear_mean > 0 && input$tau_mean_a > 0){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Linear")
          data_avail_power <<- c(data_avail_power, "Constant")
          ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                power <- PowerCalculation(size, b_input, input_a, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Constant")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_a)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_a2_power", value = FALSE)
      }
     
      
    }
    if(input$invisible_c2_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for linear proximal effect and linear expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_c = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_c[(occ*k-occ+1):(occ*k)] = replicate(occ, linear_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1
         && input$beta_linear_mean > 0 && input$tau_mean_c > 0 && min(linear_input()) >= 0 && max(linear_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Linear")
          data_avail_power <<- c(data_avail_power, "Linear")
          ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                power <- PowerCalculation(size, b_input, input_c, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Linear")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_c)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_c2_power", value = FALSE)
      }
      
    }
    if(input$invisible_e2_power)
    {### if the statement is true, we can add the result for calculating
      ### the power for linear proximal effect and quadratic expected availability
      
      days = input$days
      occ = input$occ_per_day
      Total = days*occ
      input_e = vector('numeric', length = Total)
      
      b_input = vector('numeric',Total)
      for(k in 1:days)
      {
        b_input[(occ*k-occ+1):(occ*k)] = replicate(occ, beta_linear_input()[k])
      }
      
      for(k in 1:input$days)
      {
        input_e[(occ*k-occ+1):(occ*k)] = replicate(occ, quadratic_input()[k])
      }
      
      size = input$sample_size
      if(input$days == round(input$days) && input$days > 0 && input$occ_per_day == round(input$occ_per_day) 
         && input$occ_per_day > 0 && input$sigLev >= 0 && input$sigLev <= 1 && input$sample_size > 0 
         && input$sample_size == round(input$sample_size) && min(beta_linear_input()) >= 0 && max(beta_linear_input()) <= 1
         && input$beta_linear_mean > 0 && input$tau_mean_e > 0 && min(quadratic_input()) >= 0 && max(quadratic_input()) <= 1){
        if(input$ranPro == "Constant" && input$P_intervene > 0 && input$P_intervene < 1)
        {
          delta <- input$P_intervene;
          power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
          
          if(power >= 0.5){
            HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, paste(round(power,3)*100, "%"))
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
          }
          else
          {
            HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
            data_power <<- c(data_power, "Less than 50%")
            data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
            data_ranPro_power <<- c(data_ranPro_power, input$P_intervene)
            updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
          }
          data_effect_power <<- c(data_effect_power, "Linear")
          data_avail_power <<- c(data_avail_power, "Quadratic")
          ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
          ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
          data_sig_power <<- c(data_sig_power, input$sigLev)
        }
        else if(input$ranPro == "Time-varying"){
          if(input$numbers == "re_days"){
            if(!is.null(input$file2)){
              delta <- as.vector(P_inter_days()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days){
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "Less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file2[[1]])
                  updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
            }
          }else{
            if(!is.null(input$file1)){
              delta <- as.vector(P_inter_dec()$Randomization.Probability)
              if(max(data) <= 1 && min(data) >= 0 && length(data) == input$days * input$occ_per_day){
                power <- PowerCalculation(size, b_input, input_e, delta, alpha0=input$sigLev, setup = list(days = days, occ.per.day = occ), p=2, q=3);
                if(power >= 0.5){
                  HTML(paste("<h4 style = 'color:blue';> The power is ", round(power,3)*100, "% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, paste(round(power,3)*100, "%"))
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
                }
                else
                {
                  HTML(paste("<h4 style = 'color:blue';> The power is less than 50% with sample size", size ,"when the significance level is",input$sigLev,"."))
                  data_power <<- c(data_power, "power is less than 50%")
                  data_rsize <<- c(data_rsize, as.integer(input$sample_size ))
                  data_ranPro_power <<- c(data_ranPro_power, input$file1[[1]])
                  updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
                }
                data_effect_power <<- c(data_effect_power, "Linear")
                data_avail_power <<- c(data_avail_power, "Quadratic")
                ave_effect_power <<- c(ave_effect_power, input$beta_linear_mean)
                ave_avail_power <<- c(ave_avail_power, input$tau_mean_e)
                data_sig_power <<- c(data_sig_power, input$sigLev)
              }else{
                updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
              }
            }else{
              updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
            }
          }
        }
        else{
          updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
        }
      }
      else{
        updateCheckboxInput(session, "invisible_e2_power", value = FALSE)
      }
      
    }
    ### Generate the data frame of power which contains all the vectors of information we have added through the above procedure ###
    data = data.frame("Power" = data_power, "Treatment Effect" = data_effect_power, "Average Effect" = ave_effect_power, "Availability" = data_avail_power, "Average Availability" = ave_avail_power, "Type I Error" = data_sig_power, 
                      "Sample Size" = data_rsize, "Randomization Probability" = data_ranPro_power )

  })
  

  ### Output the history table for calculating sample sizes ###
  output$result_table_size <- renderTable({
    
    data_get_size()
  })
  
  ## Output the history table for calculating powers ###
  output$result_table_power <- renderTable({
    
    data_get_power()
  })
  
  ### Download the history table of sample size as csv. or tsv. file ###
  output$downloadData_size <- downloadHandler(

    filename = function() {
      paste("Sample Size", input$filetype_size, sep = ".")
    },

    content = function(file) {
      sep <- switch(input$filetype_size,  "csv" = ",", "tsv" = "\t")
      
      # Write to a file specified by the 'file' argument
      write.table(data_get_size(), file, sep = sep,
                  row.names = FALSE)
    }
  )
  
  
  ### Download the history table of power as csv. or tsv. file ###
  output$downloadData_power <- downloadHandler(
    
    filename = function() {
      paste("Power", input$filetype_power, sep = ".")
    },
    
    content = function(file) {
      sep <- switch(input$filetype_power,  "csv" = ",", "tsv" = "\t")
      
      # Write to a file specified by the 'file' argument
      write.table(data_get_power(), file, sep = sep,
                  row.names = FALSE)
    }
  )

  
  

})
