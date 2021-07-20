#*******************************************************************************
#####HEADER#####
#*#*******************************************************************************

#purpose: calculates Rt of variants in the state of connecticut by pulling case numbers from
#john hopkins and sequencing results from the Grubaugh Lab

#initial date: 2021/07/14
#author: Tobias Koch 
#email: r.tobiaskoch#gmail.com and toby.koch@yale.edu
#initial author Jessica Rothman

#PACKAGE INSTALL
packages = c("readr", "tidyverse", "lubridate","EpiEstim", "xlsx","stats","zoo",
             "DescTools", "googledrive","googlesheets4")



## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

rm(packages)

set.seed(1234)

#*******************************************************************************
##### GLAB DATA IMPORT#####
#*#*******************************************************************************

#importing rt data from CT-Yale Variant results googlesheet
var_import <- read_sheet("12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "Rt")

#*******************************************************************************
#####GLAB DATA CLEAN#####
#*#*******************************************************************************

glab_drop_temp = c("Percent check","New_Cases","Cases Gamma","Cases Alpha",
                   "Cases Delta", "Cases Other VOC/VOI", "Cases Non-VOC/VOI")

#imported at beginning
var_data = var_import %>%
  select(-all_of(glab_drop_temp)
            ) %>%
  rename(alpha_prop = `Freq Alpha`,
         delta_prop = `Freq Delta`,
         gamma_prop = `Freq Gamma`,
         other_VOC_prop = `Freq Other VOC/VOI`,
         other_nonVOC_prop = `Freq Non-VOC/VOI`,
         ) %>% #rename to match variables in rest of code
  filter(!is.na(Date)) %>% #removes blank columns
  mutate(Epiweek = paste(year(Date), "_EW", epiweek(Date), sep = "") #creates new epiweek column to match the format from nextstrain
         )

#*******************************************************************************
##### JHOP DATA IMPORT#####
#*#*******************************************************************************
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")

#*******************************************************************************
####JHOP CASE CLEAN #####
#*#*******************************************************************************
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")

# Create CT dataset with cases and deaths (new and cumulative)
long_cases<-covid_cases %>%
  filter(Province_State == "Connecticut" & 
         Admin2 != "Out of CT" & 
         Admin2 !="Unassigned") %>%   #filters only cases in CT
  gather(Date, 
          Cumulative_Cases, 
          "1/22/20":ncol(covid_cases)) %>%  #converts from wide format to long format
  mutate(Date = mdy(Date)) %>% #converts the dates from the columnns into R date format
  select(Date, Cumulative_Cases) %>% #remove columns that are no longer needed
  group_by(Date) %>% 
  summarize(Cumulative_Cases = sum(Cumulative_Cases)) %>%
  mutate(New_Cases = Cumulative_Cases - lag(Cumulative_Cases)) #adds new cases uses the lag function


#*******************************************************************************
#####MERGE DATA#####
#*#*******************************************************************************

var_merge = var_data %>%
  left_join(long_cases, by = "Date") %>% #merges jhop data with our data and keeps only 
  mutate(alpha_cases        = New_Cases*alpha_prop,
         delta_cases        = New_Cases*delta_prop,
         gamma_cases        = New_Cases*gamma_prop,
         other_VOC_cases    = New_Cases*other_VOC_prop,
         other_nonVOC_cases = New_Cases*other_nonVOC_prop) %>%
  mutate(alpha_n        = n*alpha_prop,
         delta_n        = n*delta_prop,
         gamma_n        = n*gamma_prop,
         other_VOC_n    = n*other_VOC_prop,
         other_nonVOC_n = n*other_nonVOC_prop)

#merges by epiweek instead of by the 
var_merge_wk = var_merge %>%
  group_by(Epiweek) %>%
  summarise(weeklycases = sum(New_Cases))
  
#*******************************************************************************
#####ROLLING 7 DAY AVG FOR RT #####
#*******************************************************************************
daily_7<- var_merge %>%
  mutate(alpha_7 =         rollmean(alpha_n, k = 7, fill = NA),
         gamma_7 =         rollmean(gamma_n, k = 7, fill = NA),
         delta_7 =         rollmean(delta_n, k = 7, fill = NA),
         other_nonVOC_7 =  rollmean(other_nonVOC_n, k = 7, fill = NA),
         other_VOC_7 =     rollmean(other_VOC_n, k = 7, fill = NA),
         n_7 =             rollmean(n, k = 7, fill = NA)) %>%
  mutate(rolling_avg_alpha = alpha_7/n_7,
         rolling_avg_gamma = gamma_7/n_7,
         rolling_avg_delta = delta_7/n_7,
         rolling_avg_other_nonVOC = other_nonVOC_7/n_7,
         rolling_avg_other_VOC = other_VOC_7/n_7) %>%
  mutate(alpha_rolling_avg_daily_case = rolling_avg_alpha * New_Cases,
         gamma_rolling_avg_daily_cases = rolling_avg_gamma * New_Cases,
         delta_rolling_avg_daily_cases = rolling_avg_delta * New_Cases,
         other_nonVOC_rolling_avg_daily_cases = rolling_avg_other_nonVOC* New_Cases,
         other_VOC_rolling_avg_daily_cases = rolling_avg_other_VOC * New_Cases) %>%
  drop_na
         
#*******************************************************************************
# #arguments for the function in estimate R
#*******************************************************************************
t_start<-seq(2,nrow(daily_7)-21)
t_end<-t_start+21
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))


#*******************************************************************************
#RT FUNCTION
#*#*******************************************************************************
rt_fun <- function(v, nn, name, c){
  out = BinomCI(x=v, 
          n=nn, 
          conf.level = 0.95, 
          sides = "two.sided",
          method = "jeffreys")
  
  #adds alpha prefix to output names
  out2 = paste(name, colnames(out), sep = "_")
  colnames(out) = out2
  
  #binds binom output to daily_7
  out3 = cbind.data.frame(daily_7, out)
  
  #selects columns for only the variant being run for simplicity
  out4 = out3 %>% select(Date,
                         contains(name),
                         -ends_with("_est")) #drop est output from binom because it is identical to rolling_avg_<variant>
  
  # #Rt Calculation
  mean_Rt <-estimate_R(c,
                       method="uncertain_si",
                       config = config)
  
  smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
  smooth_spline_mean_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)
  
}

temp = rt_fun(daily_7$alpha_7, daily_7$n_7, "alpha", daily_7$alpha_rolling_avg_daily_case)




#Change the name/ path to wherever you want it saved on your computer
setwd("Output")
dir.create(as.character(Sys.Date()))
setwd(as.character(Sys.Date()))

write.xlsx(smooth_spline_alpha_mean_df,"Rt_&.15.21.xlsx",sheetName="alpha",append=TRUE)
write.xlsx(smooth_spline_gamma_mean_df,"~/Documents/Yale/Grubaugh_lab/Rt_7.15.21.xlsx",sheetName="gamma",append=TRUE)
write.xlsx(smooth_spline_delta_mean_df,"~/Documents/Yale/Grubaugh_lab/Rt_7.15.21.xlsx",sheetName="delta",append=TRUE)
write.xlsx(smooth_spline_other_nonVOC_mean_df,"~/Documents/Yale/Grubaugh_lab/Rt_7.15.21.xlsx",sheetName="other_nonVOC",append=TRUE)
write.xlsx(smooth_spline_other_VOC_mean_df,"~/Documents/Yale/Grubaugh_lab/Rt_7.15.21.xlsx",sheetName="other_VOC",append=TRUE)