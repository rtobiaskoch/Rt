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
  
#*******************************************************************************
#####ROLLING 7 DAY AVG FOR RT #####
#*******************************************************************************


daily_7<- var_merge %>%
  #7 day rolling avg for samples sequenced
  mutate(alpha_n7 =         rollmean(alpha_n, k = 7, fill = NA),
         gamma_n7 =         rollmean(gamma_n, k = 7, fill = NA),
         delta_n7 =         rollmean(delta_n, k = 7, fill = NA),
         other_nonVOC_n7 =  rollmean(other_nonVOC_n, k = 7, fill = NA),
         other_VOC_n7 =     rollmean(other_VOC_n, k = 7, fill = NA),
         n_7 =              rollmean(n, k = 7, fill = NA)) %>%
  #7 day rolling avg for frequency(prop aka proportion)
  mutate(alpha_prop7 =         rollmean(alpha_prop, k = 7, fill = NA),
         gamma_prop7 =         rollmean(gamma_prop, k = 7, fill = NA),
         delta_prop7 =         rollmean(delta_prop, k = 7, fill = NA),
         other_nonVOC_prop7 =  rollmean(other_nonVOC_prop, k = 7, fill = NA),
         other_VOC_prop7 =     rollmean(other_VOC_prop, k = 7, fill = NA))%>%
  #7 day rolling avg calculate by 7 day roll avg frequency pf variant * new cases
  mutate(alpha_cases7 = alpha_prop7 * New_Cases,
         gamma_cases7 = gamma_prop7 * New_Cases,
         delta_cases7 = delta_prop7 * New_Cases,
         other_nonVOC_cases7 = other_nonVOC_prop7 * New_Cases,
         other_VOC_cases7 = other_VOC_prop7 * New_Cases) %>%
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
#CI FUNCTION
#*#******************************************************************************
ci_fun <- function(v, nn, name, c){
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
  
  #selects columns for only the variant being run for simplicity and Rt function
  out4 = out3 %>% select(Date,
                         contains(name),
                         -ends_with("_est")) #drop est output from binom because it is identical to rolling_avg_<variant>
}

alpha_df = ci_fun(daily_7$alpha_n7, daily_7$n_7, "alpha")
gamma_df = ci_fun(daily_7$gamma_n7, daily_7$n_7, "gamma")
delta_df = ci_fun(daily_7$delta_n7, daily_7$n_7, "delta")
other_nonVOC_df = ci_fun(daily_7$other_nonVOC_n7, daily_7$n_7, "other_nonVOC")
other_VOC_df = ci_fun(daily_7$other_VOC_n7, daily_7$n_7, "other_VOC")


#*******************************************************************************
#RT FUNCTION
#*#******************************************************************************
  
#Rt Calculation
#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe

rt_fun= function(c, df, name){

  mean_Rt = estimate_R(c,
                       method="uncertain_si",
                       config = config)
  
   smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
   smooth_spline_mean_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)
   
   #renames smooth line Rt so it can merge and is more comprehensible
   smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                                  day =`smooth_spline_mean$x`,
                                  Rt = `smooth_spline_mean$y`)
   
   #merges the Rt value with the other variant data and renames Rt to have variant suffix
   merge = df %>%
     arrange(Date)%>%
     mutate(day = 1:nrow(df))%>%
     left_join(smooth_spline_mean_df) %>%
     rename_with(.fn = ~paste0(name,"_",.), .cols = Rt )
  
}

alpha_rt = rt_fun(daily_7$alpha_cases7, alpha_df, "alpha")
gamma_rt = rt_fun(daily_7$gamma_cases7, gamma_df, "gamma")
delta_rt = rt_fun(daily_7$delta_cases7, delta_df, "delta")
other_nonVOC_rt = rt_fun(daily_7$other_nonVOC_cases7, other_nonVOC_df, "other_nonVOC")
other_VOC_rt = rt_fun(daily_7$other_VOC_cases7, other_VOC_df, "other_VOC")

rt_list = list(alpha_rt,
               gamma_rt,
               delta_rt,
               other_nonVOC_rt,
               other_VOC_rt)

# WORKS UP UNTIL HERE
# bind_fun = function()
# 
# test = map_dfr(rt_list, select(Date, 
#                                ends_with("_Rt")))


#Change the name/ path to wherever you want it saved on your computer
setwd("Output")
dir.create(as.character(Sys.Date()))
setwd(as.character(Sys.Date()))

