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
                   "Cases Delta", "Cases Iota", "Cases Non-VOC/VOI")

#imported at beginning
var_data = var_import %>%
  select(-all_of(glab_drop_temp)
            ) %>%
  rename(alpha_prop = `Freq Alpha`,
         delta_prop = `Freq Delta`,
         gamma_prop = `Freq Gamma`,
         iota_prop = `Freq Iota`,
         nonVOC_prop = `Freq Non-VOC/VOI`,
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
# Create CT dataset with cases and deaths (new and cumulative)
long_cases<-covid_cases %>%
  filter(Province_State == "Connecticut" & 
         Admin2 != "Out of CT" & 
         Admin2 !="Unassigned") %>%   #filters only cases in CT
  gather(Date, 
          Cumulative_Cases, 
          "1/22/20":ncol(covid_cases)) %>%  #converts from wide format to long format 1/22/20 is the first day of covid cases in the US and shouldn't change
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
         iota_cases    = New_Cases*iota_prop,
         nonVOC_cases = New_Cases*nonVOC_prop) %>%
  mutate(alpha_n        = n*alpha_prop,
         delta_n        = n*delta_prop,
         gamma_n        = n*gamma_prop,
         iota_n    = n*iota_prop,
         nonVOC_n = n*nonVOC_prop)
  
#*******************************************************************************
#####ROLLING 7 DAY AVG FOR RT #####
#*******************************************************************************

daily_7<- var_merge %>%
  #7 day rolling avg for samples sequenced
  mutate(alpha_n7 =         rollmean(alpha_n, k = 7, fill = NA),
         gamma_n7 =         rollmean(gamma_n, k = 7, fill = NA),
         delta_n7 =         rollmean(delta_n, k = 7, fill = NA),
         nonVOC_n7 =  rollmean(nonVOC_n, k = 7, fill = NA),
         iota_n7 =     rollmean(iota_n, k = 7, fill = NA),
         n_7 =              rollmean(n, k = 7, fill = NA)) %>%
  #7 day rolling avg for frequency(prop aka proportion)
  mutate(alpha_prop7 =         rollmean(alpha_prop, k = 7, fill = NA),
         gamma_prop7 =         rollmean(gamma_prop, k = 7, fill = NA),
         delta_prop7 =         rollmean(delta_prop, k = 7, fill = NA),
         nonVOC_prop7 =  rollmean(nonVOC_prop, k = 7, fill = NA),
         iota_prop7 =     rollmean(iota_prop, k = 7, fill = NA))%>%
  #7 day rolling avg calculate by 7 day roll avg frequency pf variant * new cases
  mutate(alpha_cases7 = alpha_prop7 * New_Cases,
         gamma_cases7 = gamma_prop7 * New_Cases,
         delta_cases7 = delta_prop7 * New_Cases,
         nonVOC_cases7 = nonVOC_prop7 * New_Cases,
         iota_cases7 = iota_prop7 * New_Cases) %>%
  drop_na #drops future dates and first 3 days because of rollmean 
         

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
  cname = paste(name, colnames(out), sep = "_")
  colnames(out) = cname
  
  #binds binom output to daily_7
  out3 = cbind.data.frame(daily_7, out)
  
  #selects columns for only the variant being run for simplicity and Rt function
  out4 = out3 %>% select(Date,
                         New_Cases,
                         contains(name)
                         ) %>%
      mutate_at(vars(ends_with("est")), #searches for column with est
              funs(.*New_Cases)
               )%>% #output is <variant_. couldnt figure out how to change
     rename_at(vars(ends_with("est")),
               funs(paste("I"))
                    # use this if you need variant prefix funs(paste(name,"I",sep = "_"))
                )
    }


alpha_df = ci_fun(daily_7$alpha_n7, daily_7$n_7, "alpha")
gamma_df = ci_fun(daily_7$gamma_n7, daily_7$n_7, "gamma")
delta_df = ci_fun(daily_7$delta_n7, daily_7$n_7, "delta")
nonVOC_df = ci_fun(daily_7$nonVOC_n7, daily_7$n_7, "nonVOC")
iota_df = ci_fun(daily_7$iota_n7, daily_7$n_7, "iota")

#works till here


#*******************************************************************************
#RT FUNCTION ####
#*#******************************************************************************
  
#Rt Calculation
#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe

#run for everything -delta
rt_fun= function(df, name){
  
  non0 <- min(which(df$I > 0)) #1st day with cases of variant to start the R estimate otherise R estimate artificially high
  
  df2 = df[non0:nrow(df),] #dataframe filtered where there is the first case of variant to end of dataset
  
  #input of interval for R estimate
  t_start<-seq(2, nrow(df2)-21) 
  t_end<-t_start+21
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end)
  )
  
  
  mean_Rt = estimate_R(df2$I, #will search for column named I which was created in the ci_fun but explicitly named here
                       method="uncertain_si",
                       config = config)
  
  #adds back in days that were filtered out to match the days in the main dataframe
  mean_Rt$R$t_start = mean_Rt$R$t_start +non0 
  mean_Rt$R$t_end = mean_Rt$R$t_end + non0
  
  smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
  smooth_spline_lower_ci<-with(mean_Rt$R, smooth.spline(mean_Rt$R$t_end, mean_Rt$R$`Quantile.0.025(R)`,cv=TRUE))
  smooth_spline_upper_ci<-with(mean_Rt$R, smooth.spline(mean_Rt$R$t_end, mean_Rt$R$`Quantile.0.975(R)`,cv=TRUE))
  
  #binds them into a dataframe
  smooth_spline_mean_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y,
                                          smooth_spline_lower_ci$y,smooth_spline_upper_ci$y)
  
  # #renames smooth line Rt so it can merge and is more comprehensible
  smooth_spline_mean_df = rename(smooth_spline_mean_df,
                                 day = `smooth_spline_mean$x`,
                                 Rt = `smooth_spline_mean$y`,
                                 rtlowci = `smooth_spline_lower_ci$y`,
                                 rtupci = `smooth_spline_upper_ci$y`)
  
  #merges the Rt value with the other variant data and renames Rt to have variant suffix
  merge = df %>%
    arrange(Date)%>% #keep in date so that the day variable lines up with the first date
    mutate(day = 1:nrow(df))%>% #used to merge with the estimate_R variable output for the day
    left_join(smooth_spline_mean_df) %>%
    rename_with(.fn = ~paste0(name,"_",.), .cols = c("Rt", "rtlowci", "rtupci")) #renames the smooth_spline output to have variant prefix
}

alpha_rt = rt_fun(alpha_df, "alpha")
gamma_rt = rt_fun(gamma_df, "gamma")
nonVOC_rt = rt_fun(nonVOC_df, "nonVOC")
iota_rt = rt_fun(iota_df, "iota")
delta_rt = rt_fun(delta_df, "delta")



#*******************************************************************************
#RT MERGE FOR EXPORT ####
#*******************************************************************************
#me trying to extract the rt's in the join
rt_list = list(alpha_rt,
               gamma_rt,
               delta_rt,
               nonVOC_rt,
               iota_rt)

#new merged file that selects only the necessary variables
rt_export <- rt_list %>% 
  reduce(left_join, by = "Date") %>%
  select(Date,
         alpha_Rt, #alpha
         alpha_rtlowci,
         alpha_rtupci,
         gamma_Rt, #gamma
         gamma_rtlowci,
         gamma_rtupci,
         delta_Rt, #delta
         delta_rtlowci,
         delta_rtupci,
         iota_Rt, #iota
         iota_rtlowci,
         iota_rtupci,
         nonVOC_Rt, #nonVOC
         nonVOC_rtlowci,
         nonVOC_rtupci,
         )

#export directly back into the google sheet that we use "CT-Yale Variant results
write_sheet(rt_export, "12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "Rt_R_out")




