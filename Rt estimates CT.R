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
         VOC_prop = `Freq Other VOC/VOI`,
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
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")

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
         VOC_cases    = New_Cases*VOC_prop,
         nonVOC_cases = New_Cases*nonVOC_prop) %>%
  mutate(alpha_n        = n*alpha_prop,
         delta_n        = n*delta_prop,
         gamma_n        = n*gamma_prop,
         VOC_n    = n*VOC_prop,
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
         VOC_n7 =     rollmean(VOC_n, k = 7, fill = NA),
         n_7 =              rollmean(n, k = 7, fill = NA)) %>%
  #7 day rolling avg for frequency(prop aka proportion)
  mutate(alpha_prop7 =         rollmean(alpha_prop, k = 7, fill = NA),
         gamma_prop7 =         rollmean(gamma_prop, k = 7, fill = NA),
         delta_prop7 =         rollmean(delta_prop, k = 7, fill = NA),
         nonVOC_prop7 =  rollmean(nonVOC_prop, k = 7, fill = NA),
         VOC_prop7 =     rollmean(VOC_prop, k = 7, fill = NA))%>%
  #7 day rolling avg calculate by 7 day roll avg frequency pf variant * new cases
  mutate(alpha_cases7 = alpha_prop7 * New_Cases,
         gamma_cases7 = gamma_prop7 * New_Cases,
         delta_cases7 = delta_prop7 * New_Cases,
         nonVOC_cases7 = nonVOC_prop7 * New_Cases,
         VOC_cases7 = VOC_prop7 * New_Cases) %>%
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
              funs(.*New_Cases) #output is <variant_. couldnt figure out how to change
    )
  
}


alpha_df = ci_fun(daily_7$alpha_n7, daily_7$n_7, "alpha")
gamma_df = ci_fun(daily_7$gamma_n7, daily_7$n_7, "gamma")
delta_df = ci_fun(daily_7$delta_n7, daily_7$n_7, "delta")
nonVOC_df = ci_fun(daily_7$nonVOC_n7, daily_7$n_7, "nonVOC")
VOC_df = ci_fun(daily_7$VOC_n7, daily_7$n_7, "VOC")

#works till here


#*******************************************************************************
#RT FUNCTION
#*#******************************************************************************
  
#Rt Calculation
#generates Rt calculation, smooths the line then merges with 
#the other variant data in a dataframe

rt_fun= function(c7, df, name){
  
  non0 <- min(which(c7 > 0)) #finds the row of the first non zero value for <variant>_cases7
  t_beg = if(non0>2){
    non0-2
  }else{2}
  
  t_start<-seq(non0,length(c7)-21)
  t_end<-t_start+21
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
  
   mean_Rt = estimate_R(c7[non0:length(c7)],
                       method="uncertain_si",
                       config = config)
  
   smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
   smooth_spline_lower_ci<-with(mean_Rt_alpha$R, smooth.spline(mean_Rt_alpha$R$t_end, mean_Rt_alpha$R$`Quantile.0.025(R)`,cv=TRUE))
   smooth_spline_upper_ci<-with(mean_Rt_alpha$R, smooth.spline(mean_Rt_alpha$R$t_end, mean_Rt_alpha$R$`Quantile.0.975(R)`,cv=TRUE))
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

alpha_rt = rt_fun(alpha_df$alpha_est, alpha_df, "alpha")
gamma_rt = rt_fun(gamma_df$gamma_est, gamma_df, "gamma")
delta_rt = rt_fun(delta_df$delta_est, delta_df, "delta")
nonVOC_rt = rt_fun(nonVOC_df$nonVOC_est, nonVOC_df, "nonVOC")
VOC_rt = rt_fun(VOC_df$VOC_est, VOC_df, "VOC")

#me trying to extract the rts using an apply function
rt_list = list(alpha_rt,
               gamma_rt,
               delta_rt,
               nonVOC_rt,
               VOC_rt)

rt_names = c("Date", "alpha_Rt","gamma_Rt","delta_Rt",
             "nonVOC_Rt", "VOC_Rt")

#new
joined <- rt_list %>% 
  reduce(left_join, by = "Date")

#beginnings of an apply function
rt_export = cbind.data.frame(Date = rt_list[[1]][,"Date"],
                             alpha_rt = rt_list[[1]][,"alpha_Rt"],
                             gamma_rt = rt_list[[2]][,"gamma_Rt"],
                             delta_rt = rt_list[[3]][,"delta_Rt"],
                             nonVOC_rt = rt_list[[4]][,"nonVOC_Rt"],
                             VOC_rt = rt_list[[5]][,"VOC_Rt"]
                               )
rt_export = rt_export %>% drop_na()


#In your current working directory searches for a folder called Output and if it doesn't exist it will create one
ifelse(!dir.exists(file.path("Output")), 
       dir.create(file.path("Output")), 
       FALSE)
setwd("Output")

filename = paste(min(rt_export$Date), "to", max(rt_export$Date), "CT Rt.csv")

write.csv(rt_export, filename)




