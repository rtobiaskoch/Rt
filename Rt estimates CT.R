rm(list=ls(all=T))
graphic.off()

packages = c("readr", "tidyverse", "EpiEstim", "xlsx","stats","zoo","DescTools", "googledrive","googlesheets4")

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


# library(readr)
# library(tidyverse)
# library(EpiEstim)
# library(xlsx)
# library(stats)
# library(zoo)
# library(DescTools)


set.seed(1234)

#***********DATA IMPORT************************#
# Read in case data from JHU CSSE COVID-19 dataset
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
# Read in death data from JHU CSSE COVID-19 dataset
covid_deaths<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv")
#importing data from the googlesheet
Rt_import <- read_sheet("12xYePgxeF3pi0YiGnDCzmBnPPqZEASuobZ1DXeWZ7QA", sheet = "Rt")



# Create CT dataset with cases and deaths (new and cumulative)
long_cases<-gather(covid_cases, 
                   Date, 
                   Cumulative_Cases, 
                   '4/1/21':'7/12/21', 
                   factor_key=TRUE)
long_CT_cases<-long_cases[which(long_cases$Province_State=="Connecticut"),]
long_CT_cases<-subset(long_CT_cases,Admin2!="Out of CT"&Admin2!= "Unassigned")
long_CT_cases<-aggregate(x=long_CT_cases$Cumulative_Cases,by=list(long_CT_cases$Date),FUN=sum)
names(long_CT_cases)<-c("Date","Cumulative_Cases")
long_CT_cases<-long_CT_cases %>%
  mutate(New_Cases = Cumulative_Cases - lag(Cumulative_Cases))
long_deaths<-gather(covid_deaths, 
                    Date, 
                    Cumulative_Deaths,
                    '4/1/21':'7/12/21', 
                    factor_key=TRUE)
long_CT_deaths<-long_deaths[which(long_deaths$Province_State=="Connecticut"),]
long_CT_deaths<-subset(long_CT_deaths,Admin2!="Out of CT"&Admin2!= "Unassigned")
long_CT_deaths<-aggregate(x=long_CT_deaths$Cumulative_Deaths,by=list(long_CT_deaths$Date),FUN=sum)
names(long_CT_deaths)<-c("Date","Cumulative_Deaths")
long_CT_deaths<-long_CT_deaths%>%
  mutate(New_Deaths = Cumulative_Deaths - lag(Cumulative_Deaths))

covid_CT_wide<-cbind.data.frame(long_CT_cases$Date, long_CT_cases$Cumulative_Cases,
                                long_CT_cases$New_Cases,
                                long_CT_deaths$Cumulative_Deaths, long_CT_deaths$New_Deaths)
names(covid_CT_wide)<-c("Date", "Cumulative_Cases","New_Cases","Cumulative_Deaths","New_Deaths")
covid_CT_wide[is.na(covid_CT_wide)] <- 0
covid_CT_wide$New_Cases<-ifelse(covid_CT_wide$New_Cases<0,
                                0,
                                covid_CT_wide$New_Cases)
covid_CT_wide$New_Deaths<-ifelse(covid_CT_wide$New_Deaths<0,
                                 0,
                                 covid_CT_wide$New_Deaths)

#New Haven
covid_CT_wide$Date<-as.character(covid_CT_wide$Date)
covid_CT_wide$Date<-as.Date(covid_CT_wide$Date,"%m/%d/%y")

#imported at beginning
Gsheet_data = Rt_import %>%
  select(-c(`Percent check`) %>%
  rename(alpha_prop = `Freq Alpha`,
         delta_prop = `Freq Delta`,
         gamma_prop = `Freq Gamma`,
         other_VOC_prop = `Freq Other VOC/VOI`,
         other_nonVOC_prop = `Freq Non-VOC/VOI`) %>% #rename to match variables in rest of code
  filter(!is.na(`New_Cases`))#removes blank columns

week_int = sort(rep( 
                    seq(as.Date(min(Gsheet_data$Date))+6, 
                        as.Date(max(Gsheet_data$Date)), 
                        by = 7), #creates weekly interval dates
                    7) #repeates it 7 times for each day of the week
                ) #makes it in sequential order

week_int = week_int[1:nrow(Gsheet_data)] #filters to match the number of rows from data
  
Gsheet_data = Gsheet_data %>%
  mutate(week_int = week_int)

wide = Gsheet_data

B117_jeffreys<-B117_jeffreys
covid_CT_wide_daily<-covid_CT_wide

#change to the intervals in nextstrain
weeklycases<-rollapply(covid_CT_wide$New_Cases, width=7, FUN=sum, by=7)

daily_7<-Gsheet_data%>%
  mutate(alpha_7 = zoo::rollmean(alpha_cases, k = 7, fill = NA),
                               gamma_7 = zoo::rollmean(wide$gamma_cases, k = 7, fill = NA),
                               delta_7 = zoo::rollmean(wide$delta_cases, k = 7, fill = NA),
                               other_nonVOC_7 = zoo::rollmean(wide$other_nonVOC_cases, k = 7, fill = NA),
                               other_VOC_7 = zoo::rollmean(wide$other_VOC_cases, k = 7, fill = NA),
                               n_7 = zoo::rollmean(n, k = 7, fill = NA))
daily_7$rolling_avg_alpha<-daily_7$alpha_7/daily_7$n_7
daily_7$rolling_avg_gamma<-daily_7$gamma_7/daily_7$n_7
daily_7$rolling_avg_delta<-daily_7$delta_7/daily_7$n_7
daily_7$rolling_avg_other_nonVOC<-daily_7$other_nonVOC_7/daily_7$n_7
daily_7$rolling_avg_other_VOC<-daily_7$other_VOC_7/daily_7$n_7
rolling_avg<-cbind.data.frame(daily_7$rolling_avg_alpha,daily_7$rolling_avg_gamma,daily_7$rolling_avg_delta,
                              daily_7$rolling_avg_other_nonVOC,daily_7$rolling_avg_other_VOC)

alpha_rolling_avg_daily_cases<-daily_7$rolling_avg_alpha*daily_7$New_Cases
gamma_rolling_avg_daily_cases<-daily_7$rolling_avg_gamma*daily_7$New_Cases
delta_rolling_avg_daily_cases<-daily_7$rolling_avg_delta*daily_7$New_Cases
other_nonVOC_rolling_avg_daily_cases<-daily_7$rolling_avg_other_nonVOC*daily_7$New_Cases
other_VOC_rolling_avg_daily_cases<-daily_7$rolling_avg_other_VOC*daily_7$New_Cases
daily_cases_rolling_avg<-cbind.data.frame(alpha_rolling_avg_daily_cases,gamma_rolling_avg_daily_cases,
                                          delta_rolling_avg_daily_cases,other_nonVOC_rolling_avg_daily_cases,
                                          other_VOC_rolling_avg_daily_cases)

####ALPHA####
daily_alpha_confint<-BinomCI(x=daily_7$alpha_7[4:length(daily_7)], n=daily_7$n_7[4:length(daily_7)], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(as.Date("2021/4/1"), as.Date("2021/7/12"), by = "day")
alpha_jeffreys<-cbind.data.frame(day,daily_alpha_confint)
T <- 1 #start when there have been 12 cases
t_start <- seq(2, T-35) 
t_end <- t_start + 35
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))
alpha_jeffreys_cases_mean<-alpha_jeffreys$est*covid_CT_wide_daily$New_Cases
mean_Rt_alpha<-estimate_R(alpha_jeffreys_cases$alpha_jeffreys_cases_mean, 
                                  method="uncertain_si",
                                  config = config)
smooth_spline_alpha_mean<- with(mean_Rt_alpha$R, smooth.spline(mean_Rt_alpha$R$`t_end`, mean_Rt_alpha$R$`Mean(R)`, cv = TRUE))
smooth_spline_alpha_mean_df<-cbind.data.frame(smooth_spline_alpha_mean$x,smooth_spline_alpha_mean$y)


####GAMMA####
daily_gamma_confint<-BinomCI(x=daily_7$gamma_7[4:length(daily_7)], n=daily_7$n_7[4:length(daily_7)], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(as.Date("2021/4/1"), as.Date("2021/7/12"), by = "day")
gamma_jeffreys<-cbind.data.frame(day,daily_gamma_confint)
T <- 1 #start when there have been 12 cases
t_start <- seq(2, T-35) 
t_end <- t_start + 35
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))
gamma_jeffreys_cases_mean<-gamma_jeffreys$est*covid_CT_wide_daily$New_Cases
mean_Rt_gamma<-estimate_R(gamma_jeffreys_cases$gamma_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
smooth_spline_gamma_mean<- with(mean_Rt_gamma$R, smooth.spline(mean_Rt_gamma$R$`t_end`, mean_Rt_gamma$R$`Mean(R)`, cv = TRUE))
smooth_spline_gamma_mean_df<-cbind.data.frame(smooth_spline_gamma_mean$x,smooth_spline_gamma_mean$y)


####DELTA####
daily_delta_confint<-BinomCI(x=daily_7$delta_7[4:length(daily_7)], n=daily_7$n_7[4:length(daily_7)], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(as.Date("2021/4/1"), as.Date("2021/7/12"), by = "day")
delta_jeffreys<-cbind.data.frame(day,daily_delta_confint)
T <- 1 #start when there have been 12 cases
t_start <- seq(2, T-35) 
t_end <- t_start + 35
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))
delta_jeffreys_cases_mean<-delta_jeffreys$est*covid_CT_wide_daily$New_Cases
mean_Rt_delta<-estimate_R(delta_jeffreys_cases$delta_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
smooth_spline_delta_mean<- with(mean_Rt_delta$R, smooth.spline(mean_Rt_delta$R$`t_end`, mean_Rt_delta$R$`Mean(R)`, cv = TRUE))
smooth_spline_delta_mean_df<-cbind.data.frame(smooth_spline_delta_mean$x,smooth_spline_delta_mean$y)


####OTHER NON-VOC####
daily_other_nonVOC_confint<-BinomCI(x=daily_7$other_nonVOC_7[4:length(daily_7)], n=daily_7$n_7[4:length(daily_7)], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(as.Date("2021/4/1"), as.Date("2021/7/12"), by = "day")
other_nonVOC_jeffreys<-cbind.data.frame(day,daily_other_nonVOC_confint)
T <- 1 #start when there have been 12 cases
t_start <- seq(2, T-35) 
t_end <- t_start + 35
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))
other_nonVOC_jeffreys_cases_mean<-other_nonVOC_jeffreys$est*covid_CT_wide_daily$New_Cases
mean_Rt_other_nonVOC<-estimate_R(other_nonVOC_jeffreys_cases$other_nonVOC_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
smooth_spline_other_nonVOC_mean<- with(mean_Rt_other_nonVOC$R, smooth.spline(mean_Rt_other_nonVOC$R$`t_end`, mean_Rt_other_nonVOC$R$`Mean(R)`, cv = TRUE))
smooth_spline_other_nonVOC_mean_df<-cbind.data.frame(smooth_spline_other_nonVOC_mean$x,smooth_spline_other_nonVOC_mean$y)


####OTHER VOC####
daily_other_VOC_confint<-BinomCI(x=daily_7$other_VOC_7[4:length(daily_7)], n=daily_7$n_7[4:length(daily_7)], conf.level = 0.95, sides = "two.sided",method = "jeffreys")
day<-seq(as.Date("2021/4/1"), as.Date("2021/7/12"), by = "day")
other_VOC_jeffreys<-cbind.data.frame(day,daily_other_VOC_confint)
T <- 1 #start when there have been 12 cases
t_start <- seq(2, T-35) 
t_end <- t_start + 35
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end))
other_VOC_jeffreys_cases_mean<-other_VOC_jeffreys$est*covid_CT_wide_daily$New_Cases
mean_Rt_other_VOC<-estimate_R(other_VOC_jeffreys_cases$other_VOC_jeffreys_cases_mean, 
                          method="uncertain_si",
                          config = config)
smooth_spline_other_VOC_mean<- with(mean_Rt_other_VOC$R, smooth.spline(mean_Rt_other_VOC$R$`t_end`, mean_Rt_other_VOC$R$`Mean(R)`, cv = TRUE))
smooth_spline_other_VOC_mean_df<-cbind.data.frame(smooth_spline_other_VOC_mean$x,smooth_spline_other_VOC_mean$y)

