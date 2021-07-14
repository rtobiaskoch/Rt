library(readr)
library(tidyverse)
library(EpiEstim)
library(xlsx)
library(stats)
library(zoo)
library(DescTools)
set.seed(1234)
# Read in case data from JHU CSSE COVID-19 dataset
covid_cases<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv")
# Read in death data from JHU CSSE COVID-19 dataset
covid_deaths<-read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv")
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
  group_by(Date) %>%
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
  group_by(Date) %>%
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
end_of_week_date<-seq(as.Date("2021/4/7"), as.Date("2021/7/12"), by = "week") #change the second date to current date (last date of data)
weeklycases<-rollapply(covid_CT_wide$New_Cases, width=7, FUN=sum, by=7)
wide<-cbind.data.frame(end_of_week_date,weeklycases)


CT_Yale_Variant_results_Rt <- read_csv("CT-Yale Variant results - Rt.csv")

wide$alpha_prop<-c() #read in "alpha" proportions
wide$gamma_prop<-c() #read in "gamma" proportions
wide$delta_prop<-c() #read in "delta" proportions
wide$other_nonVOC_prop<-c() #read in "non-VOC/VOI" proportions 
wide$other_VOC_prop<-1-(wide$alpha_prop+wide$gamma_prop+wide$delta_prop+wide$other_nonVOC_prop)

wide$alpha_cases<-wide$weeklycases*wide$alpha_prop
wide$gamma_cases<-wide$weeklycases*wide$gamma_prop
wide$delta_cases<-wide$weeklycases*wide$delta_prop
wide$other_nonVOC_cases<-wide$weeklycases*wide$other_nonVOC_prop
wide$other_VOC_cases<-wide$weeklycases*wide$other_VOC_prop
wide$n<-wide$alpha_cases+wide$gamma_cases+wide$delta_cases+wide$other_nonVOC_cases+wide$other_VOC_cases

daily_7<-daily%>%dplyr::mutate(alpha_7 = zoo::rollmean(wide$alpha_cases, k = 7, fill = NA),
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

