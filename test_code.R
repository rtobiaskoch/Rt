
non0 <- min(which(delta_df$delta_est > 0))

t_beg = if(non0>2){
  non0-2
}else{2}

t_start<-seq(2,length(delta_df$delta_est)-18)
t_end<-t_start+18
config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                           std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                           n1=500,n2=50,t_start=t_start, t_end=t_end)
)

case7_R = delta_df$delta_est[t_beg:length(delta_df$delta_est)]                  

mean_Rt = estimate_R(case7_R,
                     method="uncertain_si",
                     config = config)

smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
smooth_spline_mean_alpha_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)

#renames smooth line Rt so it can merge and is more comprehensible
smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                               day =`smooth_spline_mean$x`,
                               Rt = `smooth_spline_mean$y`)

#merges the Rt value with the other variant data and renames Rt to have variant suffix
merge = delta_df %>%
  arrange(Date)%>%
  mutate(day = 1:nrow(delta_df))%>%
  left_join(smooth_spline_mean_df) %>%
  rename_with(.fn = ~paste0("alpha","_",.), .cols = Rt )


#*************************************************************************************************************************


  
  non0 <- min(which(delta_df$delta_cases7 > 0))
  
  t_beg = if(non0>2){
    non0
  }else{2}
  
  t_start<-seq(t_beg,length(delta_df$delta_cases7)-19)
  t_end<-t_start+21
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
                        
              case7_R = delta_df$delta_cases7[t_beg:length(delta_df$delta_cases7)]                  
                        
                        mean_Rt = estimate_R(case7_R,
                                             method="uncertain_si",
                                             config = config)
                        
                        smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
                        smooth_spline_mean_alpha_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)
                        
                        #renames smooth line Rt so it can merge and is more comprehensible
                        smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                                                       day =`smooth_spline_mean$x`,
                                                       Rt = `smooth_spline_mean$y`)
                        
                        #merges the Rt value with the other variant data and renames Rt to have variant suffix
                        merge = delta_df %>%
                          arrange(Date)%>%
                          mutate(day = 1:nrow(delta_df))%>%
                          left_join(smooth_spline_mean_df) %>%
                          rename_with(.fn = ~paste0("alpha","_",.), .cols = Rt )


#*************************************************************************************************************************
                        
                        non0 <- min(which(daily_7$alpha_cases7 > 0))
                        
                        t_beg = if(non0>2){
                          non0
                        }else{2}
                        
                        t_start<-seq(2,length(daily_7$alpha_cases7)-21)
                        t_end<-t_start+21
                        config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                                                   std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                                                   n1=500,n2=50,t_start=t_start, t_end=t_end)
                        )
                        
                        case7_R = daily_7$alpha_cases7[t_beg:length(daily_7$alpha_cases7)]                  
                        
                        mean_Rt = estimate_R(case7_R,
                                             method="uncertain_si",
                                             config = config)
                        
                        smooth_spline_mean<- with(mean_Rt$R, smooth.spline(mean_Rt$R$`t_end`, mean_Rt$R$`Mean(R)`, cv = TRUE))
                        smooth_spline_mean_alpha_df<-cbind.data.frame(smooth_spline_mean$x,smooth_spline_mean$y)
                        
                        #renames smooth line Rt so it can merge and is more comprehensible
                        smooth_spline_mean_df = rename(smooth_spline_mean_df, 
                                                       day =`smooth_spline_mean$x`,
                                                       Rt = `smooth_spline_mean$y`)
                        
                        #merges the Rt value with the other variant data and renames Rt to have variant suffix
                        merge = alpha_df %>%
                          arrange(Date)%>%
                          mutate(day = 1:nrow(alpha_df))%>%
                          left_join(smooth_spline_mean_df) %>%
                          rename_with(.fn = ~paste0("alpha","_",.), .cols = Rt )
                        