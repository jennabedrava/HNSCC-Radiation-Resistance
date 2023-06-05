#########################################################################################################################
## Project: TRIP13 Radiation Resistance In Vivo
## Objectives: Investigate the effect of phosphorylation of Y56 on radiation resistance in vivo
## Date: Last updated 03/18/2022
## Author: Jenna Bedrava
#########################################################################################################################

library(tidyverse)
library(nlme)
library(sjmisc)
library(Hmisc)
library(predictmeans)
library(scales)
library(ggpubr)

#### CLEANING DATA FOR MODELS - EACH ROW IS ONE MOUSE ####
data<- read.csv("/Users/jennabedrava/Library/Mobile Documents/com~apple~CloudDocs/UMICH/BIOS 699/project2/Copy of MouseTum_Y56F_Clean.csv")

#creating group and tumor_position variables
data <- data %>% 
  mutate(Group = ifelse(substr(Mouse,9,12)=="TRIP",0,
                ifelse(substr(Mouse,9,12)=="Y56F",1,999)),
          Tumor_Position = ifelse(Tumor_Position== "Left (NoIR)","Left",
                           ifelse(Tumor_Position == "Right(IR)","Right","9999")))
#creating time variable
dates <- data %>% select(Date) %>% distinct(Date) %>%
  mutate(Time=row_number())
data <- merge(data, dates, by="Date")

#converting to wide format
data_wide <- reshape(as.data.frame(data),
                         idvar=c("Mouse","Time","Date","Group"),
                         v.names=c("Volume","Breadth","Length"),
                         timevar="Tumor_Position",
                         direction="wide")

#compute log of ratio - note that we scale by 0.01 in order to avoid problem of dividing by zero
data_wide$Volume_Ratio = (data_wide$Volume.Left+0.01)/(data_wide$Volume.Right+0.01)
data_wide$Log_Volume_Ratio = log(data_wide$Volume_Ratio)

#create ID variable that represents each mouse
ids <- data_wide %>% select(Mouse) %>% distinct(Mouse) %>%
  mutate(id=row_number())

data <- merge(data_wide, ids,by="Mouse") %>%
  arrange(id,Time)

#making day variable
data <- mutate(data, day = case_when(Time == 1 ~ 10, 
                                     Time == 2 ~ 12,
                                     Time == 3 ~ 14,
                                     Time == 4 ~ 17,
                                     Time == 5 ~ 18,
                                     Time == 6 ~ 19,
                                     Time == 7 ~ 20,
                                     Time == 8 ~ 21,
                                     Time == 9 ~ 24,
                                     Time == 10 ~ 25,
                                     Time == 11 ~ 26,
                                     Time == 12 ~ 27,
                                     Time == 13 ~ 28,
                                     Time == 14 ~ 31,
                                     Time == 15 ~ 33,
                                     Time == 16 ~ 35))

#adding knot variables to data, fixing Time variable
data <- data %>%
  mutate(Time=day-9)%>%
  mutate(Time_y56_knot = ifelse(Group==1 & Time>=15,Time-15,0), 
         Time_trip_knot = ifelse(Group==0 & Time>=8,Time-8,0))

#plotting log of ratio by group
y56.labs <- c("0" = "TRIP 13", "1" = "Y56")
ggplot(data = data) + aes(x = day, y = Log_Volume_Ratio) + 
  geom_line(aes(group = Mouse)) +
  geom_smooth(aes(group = 1), method="loess", formula=y~x) +
  facet_grid(. ~Group, labeller = labeller(Group = y56.labs)) + 
  labs(x = "Day", y = "Log Volume Ratio", 
       title = "Ratio of No Irradiation to Irradiation Tumor Growth Over Time",
       subtitle = "Each line represents one mouse.") +
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size =22), axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        plot.title = element_text(size=22, face = "bold", hjust=0.5),
        plot.subtitle = element_text(size=16, hjust=0.5),
        legend.text = element_text(size = 20), legend.title = element_text(size = 20))


#### CLEANING DATA FOR SPAGHETTI PLOT - EACH ROW IS HALF A MOUSE ####
data_half <- read.csv("/Users/jennabedrava/Library/Mobile Documents/com~apple~CloudDocs/UMICH/BIOS 699/cleaned_data.csv")

#converting data_half to long 
data_half <- pivot_longer(data_half, cols = X7.19.21:X8.13.21, names_to = "Date", values_to = "Volume")
data_half <- rename(data_half, Tumor_Position = Tumor.Position)

#adding day variable
data_half <- mutate(data_half, 
                    Day = rep(c(10, 12, 14, 17, 18, 19, 20, 21, 24, 25, 26, 27, 28, 31, 33, 35), 18))

#adding time variable (1 through 16)
data_half <- mutate(data_half, time = rep(1:16, 18))

#adding group variable
data_half <- mutate(data_half, group = ifelse(substr(Mouse,4,7)=="TRIP",0,
                                            ifelse(substr(Mouse,4,7)=="Y56F",1,999)))

#adding treatment variable
data_half <- mutate(data_half, treatment = ifelse(substr(Mouse,9,10)=="IR",1,
                                              ifelse(substr(Mouse,9,12)=="NoIR",0,999)))

#making spaghetti plot with both TRIP13 and Y56 group faceted
data_half <- mutate(data_half, Time = Day - 9)

y56.labs <- c("0" = "TRIP13", "1" = "Y56F")
ggplot(data = na.omit(data_half)) + aes(x = Time, y = Volume) + 
  geom_line(aes(group = Mouse, color = as.factor(group_x_treatment))) +
  facet_grid(. ~group, labeller = labeller(group = y56.labs)) + 
  labs(x = "Day", y = expression(Volume ~(mm^3)), title = "Tumor Growth Over Time", color = "Treatment x Group",
       subtitle = "Each line represents one half of a mouse.") +
  scale_color_manual(values = c("#FF7F50", "#619CFF", "red4", "darkblue"), 
                     labels = c("TRIP13-NoIR", "Y56F-NoIR", "TRIP13-IR", "Y56F-IR")) +
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size =12), axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        plot.subtitle = element_text(size=12, hjust=0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 12))

#making plot of 4 loess curves - one for each group*treatment combination
data_half <- mutate(data_half, group_x_treatment = ifelse(group == 0 & treatment == 0, 1,
                                                          ifelse(group == 1 & treatment == 0, 2,
                                                                ifelse(group == 0 & treatment == 1, 3,
                                                                       ifelse(group == 1 & treatment == 1, 4, NA)))))
#1 = TRIP, NO IR
#2 = Y56, NO IR
#3 = TRIP, IR
#4 = Y56, IR

ggplot(data = data_half) + aes(x = Time, y = Volume) + 
  geom_smooth(aes(color = as.factor(group_x_treatment)), method="loess") +
  labs(x = "Day", y = expression(Volume ~(mm^3)), title = "Tumor Growth Over Time", color = "Treatment x Group") +
  scale_color_manual(values = c("#FF7F50", "#619CFF", "red4", "darkblue"), 
                     labels = c("TRIP13-NoIR", "Y56F-NoIR", "TRIP13-IR", "Y56F-IR")) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "dashed")) +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13", y=850), colour="#FF7F50", angle=90, size=6) + 
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") +
  geom_text(aes(x=15, label="\nY56F", y=870), colour="#619CFF", angle=90, size=6) + 
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size =12), axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        plot.subtitle = element_text(size=12, hjust=0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 12))

#loess curves of log ratio 
all <- ggplot(data = data) + aes(x = Time, y = (Log_Volume_Ratio)) + 
  geom_smooth(aes(color = as.factor(Group)), method="loess") +
  #geom_segment(x = 14, xend = 24, y = (-0.956696 + 0.076481*14), yend = (-0.956696 + 0.076481*24), 
  #color = "#619CFF", linetype = "dashed") +
  labs(x = "Day", y = "Log (NoIR Volume / IR Volume)", title = "Log of Volume Ratio Over Time", color = "Group") +
  scale_color_manual(values = c("#FF7F50", "#619CFF"), 
                     labels = c("TRIP13", "Y56F" )) +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13", y=4.5), colour="#FF7F50", angle=90, size=6) + 
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") +
  geom_text(aes(x=15, label="\nY56F", y=4.5), colour="#619CFF", angle=90, size=6) + 
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size =12), axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        plot.subtitle = element_text(size=12, hjust=0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 12))


#loess curves of log ratio - excluding mouse 7 and 9
no7and9 <- ggplot(data = filter(data, id != 7, id != 9)) + aes(x = Time, y = (Log_Volume_Ratio)) + 
  geom_smooth(aes(color = as.factor(Group)), method="loess") +
  #geom_segment(x = 14, xend = 24, y = (-0.956696 + 0.076481*14), yend = (-0.956696 + 0.076481*24), 
               #color = "#619CFF", linetype = "dashed") +
  labs(x = "Day", y = "Log (NoIR Volume / IR Volume)", title = "Log of Volume Ratio Over Time \n Excluding Mouse 7 and 9", 
       color = "Group") +
  scale_color_manual(values = c("#FF7F50", "#619CFF"), 
                     labels = c("TRIP13", "Y56F" )) +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13", y=2.5), colour="#FF7F50", angle=90, size=6) + 
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") +
  geom_text(aes(x=15, label="\nY56F", y=2.5), colour="#619CFF", angle=90, size=6) + 
  geom_hline(yintercept = 0, color = "black") +
  theme_bw() + 
  theme(strip.background = element_rect(colour="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size =12), axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        plot.subtitle = element_text(size=12, hjust=0.5),
        legend.text = element_text(size = 10), legend.title = element_text(size = 12))

ggarrange(all, no7and9)

#### LINEAR MIXED MODELS ####
#linear mixed effects with random intercept only
model1 = lme(Log_Volume_Ratio ~ day + as.factor(Group) + day*as.factor(Group),  
             random = ~ 1 | id, 
             data = data, method = "ML", na.action = na.exclude)
summary(model1)

#linear mixed effects with random intercept and slope
model2 = lme(Log_Volume_Ratio ~ day + as.factor(Group) + day*as.factor(Group), 
             random = ~ day | id, 
             data = data, method = "ML", na.action = na.exclude)
summary(model2)
anova(model1, model2)

#linear mixed effects with spline terms (knots at time irradiation begins)
model3 = lme(Log_Volume_Ratio ~ Time + as.factor(Group) + Time*as.factor(Group) + 
               Time_trip_knot + Time_y56_knot, random = ~ Time | id, 
             data = data, method = "REML", na.action = na.exclude)
summary(model3) #FINAL MODEL
model3b <-update(model3,~.,corr=corCAR1(form=~ Time | id)) #can't specify autoregressive covariance structure because convergence

t.test(Time ~ Group, data = data) #significant at alpha = 0.10, not at alpha = 0.05

#linear mixed effects model without main group effect and without day*group interaction
model4 = lme(Log_Volume_Ratio ~ Time + Time_trip_knot + Time_y56_knot, 
             random = ~ Time | id, 
             data = data, method = "ML", na.action = na.exclude)
summary(model4) 

#comparing models
anova(model1, model3) 
anova(model2, model3) #model3 is best model
anova(model4, model3) #model3 is best model 

#making nice model output
sjPlot::tab_model(model3, 
                  show.intercept = FALSE,
                  show.re.var= TRUE,
                  pred.labels =c("Day", "Y56F Group", "TRIP13 Knot",
                                 "Y56F Knot", "Day * Y56F Group"),
                  dv.labels= "Log Volume Ratio")

#exponentiated parameters
exp(-0.2683625)
exp(0.0764815)
exp(2.4462231)
exp(-0.1099129)
exp(0.3414225)
exp(-0.3301048)

#exponentiated parameter confidence intervals
intervals(model3)
#intercept
exp(-1.7331991)
exp(1.19647413)
#time
exp(-0.1710969)
exp(0.32405988)
#group1
exp(-0.1987080)
exp(5.09115416)
#trip knot
exp(-0.3765677)
exp(0.15674190)
#y56 knot
exp(0.1269005)
exp(0.55594454)
#time x group1
exp(-0.6416366)
exp(-0.01857297)

#slope of y56 after radiation
0.0764815 + 0.3414225

#intercept and slope of trip13 before radiation
-0.2683625 + 2.4462231
0.0764815 - 0.3301048
#intercept and slope of trip13 after radiation
-0.2683625 + 2.4462231
0.0764815 - 0.1099129 - 0.3301048

exp(0.417904)
exp(-0.2536233)
exp(-0.3635362) 

#### MODEL DIAGNOSTICS PLOTS ####
#residual plots
predictmeans::residplot(model3)

#Pearson's residuals vs fitted values plot
resid_vs_fitted <- ggplot(data = data.frame(fitted_values = fitted(model3), 
                         residuals =  residuals(model3, type = c("pearson")),
                         group = data$Group, id = data$id)) + 
  geom_point(aes(x = fitted_values, y = residuals, color = as.factor(id))) + 
  labs(x = "Fitted Values", y = "Pearson Residuals", 
       title = "Residuals vs. Fitted Values", color = "Mouse") +
  scale_color_brewer(palette="Set1") +
  theme_bw() + 
  theme(axis.title=element_text(size=14),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#Pearson's residuals vs fitted values plot WITHOUT MOUSE 9
resid_vs_fitted_no9 <- ggplot(data = filter(data.frame(fitted_values = fitted(model3), 
                         residuals =  residuals(model3, type = c("pearson")),
                         group = data$Group, id = data$id), id != 9)) + 
  geom_point(aes(x = fitted_values, y = residuals, color = as.factor(id))) + 
  labs(x = "Fitted Values", y = "Pearson Residuals", 
       title = "Residuals vs. Fitted Values Excluding Mouse 9", color = "Mouse") +
  scale_color_brewer(palette="Set1") +
  theme_bw() + 
  theme(axis.title=element_text(size = 14),
        plot.title = element_text(size=14, face = "bold", hjust=0.5),
        legend.text = element_text(size = 12), legend.title = element_text(size = 12))

#combining residuals vs fitted plots
ggarrange(resid_vs_fitted, resid_vs_fitted_no9)

#checking normality of residuals
qq_resid <- ggplot(data.frame(x = residuals(model3, type = "pearson"))) +
  stat_qq(aes(sample = x), size=2) +
  stat_qq_line(aes(sample = x)) +
  theme_bw() +
  labs(title = "Q-Q Plot of Residuals") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab("Theoretical Quantiles") +
  theme(plot.title = element_text(size=14, face = "bold", hjust=0.5),
         text = element_text(size = 12),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14))

#checking normality of residuals WITHOUT MOUSE 9
model3_no9 = lme(Log_Volume_Ratio ~ Time + as.factor(Group) + Time*as.factor(Group) + 
               Time_trip_knot + Time_y56_knot, random = ~ Time | id, 
             data = filter(data, id != 9), method = "REML", na.action = na.exclude)
summary(model3_no9) #FINAL MODEL WITHOUT MOUSE 9

#making nice model output for model without mouse 9
sjPlot::tab_model(model3_no9, 
                  show.intercept = FALSE,
                  show.re.var= TRUE,
                  pred.labels =c("Day", "Y56F Group", "TRIP13 Knot",
                                 "Y56F Knot", "Day * Y56F Group"),
                  dv.labels= "Log Volume Ratio")

qq_resid_no9 <- ggplot(data.frame(x = residuals(model3_no9, type = "pearson"))) +
  stat_qq(aes(sample = x), size=2) +
  stat_qq_line(aes(sample = x)) +
  theme_bw() +
  labs(title = "Q-Q Plot of Residuals Excluding Mouse 9") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab("Theoretical Quantiles") +
  theme(plot.title = element_text(size=14, face = "bold", hjust=0.5),
        text = element_text(size = 12),
        axis.text.x= element_text(size=14),
        axis.text.y= element_text(size=14)) + ylim(-4, 4)

ggarrange(qq_resid, qq_resid_no9)

#qqplot of random intercepts
intercepts <- ggplot(data.frame(ranef(model3), 
                  Mouse = c("M1", "M2","M3","M4","M5","M6","M7", "M8","M9"))) +
  stat_qq(aes(sample = X.Intercept.), position = "identity", size=2) +
  stat_qq_line(aes(sample = X.Intercept.)) +
  theme_bw() +
  labs (title = "Q-Q Plot of Random Intercept") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab ("Theoretical Quantiles") +
  theme (plot.title = element_text (size = 14, face = "bold", hjust=0.5),
         text = element_text (size = 12),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14))

#qqplot of random slopes
slopes <- ggplot(data.frame(ranef(model3), 
                                Mouse = c("M1", "M2","M3","M4","M5","M6","M7", "M8","M9"))) +
  stat_qq(aes(sample = Time), size=2) +
  stat_qq_line(aes(sample = Time)) +
  theme_bw() +
  labs (title = "Q-Q Plot of Random Slope") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab ("Theoretical Quantiles") +
  theme (plot.title = element_text (size = 14, face = "bold", hjust=0.5),
         text = element_text (size = 12),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14))

#combining qqplots of random intercepts and random slopes
ggarrange(intercepts, slopes)

#qqplot of random intercepts WITHOUT MOUSE 9 
intercepts_no9 <- ggplot(data.frame(ranef(model3_no9), 
                                Mouse = c("M1", "M2","M3","M4","M5","M6","M7", "M8"))) +
  stat_qq(aes(sample = X.Intercept.), position = "identity", size=2) +
  stat_qq_line(aes(sample = X.Intercept.)) +
  theme_bw() +
  labs (title = "Q-Q Plot of Random Intercept \n Excluding Mouse 9") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab ("Theoretical Quantiles") +
  theme (plot.title = element_text (size = 14, face = "bold", hjust=0.5),
         text = element_text (size = 12),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14))

#qqplot of random slopes
slopes_no9 <- ggplot(data.frame(ranef(model3_no9), 
                            Mouse = c("M1", "M2","M3","M4","M5","M6","M7", "M8"))) +
  stat_qq(aes(sample = Time), size=2) +
  stat_qq_line(aes(sample = Time)) +
  theme_bw() +
  labs (title = "Q-Q Plot of Random Slope \n Excluding Mouse 9") +
  theme_bw() +
  ylab("Sample Quantiles") +
  xlab ("Theoretical Quantiles") +
  theme (plot.title = element_text (size = 14, face = "bold", hjust=0.5),
         text = element_text (size = 12),
         axis.text.x= element_text(size=14),
         axis.text.y= element_text(size=14))

#combining qqplots of random intercepts and random slopes
ggarrange(intercepts_no9, slopes_no9, widths=c(0.95,1))

#combining all qqplots of random intercepts and random slopes
ggarrange(intercepts, slopes, intercepts_no9, slopes_no9)

#BLUP plots
#loaded in mousetumor_predictions from class
blup1 <- ggplot(data = mousetumor_predictions %>% filter(id==1),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 1: BLUP vs. Observed") +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13 IR", y=0), colour="#FF7F50", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#FF7F50") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#FF7F50") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup2 <- ggplot(data = mousetumor_predictions %>% filter(id==2),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 2: BLUP vs. Observed") +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13 IR", y=0), colour="#FF7F50", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#FF7F50") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#FF7F50") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup3 <- ggplot(data = mousetumor_predictions %>% filter(id==3),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 3: BLUP vs. Observed") +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13 IR", y=-0.8), colour="#FF7F50", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#FF7F50") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#FF7F50") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup4 <- ggplot(data = mousetumor_predictions %>% filter(id==4),
                aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 4: BLUP vs. Observed") +
  geom_vline(xintercept = 8, linetype = "dashed", color = "#FF7F50") + 
  geom_text(aes(x=8, label="\nTRIP13 IR", y=0.75), colour="#FF7F50", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#FF7F50") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#FF7F50") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup5 <- ggplot(data = mousetumor_predictions %>% filter(id==5),
                aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) +
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 5: BLUP vs. Observed") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") + 
  geom_text(aes(x=15, label="\nY56F IR", y=2.4), colour="#619CFF", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#619CFF") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#619CFF") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup6 <- ggplot(data = mousetumor_predictions %>% filter(id==6),
                aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 6: BLUP vs. Observed") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") + 
  geom_text(aes(x=15, label="\nY56F IR", y=0.4), colour="#619CFF", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#619CFF") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#619CFF") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup7<-ggplot(data = mousetumor_predictions %>% filter(id==7),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 7: BLUP vs. Observed") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") + 
  geom_text(aes(x=15, label="\nY56F IR", y=5.75), colour="#619CFF", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#619CFF") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#619CFF") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup8<-ggplot(data = mousetumor_predictions %>% filter(id==8),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 8: BLUP vs. Observed") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") + 
  geom_text(aes(x=15, label="\nY56F IR", y=1.5), colour="#619CFF", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#619CFF") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#619CFF") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

blup9<-ggplot(data = mousetumor_predictions %>% filter(id==9),
              aes(x = Time, y = Prediction, group = as.character(id), color = as.character(id), label = as.character(id))) + 
  labs(x = "Day", y = expression(log(VolumeRatio ~(mm^3))), title = "Mouse 9: BLUP vs. Observed") +
  geom_vline(xintercept = 15, linetype = "dashed", color = "#619CFF") + 
  geom_text(aes(x=15, label="\nY56F IR", y=-1.2), colour="#619CFF", angle=90, text=element_text(size=10)) +
  theme_bw() + 
  geom_line(aes(y = Prediction), 
            size = 1.5,
            linetype = "solid", color="#619CFF") + 
  geom_line(aes(y = Log_Volume_Ratio), 
            size = 1.5,
            linetype = "dashed", color="#619CFF") + 
  scale_color_brewer(palette = "Dark2") + 
  guides(color = "none") + 
  theme(plot.title = element_text(size = 12, hjust=0.5))

ggarrange(blup1, blup2, blup3, blup4, blup5, blup6, blup7, blup8, blup9)
annotate_figure(ggarrange(blup1, blup2, blup3, blup4, blup5, blup6, blup7, blup8, blup9), 
                top = text_grob("Predictive vs. Observed Plots for Each Mouse", face = "bold", size = 14))
