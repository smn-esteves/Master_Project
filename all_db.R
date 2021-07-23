library(readxl)
DadosTB_IDTC_efetivo_2019 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2019.xlsx")

DadosTB_IDTC_efetivo_2018 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2018.xlsx")

DadosTB_IDTC_efetivo_2017 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2017.xlsx")

DadosTB_IDTC_efetivo_2016 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2016.xlsx")

DadosTB_IDTC_efetivo_2015 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2015.xlsx")

DadosTB_IDTC_efetivo_2014 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2014.xlsx")

DadosTB_IDTC_efetivo_2013 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2013.xlsx")

DadosTB_IDTC_efetivo_2012 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2012.xlsx")

DadosTB_IDTC_efetivo_2011 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2011.xlsx")

DadosTB_IDTC_efetivo_2010 <- read_excel("D:/MESTRADO/tese/data & code/DadosTB_IDTC_efetivo_2010.xlsx")

db_bTB<- rbind(DadosTB_IDTC_efetivo_2017, DadosTB_IDTC_efetivo_2018, DadosTB_IDTC_efetivo_2019, DadosTB_IDTC_efetivo_2016, DadosTB_IDTC_efetivo_2015, DadosTB_IDTC_efetivo_2014, DadosTB_IDTC_efetivo_2013, DadosTB_IDTC_efetivo_2012, DadosTB_IDTC_efetivo_2011, DadosTB_IDTC_efetivo_2010)

db_bTB <- db_bTB[!db_bTB$DiCo=="0801",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0802",] 
db_bTB <- db_bTB[!db_bTB$DiCo=="0803",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0804",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0805",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0806",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0807",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0808",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0809",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0810",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0811",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0812",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0813",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0814",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0815",]
db_bTB <- db_bTB[!db_bTB$DiCo=="0816",]

#db_bTB <- db_bTB[!db_bTB$Total_number_animals==1,]

library(dplyr)
db_bTB_PS<- db_bTB %>% 
  group_by(time=Year) %>%
  summarise(I= sum(Positive))
  #mutate(PercentagePs = ps/sum(ps)*100)
db_bTB_PS$I <- round(db_bTB_PS$I)

library(ggthemes)
library(ggplot2)

db_bTB_PS%>% 
  ggplot(aes(x= time, y = I)) + 
  geom_point() + 
  labs(title= "Number of positive animals", x = "Year", y = "Number of positive animals")

x<-db_bTB%>% 
  group_by(time=Year, Herd) %>%
  summarise(prevalence = round(((Positive))/((Total_tests))*100)) #percentage of positive animals per year

y<- x%>%
  group_by(time) %>%
  summarise(prevalence1=round(mean(prevalence),digits=3))

mean(y$prevalence)

x1<-db_bTB%>% 
  group_by(time=Year, Herd) %>%
  summarise(testing_rate = round((sum(Total_tests))/(max(Total_number_animals))*100 )) #number of tests per year 

y1<- x1%>%
  group_by(time) %>%
  summarise(testing_rate=round(mean(testing_rate)))

mean(y1$testing_rate)

e<-db_bTB%>% 
  group_by(Year, Herd) %>%
  summarise(max(Total_number_animals)) 

i<-db_bTB%>% 
  group_by(Year, Herd) %>%
  summarise(round(mean(Total_number_animals)))

herdsize <- db_bTB %>% 
  group_by(Herd, Year) %>%
  summarise(herd_size=(max(Total_number_animals)))

mean_hs_year<- herdsize %>%
  group_by(Year) %>%
  summarise(herd_size=round(mean(herd_size)))


summary(herdsize)

ps<-db_bTB%>% 
  group_by(time=Year, Herd) %>%
  summarise(cullingrate=(Positive/Total_tests)*100) #number of positive/ killed animals per tests per year 
ps1<- ps%>%
  group_by(time) %>%
  summarise(culling_rate=round(mean(cullingrate),digits=3))

mean(ps1$culling_rate)

ps_herd<- db_bTB%>% 
  group_by(Year,Herd)%>%
  summarise(Total_number_animals, Positive)

plot(ps_herd$Total_number_animals, ps_herd$Positive)

i<- db_bTB %>% 
  group_by(Herd) %>%
  summarise(herd_size=round(mean(Total_number_animals)), Positive) %>%
  mutate(herd_size_range1 = case_when(herd_size <= 5 ~ '0-5',
                                      herd_size <= 10 ~ '6-10',
                                      herd_size <= 15 ~ '11-15',
                                      herd_size <= 20 ~ '16-20',
                                      herd_size <= 25 ~ '21-25',
                                      herd_size <= 30 ~ '26-30',
                                      herd_size <= 35 ~ '31-35',
                                      herd_size <= 40 ~ '36-40',
                                      herd_size <= 45 ~ '41-45',
                                      herd_size <= 50 ~ '46-50',
                                      herd_size <= 55 ~ '51-55',
                                      herd_size <= 60 ~ '56-60',
                                      herd_size <= 65 ~ '61-65',
                                      herd_size <= 70 ~ '66-70',
                                      herd_size <= 75 ~ '71-75',
                                      herd_size <= 80 ~ '76-80',
                                      herd_size <= 85 ~ '81-85',
                                      herd_size <= 90 ~ '86-90',
                                      herd_size <= 95 ~ '91-95',
                                      herd_size <= 100 ~ '96-100'))
i %>%
  ggplot(aes(x = herd_size_range1, y = Positive)) +
  geom_point()

t<-db_bTB%>% 
  filter(Year == 2010, Positive >= 1)

sum(t$Positive)

