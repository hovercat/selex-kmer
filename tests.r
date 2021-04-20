library(dplyr)
library(tidyr)
library(readr)

tags_ef07 <- read.table("../tagfile.csv", sep=";", header=TRUE)
r0_r2 <- read.table("scores_R0_R2.csv", sep=";", header=TRUE)

df <- list.files(path=".", full.names=TRUE, pattern="scores_*") %>%
  lapply(read_delim, delim=";") %>%
  bind_rows()

df <- read_delim("scores_R4_R5.csv", delim=";")

dff <- df %>%
  inner_join(tags_ef07)


View(df %>% 
       inner_join(tags_ef07) %>% 
       filter((name1 == "R0" & name2 =="R2") | (name1 == "R5" & name2 == "R6"))%>% 
       group_by(sequence, tag) %>% 
       summarize(score = sum(score)))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter((name1 == "R2" & name2 =="R3") | (name1 == "R6" & name2 == "R7"))%>% 
       group_by(sequence, tag) %>% 
       summarize(score = sum(score)))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter((name1 == "R3" & name2 =="R4") | (name1 == "R7" & name2 == "R8"))%>% 
       group_by(sequence, tag) %>% 
       summarize(score = sum(score)))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter((name1 == "R4" & name2 =="R5") | (name1 == "R8" & name2 == "R9"))%>% 
       group_by(sequence, tag) %>%
       summarize(score = sum(score)))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R4" & name2 =="R5"))
View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R8" & name2 =="R9"))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R0" & name2 =="R2"))


View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R4" & name2 =="R5"))
View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R8" & name2 =="R9"))


View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R0" & name2 =="R2"))
View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R5" & name2 =="R6"))



View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R2" & name2 =="R3"))
View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name1 == "R6" & name2 =="R7"))


View(df %>% 
       inner_join(tags_ef07) %>% 
       group_by(sequence, tag) %>% 
       summarize(score = max(score)))


View(df %>% 
       inner_join(tags_ef07) %>% 
       filter((name2 == "R8") | name2 == "R4")%>% 
       group_by(sequence, tag) %>% 
       summarize(score = max(score)))

View(df %>% 
       inner_join(tags_ef07) %>% 
       filter(name2 == "R5" | name2 == "R9")%>% 
       group_by(sequence, tag) %>% 
       summarize(score = mean(score)))
