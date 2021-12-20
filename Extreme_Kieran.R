library('tibble')
Project2021<-read.table(file="Project2021.txt", sep="")
d = Project2021$Autos
d1 = Project2021$Oil


plot(subset(Project2021,Oil< -1.5)$Oil)
sum(Project2021[1,2:30])


original <-read.table(file="original.txt", sep="",nrows = 25000)
original["date"] = row.names(original)

data = subset(original, date<=20151231 & date>=19500103)
a = data$Food

summary(Project2021)

