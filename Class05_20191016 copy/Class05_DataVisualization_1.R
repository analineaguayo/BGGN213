#' ---
#' title: "Class 5 Data Visualization"
#' output: github_document
#' ---


#Class 5 Data Visualization
#20191016 Wednesday
#Week 3
#to find out what a function does, you can put the function into the console with help before or a question mark then put the function in the parenthesese
#help(rnorm)
#?rnorm

x <- rnorm(1000)
#when you want to copy and past a script to the console, put your cursor on the line that you want to copy , and then click "command, return"
#you can also just click the RUN button up ahead
# "is.vector(x)" in the console will tell you if x is a vector or not with "TRUE" or "FALSE"
# "length(x)" will tell you how many things are in x
# if you just put "x" and then enter, all of the x values will pop up

#how many things are in x?
length(x)

mean(x)
sd(x)

summary(x)
#"summary(x)" will give you a summary of the data.. this will include: min, 1st quartile, median, mean, 3rd quartile, and max

#"mean(x) will give you the average, in this case, we had it randomized

# "boxblot(x)" will make a boxplot of our x
boxplot(x)

#"hist(x)" will give you a histogram of your x 
hist(x)

#"rug(x) will give you a distribution under your histogram of how many things are in each of the bars
#will show you where the actual points are
rug(x)

#From here down, I will be doing the In-Class Class 5 work from the website.
#This includes downloading the rscript from the website

#I have looked at the "Weight_chart" data from the uploaded "bimm143" folder
#data has been inputed with header columns and spaces between the columns

#look at the "read.table()" function using the help page to answer some of the questions
read.table(weight.chart.txt)
#oops that didn't work
help(read.table)
#this command has opened up the help topic so that I can see what kind of ways I can import the data so that r can "read" the data in a specific way
#the options are as a: table, csv, csv2, delim, delim2

#you can look up each of these options by inputting "help(read.csv)" etc.

#I am going to try to open this file with the "read.table" option, but you need to input the file name
read.table(file="bimm143_05_rstats/weight_chart.txt")
#you should be able to start putting in the subfolder of this project and then tab to help you out
#however, I was not able to do that, so I manually inputted it
#use "/" to open up the specific doc in the folder you've opened
#I opened up the data using the script above, however it is not exactly how it was when I opened the file as a preview
# in other words, there is a line with "V1" and "V2" at the top of the date (in the console), but that doesn't exist in the actual file
#this means that I have to change the way to input the data table
# we went back to the help for the "read.table" function to look at the header settings
#since there is a header in our original file, we need to right ",header=TRUE" to the end of our old script

#let's try
read.table(file="bimm143_05_rstats/weight_chart.txt",header=TRUE)

#cool, it worked! Now we do not have the "V1" and "V2" above our data in the console!
#it looks exactly the way it does in our original file (like in the preview tab)

#now let us make this data set called "weight"
weight <- read.table(file="bimm143_05_rstats/weight_chart.txt",header=TRUE)
#cool, that worked
#this means that you can input "weight" and then the data will pop up
#we have officially put in the data set to rscript!

#now to plot!

help("plot")
#I know nothing about the plot function so we should look at the help topic using the script up above
#following what the instructions say on the in-class asignment, I will input "plot(weight$Age, weight$Weight)
plot(weight$Age, weight$Weight)
#when I did this I got a scatter plot with points demarkated with a circle, there is no line connecting them
#that means that we need to add an argument after that script using "type=_" the options here are "l", "p", "b", "o"
#try them out to see what they all mean

plot(weight$Age, weight$Weight, type="l")
#this changed all of the points to a single line

plot(weight$Age, weight$Weight, type="p")
#this changed to just the dots again
#probably for POINTS

plot(weight$Age, weight$Weight, type="b")
#this resulted in POINTS with LINES inbetween

plot(weight$Age, weight$Weight, type="o")
#this put the line through the open points. this may mean OVERLAY both LINE and POINT

#now I am going to finagle with the point displaying chart to answer the rest of the questions
#what argument changes the point character to be a filled squared
#answer: "pch=15"
plot(weight$Age, weight$Weight, type="p", pch=15)
#it worked, the points have become squares
##you can try the other options: lwd=2, lwd=0.5, pch=4, pch=15, cex=1.5, cex=0.5
#idk what all of these do, but this should be something to try to do later!

#what argument changes the plot point size to be 1.5x normal size?
#answer: "cex=1.5"

#what argument changes the line width thickness to be twice the default size?
#answer: "lwd=2"

#what argument changes the y-axis limits scale between 2 and 10kg?
#answer: "ylim=c(2,10)"
#this question had other options: "main", "xlab", "ylab"
#still need to try out these other answers

plot(weight$Age, weight$Weight, typ="o", 
     pch=15, cex=1.5, lwd=2, ylim=c(2,10), 
     xlab="Age (months)", ylab="Weight (kg)", 
     main="Baby weight with age", col="blue")

plot(weight$Age, weight$Weight, typ="o", 
     pch=15, cex=1.5, lwd=2, ylim=c(2,10), 
     xlab="Age (months)", ylab="Weight (kg)", 
     main="Baby weight with age", col="blue") 




#2B Exercise

read.csv(file="bimm143_05_rstats/feature_counts.txt", header=TRUE, sep = "\t")
mouse <- read.csv(file="bimm143_05_rstats/feature_counts.txt", header=TRUE, sep = "\t")
mouse

par(mar=c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, names.arg=mouse$Feature, 
        horiz=TRUE, ylab="", 
        main="Number of features in the mouse GRCm38 genome", 
        las=1, xlim=c(0,80000))





