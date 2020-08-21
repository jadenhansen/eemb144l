##### TITLE #####
# EEMB 144L LAB 1
# Nicholas Baetge
# 8/21/20


##### Intro to RStudio Environment #####

### A - Open RStudio
### B - The environment (windows, scripts, files, environment, etc.)
### C - Working in the console (why it's quick and dirty, why it's a nightmare for research and collaboration)

##### Projects #####

### Create a new project with File > New Project. 

### When you open the project, R will recognize that as your working directory. Notice that the pathway shows up in the top bar of RStudio, and any files that you save/add are automatically added to that working directory. A folder for the entire project is also created.

### Now we have a project created, and it's waiting for us to actually contribute something to it. We're going to work today in something called a SCRIPT. It's like a text editor (as opposed to the always-active console) and only runs lines of code when you ask it to. It also allows you to keep a clear record of the work you have done. 

##### Scripts #####
#  Good habits, organization and other things

### A - How scripts work (+ comments)
### B - How to run code in a script (+ shortcuts)
### C - Creating own script using great organization
### D - Create a couple of variables - note how they show up in the Environment window when stored (saving a script means they can quickly retrieve variables)
### E - Save the script (the importance of organized folders!!!)

### *** wrap text: Tools > Global Options > Code > Soft Wrap R Source Files > Apply ***

##### Loading the package 'tidyverse' #####

### A - What are packages? Why they don't all load automatically?
### B - Load package tidyverse and open

library(tidyverse)

##### Loading data from .csv or xlsx files #####

### Loading files into RStudio is EASY if you're working in a project. All you have to do is drag and drop the file into your project folder.

# Drag and drop into the project folder, and note that it appears in the 'Files' tab in RStudio automatically. That's because R knows that it's been added to our working directory. 

# Once it's in the Project, then it's easy to load using read_csv() or read_excel().

#install.packages("readxl")
library(readxl)

# unlike .csv files, .xlsx files can have multiple sheets
excel_sheets("Input_Data/2013_2019_CalFire_Redbook.xlsx") #let's see what the excel sheets are called

calfire.data <- read_excel("Input_Data/2013_2019_CalFire_Redbook.xlsx", sheet = "Data") # we can store each sheet separately as a dataframe! this one is the data

calfire.metadata <- read_excel("Input_Data/2013_2019_CalFire_Redbook.xlsx", sheet = "Metadata")  # this one is the metadata

##### Initial data exploration #####

## ALWAYS. EXPLORE. YOUR. DATA. ##

names(calfire.data) # Shows variable (column) names
dim(calfire.data) # Dimensions of dataset
class(calfire.data) # Data class
head(calfire.data) # Shows first 6 lines of dataset
tail(calfire.data) # Shows last 6 lines of dataset

# Want to know how a function works?

?names # Single question mark brings up R documentation for that function
??name # Will bring up every function that may contain this term

# Single columns can be referred to using a '$'
county <- calfire.data$County_Unit # A vector only containing information in column 'county' from calfire.data

max_acres <- max(calfire.data$Total_Acres_Burned, na.rm = T) # Finds the maximum value from the 'Total_Acres_Burned' column in calfire.data
max_str_des <- max(calfire.data$Structures_Destroyed)
max_str_des <- max(calfire.data$Structures_Destroyed, na.rm = T) # this is a good way to show what can happen  NA 'values' are present)

##### Basic data wrangling (dplyr functions) #####

df1 <- select(calfire.data, County_Unit:Controlled_Date, Total_Acres_Burned:Civil_Fatalities) # Only select columns from CountyUnit through Controlled_Date and Total_Acres_Burnedthrough Civil_Fatalities in calfire.data, store as a new data frame 'df1'
View(df1) # Remember to look at it!


df2 <- filter(df1, County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") &  Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") # Filters so that only rows are retained if the County_Unit column is 'SANTA BARBARA', 'VENTURA', 'LOS ANGELES' or 'SAN LUIS OBISPO') and the Total_Acres_Burned is 500 or greater OR if Fire_Name is 'THOMAS', stores as new subset data frame 'df2'
View(df2) # Always look!

df3 <- arrange(df2, desc(Start_Date), Total_Acres_Burned) # Arranges data first by dates in Start_Date column in descending order, THEN by ascending value in Total_Acres_Burned column
View(df3)

df4 <- mutate_at(df3, vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0 ) # we can modify multiple existing columns. Here we can replace all the NAs in the columns from Structures_Destroyed to Civil Fatalities with 0
View(df4)

df5 <- mutate(df4, Fatalities = Fire_Fatalities + Civil_Fatalities) # Adds a NEW column called 'Fatalities' to 'df3' that is the sum of the fatalities for fire service personnel and civilians
View(df5)

#mess with time! 
library(lubridate)

df6 <- mutate(df5, 
              interv = interval(Start_Date, Controlled_Date),
              dur = as.duration(interv),
              days = as.numeric(dur, "days"))
View(df6)

### We used 15 lines (excl. comments) to do all of that! And now we have 5 dataframes! This seems a little inefficient and too much to keep track of. No fear, there is a better way. It's called "piping", and it allows you to write more streamlined code for sequential operations on a data frame. 

##### Introduction to piping #####

# In the example above, we wrote a separate line of code for each operation to manipulate the CalFire data. That's good in some ways, but tedious in others. 

# What if we wanted to restrict our data to the southern CA coast, exclude fires that burned less that 500 acres, add a column that sums the human fatalities after changing NAs to 0s, add another column that computes the duration of each fire in days, arrange the data from most recent to most distant fires AND by total acreage burned, and then designate the Thomas Fire as a Ventura County fire? We could follow the same process as above, or we can use piping. 

# The MAGICAL pipe operator %>%  (command + shift + m on a mac or control + shift + m on windows)

# Think of this as a code way of saying "and then..."

socal_fires <- calfire.data %>% 
  filter(County_Unit %in% c("SANTA BARBARA", "VENTURA", "LOS ANGELES", "ORANGE", "SAN DIEGO") &  Total_Acres_Burned >= 500 | Fire_Name == "THOMAS") %>% 
  mutate_at(vars(Structures_Destroyed:Civil_Fatalities), replace_na, 0 ) %>% 
  mutate(Fatalities = Fire_Fatalities + Civil_Fatalities,
         interv = interval(Start_Date, Controlled_Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         County_Unit = ifelse(County_Unit == "VENTURA/SANTA BARBARA", "VENTURA", County_Unit)) %>% 
  arrange(desc(Start_Date), Total_Acres_Burned) 

view(socal_fires) #boom 10 lines, 1 dataframe
  

##### Our first graphs in ggplot #####

# We're going to make a graph of acres burned in the South Coast from 2013 - 2018, with the color dependent on which county we're showing. 

# Three things you must tell R to make a graph in ggplot: 
# (1) That you're using ggplot
# (2) What data you're using (including what should be x and what should be y)
# (3) What type of graph you want to create
# Everything after that is extra to make it beautiful

socal.plot <- socal_fires %>%
  rename(start = Start_Date, 
         acres = Total_Acres_Burned) %>% 
  ggplot(aes(x = start, y = acres)) +
  geom_point()

socal.plot 

# So that's kind of cool for a first shot, but it is definitely not a finalized graph. Let's start adding levels, and changing the line type, so we can get a graph that's a little nicer. 

socal.plot  <- socal_fires %>%
  rename(start = Start_Date, 
         acres = Total_Acres_Burned,
         county = County_Unit) %>% 
  ggplot(aes(x = start, y = acres)) + # Tells R what data to use
  geom_point(aes(colour = county)) + # Tells ggplot to make a scatterplot graph
  ggtitle("California South Coast Major Fires\n2013 - 2018") + # Gives it a title
  xlab("\nDate") + # Note that \n just adds a blank line before this label
  ylab("Acres\n") + # Again, but after the label
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) #These things all 'clean up' graphic area - remove gridlines, background, etc.

socal.plot # ok ok, better, but not super informative! what if we separated the counties out?

socal.plot + facet_grid(~county) # cool, the biggest from this is that man, that Thomas fire was HUGE!

#there are so many other ways we can explore this dataset...

#what if we plotted the total number of fires?

incidents <- socal_fires %>% 
  rename(county = County_Unit) %>% 
  mutate(year = year(Start_Date),
         county = ifelse(county == "VENTURA/SANTA BARBARA", "VENTURA", county)) %>% 
  group_by(county, year) %>%
  tally() %>% 
  ungroup() %>% 
  rename(incidents = n)


incidents.plot  <-  incidents %>%
  ggplot(aes(x = year, y = incidents)) + # Tells R what data to use
  geom_point(aes(colour = county)) + # Tells ggplot to make a scatterplot graph
  geom_line(aes(colour = county)) + #...then add a line
  ggtitle("California South Coast Major Fire Incidents\n2013 - 2018") + # Gives it a title
  xlab("\nDate") + # Note that \n just adds a blank line before this label
  ylab("Number of Incidents\n") + # Again, but after the label
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(rows = "county")

incidents.plot 


# what about southern ca as a whole?

all_incidents <- socal_fires %>% 
  mutate(year = year(Start_Date)) %>% 
  group_by(year) %>% # you can break your dataframe up into multiple units and any following functions you use will apply to those units
  tally() %>% #in this case we're just summing up the number of rows
  ungroup() %>% #always ungroup after your done!!
  rename(incidents = n)


all_incidents.plot  <-  all_incidents %>%
  ggplot(aes(x = year, y = incidents)) + # Tells R what data to use
  geom_point() + # Tells ggplot to make a scatterplot graph
  geom_line() + #...then add a line
  ggtitle("California South Coast Major Fire Incidents\n2013 - 2018") + # Gives it a title
  xlab("\nYear") + # Note that \n just adds a blank line before this label
  ylab("Number of Incidents\n") + # Again, but after the label
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

all_incidents.plot 


#do the total number of acres burned follow the pattern of incidents? you try!!
#hint: you'll need to play with the data a bit and use group_by(), ungroup(), mutate(), and sum()
