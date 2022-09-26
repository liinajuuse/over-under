# R SCRIPT - HANDLING EEG DATA
# JUUSE & FOERSTER 2022
# stenio.foerster@ut.ee
# Last update: 26 September 2022 by Liina Juuse
# modifying script to fit over-under data

# Libraries and data ------------------------------------------------------

library(tidyverse)

#setwd("~/GitHub/over-under/over-under")
#folder <- c("input","output","figures", "data", "scripts") #create subfolders
#for (j in seq_along(folder)){
#  if (file.exists(folder[j])) {
#    cat("The folder already exists \n")
#  } else {
#    dir.create(folder[j])}}

setwd("~/GitHub/over-under/over-under/input")

rw <- list.files(pattern = '\\.dat')

# File inspection ---------------------------------------------------------

# Here, we are going to check if the dat files have the expected dimensions 
# in terms of number of rows and columns. In this case, a correct dat file
# must have 22116 rows and 64 columns.

SamplesTrial <- 3686 # enter the number of samples per trial // rows in file divided by nr of trials

level_key = c("019" = "Anger", "115" = "Anger", 
              "031" = "Happy", "127" = "Happy",
              "015" = "Sad", "111" = "Sad",
              "027" = "Neutral", "123" = "Neutral",
              "007" = "Disgust", "103" = "Disgust",
              "023" = "Fear", "119" = "Fear",
              "011" = "Surprise", "107" = "Surprise")

data_inspec <- data.frame()

for (i in 1:length(rw)) {
  
  nr <- nrow(read.table(file = rw[i], header = T, sep = ';'))
  nc <- ncol(read.table(file = rw[i], header = T, sep = ';'))
  di <- (nr * nc)
  tr <- (nr/SamplesTrial)
  fo <- ifelse(di == 1161090, 'Correct', 'Incorrect') #change 'di' values based on rows x columns value
  data_inspec <- rbind(data_inspec, data.frame('File' = rw[i], 'Rows' = nr, 'Columns' = nc, 'Cells' = di, 'Trials' = tr, 'Format' = fo))
  print(paste(rw[i], 'Done!', sep = ' '))
}

#setwd("~/GitHub/over-under/over-under/output")

#write.csv(data_inspec, 'data_inspect.csv', row.names = F)

fi <- data_inspec[data_inspec$Format == 'Correct', ]$File # dat files with correct dimensions. This object will be used in the following steps

fa <- data_inspec[data_inspec$Format == 'Incorrect', ]$File # dat files with incorrect dimensions. These need to be analysed separately


# Two STDEV rule ----------------------------------------------------------

nTrials <- 6 # enter the number of trials

data_2stdev <- data.frame()
miss_data_2stdev <- data.frame()

t0 <- Sys.time()

for (j in 1:length(fi)) {
  
  raw.data <- read.table(file = fi[j], header = T, sep = ';')
  
  # Attributing some info
  nSamplesTrial <- (nrow(raw.data)/nTrials) # this calculates the number of samples in each trial based on the number of rows present in the raw data
  vTrials <- unlist(lapply((1:nTrials), FUN = rep, nSamplesTrial)) # this will generate a vector containing the trials. Each trial is repeated "nSampleTrial" times. The vector length is equal to the number of rows in the raw data
  ElectNames <- colnames(raw.data) # store the electrode names. This will be used to assembly the final table
  
  data_01 <- data.frame('Trial' = vTrials, raw.data) # raw data with assigned trials
  
  # Step 1: Wide table
  # In this step, we will generate a wide table that will be used to average the 
  # electrodes respecting the row order in each trial. This is just an intermediary 
  # step, you don't need to pay attention to the structure of the 
  # objects produced in this step.
  
  u <- unique(vTrials) # a simple vector containing the trial ids
  d <- matrix(nrow = nSamplesTrial) # this is just a backbone matrix used to concatenate the following matrices in a single object (don't pay attention to this object)
  
  for (i in 1:length(u)) {
    
    a <- data_01[data_01$Trial == u[i], ]
    
    d <- cbind(d, a)
  }
  
  # Step 2: Averaging trials
  # Averaging the trials, respecting the row order. For example, if there are two trials with two rows each, 
  # then the code below will calculate the the average between: 1) the values present in the first cell of
  # each trial; and 2) the values present in the second cell of each trial. 
  
  v <- paste(paste('\\b', ElectNames, sep = ''), '\\b', sep = '') # vector with electrode names in a regular expression format. Don't change this. This will make sure that the electrode names will be correctly called during the loop.
  f <- matrix(nrow = nSamplesTrial)
  u <- unique(vTrials) # a simple vector containing the trial ids
  
  for (i in 1:length(v)) {
    
    a <- d[, grep(colnames(d), pattern = v[i], ignore.case = F)]
    s <- (apply(a, MARGIN = 2, FUN = sd) * 2) + abs(apply(a, MARGIN = 2, FUN = mean)) # two standard deviations of the mean
    
    for (k in 1:ncol(a)) {
      
      a[, k][a[, k] > s[k]] <- NA # assigning NA to values above the mean plus 2 stdev
      a[, k][a[, k] < -s[k]] <- NA # assigning NA to values below (more negatives than) the mean plus 2 stdev
      a[, k][a[, k] == 0] <- NA # assigning NA to zero values
      
      # storing missing data
      m <- length(a[is.na(a)]) # amount of NAs attributed to each electrode for each participant across the trials
      p <- substr(fi[j], 6, 8)
      e <- ElectNames[i]
      g <- fi[j]
      
      # missing data table
      miss_data_2stdev <- rbind(miss_data_2stdev, data.frame('Electrode' = e, 'Participant' = p,  'File' = g, 'Trial' = u[k], 'Count_data' = nrow(a)*ncol(a), 'Missing_count' = m, 'Missing_perct' = round(m/(nrow(a)*ncol(a))*100, digits = 4)))
    }
    
    # Electrode average table (average among trials, respecting the row order, of course)
    b <- data.frame(apply(a, MARGIN = 1, FUN = mean, na.rm = T))
    colnames(b) <- ElectNames[i]
    f <- cbind(f, b)
  }
  
  data_avgTrials <- f[, -1] # this is the table with the electrodes averaged by trials (averaged respecting the row order)
  #head(data_avgTrials)
  #dim(data_avgTrials)
  
  
  # Step 3: Flip the table
  
  data_02 <- data_avgTrials[1:615, ] # subset rows (-199 to 1000)
  
  data_02 <- as.data.frame(t(data_02)) # flip the table
  
  column.names <- paste(rep('T_', ncol(data_02)), seq(1, ncol(data_02), by = 1), sep = '')
  
  colnames(data_02) <- column.names
  
  data_02 <- data.frame('Electrode' = rownames(data_02), data_02)
  
  rownames(data_02) <- NULL
  
  # Step 4: Adding metadata
  
  
  # Extracting the meta info from the name of the file and add the info into the data frame
  #tmp$stdv <- substr(fi[j],nchar(fi[j])-8,nchar(fi[j])-7)
  data_02$Subj <- substr(fi[j], 6, 8)
  data_02$Group <- substr(fi[j], 12, 13)
  data_02$Age <- substr(fi[j], 10, 11)
  data_02$Sex <- substr(fi[j],9 ,9)
  data_02$Marker <- substr(fi[j], nchar(fi[j])-6, nchar(fi[j])-4)
  data_02$Condition <- ifelse(data_02$Marker %in% c("007", "011", "015", "019", "023", "027", "031"), "Ekman", "Under")
  data_02$Emotion <- data_02$Marker
  data_02$Emotion <- recode(data_02$Emotion, !!!level_key)
  #tmp$order <- substr(fi[j], 12, 12)
  #tmp$series <- substr(fi[j], 14, nchar(fi[j])-29)
  #tmp$eventType <- substr(fi[j], nchar(fi[j])-8, nchar(fi[j])-4)
  
  data_02 <- data_02[, c('Electrode', 'Subj', 'Group', 'Age', 'Sex', 'Marker', 'Condition', 'Emotion', column.names)]
  
  data_2stdev <- rbind(data_2stdev, data_02)
  print(paste(fi[j], 'Done!', sep = ' '))
  
}

(Sys.time()-t0)

# Plotting the missing data
miss_data_2stdev %>% group_by(Electrode) %>% summarise('elec_avgper' = mean(Missing_perct, na.rm = T)) -> md_2stdev

p <- ggplot(data = md_2stdev, mapping = aes(x = reorder(Electrode, elec_avgper), y = elec_avgper)) + geom_bar(stat = 'identity', width = 0.7, fill = '#85d0ff') + coord_flip()
p <- p + geom_hline(yintercept = mean(md_2stdev$elec_avgper), linetype = 'dashed', size = 0.3, colour = 'red') # mean percentage of missing data based on the electrodes displayed on the graph
p <- p + xlab('Electrode') + ylab('Mean percentage of missing data') + ggtitle(label = 'Electrodes with missing data')
p <- p + theme_bw()
p <- p + theme(axis.title = element_text(face = 'bold', size = 12)
               ,axis.text = element_text(size = 9, colour = 'black')
               ,title = element_text(face = 'bold', size = 12)
               ,panel.grid = element_line(size = 0.3)); p


# Threshold rule ----------------------------------------------------------

# Applying the threshold rule to clean the data

# fi <- list.files(pattern = '\\.dat')

Threshold <- 100 # enter the threshold (micro volts) used to clean the data

nTrials <- 6 # enter the number of trials

data_threshold <- data.frame()
miss_data_thresh <- data.frame()

t0 <- Sys.time()

for (j in 1:length(fi)) {
  
  raw.data <- read.table(file = fi[j], header = T, sep = ';')
  
  # Attributing some info
  nSamplesTrial <- (nrow(raw.data)/nTrials) # this calculates the number of samples in each trial based on the number of rows present in the raw data
  vTrials <- unlist(lapply((1:nTrials), FUN = rep, nSamplesTrial)) # this will generate a vector containing the trials. Each trial is repeated "nSampleTrial" times. The vector length is equal to the number of rows in the raw data
  ElectNames <- colnames(raw.data) # store the electrode names. This will be used to assembly the final table
  
  data_01 <- data.frame('Trial' = vTrials, raw.data) # raw data with assigned trials
  
  # Step 1: Wide table
  # In this step, we will generate a wide table that will be used to average the 
  # electrodes respecting the row order in each trial. This is just an intermediary 
  # step, you don't need to pay attention to the structure of the 
  # objects produced in this step.
  
  u <- unique(vTrials) # a simple vector containing the trial ids
  d <- matrix(nrow = nSamplesTrial) # this is just a backbone matrix used to concatenate the following matrices in a single object (don't pay attention to this object)
  
  for (i in 1:length(u)) {
    
    a <- data_01[data_01$Trial == u[i], ]
    
    d <- cbind(d, a)
  }
  
  # Step 2: Averaging trials
  # Averaging the trials, respecting the row order. For example, if there are two trials with two rows each, 
  # then the code below will calculate the the average between: 1) the values present in the first cell of
  # each trial; and 2) the values present in the second cell of each trial. 
  
  v <- paste(paste('\\b', ElectNames, sep = ''), '\\b', sep = '') # vector with electrode names in a regular expression format. Don't change this. This will make sure that the electrode names will be correctly called during the loop.
  f <- matrix(nrow = nSamplesTrial)
  u <- unique(vTrials) # a simple vector containing the trial ids
  
  for (i in 1:length(v)) {
    
    a <- d[, grep(colnames(d), pattern = v[i], ignore.case = F)]
    #s <- (apply(a, MARGIN = 2, FUN = sd) * 2) + abs(apply(a, MARGIN = 2, FUN = mean)) # two standard deviations of the mean
    
    for (k in 1:ncol(a)) {
      
      a[, k][a[, k] > Threshold] <- NA # assigning NA to values above the user-specified threshold
      a[, k][a[, k] < -Threshold] <- NA # assigning NA to values below (more negatives than) the user-specified threshold
      a[, k][a[, k] == 0] <- NA # assigning NA to zero values
      
      # storing missing data
      m <- length(a[is.na(a)]) # amount of NAs attributed to each electrode for each participant across the trials
      p <- substr(fi[j], 6, 8)
      e <- ElectNames[i]
      
      # missing data table
      miss_data_thresh <- rbind(miss_data_thresh, data.frame('Electrode' = e, 'Participant' = p, 'Trial' = u[k], 'Sel_Threshold' = Threshold, 'Count_data' = nrow(a)*ncol(a), 'Missing_count' = m, 'Missing_perct' = round(m/(nrow(a)*ncol(a))*100, digits = 4)))
    }
    
    # Electrode average table (average among trials, respecting the row order, of course)
    b <- data.frame(apply(a, MARGIN = 1, FUN = mean, na.rm = T))
    colnames(b) <- ElectNames[i]
    f <- cbind(f, b)
  }
  
  data_avgTrials <- f[, -1] # this is the table with the electrodes averaged by trials (averaged respecting the row order)
  #head(data_avgTrials)
  #dim(data_avgTrials)
  
  
  # Step 3: Flip the table
  
  data_02 <- data_avgTrials[1:615, ] # subset rows (-199 to 1000)
  
  data_02 <- as.data.frame(t(data_02)) # flip the table
  
  column.names <- paste(rep('T_', ncol(data_02)), seq(1, ncol(data_02), by = 1), sep = '')
  
  colnames(data_02) <- column.names
  
  data_02 <- data.frame('Electrode' = rownames(data_02), data_02)
  
  rownames(data_02) <- NULL
  
  # Step 4: Adding metadata
  # Extracting the meta info from the name of the file and add the info into the data frame
  #tmp$stdv <- substr(fi[j],nchar(fi[j])-8,nchar(fi[j])-7)
  data_02$Subj <- substr(fi[j], 6, 8)
  data_02$Group <- substr(fi[j], 12, 13)
  data_02$Age <- substr(fi[j], 10, 11)
  data_02$Sex <- substr(fi[j],9 ,9)
  data_02$Marker <- substr(fi[j], nchar(fi[j])-6, nchar(fi[j])-4)
  data_02$Condition <- ifelse(data_02$Marker %in% c("007", "011", "015", "019", "023", "027", "031"), "Ekman", "Under")
  data_02$Emotion <- data_02$Marker
  data_02$Emotion <- recode(data_02$Emotion, !!!level_key)
  
  data_02 <- data_02[, c('Electrode', 'Subj', 'Group', 'Age', 'Sex', 'Marker', 'Condition', 'Emotion', column.names)]
  
  data_threshold <- rbind(data_threshold, data_02)
  print(paste(fi[j], 'Done!', sep = ' '))
  
}

(Sys.time()-t0)


# Plotting missing data ---------------------------------------------------


# Extra: plotting the missing data

miss_data_thresh %>% group_by(Electrode) %>% summarise('mperc' = mean(Missing_perct)) %>% filter(mperc > 0) -> h

p <- ggplot(data = h, mapping = aes(x = reorder(Electrode, mperc), y = mperc)) + geom_bar(stat = 'identity', width = 0.7, fill = '#85d0ff') + coord_flip()
p <- p + geom_hline(yintercept = mean(h$mperc), linetype = 'dashed', size = 0.3, colour = 'red') # average percentage of missing data
p <- p + xlab('Electrode') + ylab('Percentage of missing data') + ggtitle(label = 'Electrodes with missing data')
p <- p + theme_bw()
p <- p + theme(axis.title = element_text(face = 'bold', size = 12)
               ,axis.text = element_text(size = 9, colour = 'black')
               ,title = element_text(face = 'bold', size = 12)
               ,panel.grid = element_line(size = 0.3)); p


# Saving elements ---------------------------------------------------------


# Saving threshold loop file
data_threshold = na.omit(data_threshold)
save(data_threshold, file="data-threshold.RData")

#save miss_data_thresh file
save(miss_data_thresh, file="data-miss-thresh.RData")

# removing NAs
data_thres = data_threshold
data_thres = na.omit(data_thres)

# Subsetting area of interest ---------------------------------------------

oz = data_thres %>% filter(Electrode=="T7"|Electrode=="T8"|Electrode=="C5"|
                             Electrode=="C6"|Electrode=="F7"|Electrode=="F8"|
                             Electrode=="CP3"|Electrode=="CP4"|
                             Electrode=="O1"|Electrode=="Oz"|Electrode=="O2"|
                             Electrode=="P3"|
                             Electrode=="Pz"|Electrode=="P4"|Electrode=="PO3"|
                             Electrode=="PO4")

oz = oz %>% mutate(ROI = recode(Electrode,
                                "O1" = "Visual",
                                "O2" = "Visual",
                                "Oz" = "Visual",
                                "P3" = "Visual",
                                "Pz" = "Visual",
                                "P4" = "Visual",
                                "PO4" = "Visual",
                                "PO3" = "Visual",
                                .default = "Verbal"
)) 
oz = oz[, c(1:8, 624, 9:623)]
