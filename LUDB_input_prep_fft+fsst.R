library(shiva)
library(ggplot2)
library(seewave)
library(fftw)
library(SynchWave)
# Set WFDB Path:
options(wfdb_path = 'wsl /usr/local/bin') 

# Load LUDB CSV file:
sample_info <-
  read.csv('C:/Users/shsan/Documents/Medical School/M2 Year/Darbar 2023-2024/PhysioNet Projects/lobachevsky-university-electrocardiography-database-1.0.1/ludb.csv')

# Include only those with normal sinus rhythm / brady / tachy
sinus_samples_index <- sample_info$ID[which(sample_info$Rhythms == "Sinus rhythm" | sample_info$Rhythms == "Sinus bradycardia" | sample_info$Rhythms == "Sinus tachycardia")] 

#build empty sample/solutions matrix
samples <- array(0, c(length(sinus_samples_index), 5000, 2))

# samples x 5000 x 4 signal types
solutions <- array(0,c(length(sinus_samples_index), 5000, 4)) #REMOVE!
  

#for loop predefine
rec_dir <- fs::path("C:/Users/shsan/Documents/Medical School/M2 Year/Darbar 2023-2024/PhysioNet Projects/lobachevsky-university-electrocardiography-database-1.0.1/data")
ann = 'i'

#Expand for all leads?

for (rec in 1:length(sinus_samples_index)) {
  test <- read_wfdb(record = sinus_samples_index[rec],
                    record_dir = rec_dir,
                    annotator = ann) #read file

  
  #build solutions array:
  dimmensions <- dim(test$signal)
  markers <- matrix(0, dimmensions[1])
  
  pwaves <- which(test$annotation$type == 'p')
  if (length(pwaves) > 0) {
    for (i in 1:length(pwaves)) {
      bounds <- test$annotation$sample[c(pwaves[i] - 1, pwaves[i] + 1)]
      
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 1 #'p' #remove "+1"??
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 1] <- 1
    }
  } else{
    print(paste("No P-Wave on Sample", sinus_samples_index[rec]))
  }
  
  qrswaves <- which(test$annotation$type == 'N')
  if (length(qrswaves) > 0) {
    for (i in 1:length(qrswaves)) {
      bounds <- test$annotation$sample[c(qrswaves[i] - 1, qrswaves[i] + 1)]
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 2 #'N'
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 2] <- 1
    }
  } else{
    print(paste("No QRS-Complex on Sample", sinus_samples_index[rec]))
  }
  
  
  twaves <- which(test$annotation$type == 't')
  if (length(twaves) > 0) {
    for (i in 1:length(twaves)) {
      bounds <- test$annotation$sample[c(twaves[i] - 1, twaves[i] + 1)]
      markers[(bounds[1] + 1):(bounds[2] + 1), ] <- 3 #'t'
      solutions[rec,(bounds[1] + 1):(bounds[2] + 1), 3] <- 1
      
    }
    } else{
      print(paste("No T-Wave on Sample", sinus_samples_index[rec]))
    }
  
  samples[rec, , 1] <- test$signal$i #can change to 1 thru 12 to be all leads
  samples[rec, , 2] <- markers
}


# Padding signal prior to first and last annotation with zeros:
for (i in 1:dim(samples)[[1]]){
  samples[i,0:(which(samples[i,,2] != 0 )[[1]]-1),1] <- 0 # can switch to other value
  samples[i, (tail(which(samples[i,,2] !=0), n = 1) + 1) :dim(samples)[[2]],1] <- 0 # can switch to other value
}

#FFT.test-------------------------------------------------------
fft_samp <- array(0, c(length(sinus_samples_index), 5000, 2))
fft_result <- apply(samples[, , 1], 1, function(row) fft(row))
fft_result_real <- Re(fft_result/5000) # Re pulls real numbers, used normalization
fft_result_imag <- Im(fft_result/5000)
fft_samp[, ,1] <- t(fft_result_real)
fft_samp[, ,2] <- t(fft_result_imag)

#Calculate the frequencies corresponding to the transformed values by rec.
n <- length(fft_samp[rec, , 1])
sampling_rate <- 500 

# Plot the magnitude spectrum using real components.
plot(frequencies, Mod(fft_samp[rec, ,1]), type = "b", xlab = "Frequency (Hz)", ylab = "Magnitude", main = "FFT ")

#FSST-----------------------------------------------------------------
fsst_dat <- data.frame(Signal = samples[rec, , 1],Time = 1:5000 / 500)
par(mar = c(2.5, 0.5, 2.25, 1.5))
fsst_dat<- fsst_dat[, c(2,1)]
plot(fsst_dat)
fsst_dat2 <- synsq_cwt_fw(tt = fsst_dat[[1]], x = fsst_dat[[2]], nv=16, opt=NULL) #synchrosqueezed output
View(fsst_dat2$Tx)


# Graph Verification ------------------------------------------------------
test_sample = 170

test_frame <- data.frame(Time = 1:5000 / 500, Signal = samples[test_sample, , 1])
#color code
colors <- samples[test_sample, , 2]
colors[colors == 1] <- 'p'
colors[colors == 2] <- 'N'
colors[colors == 3] <- 't'

test_plot <-
  ggplot(test_frame, aes(Time, Signal, color = colors)) + geom_path(linewidth =
                                                                      1, aes(group = 1)) + geom_point() + scale_x_continuous(breaks = seq(0, 10, 1))
test_plot


# Train Model: ---------------------------------------------------------

# Split into Train / Test data
split <- 0.7
train_size <- round(nrow(samples)*0.7)
train_rows <- sample(nrow(samples), train_size)

train_samples <- samples[train_rows,,]
test_samples <- samples[-train_rows,,]

