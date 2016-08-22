#' A Reference Class to represent a bank account.
#'
#' @field dat The ECG signal data.
#' @field window Default ECG window when not provided.
#' @field P,Q,R,S,T Peak locations and detection quality.
#'
#' @import methods
#' @export
eek <-
  setRefClass("eek",
              fields = c(
                "dat",
                "window",
                "P",
                "Q",
                "R",
                "S",
                "T"
              )
  )

eek$methods(initialize = function(file, channel = 1){

  "Initialize eek object using a single ECG channel."

  # Read ECG file and select i-th channel
  data <- read.delim(file, skip = 1, header = TRUE,
                     stringsAsFactors = FALSE)
  data <- data.frame("Time" = data[, 1], "ECG" = data[, channel + 1],
                     stringsAsFactors = FALSE)

  # Convert min:sec time to decimal seconds
  if(!is.numeric(data$Time)){

    minsec <- lapply(strsplit(data$Time, ":"), as.numeric)
    minsec <- unlist(lapply(minsec, function(x) x[1] * 60 + x[2]))
    data$Time <- minsec
  }

  # Zero-scale time
  data$Time <- data$Time - data$Time[1]
  window <<- 1:nrow(data)
  dat <<- data
})

eek$methods(filter = function(n = 1, lo = 1/40, hi = 1/2){

  "Filter ECG signal using a bandpass filter."

  bf <- signal::butter(n, c(lo, hi), type = "pass")
  dat$Filter <<- signal::filtfilt(bf, dat$ECG)
})

eek$methods(getR = function(minheight = .5, maxrate = 300){

  "Locate all R peaks and assign quality score."

  if(is.null(dat$Filter)) stop("Call eek$filter() before peak finding.")

  # Calculate minpeakdistance from maxrate
  totaltime <- floor(max(dat$Time)) - min(dat$Time)
  totalstep <- which.max(dat$Time >= floor(max(dat$Time)))
  persec <- totalstep / totaltime
  peakd <- ceiling(60 / maxrate * persec)

  # Find putative R peaks
  peaks <- pracma::findpeaks(
    x = dat$Filter,
    minpeakheight = minheight,
    minpeakdistance = peakd
  )

  # Use minimum peak height as a threshold to calculate a quality score
  peakind <- order(peaks[, 2])
  quality <- vector("numeric", length(peakind))
  quality[1] <- 0 # First quality score always zero
  lastRend <- peaks[peakind[1], 4] # Check quality from first peak onward
  for(i in 2:(length(peakind) - 1)){

    # Quality score for range from (i-1)th peak to (i+1)th peak
    nextRstart <- peaks[peakind[i + 1], 3]
    currentRall <- peaks[peakind[i], 3]:peaks[peakind[i], 4]
    qualityRange <- setdiff(lastRend:nextRstart, currentRall)
    quality[i] <- 1 - sum(dat[qualityRange, "Filter"] > (minheight / 2)) /
      (nextRstart - lastRend + 1)
    lastRend <- peaks[peakind[i], 4]
  }
  quality[length(peakind)] <- 0 # Last quality score always zero

  # Document putative R peaks
  final <- data.frame(peaks[peakind, c(3, 2, 4, 1)], quality^4)
  colnames(final) <- c("start", "peak", "end", "height", "quality")
  R <<- final
})

eek$methods(getPT = function(minheight = 0, maxrate = 300){

  "Locate all P and T peaks and assign quality score."

  if(!is.data.frame(R)) stop("Call eek$getR() before baseline correction.")

  # Calculate minpeakdistance from maxrate
  totaltime <- floor(max(dat$Time)) - min(dat$Time)
  totalstep <- which.max(dat$Time >= floor(max(dat$Time)))
  persec <- totalstep / totaltime
  peakd <- ceiling(60 / (maxrate*2) * persec)

  # Prepare output containers
  t.peaks <- vector("list", nrow(R)-1)
  p.peaks <- vector("list", nrow(R)-1)

  # Divide ECG into a series of R-R windows
  for(i in 1:(nrow(R) - 1)){

    Rstart <- R[i, "peak"]
    Rend <- R[i+1, "peak"]

    # Find putative PT peaks
    peaks <- pracma::findpeaks(
      x = dat$Filter[Rstart:Rend],
      minpeakheight = minheight,
      minpeakdistance = peakd
    )

    # Shift peak locations based on R start
    peaks[, c(2, 3, 4)] <- Rstart + peaks[, c(2, 3, 4)] - 1

    # Make sure T-wave is always index 1
    if(peaks[2, 3] < peaks[1, 3]){ # If T-wave is SECOND...

      peaks[c(1, 2), ] <- peaks[c(2, 1), ]
    }

    # Call first peak "T" and assign quality score
    df <- as.data.frame(t(c(peaks[1, c(3, 2, 4, 1)], R[i, "quality"])))
    colnames(df) <- c("start", "peak", "end", "height", "quality")
    rownames(df) <- NULL
    t.peaks[[i]] <- df

    # Call last peak "P" and assign quality score
    df <- as.data.frame(t(c(peaks[2, c(3, 2, 4, 1)], R[i+1, "quality"])))
    colnames(df) <- c("start", "peak", "end", "height", "quality")
    rownames(df) <- NULL
    p.peaks[[i]] <- df
  }

  # Save peak data in single data.frame
  T <<- do.call("rbind", t.peaks)
  P <<- do.call("rbind", p.peaks)
})

eek$methods(getQS = function(minheight = 0, maxrate = 300){

  "Locate all Q and S peaks and assign quality score."

  if(!is.data.frame(R)) stop("Call eek$getR() before baseline correction.")

  # Calculate minpeakdistance from maxrate
  totaltime <- floor(max(dat$Time)) - min(dat$Time)
  totalstep <- which.max(dat$Time >= floor(max(dat$Time)))
  persec <- totalstep / totaltime
  peakd <- ceiling(60 / (maxrate*2) * persec)

  # Prepare output containers
  s.peaks <- vector("list", nrow(R)-1)
  q.peaks <- vector("list", nrow(R)-1)

  # Divide ECG into a series of R-R windows
  for(i in 1:(nrow(R) - 1)){

    Rstart <- R[i, "peak"]
    Rend <- R[i+1, "peak"]

    # Find putative QS peaks
    peaks <- pracma::findpeaks(
      x = -1 * dat$Filter[Rstart:Rend],
      minpeakheight = minheight,
      minpeakdistance = peakd
    )

    # Shift peak locations based on R start
    peaks[, c(2, 3, 4)] <- Rstart + peaks[, c(2, 3, 4)] - 1

    # Make sure S-wave is always index 1
    if(peaks[2, 3] < peaks[1, 3]){ # If S-wave is SECOND...

      peaks[c(1, 2), ] <- peaks[c(2, 1), ]
    }

    # Call first peak "S" and assign quality score
    df <- as.data.frame(t(c(peaks[1, c(3, 2, 4, 1)], R[i, "quality"])))
    colnames(df) <- c("start", "peak", "end", "height", "quality")
    rownames(df) <- NULL
    s.peaks[[i]] <- df

    # Call last peak "Q" and assign quality score
    df <- as.data.frame(t(c(peaks[2, c(3, 2, 4, 1)], R[i+1, "quality"])))
    colnames(df) <- c("start", "peak", "end", "height", "quality")
    rownames(df) <- NULL
    q.peaks[[i]] <- df
  }

  # Save peak data in single data.frame
  S <<- do.call("rbind", s.peaks)
  Q <<- do.call("rbind", q.peaks)
})

eek$methods(smooth = function(l = 65, sd = .25){

  "Smooth ECG signal by convolving with Gaussian window."

  if(is.null(dat$Filter)) stop("Call eek$filter() before zoning.")

  # Convole the Gaussian window
  gw <- signal::gausswin(l, w = 1/sd)
  gw <- gw / sum(gw)
  smoothed <- stats::convolve(dat$Filter, gw, type = "open")

  # Align convolution with raw data
  offset <- (l - 1) / 2
  subset <- (1 + offset):(length(dat$Filter) + offset)
  dat$Gaus <<- smoothed[subset]
})

eek$methods(qplot = function(view){

  "Visualize an ECG window including peak locations."

  if(missing(view)) view <- window

  plot(dat$Time[view], dat$ECG[view], type = "l",
       xlab = "Time (sec)", ylab = "ECG (mV)",
       ylim = c(min(dat$ECG[view]) - 1,
                max(dat$ECG[view]) + 1))

  if(!is.null(dat$Filter)){

    points(dat$Time[view], dat$Filter[view], col = "orange", type = "l")
  }

  if(!is.null(dat$Gaus)){

    points(dat$Time[view], dat$Gaus[view], col = "green", type = "l")
  }

  for(wave in list(P, Q, R, S, T)){

    if(is.data.frame(wave)){

      range <- wave$peak %in% view
      points(dat$Time[wave$peak[range]], dat$ECG[wave$peak[range]], col = "red")
    }
  }
})

eek$methods(export = function(file = paste0(getwd(), "/eek-peaks.txt")){

  "Export ECG annotations including peak locations."

  # Combine P, Q, R, S, and T peak data
  cleaned <- lapply(list("P", "Q", "R", "S", "T"),
                    function(wave){

                      tryCatch(data.frame("ID" = wave, get(wave)),
                               error = function(e) data.frame())
                    })
  out <- do.call("rbind", cleaned)

  # Combine time, peak location, and peak label (see: PhysioBank)
  out <- out[order(out$peak), c("peak", "ID", "quality")]
  final <- data.frame(dat$Time[out$peak], out)

  # Save output
  write.table(final, file = file, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)

  return(TRUE)
})
