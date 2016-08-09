#' A Reference Class to represent a bank account.
#'
#' @field balance A length-one numeric vector.
eek <-
  setRefClass("eek",
              fields = c(
                "dat",
                "window",
                "bounds",
                "startzones",
                "endzones",
                "R"
              )
  )

# Methods describe as doc string inside method definition...
eek$methods(read = function(file){

  # "DESCRIBE METHOD HERE LIKE THIS"

  data <- read.delim(file, skip = 1, header = TRUE, stringsAsFactors = FALSE)
  colnames(data) <- c("Time", "ECG")

  minsec <- lapply(strsplit(data$Time, ":"), as.numeric)
  minsec <- unlist(lapply(minsec, function(x) x[1] * 60 + x[2]))
  data$Time <- minsec

  dat <<- data
  window <<- 1:nrow(dat)
})

eek$methods(filter = function(n = 1, lo = 1/40, hi = 1/2){

  bf <- signal::butter(n, c(lo, hi), type = "pass")
  dat$Filter <<- signal::filtfilt(bf, dat$ECG)
})

eek$methods(zone = function(l = 65, sd = .25){

  if(is.null(dat$Filter)) stop("Call eek$filter() before zoning.")

  # Convole the Gaussian window
  gw <- signal::gausswin(l, w = 1/sd)
  gw <- gw / sum(gw)
  smoothed <- stats::convolve(dat$Filter, gw, type = "open")

  # Align convolution with raw data
  offset <- (l - 1) / 2
  subset <- (1 + offset):(length(dat$Filter) + offset)
  dat$Gaus <<- smoothed[subset]

  # Find all of the zero-crossings
  res <- EMD::extrema(dat$Gaus)
  bounds <<- res$cross[, 1]
})

eek$methods(getR = function(minheight = .5, maxrate = 300){

  if(is.null(dat$Filter)) stop("Call eek$filter() before peak finding.")

  # Index putative R peaks
  totaltime <- floor(max(dat$Time)) - min(dat$Time)
  totalstep <- which.max(dat$Time >= floor(max(dat$Time)))
  persec <- totalstep / totaltime
  peakd <- ceiling(60 / maxrate * persec)

  peaks <- pracma::findpeaks(
    x = dat$Filter,
    minpeakheight = minheight,
    minpeakdistance = peakd
  )

  # Use minimum peak height to calculate quality score
  peakind <- order(peaks[, 2])
  quality <- vector("numeric", length(peakind))
  quality[1] <- 0 # First quality score always zero
  lastRend <- peaks[peakind[1], 4] # Check quality from first peak onward
  for(i in 2:length(peakind)){

    # Quality score for range from (i-1)th peak to (i+1)th peak
    nextRstart <- peaks[peakind[i + 1], 3]
    currentRall <- peaks[peakind[i], 3]:peaks[peakind[i], 4]
    qualityRange <- setdiff(lastRend:nextRstart, currentRall)
    quality[i] <- 1 - sum(dat[qualityRange, "Filter"] > (minheight / 2)) /
      (nextRstart - lastRend + 1)
    lastRend <- peaks[peakind[i], 4]
  }

  final <- data.frame(peaks[peakind, c(3, 2, 4, 1)], quality^4)
  colnames(final) <- c("start", "peak", "end", "height", "quality")
  R <<- final
})

eek$methods(importance = function(threshold){

  # Find importance zones (area between bounds)
  startzones <<- vector("numeric")
  endzones <<- vector("numeric")
  s <- 1
  i <- 1
  for(e in bounds){

    # First: Make sure zero-crossings flank an extrema
    if(any(res$minindex[, 1] %in% s:e) | any(res$maxindex[, 1] %in% s:e)){

      startzones[i] <<- s
      endzones[i] <<- e
      i <- i + 1
    }

    s <- e + 1
  }
})

eek$methods(qplot = function(view){

  if(missing(view)){

    view <- window
  }

  plot(dat$Time[view], dat$ECG[view], type = "l",
       xlab = "Time (sec)", ylab = "ECG (mV)")

  if(!is.null(dat$Filter)){

    points(dat$Time[view], dat$Filter[view], col = "orange", type = "l")
  }

  if(!is.null(dat$Gaus)){

    points(dat$Time[view], dat$Gaus[view], col = "green", type = "l")
  }

  if(is.numeric(bounds)){

    for(mark in bounds[bounds > min(view) & bounds < max(view)]){

      abline(v = dat$Time[mark], col = "pink")
    }
  }

  if(is.numeric(startzones) & is.numeric(endzones)){

    range <- startzones >= min(view) & endzones <= max(view)
    st <- startzones[range]
    en <- endzones[range]
    for(i in 1:length(st)){

      polygon(x = c(rep(dat$Time[st[i]], 2),
                    rep(dat$Time[en[i]], 2)),
              y = c(-100, 100, 100, -100),
              col = "pink", angle = 0,
              density = 5)
    }
  }
})
