peakw = function(cs, freq, req = 2) {
  for (index in 1:length(freq)) {
    start = index - req
    end = index + req
    if (start < 1 | end > length(freq)) {
      next
    }
    if (all(freq[start:end] <= freq[index])) {
      return(index)
    }
  }
  warning("No peak identified!")
  return(NA)
}

peakglobal = function(freq) {
  max(which(freq == max(freq)))
  }
