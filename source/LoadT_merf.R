

clinc_merge <- function(indat,inclinc,intype){
  dat = as.data.frame(indat)
  clinc = as.data.frame(inclinc)
  inset_samples <- intersect(colnames(dat),clinc$submitter)
  dat <- dat[,colnames(dat) %in% inset_samples]
  clinc <- clinc[clinc$submitter %in% inset_samples,]
  
  clinc <- clinc[order(clinc[,intype]),]
  return(clinc)
}

dat_merge <- function(indat,inclinc,intype){
  dat = as.data.frame(indat)
  clinc = as.data.frame(inclinc)
  inset_samples <- intersect(colnames(dat),clinc$submitter)
  dat <- dat[,colnames(dat) %in% inset_samples]
  clinc <- clinc[clinc$submitter %in% inset_samples,]
  
  clinc <- clinc[order(clinc[,intype]),]
  dat <- dat[,clinc$submitter]
  return(dat)
}