

Interf <- function(indat){
  len <- length(lengths(indat))
  nam <- names(indat)
  inall <- c()
  for (i in 1:len) {
    assign(paste0("len",i),indat[[i]])
    inall <- unique(c(inall,get(paste0("len",i))))
  }
  medat <- data.frame(all = inall)
  rownames(medat) <- medat$all
  for (k in 1:length(nam)) {
    # nam2 <- nam[k]
    gene2 <- get(paste0("len",k))
    medat <- within(medat,{
      nam2 <- ifelse(inall %in% gene2,1,0)
    })
    names(medat)[k+1] <- nam[k]
  }
  medat <- medat[,-1]
  medat <- cbind(medat, SUM=rowSums(medat))
  
  return(medat)
}



##  字符串转为数值
char2numf <- function(indat){
  NumDat <- indat
  rname <- rownames(indat)
  NumDat <- apply(as.matrix(NumDat),2,function(x) as.numeric(x))
  NumDat <- as.data.frame(NumDat)
  rownames(NumDat) <- rname
  
  return(NumDat)
}


##  
NA2Numf <- function(indat){
  rname <- rownames(indat)
  NADat <- as.data.frame(indat)
  NADat <- apply(NADat, 2, function(x){
    x[is.na(x)] = mean(x,na.rm = T)
  })
  
  rownames(NADat) <- rname
  return(NADat)
}

##  
row4 <- function(indat, n, num){
  vars = names(table(indat[,n]))
  alldat <- data.frame()
  for (i in 1:length(vars)) {
    dat <- indat[indat[,n] == vars[i],]
    if (length(rownames(dat)) > num) {
      dat <- dat[1:num,]
      assign(paste0("dat",i),dat)
      alldat <- rbind(alldat,dat)
    }else{
      assign(paste0("dat",i),dat)
      alldat <- rbind(alldat,dat)
    }
  }
  return(alldat)
}