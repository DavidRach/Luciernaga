unmix_ff <- function(fs, control, unmixMethod, multiplier, outpath) {
  expresionData<-exprs(fs)
  expresionData<-as.data.frame(expresionData)
  expresionData<-expresionData[,-grep("SC|SS|FS", names(expresionData))]
  expresionData<-expresionData[,grep("-A", names(expresionData))]
  control<-control[names(expresionData)] #Reorder data to match the expression data - vital if using exported data from FlowJo

  if(unmixMethod=="lsfit") {
    ls_corr <- lsfit(x = t(control), y = t(expresionData), intercept = FALSE)
    unmixResult <- t(ls_corr$coefficients)
    unmixResult<-unmixResult*multiplier
  }
  else if(unmixMethod=="ginv"){
    pseudoinverse<-ginv(as.matrix(t(control)))
    colnames(pseudoinverse)<-colnames(control)
    rownames(pseudoinverse)<-rownames(control)
    unmixResult<-t(apply(expresionData,1,function(x)colSums(t(pseudoinverse)*x)))
    unmixResult<-unmixResult*multiplier
  }
  else if(unmixMethod=="qr.solve"){
    qrs<-qr.solve(as.matrix(t(control)),as.matrix(t(expresionData)))
    unmixResult<-as.data.frame(t(qrs))
    unmixResult<-unmixResult*multiplier
    unmixResult<-as.matrix(unmixResult)
  }
  else if(unmixMethod=="lm.fit"){
    lmfit<-.lm.fit(as.matrix(t(control)),as.matrix(t(expresionData)))
    unmixResult <- t(lmfit$coefficients)
    unmixResult<-unmixResult*multiplier
    unmixResult<-data.frame(unmixResult)
    colnames(unmixResult)<-rownames(control)
    unmixResult<-as.matrix(unmixResult)
  }
  else if(unmixMethod=="crosspod"){
    crspod<-solve(crossprod(as.matrix(t(control)))) %*% crossprod(as.matrix(t(control)),as.matrix(t(expresionData)))
    crspod<-crspod*multiplier
    crspod<-data.frame(t(crspod))
    colnames(crspod)<-rownames(control)
    unmixResult<-as.matrix(crspod)
  }
  else if(unmixMethod=="nnls"){
    expresionData2<-as.matrix(expresionData)
    system.time({nnlsunmixed<-apply(expresionData2, 1, function(x)nnls(as.matrix(t(control)), x)$x)})
    nnlsunmixed<-nnlsunmixed*multiplier
    nnlsunmixed<-data.frame(t(nnlsunmixed))
    colnames(nnlsunmixed)<-rownames(control)
    unmixResult<-as.matrix(nnlsunmixed)
  }
  else if(unmixMethod=="baselm"){
    system.time({lmbase<-apply(expresionData, 1, function(x)lm (x ~ t(control))$coefficients)})
    lmbase<-lmbase*multiplier
    lmbase<-data.frame(t(lmbase))
    ncol(lmbase)
    lmbase<-lmbase[,2:ncol(lmbase)]
    colnames(lmbase)<-rownames(control)
    unmixResult<-as.matrix(lmbase)
  }

  new_fr<-fr_append_cols(fs, unmixResult)
  #thecolnames <- colnames(new_fr)[c(1, 18:20, 37:42, 65:93)]
  #panel <- data.frame(fcs_colname = factor(c(thecolnames)))
  #new_fr1 <- new_fr[, levels(panel$fcs_colname)]
  #exprs(new_fr1)
  #colnames(new_fr1)

  keyword(new_fr)$`FIL`<-paste(keyword(new_fr)$`FIL`,unmixMethod)
  keyword(new_fr)$GUID<-paste(keyword(new_fr)$GUID,unmixMethod)
  keyword(new_fr)$TUBENAME<-paste(keyword(new_fr)$TUBENAME,unmixMethod, sep = "_")

  location <- paste0(outpath, "/", keyword(new_fr)$`TUBENAME`)
  write.FCS(new_fr, filename=paste0(location,".fcs"))
}
