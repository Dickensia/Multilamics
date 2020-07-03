cptac<-cbioportal$...4
cptac<- levels(as.factor(cptac))
cb1 <- levels(as.factor(cbioportal$...1))
cb2 <- levels(as.factor(cbioportal$...2))
cb3 <- levels(as.factor(cbioportal$...3))
xenarna <- levels(as.factor(cbioportal$...5))
for (i in 1:length(cb1)) {
  cb1[i] <- substr(cb1[i],9,nchar(cb1[i])-3)
}
for (i in 1:length(cb2)) {
  cb2[i] <- substr(cb2[i],28,nchar(cb2[i])-3)
}
for (i in 1:length(cb3)) {
  cb3[i] <- substr(cb3[i],13,nchar(cb3[i])-3)
}
for (i in 1:length(xenarna)) {
  xenarna[i] <- substr(xenarna[i],1,nchar(xenarna[i])-4)
}

which(cptac %in% xenarna == FALSE)
















