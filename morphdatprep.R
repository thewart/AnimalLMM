morphdat <- fread("~/Dropbox/cayomorphodata.csv")
#morphdat <- morphdat[!is.na(HeadCircum)]
#aaargh wtf 
morphdat <- morphdat[SkinFoldAbUp.Abdomen.ABOVE!=0 & SkinFoldAbDown.Abdomen.BELOW!=0]
setkey(morphdat,"ID","Trapyear")
rep <- names(morphdat[,table(ID)])[morphdat[,table(ID)]>1]
for (i in rep) {
  moo <- which(morphdat$ID == i)
  morphdat <- morphdat[-moo[-1]]
}

if (geno=="snp") morphdat <- morphdat[ID %in% G$ID]

X <- model.matrix( ~ 1 + as.factor(Trapyear) + sex*poly(Age,2),data=morphdat)[,-1]
X[,5:6] <- apply(X[,5:6],2,function(x) (x-mean(x))/sd(x))
X[,7:8] <- X[,5:6] * X[,4]

measurer <- morphdat$Measurer
misc <- names(table(measurer))[table(measurer)<10]
measurer[(measurer %in% misc) | is.na(measurer)] <- "MISC"

Z <- data.table(measurer,morphdat[,.(group)]) %>% 
  lapply(function(x) as.factor(x) %>% as.numeric) %>% 
  unlist %>% matrix(ncol=2)

Y <- morphdat[,log(cbind(WeightKg,CrownRumpcm,UpArmCircum,#HeadCircum,
                         SkinFoldAbUp.Abdomen.ABOVE+SkinFoldAbDown.Abdomen.BELOW,SkinFoldSubscap))] %>% apply(2,function(x) (x-mean(x))/sd(x))

standat <- list(N=nrow(Y),D=ncol(Y),P=ncol(X),V=ncol(Z),K=apply(Z,2,function(x) length(unique(x))),
                Y=Y,X=X,Z=t(Z))