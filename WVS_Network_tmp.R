#Author=ORCID:0000-0002-3324-2255
#Task= Work in progress degree distribution of combinatorical value selection of WVS

#1.DATA PREPARATION
#1.DATA PREPARATION
#1A. Loading Libraries and defining functions
rm(list = ls()) #deleting memory
library(psych) #for KMO
library(pvclust) #for hier clustering
library(kohonen) #for SOM clustering
library(DescTools) #for CramersV and TT
library(NbClust) #for clustering numbers
library(igraph) #fot graph analyses

#there is no predict.kmeans
predict.kmeans <- function(object,
                           newdata,
                           method = c("centers", "classes")) {
  method <- match.arg(method)
  
  centers <- object$centers
  ss_by_center <- apply(centers, 1, function(x) {
    colSums((t(newdata) - x) ^ 2)
  })
  best_clusters <- apply(ss_by_center, 1, which.min)
  
  if (method == "centers") {
    centers[best_clusters, ]
  } else {
    best_clusters
  }
}

#there is no within sum of squares 
wcss <- function(data,clClassification,verbose=TRUE,... ) {
  if(is.matrix(data)){
    mat <-data
    df <- as.data.frame(mat)
  }else if(is.data.frame(data)){
    df <- data
    mat<-as.matrix(data)
  }else{
    if(verbose){
      cat("The data has to be a matrix or a dataframe\n",sep="")
    }
    stop()
  }
  if(nrow(mat)!=length(clClassification)){
    if(verbose){
      cat("The classification vector length (",length(clClassification),
          ") differs from the number of observation in the dataset (",nrow(mat),") \n",sep="")
    }
    stop()
  }else{
    sizeCl <- summary(as.factor(clClassification))
    WCSSByClByVar <- tapply(mat, list(rep(clClassification,ncol(mat)),col(mat)),
                            function(x) var(x, na.rm=TRUE)*(length(x)-1))
    WCSSByClByVar <- as.data.frame(WCSSByClByVar)
    WCSSByClByVar <- setNames(WCSSByClByVar, names(df))
    WCSSByCl <- rowSums(WCSSByClByVar)
    WCSS <- sum(WCSSByCl)
  }
  nClusters = nlevels(as.factor(clClassification))
  list(
    WCSSByClByVar=WCSSByClByVar,
    WCSSByCl = WCSSByCl,
    WCSS = WCSS,
    nClusters=nClusters,
    sizeCl=sizeCl
  )
}


#1B. Reading the EVW WVS file
#EVS/WVS (2022). European Values Study and World Values Survey: Joint EVS/WVS 2017-2022 Dataset (Joint EVS/WVS). 
#JD Systems Institute & WVSA. Dataset Version 3.0.0, doi:10.14281/18241.19
wvs <- readRDS("~/WVS_Huntigton/raw_data/evswvs.rds") 

#reading the Huntington classification
#hciv <- read.csv(file="https://raw.githubusercontent.com/yudhanjaya/HuntingtonsCiv/master/Huntingtons.csv", header = TRUE, 
#               sep = ",", encoding="UTF-8", stringsAsFactors=FALSE) #HuntingtonCivilizationDataset
#saveRDS(hciv, file="~/WVS_Huntigton/raw_data/hciv.rds") saving the Huntington file

hciv <- readRDS(file="~/WVS_Huntigton/raw_data/hciv.rds")
hcivs <- hciv[,c(1:7)]
hcivs[12,7] <- "Antarctica" #completing the only missing value

#Adding to the dataset
wvs_civ <- merge.data.frame(wvs, hcivs, by.x = 'cntry', by.y = 'country.code', all=FALSE, all.x = TRUE)

#Checking the merge:
table(paste(wvs_civ$cntry_AN, wvs_civ$name))
dim(wvs)[1]==dim(wvs_civ)[1]

# 1C. Countries Selection
#how many countries per civilizations?
sum(table(wvs_civ$cntry_AN, wvs_civ$Huntingtons_class)>0) #surveyed countries sum
sum(table(hcivs$Huntingtons_class)) #sum countries

summary(table(wvs_civ$cntry_AN, wvs_civ$Huntingtons_class)>0) #surveyed countries per civilization 
table(hcivs$Huntingtons_class) #sum countries per civilization

#RESUKTS
#African(3/41), Antartica(0/1), Ethiopia (1/1), Japan (1/1), 
#Latin American(12/45, 26%),Buddhist(3/7, 42.8%), Sinic (7/8 87.5%)
#Islamic(19/48, 39.5%), Orthodox(16/17 94.1%), Western (26/72 36.1%)


# 1D. Variable Selection (Theoretical, Design driven)
wwse <- wvs_civ[ ,c(
               'A001','A002', 'A003', 'A004', 'A005', "A006",
#               'A027','A029','A030','A032','A034','A035','A038','A039','A040','A041','A042',
#               'A124_02','A124_03','A124_06','A124_08','A124_09', 
#               'C001_01', 'D081','D026_03','D026_05','C038','C039','C041','D054','D059','D060','D061','D078',
#               'F114A','F115','F116','F117','F118','F119','F120','F121','F122','F123','F132','F144_02','E290',
#                'A165','D001_B','G007_18_B','G007_33_B','G007_34_B','G007_35_B','G007_36_B',
                'Huntingtons_class', 'name', 'region' ,
                'size_5c', 'X001', 'X003', 'G027A', 'X007', 'x026_01', 'X011', 'X013', 'X025R'
               )
              ]

wwse$Huntingtons_class <- as.factor(wwse$Huntingtons_class)
wwse$name <- as.factor(wwse$name)
wwse$region <- as.factor(wwse$region)

#removing NAs
wwse <- na.omit(wwse[,c(1:18)])

#distribution                
round((table(wwse$A001) / sum(table(wwse$A001))) * 100,1)
round((table(wwse$A002) / sum(table(wwse$A001))) * 100,1)
round((table(wwse$A003) / sum(table(wwse$A001))) * 100,1)
round((table(wwse$A004) / sum(table(wwse$A001))) * 100,1)
round((table(wwse$A005) / sum(table(wwse$A001))) * 100,1)
round((table(wwse$A006) / sum(table(wwse$A001))) * 100,1)

#combinatorical values
dsum <- 1000000*wwse$A001 + 100000*wwse$A002 + 10000*wwse$A003 + 1000*wwse$A004 + 100*wwse$A004 + 10*wwse$A005 + wwse$A006
length(dsum)
dsumt <- as.data.frame(table(dsum))


dsumts <- dsumt[order(-dsumt$Freq),] #reordering into decresing order by freq, however it keeps rowmanes the original order
rownames(dsumts)<-NULL

sumdstf <- sum(dsumts$Freq)
dsumts$p1 <- dsumts$Freq / sumdstf
sum(dsumts$p1)

head(dsumts)

#seeing the distribution of value-orders:
plot(dsumts$p1~rownames(dsumts))
plot(dsumts$p1~rownames(dsumts), log="xy", main="Lehetséges értékkombinációkból a választások valószínűsége",
     xlab="x=kombináció azonosítója", ylab="y=p(x)", sub="--piros szagatott az egyenletes eloszlás értéke--" )
abline(h=1/(4*4*4*4*4*4), col="red", lwd=3, lty=2)

#forditva van kodolva.

ta01 <- (4-wwse$A001) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))
ta02 <- (4-wwse$A002) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))
ta03 <- (4-wwse$A003) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))
ta04 <- (4-wwse$A004) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))
ta05 <- (4-wwse$A005) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))
ta06 <- (4-wwse$A006) / ((4-wwse$A001) + (4-wwse$A002) + (4-wwse$A003) + (4-wwse$A004) + (4-wwse$A005) + (4-wwse$A006))

mean(ta01, na.rm = TRUE) #family
mean(ta02, na.rm = TRUE) #friends
mean(ta03, na.rm = TRUE) #leisure time
mean(ta04, na.rm = TRUE) #politics
mean(ta05, na.rm = TRUE) #work
mean(ta06, na.rm = TRUE) #religion

hist(ta01)
hist(ta02)
hist(ta03)
hist(ta04)
hist(ta05)
hist(ta06)

tadf <- data.frame(ta01, ta02, ta03, ta04, ta05, ta06)

kmeans(na.omit(tadf), centers = 2)$tot.withinss
kmeans(na.omit(tadf), centers = 3)$tot.withinss
kmeans(na.omit(tadf), centers = 4)$tot.withinss
kmeans(na.omit(tadf), centers = 5)$tot.withinss
kmeans(na.omit(tadf), centers = 6)$tot.withinss
kmeans(na.omit(tadf), centers = 7)$tot.withinss
kmeans(na.omit(tadf), centers = 8)$tot.withinss
kmeans(na.omit(tadf), centers = 9)$tot.withinss
kmeans(na.omit(tadf), centers = 10)$tot.withinss

summary(tadf)
psych::describe(tadf)
boxplot(tadf)


#15 edges (dynamics network)
d12<-(wwse$A001 - wwse$A002)
d13<-(wwse$A001 - wwse$A003)
d14<-(wwse$A001 - wwse$A004)
d15<-(wwse$A001 - wwse$A005)
d16<-(wwse$A001 - wwse$A006)

d23<-(wwse$A002 - wwse$A003)
d24<-(wwse$A002 - wwse$A004)
d25<-(wwse$A002 - wwse$A005)
d26<-(wwse$A002 - wwse$A006)

d34<-(wwse$A003 - wwse$A004)
d35<-(wwse$A003 - wwse$A005)
d36<-(wwse$A003 - wwse$A006)

d45<-(wwse$A003 - wwse$A004)
d46<-(wwse$A003 - wwse$A005)

d56<-(wwse$A003 - wwse$A004)

#15 edges other direction
d21<- -(wwse$A001 - wwse$A002)
d31<- -(wwse$A001 - wwse$A003)
d41<- -(wwse$A001 - wwse$A004)
d51<- -(wwse$A001 - wwse$A005)
d61<- -(wwse$A001 - wwse$A006)

d32<- -(wwse$A002 - wwse$A003)
d42<- -(wwse$A002 - wwse$A004)
d52<- -(wwse$A002 - wwse$A005)
d62<- -(wwse$A002 - wwse$A006)

d43<- -(wwse$A003 - wwse$A004)
d53<- -(wwse$A003 - wwse$A005)
d63<- -(wwse$A003 - wwse$A006)

d54<- -(wwse$A003 - wwse$A004)
d64<- -(wwse$A003 - wwse$A005)

d65<- -(wwse$A003 - wwse$A004)


df <- data.frame(d12,d13,d14,d15,d16,d23,d24,d25,d26,d34,d35,d36,d45,d46,d56, 
                 d21,d31,d41,d51,d61,d32,d42,d52,d62,d43,d53,d63,d54,d64,d65)

#how many puts family above all:
famfirst <- (d12<0) & (d13<0) & (d14<0) & (d15<0) & (d16<0)
sum(famfirst) / length(famfirst)

#how many puts friends above all:
frfirst <- (d21<0) & (d23<0) & (d24<0) & (d25<0) & (d26<0)
sum(frfirst) / length(frfirst)

#how many puts leisure time above all:
lefirst <- (d31<0) & (d32<0) & (d34<0) & (d35<0) & (d36<0)
sum(lefirst) / length(lefirst)

#how many puts politics time above all:
pofirst <- (d41<0) & (d42<0) & (d43<0) & (d45<0) & (d46<0)
sum(pofirst) / length(pofirst)

#how many puts work above all:
wofirst <- (d51<0) & (d52<0) & (d53<0) & (d54<0) & (d56<0)
sum(wofirst) / length(wofirst)

#how many puts religion above all:
refirst <- (d61<0) & (d62<0) & (d63<0) & (d64<0) & (d65<0)
sum(refirst) / length(refirst)


summary(df)
psych::describe(df)
boxplot(df)

#GRAPH
dt1 <- psych::describe(df)
g1<-data.frame( V1=c(1,1,1,1,1,2,2,2,2,3,3,3,4,4,5,2,3,4,5,6,3,4,5,6,4,5,6,5,6,6) , 
                V2=c(2,3,4,5,6,3,4,5,6,4,5,6,5,6,6,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5) , 
                rowname1 = c('d12','d13','d14','d15','d16','d23','d24','d25','d26','d34','d35','d36','d45','d46','d56', 
                'd21','d31','d41','d51','d61','d32','d42','d52','d62','d43','d53','d63','d54','d64','d65'),
                   mean1= dt1$mean , 
                   sd1= dt1$sd , 
                   var1= dt1$sd^2 ,
                   median1 = dt1$median
              )

graph1 <-graph_from_data_frame(g1[c(1:2)],directed = TRUE)
graph1w <- set_edge_attr(graph1, name="weight", value=g1$mean1)
is.weighted(graph1w)

layout2 <- layout_in_circle(graph1w)
plot(graph1w, edge.label=round(E(graph1w)$weight,1), layout=layout2)

pgraph1<-subgraph.edges(graph1w, eids=which(E(graph1w)$weight>0))
is.weighted(pgraph1)

layout1 <- layout.fruchterman.reingold(pgraph1, weights = (2/(E(pgraph1)$weight)) )
plot(pgraph1, edge.label=round(E(pgraph1)$weight,1), layout=layout1)

#let!s take a look at variances of differences of espoused values:

g1[order(-g1$var1),c('rowname1', 'mean1', 'var1')]
plot(g1[order(-g1$var1),c('rowname1', 'mean1', 'var1')]$var1)

#ez alapjan, 3,6,4,5,1,2 a sorrend, vegyuk az elso 4-et

pgraph1n3645<-subgraph(pgraph1, vids=c(3,6,4,5))
is.weighted(pgraph1n3645)

layout1 <- layout.fruchterman.reingold(pgraph1n3645, weights = (1/(E(pgraph1n3645)$weight)) )
plot(pgraph1n3645, edge.label=round(E(pgraph1n3645)$weight,1), layout=layout1)

layout1 <- layout.fruchterman.reingold(pgraph1n3645, weights = (1/(E(pgraph1n3645)$weight)) , dim=3)


grafcs <- kmeans(df, 5)

wwsedf <- cbind(wwse, df, csop=grafcs$cluster) #5 klaszter

cs1ms<-psych::describe(wwsedf[wwsedf$csop==1,c(22:36)])$mean
cs2ms<-psych::describe(wwsedf[wwsedf$csop==2,c(22:36)])$mean
cs3ms<-psych::describe(wwsedf[wwsedf$csop==3,c(22:36)])$mean
cs4ms<-psych::describe(wwsedf[wwsedf$csop==4,c(22:36)])$mean
cs5ms<-psych::describe(wwsedf[wwsedf$csop==5,c(22:36)])$mean


cs5df <- data.frame(cs1ms, cs2ms, cs3ms, cs4ms, cs5ms)

plot(as.numeric(cs5df[1,]))
plot(as.numeric(cs5df[2,]))
plot(as.numeric(cs5df[3,]))
plot(as.numeric(cs5df[4,]))
plot(as.numeric(cs5df[5,]))
plot(as.numeric(cs5df[6,]))
plot(as.numeric(cs5df[7,]))
plot(as.numeric(cs5df[8,]))
plot(as.numeric(cs5df[9,]))
plot(as.numeric(cs5df[10,]))
plot(as.numeric(cs5df[11,]))
plot(as.numeric(cs5df[12,]))
plot(as.numeric(cs5df[13,]))
plot(as.numeric(cs5df[14,]))
plot(as.numeric(cs5df[15,]))
