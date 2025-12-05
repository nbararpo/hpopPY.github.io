
library(readxl)
library(sva)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(variancePartition)
library('variancePartition')
library('edgeR')
library('BiocParallel')
library(DESeq2)
library(tidyverse)
```

#load data

rm(list=ls())


meta<- read.csv("/Users/nasimbararpour/Library/CloudStorage/Box-Box/NasimBararpour/hPOP/MetaData/FinalVersion/Meta/metadata_HPOP_FV5.csv")

rownames(meta)<- meta$Finalcode
# Reformat metadata
meta <- meta[order(meta$Finalcode),]
meta$Height <- as.numeric(meta$Height)
meta$Height[meta$Height < 50 &!is.na(meta$Height)] <- meta$Height[meta$Height < 50 & !is.na(meta$Height)]*100
meta$Weight <- as.numeric(meta$Weight)
meta$bmi <- as.numeric(meta$Weight)/(as.numeric(meta$Height)/100)^2
meta$bmi_groups <- cut(as.numeric(meta$bmi), breaks = c(0,25,30,70),labels=c('Normal','OverWeight','Obese'))
meta$Pt_Age_Range[meta$Pt_Age_Range=='>70'] <- '>60'
meta$Pt_Age_Range[meta$Pt_Age_Range=='60 - 70'] <- '>60'
meta$Pt_Age_Range <- factor(meta$Pt_Age_Range,
                            levels=c("20 - 30",'30 - 40', '40 - 50', '50 - 60','>60'))# '60 - 70', '>70'))

metdata<-meta




### merge metadata with urine


data<- merge(metdata,df, by.x= "Finalcode", by.y= "name")
rownames(data)<- data$Finalcode

#dim(data)#only 123 samples are left
#rownames(data)<- data$Finalcode#take sample as row names
#dx<- data[-1]# remove name of samples

```


#separate metadata from expression data
#For merged data
```{r}
metdafv<-data[1:210] ### metadata
dExp<- data[211:575] #expression Data ( metabolome and lipidome data)
```



#make non-ordered Ordinal factor
```{r}
metdafv<- metdata
unique(metdafv$age)

metdafv$age<- factor(metdafv$Pt_Age_Range)
unique(metdafv$age)

#Age
metdafv$Pt_Age_Range[metdafv$Pt_Age_Range == "20 - 30"]<- "A"
#metdafv$Pt_Age_Range[metdafv$Pt_Age_Range == "20 - 30" ] <- "A"
metdafv$Pt_Age_Range[metdafv$Pt_Age_Range== "30 - 40" ]<- "B"
metdafv$Pt_Age_Range[metdafv$Pt_Age_Range== "40 - 50" ]<- "C"
metdafv$Pt_Age_Range[metdafv$Pt_Age_Range== "50 - 60" ]<- "D"
metdafv$Pt_Age_Range[metdafv$Pt_Age_Range== ">60" ]<- "E"
metdafv$Pt_Age_Range

###Eth
metdafv$Eth[metdafv$Eth=="Asian"]<-"A"
metdafv$Eth[metdafv$Eth=="Caucasian"]<-"C"
metdafv$Eth[metdafv$Eth=="Indian"]<-"I"
metdafv$Eth[metdafv$Eth=="Other"]<-"O"
metdafv$Eth

##bmi_group
metdafv$bmi_category[metdafv$bmi_category=="normal"]<- "A"
metdafv$bmi_category[metdafv$bmi_category=="overweight"]<- "B"
metdafv$bmi_category[metdafv$bmi_category=="obese"]<- "C"
#Missing bmi_group to NN variable
metdafv$bmi_category[metdafv$bmi_category=="no"]<- "NN"

metdafv$bmi_category


#Site
metdafv$Conf_Site_Code[metdafv$Conf_Site_Code=="1"]<- "A.B"
metdafv$Conf_Site_Code[metdafv$Conf_Site_Code=="2"]<- "B.T"
metdafv$Conf_Site_Code[metdafv$Conf_Site_Code=="3"]<- "C.D"
metdafv$Conf_Site_Code[metdafv$Conf_Site_Code=="4"]<- "D.O"
metdafv$Conf_Site_Code[metdafv$Conf_Site_Code=="5"]<- "E.C"

metdafv$Conf_Site_Code


#gender
#missing gender considered as Female
metdafv$Sex[metdafv$Sex=="no"]<- "Female"

metdafv$Sex
##
```
#Remove samples that is not present

```{r}
#Remove samples that
#meta <- metdafv[metdafv$Finalcode %in% rownames(horm_data), ]

#head(meta, n=10)
#head(horm_data, n=10)
#dim(meta)

```

#Data Preparation
```{r}
#rownames(dExp)<- dExp$Finalcode
#Pro_data<-horm_data[-c(1:2)]
geneEx<- data.frame(dExp)
#save column Name
colnames<- colnames(geneEx)
```

#start model
```{r}
geneEx[geneEx == 0] <- NA
#colnames(geneEx)<- 1:198 ##to be handled by DESEQ
geneExpr01<-t(geneEx) ###to keep sample in column for DERAM

geneExpr01[is.na(geneExpr01)] <- 0
```
#Dream from DESEQ2

```{r}
geneExpr = DGEList(geneExpr01)

geneExpr = calcNormFactors(geneExpr )

param = SnowParam(4, "SOCK", progressbar=TRUE)
```

#Before appling model
```{r}
library(ggplot2)
meta<-metdafv # to keep the same file name
meta$Pt_Age_Range <- gsub(" ","", meta$Pt_Age_Range)
```
#Model 1: Age as fix effect and all other factor as randome effect
#define model
```{r}
form <- ~ Pt_Age_Range +  (1|Conf_Site_Code) +(1|Sex)+(1|Eth)+ (1|bmi_category)

# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, meta, BPPARAM=param, plot = TRUE )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream(vobjDream,form, meta )
fitmm = eBayes(fitmm)

head(fitmm$design, 4)
```


# Get results of hypothesis test on coefficients of interest
#Section Save DATA
#Save data
#Model 1:Age
```{r}
#baseline is Age-RangeA
TableBA01<-topTable(fitmm, coef='Pt_Age_RangeB', number=10000 )#compared with Age-RangeA
TableCA01<-topTable(fitmm, coef='Pt_Age_RangeC', number=10000 )
TableDA01<-topTable(fitmm, coef='Pt_Age_RangeD', number=10000 )
TableEA01<-topTable(fitmm, coef='Pt_Age_RangeE', number=10000 )

#Take gene names
Id<- data.frame(colnames)


TableBA01$model<- "AgeB-AgeA"
TableCA01$model<- "AgeC-AgeA"
TableDA01$model<- "AgeD-AgeA"
TableEA01$model<- "AgeE-AgeA"


TableFV01<-rbind(TableBA01,TableCA01,TableDA01, TableEA01)
getwd() #Users/nasimbararpour/Library/CloudStorage/Box-Box/NasimBararpour/hPOP/HPOP Call/SlenoLab/Plasma_PFP
```


#Save data
```{r}
#Log T

write.csv(TableFV01, "DREAM_DESEQ/Age/PlsPFP_baseLineContrast.csv")
```
#Contrast

#METHOD 2
#Contrast matrix
#Model 1: Age
```{r}
form <- ~ 0 + Pt_Age_Range +  (1|Conf_Site_Code) +(1|Sex)+(1|Eth)+ (1|bmi_category)

L = makeContrastsDream(form, meta,
                       contrasts = c("Pt_Age_RangeB - Pt_Age_RangeA",
                                     "Pt_Age_RangeC - Pt_Age_RangeB",
                                     "Pt_Age_RangeD - Pt_Age_RangeC",
                                     "Pt_Age_RangeE - Pt_Age_RangeD"))

# Visualize contrast matrix
plotContrasts(L)
```
# fit dream model with contrasts
```{r}
fit = dream(vobjDream, form, meta, L)
fit = eBayes(fit)
# extract results from first contrast
topTable( fit, coef= "Pt_Age_RangeB - Pt_Age_RangeA", number=3)
head(fit$coefficients, n=4)
```
#extract variable
```{r}
TableContrast1<- topTable(fit, coef= "Pt_Age_RangeB - Pt_Age_RangeA", number=10000 )
TableContrast2<- topTable(fit, coef= "Pt_Age_RangeC - Pt_Age_RangeB", number=10000 )
TableContrast3<- topTable(fit, coef= "Pt_Age_RangeD - Pt_Age_RangeC", number=10000 )
TableContrast4<- topTable(fit, coef= "Pt_Age_RangeE - Pt_Age_RangeD", number=10000 )
```

#Merge row names with Genome ID and save
```{r}
#Take gene names
Id<- data.frame(colnames)


TableContrast1$model<- "Pt_Age_RangeB - Pt_Age_RangeA"
TableContrast2$model<- "Pt_Age_RangeC - Pt_Age_RangeB"
TableContrast3$model<- "Pt_Age_RangeD - Pt_Age_RangeC"
TableContrast4$model<- "Pt_Age_RangeE - Pt_Age_RangeD"


TableContrastFV<- rbind(TableContrast1, TableContrast2,
                        TableContrast3,TableContrast4)
head(TableContrastFV, n=5)
```

#save data
```{r}
#contrast
write.csv(TableContrastFV,"DREAM_DESEQ/Age/PlsPFP_TableContrastFV.csv")

```

#variancePartition plot
A variancePartition analysis can indicate important variables that should be included as fixed or random effects in the dream analysis.
```{r}
# Note: this could be run with either vobj from voom()
# or vobjDream from voomWithDreamWeights()
# The resuylts are similar
form =  ~ (1|Pt_Age_Range) +  (1|Conf_Site_Code) +(1|Sex)+(1|Eth)+ (1|bmi_category)

vp = fitExtractVarPartModel(vobjDream, form, meta)

plotVarPart( sortCols(vp))
```

#correlation of factors
#Standard application
```{r}
form =  ~ Pt_Age_Range +  Conf_Site_Code +Sex+ Eth+ bmi_category
# Compute Canonical Correlation Analysis (CCA)
# between all pairs of variables
# returns absolute correlation value
C = canCorPairs(form, meta)
# Plot correlation matrix
plotCorrMatrix( C )
```
#Model 2: site as fix effect and all other factor as randome effect
```{r}
library(ggplot2)
form <- ~ Conf_Site_Code+ (1|Pt_Age_Range) +(1|Sex)+(1|Eth)+ (1|bmi_category)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, meta, BPPARAM=param, plot = TRUE )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream(vobjDream, form, meta )
fitmm = eBayes(fitmm)

topTable( fitmm, coef= "Conf_Site_CodeB.T", number=3)

topTable( fitmm, coef= "Conf_Site_CodeC.D", number=3)

head(fitmm$design, 3)
```
#Save Data
#Model 2: Conf_Site_Code
```{r}
Table1<-topTable(fitmm, coef='Conf_Site_CodeB.T', number=10000 )#contast with baseline which is first conference site
Table2<-topTable(fitmm, coef='Conf_Site_CodeC.D', number=10000 )
Table3<-topTable(fitmm, coef='Conf_Site_CodeD.O', number=10000 )
Table4<-topTable(fitmm, coef='Conf_Site_CodeE.C', number=10000 )


Table1$model<- "SiteT.2-SiteB.1"
Table2$model<- "SiteD.3-SiteB.1"
Table3$model<- "SiteO.4-SiteB.1"
Table4$model<- "SiteC.5-SiteB.1"


TableFV<-rbind(Table1,Table2,Table3,Table4)
head(TableFV, n=10)
```


#save data
```{r}

#Log T
write.csv(TableFV, "DREAM_DESEQ/Site/PlsPFP_ConSitebaseLineContrast.csv")


```
####Site
#Contrast
#Model 2: Conf_Site_Code
```{r}
form <- ~ 0+ Conf_Site_Code+  (1|Pt_Age_Range) +(1|Sex)+(1|Eth)+ (1|bmi_category)

L = makeContrastsDream( form, meta,
                        contrasts = c("Conf_Site_CodeB.T - Conf_Site_CodeA.B",
                                      "Conf_Site_CodeC.D - Conf_Site_CodeB.T",
                                      "Conf_Site_CodeD.O - Conf_Site_CodeC.D",
                                      "Conf_Site_CodeE.C - Conf_Site_CodeD.O"))
# Visualize contrast matrix
plotContrasts(L)
```



# fit dream model with contrasts
```{r}
fit = dream( vobjDream, form,meta, L)
fit = eBayes(fit)
# extract results from first contrast
topTable( fit, coef= "Conf_Site_CodeB.T - Conf_Site_CodeA.B", number=3)
head(fit$coefficients, n=4)
```
#Extract variable
```{r}
TableContrast1<- topTable( fit, coef= "Conf_Site_CodeB.T - Conf_Site_CodeA.B", number=10000 )
TableContrast2<- topTable( fit, coef= "Conf_Site_CodeC.D - Conf_Site_CodeB.T", number=10000 )
TableContrast3<- topTable( fit, coef= "Conf_Site_CodeD.O - Conf_Site_CodeC.D", number=10000 )
TableContrast4<- topTable( fit, coef= "Conf_Site_CodeE.C - Conf_Site_CodeD.O", number=10000 )
```

#Merge row names with Genome ID and save
```{r}


TableContrast1$model<- "siteT.2-siteB.1"
TableContrast2$model<- "siteD.3-siteT.2"
TableContrast3$model<- "siteO.4-siteD.3"
TableContrast4$model<- "siteC.5-siteO.4"

TableContrastFV<- rbind(TableContrast1, TableContrast2,
                        TableContrast3,TableContrast4)
```

#save data
```{r}

write.csv(TableContrastFV,"DREAM_DESEQ/Site/PlsPFP_TableContrst.csv")


```

#Model 3: Bmi as fix effect and all other factor as randome effect
```{r}
library(ggplot2)
form <- ~ bmi_category+ (1|Pt_Age_Range) +(1|Sex)+(1|Eth)+ (1|Conf_Site_Code)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, meta, BPPARAM=param, plot = TRUE )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream(vobjDream, form, meta )
fitmm = eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 10)
```
#save Data
#Model 3: bmi-group
```{r}
TableBA01<-topTable(fitmm, coef='bmi_categoryB', number=10000 )# baseline is normal A
TableCA01<-topTable(fitmm, coef='bmi_categoryC', number=10000 )
TableDA01<-topTable(fitmm, coef='bmi_categoryNN', number=10000 )

TableBA01$model<- "bmi_groupsB-bmi_groupsA"
TableCA01$model<- "bmi_groupsC-bmi_groupsA"
TableDA01$model<- "bmi_groupsNN-bmi_groupsA"


TableFVBMI<-rbind(TableBA01,TableCA01,TableDA01)
```

#save data
```{r}
#Log T

write.csv(TableFVBMI, "DREAM_DESEQ/BMI/PlsPFP_BMIbaseLineContrast.csv")

```


#Contrast
#Model 3:bmi-Group
```{r}
form <- ~ 0+ bmi_category+ (1|Conf_Site_Code) +(1|Sex)+(1|Eth)+ (1|Pt_Age_Range)

L = makeContrastsDream(form, meta,
                       contrasts = c("bmi_categoryB-bmi_categoryA",
                                     "bmi_categoryC - bmi_categoryB",
                                     "bmi_categoryNN - bmi_categoryC"
                       ))

# Visualize contrast matrix
plotContrasts(L)
# fit dream model with contrasts
```

#fit model
```{r}
fit = dream(vobjDream, form, meta, L)
fit = eBayes(fit)
head(fit$coefficients, n=3)
# extract results from first contrast
topTable( fit, coef= "bmi_categoryB - bmi_categoryA", number=3)
```

#Extract variable: Bmi-Group
```{r}
TableContrast1<- topTable( fit, coef= "bmi_categoryB - bmi_categoryA", number=10000 )
TableContrast2<- topTable( fit, coef= "bmi_categoryC - bmi_categoryB", number=10000 )
TableContrast3<- topTable( fit, coef= "bmi_categoryNN - bmi_categoryC", number=10000 )

TableContrast1$model<- "BmiB-BmiA"
TableContrast2$model<- "BmiC-BmiB"
TableContrast3$model<- "BmiNN-BmiC"

TableContrastFVBMI<- rbind(TableContrast1, TableContrast2,
                           TableContrast3)
```


#save data
```{r}
#POS
write.csv(TableContrastFVBMI,"DREAM_DESEQ/BMI/PlsPFPP_BMIContrstFVLogT2.csv")

```

#Model 4: Ethnicity
```{r}
form <- ~Eth +  (1|Pt_Age_Range) +(1|Sex)+(1|bmi_category)+ (1|Conf_Site_Code)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, meta, BPPARAM=param, plot = TRUE )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream(vobjDream, form, meta )
fitmm = eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 3)
```
#Save Data
#Model 4:Ethnecity
```{r}
TableBA01<-topTable(fitmm, coef='EthC', number=10000 )#Asian are baseline Ethnecity
TableCA01<-topTable(fitmm, coef='EthI', number=10000 )
TableDA01<-topTable(fitmm, coef='EthO', number=10000 )


TableBA01$model<- "EthCau-EthAsian"
TableCA01$model<- "EthIn-EthAsian"
TableDA01$model<- "EthOther-EthAsian"


TableFVEth<-rbind(TableBA01,TableCA01,TableDA01)

head(TableFVEth, n=6)
```




#save data
```{r}
#Log T
write.csv(TableFVEth, "DREAM_DESEQ/Ethnic/PlsPFPLogEthTableFV.csv")

```
#Contrast
#Method 4: Ethnicity
```{r}
form <- ~ 0+ Eth+(1|Conf_Site_Code) +(1|Sex)+(1|bmi_category)+ (1|Pt_Age_Range)

L = makeContrastsDream( form, meta,
                        contrasts = c("EthC-EthA",
                                      "EthI-EthC",
                                      "EthO-EthI" ))

# Visualize contrast matrix
plotContrasts(L)
# fit dream model with contrasts
```

#fit model
```{r}
fit = dream(vobjDream, form, meta, L)
fit = eBayes(fit)
head(fit$coefficients, n=4)
# extract results from first contrast
topTable( fit, coef= "EthI - EthC", number=3)

TableContrast1<- topTable( fit, coef= "EthC - EthA", number=10000 )
TableContrast2<- topTable( fit, coef= "EthI - EthC", number=10000 )
TableContrast3<- topTable( fit, coef= "EthO - EthI", number=10000 )

```
#Merge row names with Genome ID and save
```{r}

TableContrast1$model<- "EthC - EthA"
TableContrast2$model<- "EthI - EthC"
TableContrast3$model<- "EthO - EthI"

TableContrastFVEth<- rbind(TableContrast1, TableContrast2,
                           TableContrast3)

# Get results of hypothesis test on coefficients of interest
```

#save data
```{r}
##LogT
#POS
write.csv(TableContrastFVEth,"DREAM_DESEQ/Ethnic/PlsPFP_EthContrstFVLogT.csv")

```

#Model 5: Sex
```{r}
form <- ~Sex + (1|Pt_Age_Range) +(1|Eth)+(1|bmi_category)+ (1|Conf_Site_Code)
# estimate weights using linear mixed model of dream
vobjDream = voomWithDreamWeights(geneExpr, form, meta, BPPARAM=param, plot = TRUE )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test
fitmm = dream(vobjDream, form, meta )
fitmm = eBayes(fitmm)

# Examine design matrix
head(fitmm$design, 3)
```
#Save Data
#Model 4:Gender

```{r}
TableBA01<-topTable(fitmm, coef='SexMale', number=10000 )



TableBA01$model<- "Male -Female"
head(TableBA01, n=10)
```


```{r}
#Log T

#POS
write.csv(TableBA01, "DREAM_DESEQ/Sex/PlsPFPLogTGenderTableFV.csv")

```


##Summerize data
#1 Table POS
```{r}
##Summarized data within each factor data compared with baseline

rm(list=ls())

```
#setwd("~/Library/CloudStorage/Box-Box/NasimBararpour/hPOP/HPOP Call/SlenoLab/Plasma_PFP/DREAM_DESEQ")
#Age
```{r}
dAgebase<- read.csv("Age/PlsPFP_baseLineContrast.csv")
#Age<- TableFV.Age %>% pivot_wider(names_from = model, values_from = P.Value)

dAgecon<- read.csv( "Age/PlsPFP_baseLineContrast.csv")


TableFV.Age<- rbind(dAgebase, dAgecon)
write.csv(TableFV.Age, "Age/PlsAge_baseLineContrast.csv")



#BMI
#baseline
dbmibase<- read.csv("BMI/PlsPFP_BMIbaseLineContrast.csv")

#contrast
dbmicon<- read.csv("BMI/PlsPFPP_BMIContrstFVLogT2.csv")


TableFV.Bmi<- rbind(dbmibase, dbmicon)
write.csv(TableFV.Bmi, "BMI/PlsPFPBMI_contrastbaseTable.csv")


```

#Site
```{r}
#baseline
dSitebase<- read.csv( "Site/PlsPFP_ConSitebaseLineContrast.csv")


#contrast

dSitecon<- read.csv( "Site/PlsPFP_TableContrst.csv")

TableFV.Site<- rbind(dSitePOS, dSiteNEG)
write.csv(TableFV.Site,"Conf_Site/PlsSite_ContrastTable.csv")
```

#Eth
```{r}
#baseline
dEthFVPos<- read.csv( "Ethnicity/PlsPos_EthLogEthTableFV.csv")

dEthFVNEG<- read.csv( "Ethnicity/PlsNEG_EthLogEthTableFV.csv")

dEthFVPos$mode<- "pos"
dEthFVNEG$mode<- "neg"

TableFV.Eth<- rbind(dEthFVPos, dEthFVNEG)
write.csv(TableFV.Eth,"Ethnicity/PlsEth_ baseLineContrast.csv")

#Contrast
dEthFVPos<- read.csv( "Ethnicity/PlsPos_EthContrstFVLogT.csv")

dEthFVNEG<- read.csv( "Ethnicity/PlsNEG_EthContrstFVLogT.csv")

dEthFVPos$mode<- "pos"
dEthFVNEG$mode<- "neg"

TableFV.Eth<- rbind(dEthFVPos, dEthFVNEG)
write.csv(TableFV.Eth,"Ethnicity/lsEth_ContrastTable")



```

#GENDER
```{r}
dSexFVPos<- read.csv("Sex/PlsPos_GenderTableFV.csv")
dSexFVNEG<- read.csv("Sex/PlsNeg_LogTGenderTableFV.csv")

dSexFVPos$mode<- "pos"
dSexFVNEG$mode<- "neg"

TableFV.Sex<- rbind(dSexFVPos, dSexFVNEG)
write.csv(TableFV.Sex,"Sex/PlsSex_ baseLineContrast.csv")
```

```{r}
#read data
#age
rm(list = ls())
agebaseline<- read.csv("Age/PlsAge_baseLineContrast.csv")
agecontrast<- read.csv("Age/PlsAge_ContrastTable.csv")

#bmi
bmibaseline<- read.csv("BMI/PlsBMI_baseLineContrast.csv")
bmicontrast<- read.csv("BMI/PlsBMI_ContrastTable.csv")

#Eth
Ethbaseline<- read.csv("Ethnicity/PlsEth_ baseLineContrast.csv")
Ethcontrast<- read.csv("Ethnicity/PlsEth_ContrastTable.csv")

#Conf_Site
sitebaseline<- read.csv("Conf_Site/PlsSite_ baseLineContrast.csv")
sitecontrast<- read.csv("Conf_Site/PlsSite_ContrastTable.csv")

#SEX
sexbaseline<- read.csv("Sex/PlsSex_ baseLineContrast.csv")

```



#contrast stairpwise measurements
```{r}
rm(list = ls())
```

#Significant
```{r}
#baseline
Sig.Age<- agebaseline %>% pivot_wider(names_from = model, values_from = P.Value)

Sig.Bmi<-bmibaseline %>% pivot_wider(names_from = model, values_from = P.Value)

Sig.Site<-sitebaseline %>% pivot_wider(names_from = model, values_from = P.Value)
Sig.Eth<- Ethbaseline %>% pivot_wider(names_from = model, values_from = P.Value)

#contrast
Sig.Age<- agecontrast %>% pivot_wider(names_from = model, values_from = P.Value)

Sig.Bmi<-bmicontrast %>% pivot_wider(names_from = model, values_from = P.Value)
Sig.Site<-sitecontrast %>% pivot_wider(names_from = model, values_from = P.Value)
Sig.Eth<- Ethcontrast %>% pivot_wider(names_from = model, values_from = P.Value)
```

#Upset Plot
```{r}
library(UpSetR)
## cast multiple value columns
```

#age
```{r}
Sig.Age[is.na(Sig.Age)] <- 0.9
#for baseline
Age<-Sig.Age[11:14]#select columns in Sig.Age
#for contrast
Age<-Sig.Age[11:14]


test_up <- Age
test_up[Age > 0.05] <- 0
test_up[Age <= 0.05] <- 1
```
#Plot
#upset plot
```{r}
#dataframe
df<- data.frame(test_up)
upset(df, main.bar.color = "black")
```
#Bmi
```{r}
Sig.Bmi[is.na(Sig.Bmi)] <- 0.9
#bmi for baseline
Bmi<-Sig.Bmi[11:13]
#bmi for contrast
Bmi<-Sig.Bmi[10:12]

test_up <- Bmi
test_up[Bmi > 0.05] <- 0
test_up[Bmi <= 0.05] <- 1
```
#Plot Bmi
#upset plot
```{r}
#dataframe
df<- data.frame(test_up)

upset(df, main.bar.color = "black")
```

#Eth
```{r}
Sig.Eth[is.na(Sig.Eth)] <- 0.9
#for baseline
Eth<-Sig.Eth[11:13]
#for contrast
Eth<-Sig.Eth[11:13]


test_up <- Eth
test_up[Eth > 0.05] <- 0
test_up[Eth <= 0.05] <- 1
```

#Plot Eth
#upset plot
```{r}
#dataframe
df<- data.frame(test_up)

upset(df, main.bar.color = "black")
```
#Site
```{r}
Sig.Site[is.na(Sig.Site)] <- 0.9
#for baseline
Site<-Sig.Site[11:14]
#for contrast
Site<-Sig.Site[11:14]

test_up <- Site
test_up[Site > 0.05] <- 0
test_up[Site <= 0.05] <- 1

#colnames(test_up)<- c("site2-site1", "site3-site2", "site4-site3", "site5-site4")

```
#Plot Site
#upset plot
```{r}
#dataframe
df<- data.frame(test_up)

upset(df, main.bar.color = "black")
```

# sex

#Sex factor visualize as barplot
```{r}
sexFCA.increasing <- sexbaseline[intersect(which(sexbaseline$logFC > 0), which(sexbaseline$P.Value <  0.05)),]  # 21

sexFCA.dec <- sexbaseline[intersect(which(sexbaseline$logFC < 0), which(sexbaseline$P.Value <  0.05)),]  # 5


sexFCA.increasing$Log2FC<-1
sexFCA.increasing$status<- sexFCA.increasing$Log2FC
sexFCA.increasing$status<- "inc"

sexFCA.dec$Log2FC<--1
sexFCA.dec$status<- sexFCA.dec$Log2FC
sexFCA.dec$status<- "dec"

df<- rbind(sexFCA.dec, sexFCA.increasing)
sexFCA.dec$logFC
```

#plot
```{r}
library(ggpubr)


A<- ggplot(df, aes(x=status,y=Log2FC, fill=as.factor(status) )) +
              geom_bar(stat="identity") +
              scale_fill_manual(values =c( "blue","red4")) +
              theme(legend.position="right")+ labs(x = "Status in Male")


##In case separately done
B<-  ggplot(sexFCA.dec, aes(x=status, fill=as.factor(status) )) +
              geom_bar( ) +
              scale_fill_manual(values = "blue") +
              theme(legend.position="right")+ labs(x = "Status in Male")

bar.dec03<- B + scale_y_reverse() #Oposit direction
```


```{r}
library(patchwork)

plot =
              A + patchwork::plot_layout(ncol = 3, nrow= 1) #widths = c(0,0))#next together
plot
```


```{r}
p.bar<-ggarrange(A, A,bar.dec03, bar.dec03,
                 labels =c("","","",""),
                 ncol = 3, nrow = 2)%>%
              ggexport(filename = "test.pdf")
p.bar

getwd()

# Export to pdf
ggarrange(bxp, dp, dens, ncol = 2) %>%
              ggexport(filename = "barplotSex.pdf")
```


