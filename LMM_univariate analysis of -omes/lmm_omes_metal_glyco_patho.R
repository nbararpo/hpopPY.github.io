
# use variancePartition
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(readxl)
require(lattice)
require(car)
require(dplyr)
require(variancePartition)
require(limma)
require(edgeR)
require(purrr)
param=SnowParam(4, "SOCK", progressbar=TRUE)
#
comp="/Users/yuewu/"
pardir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/hpop/mixed_linear_model/")
resdir=paste0(pardir,"res/")
datadir=paste0(pardir,"data/")
set.seed(1)
setwd(resdir)
#
subfeatureslist=c("Finalcode","Conf_Site","Pt_Age_Range","Sex","bmi_category","Ethnic_1","h.POP.ID.")
featurenames=c("Age","Sex","bmi_category","Eth","Conf_Site")
#
metadf=read_excel(paste0(datadir,"metadata_HPOP_FV4.xlsx"),na=c("","NA"))
metadf[["Finalcode"]]=str_replace_all(metadf[["Finalcode"]],pattern=fixed("-"),replacement=".")
metadf=metadf[rowSums(is.na(metadf[,subfeatureslist]))==0,]#rem na
# feature factorize
f_fac_lev_list=list("Age"=list("A"="20 - 30","B"="30 - 40","C"="40 - 50","D"="50 - 60","E"=">60"),
                    "Sex"=list("A"="Female","B"="Male"),
                    "bmi_category"=list("A"="underweight","B"="normal","C"="overweight","D"="obese"),
                    "Eth"=list("A"="Asian","B"="Indian","C"="Caucasian","D"="Other"))
metadf$Age=metadf$Pt_Age_Range
metadf$Eth=metadf$Ethnic_1
metadf$Eth[metadf$Ethnic_2=="Indian"]="Indian"
for(classname in names(f_fac_lev_list)){
              labellist=f_fac_lev_list[[classname]]
              for(labelname in names(labellist)){
                            metadf[metadf[,classname]==labellist[[labelname]],classname]<-labelname
              }
}
#
omics_files=list(metal="MetalomeData.csv",glyco="HPOP_glycomics_filtered.csv",patho="HPOP_MutiModalDataEx_Vsep_pathupd.csv")
omics_list=list()
for(omicfile_name in names(omics_files)){
              omicfile=omics_files[[omicfile_name]]
              # the omics data
              omics_tab=read.table(paste0(datadir,omicfile),sep=",",header=TRUE,row.names=1)
              # select intersect samples between the two tables
              idoverlap=intersect(colnames(omics_tab),metadf[["Finalcode"]])
              omics_tab=omics_tab[,idoverlap]
              reordind=sapply(idoverlap,function(x) which(metadf[["Finalcode"]]==x))
              metadf2=metadf[reordind,]
              metadf2=as.data.frame(metadf2)
              rownames(metadf2)=metadf2[["Finalcode"]]
              # transformation of the omics
              # matExpr=DGEList(omics_tab)
              # matExpr=calcNormFactors(matExpr)
              omics_tab=as.matrix(omics_tab)
              if(omicfile_name=="metal"){
                            omics_tab=log10(omics_tab)
              }
              for(feati in seq(dim(omics_tab)[1])){
                            omics_tab[feati,]=(omics_tab[feati,]-mean(omics_tab[feati,],na.rm=TRUE))/sd(omics_tab[feati,],na.rm=TRUE)
              }
              omics_list[[omicfile_name]]=omics_tab
              #
              table_baseline_list=vector(mode="list")
              table_contr_list=vector(mode="list")
              for(featname in featurenames[-length(featurenames)]){#no conference site in the results
                            formula_mlm_str=paste0("~")
                            for(equterm in featurenames){
                                          if(equterm!=featurenames[1]){
                                                        formula_mlm_str=paste0(formula_mlm_str,"+")
                                          }
                                          if(equterm==featname){
                                                        formula_mlm_str=paste0(formula_mlm_str,equterm)
                                          }else{
                                                        formula_mlm_str=paste0(formula_mlm_str,"(1|",equterm,")")
                                          }
                            }
                            formula_mlm=as.formula(formula_mlm_str)
                            # vobjDream=voomWithDreamWeights(matExpr,formula_mlm,metadf,BPPARAM=param)#transformation
                            namesord=rownames(omics_tab)
                            # baseline
                            fitmm=dream(omics_tab,formula_mlm,metadf2)#regression
                            fitmm=eBayes(fitmm)
                            alllevels=sort(unique(metadf2[,featname]))
                            #
                            listtab=list()
                            levels_comp=alllevels[-1]
                            for(levelele in levels_comp){
                                          table1<-topTable(fitmm,coef=paste0(featname,levelele),number=1000)
                                          identname=paste0(featname,levelele," ",f_fac_lev_list[[featname]][[levelele]]," vs A ",f_fac_lev_list[[featname]][["A"]])
                                          colnames(table1)=paste0(colnames(table1)," ",identname)
                                          listtab[[levelele]]=table1[namesord,]
                            }
                            table_baseline_list[[featname]]=Reduce(cbind,listtab)
                            # contrast
                            # build the contrast matrix
                            cont_mat=c()
                            cont_mat_name=c()
                            for(levelele_i in seq(from=2,to=length(alllevels))){
                                          addstr=paste0(featname,alllevels[levelele_i]," - ",featname,alllevels[levelele_i-1])
                                          cont_mat=c(cont_mat,addstr)
                                          cont_mat_name=c(cont_mat_name,paste0(featname,".",as.character(levelele_i-1)))
                            }
                            names(cont_mat)=cont_mat_name
                            formula_mlm_str %>% str_replace_all(string=.,pattern="~",replacement="~ 0 + ") %>% as.formula() -> formula_mlm
                            L=makeContrastsDream(formula_mlm,metadf2,contrasts=cont_mat)
                            # plotContrasts(L)
                            fitmm=dream(omics_tab,formula_mlm,metadf2,L)
                            fitmm=eBayes(fitmm)
                            listtab=list()
                            for(levelele_i in seq(from=2,to=length(alllevels))){
                                          levelele=alllevels[levelele_i]
                                          table1<-topTable(fitmm,coef=paste0(featname,".",as.character(levelele_i-1)),number=1000)
                                          identname=paste0(featname," ",f_fac_lev_list[[featname]][[levelele_i]]," vs ",f_fac_lev_list[[featname]][[levelele_i-1]])
                                          colnames(table1)=paste0(colnames(table1)," ",identname)
                                          listtab[[levelele]]=table1[namesord,]
                            }
                            table_contr_list[[featname]]=Reduce(cbind,listtab)
              }
              table_baseline=Reduce(cbind,table_baseline_list)
              table_contr=Reduce(cbind,table_contr_list)
              save(table_baseline,file=paste0(omicfile_name,"_baseline.RData"))
              save(table_contr,file=paste0(omicfile_name,"_contrast.RData"))
              save(metadf2,omics_tab,file=paste0(omicfile_name,"_data.RData"))
              #variance sparation plot
              formula_mlm<- ~ (1|Age)+(1|Conf_Site)+(1|Sex)+(1|Eth)+(1|bmi_category)###1
              vp=fitExtractVarPartModel(omics_tab,formula_mlm,metadf2)
              p=plotVarPart(sortCols(vp))##violinplot
              ggsave(filename=paste0(omicfile_name,"varpartplot.pdf"),plot=p)
}
# save the omics matrix
omics_list2=omics_list
for(omic in names(omics_list2)){
              arra=omics_list2[[omic]]
              colnames(arra)<-str_replace_all(string=colnames(arra),pattern="\\.",replacement="-")
              omics_list2[[omic]]=arra
}
save(omics_list2,file="metal_glyco_patho_omics_arra.RData")
# kmeans clustering of samples in metalome
res_kmeans_dir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/hpop/kmeans/")
setwd(res_kmeans_dir)
wss<-function(k){
              return(kmeans(loctab,k,nstart=20)$tot.withinss)
}
for(omicsgr in names(omics_list)){
              loctab=t(omics_list[[omicsgr]])
              # k2<-kmeans(loctab,centers=2,nstart=10)
              pca_res=prcomp(loctab,center=FALSE)
              scores=pca_res$x
              # summary(pca_res)
              pdf(paste0(omicsgr,"_pca.pdf"))
              plot(scores[,c(1,2)])#,col=k2$cluster
              dev.off()
              #
              k.values<-1:30
              wss_values<-map_dbl(k.values,wss)
              pdf(paste0(omicsgr,"_k_seq.pdf"))
              plot(k.values,wss_values,xlab="Number of clusters K",ylab="Total within-clusters sum of squares")
              dev.off()
}
# plot with the selected features
# plot metalome colored by Sex, Age, Ethnicity,
# plot glycome colored by Ethnicity, age, bmi
colorgrouplist=list(metal=c("Sex","Age","Eth"),glyco=c("Eth","Age","bmi_category"),patho=c("Sex","Age","Eth"))
for(omicsgr in names(omics_list)){
              loctab=t(omics_list[[omicsgr]])
              pca_res=prcomp(loctab,center=FALSE)
              scores=pca_res$x
              colorgroups=colorgrouplist[[omicsgr]]
              for(grp in colorgroups){
                            coldf=as.data.frame(metadf[,c(grp,"Finalcode")])
                            rownames(coldf)=coldf[,"Finalcode"]
                            colvec=coldf[rownames(scores),grp]
                            plotdf=data.frame(pc1=scores[,1],pc2=scores[,2],classes=as.factor(unlist(f_fac_lev_list[[grp]][colvec])))
                            p<-ggplot(plotdf,aes(x=pc1,y=pc2,color=classes))+geom_point()+ggtitle("pca scatter plot")+theme(strip.background=element_blank(),strip.text.x=element_blank())#
                            ggsave(paste0(omicsgr,grp,"_pca.pdf"))
              }
}
# table_baseline_pval=table_baseline[,str_which(string=colnames(table_baseline),pattern=fixed("adj.P.Val"))]
# ind=which(table_baseline_pval<0.05,arr.ind=TRUE)
# colnames(table_baseline_pval)[ind[,2]]
# table_baseline_pval[ind]
setwd(resdir)
# reformat mixed linear model results into tables
for(omic in names(omics_files)){
              load(paste0(omic,"_contrast.RData"))
              allcols=colnames(table_contr)
              resdf=c()
              for(feat in names(f_fac_lev_list)){
                            alllevels=f_fac_lev_list[[feat]]
                            for(levelele_i in seq(from=2,to=length(alllevels))){
                                          contrlev=paste0(alllevels[[levelele_i]]," vs ",alllevels[[levelele_i-1]])
                                          selecols=allcols[str_which(string=allcols,pattern=paste0(" ",feat," ",contrlev))]
                                          loctab=table_contr[,selecols]
                                          colnames(loctab)<-colnames(loctab)%>%str_extract_all(string=.,pattern="^[^\\s]+\\s+")%>%str_replace_all(string=.,pattern="\\s+",replacement="")
                                          loctab=loctab[loctab[,"adj.P.Val"]<0.05,]
                                          nele=dim(loctab)[1]
                                          if(nele==0){
                                                        next
                                          }
                                          featnames=rownames(loctab)
                                          tempdf=data.frame(Feature=featnames,"Feature identifier"=featnames,logFC=loctab[,"logFC"],AveExpr=loctab[,"AveExpr"],t=loctab[,"t"],"P.Value"=loctab[,"P.Value"],"adj.P.Val"=loctab[,"adj.P.Val"],B=loctab[,"B"],"z.std"=loctab[,"z.std"],Contrst=rep(feat,times=nele),model=rep(contrlev,times=nele))
                                          resdf=rbind(resdf,tempdf)
                            }
              }
              write.table(resdf,file=paste0(omic,"res_tab.txt"))
}
