# mediation analysis
# treatment -> mediation -> outcome
# microbiome (dfmetaGen) -> metabolomics (dfBioC01) -> RNA (trans) (select ethnicity different ones)
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(tidyr)
require(igraph)
require(Hmisc)
require(npreg)
require(readxl)
require(clusterProfiler)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(circlize)
require(ComplexHeatmap)
require(mediation)
#
require(foreach)
require(doSNOW)
cl<-makeSOCKcluster(8)
registerDoSNOW(cl)
#
set.seed(1)
folddir="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/hpop/"
load(paste0(folddir,"data/HPOP_MutiModalDataEx_Vsep.RData"))
omicslist=ls()
omicslist=setdiff(omicslist,c("folddir","cl"))
resdir=paste0(folddir,"mediation_ana/res/");
setwd(resdir)
# omics list
listomics=list()
for(omics in omicslist){
              listomics[[omics]]=t(get(omics))
}
#
tab_clinic=read.table(paste0(folddir,"data/PhenoAgeDiff.csv"),header=TRUE,sep=",")
clcfeatures=c("bmi","ageY","albumin","fasting_creatinine","fasting_glucose","hscrp","alk_ptase_total","phenoage","deltaAge","Sex")
tab_clinic=tab_clinic[,c("Finalcode",clcfeatures)]
tab_clinic$sex=ifelse(tab_clinic[,"Sex"]=="Male",1,0)
outcome_flist=c("albumin","fasting_creatinine","fasting_glucose","hscrp","alk_ptase_total","phenoage","deltaAge")
covar_cols=c("sex","bmi","ageY")
temptab_cl=tab_clinic[,c("Finalcode",outcome_flist,covar_cols)]
rownames(temptab_cl)=temptab_cl[,1]
temptab_cl=as.matrix(temptab_cl[,-1])
listomics[["clincal"]]=temptab_cl[,outcome_flist]
#
temptab=listomics[["dfmetaGen"]]
listomics[["dfmetaGen"]]=temptab[,colnames(temptab)!=""]
# listomics[["dfmetaGen"]]=scale(listomics[["dfmetaGen"]])
# target id list
id_interest=c("CDH13","TERF1","TERF2","TERF2IP","CYP27A1","SLC27A1","ST6GALNAC1","GALNT10","CD6","PTGDS","L18R1","DNER","CCL8","F12","CD6","TNFRSF9","SHBG","CTRC","HPR","CCL20","GP1BA","SERPINA12","ADM","APOE","CDCP1","CCL4","HAO1","HAVCR1","AFM","MST1","FCGR2B","TNFRSF10B","MBL2")
idens=mapIds(org.Hs.eg.db,keys=id_interest,column="ENSEMBL",keytype="SYMBOL",multiVals="first")
idens=unique(idens[!is.na(idens)])
#
trans_sele=idens
temptab=listomics[["trans"]]
matchmask=sapply(colnames(temptab),function(x) any(str_detect(string=x,pattern=fixed(trans_sele))))
temptab=temptab[,matchmask]
listomics[["trans"]]=temptab
#
omicschtab=data.frame(treatment=c("dfmetaGen"),mediation=c("dfBioC01"),outcome=c("trans"))
# correlation based filter of tri
progress<-function(n) cat(sprintf("task %d is complete\n", n))
opts<-list(progress=progress)
mergtab_sele=c()
mergtablist=list()
pvalcut_trans=0.01
for(rowi in seq(nrow(omicschtab))){
              locvec=unlist(omicschtab[rowi,])
              pairs=list("1"=c(locvec[1:2]),"2"=c(locvec[2:3]),"3"=c(locvec[c(1,3)]))
              asstablist=list()
              for(asslev in names(pairs)){
                            omicpair=pairs[[asslev]]
                            omics1=listomics[[omicpair[1]]]
                            omics2=listomics[[omicpair[2]]]
                            mergtab=merge(omics1,omics2,by=0)
                            tabstat_cor=c()
                            omic2names=colnames(omics2)
                            omic1names=colnames(omics1)
                            tabstat_cor_list<-foreach(ind=seq(length(omic1names)),.packages=c("ppcor"),.options.snow=opts)%dopar%{
                                          set.seed(ind)
                                          omics1ele=omic1names[ind]
                                          corvec=rep(NA,times=dim(omics2)[2])
                                          pvec=corvec
                                          for(omics2ele_i in seq(length(omic2names))){
                                                        omics2ele=omic2names[omics2ele_i]
                                                        xvec=mergtab[,omics1ele]
                                                        yvec=mergtab[,omics2ele]
                                                        nonaind=(!is.na(xvec))&(!is.na(yvec))
                                                        if(length(unique(yvec[nonaind]))<=1){
                                                                      next
                                                        }
                                                        corres=cor.test(x=xvec[nonaind],y=yvec[nonaind],method="spearman")
                                                        pvec[omics2ele_i]=corres$p.value
                                                        corvec[omics2ele_i]=corres$estimate
                                          }
                                          data.frame(feat1=rep(omics1ele,times=length(corvec)),feat2=omic2names,corr=corvec,pval=pvec)
                            }
                            tabstat_cor=Reduce("rbind",tabstat_cor_list)
                            tabstat_cor$padj=p.adjust(tabstat_cor$pval,method="fdr")
                            if("trans" %in% omicpair){
                                          pcut=pvalcut_trans
                            }else{
                                          pcut=0.05
                            }
                            seleind=which(c(tabstat_cor$pval<pcut))
                            seleind=seleind[!is.na(seleind)]
                            subtab=tabstat_cor[seleind,c("feat1","feat2","padj")]
                            colnames(subtab)=c(omicpair,paste0("padj",asslev))
                            asstablist[[asslev]]=subtab
              }
              mergtab=merge(merge(asstablist[[1]],asstablist[[2]],by=locvec[2]),asstablist[[3]],by=locvec[c(1,3)])
              mergtab=mergtab[,locvec]
              colnames(mergtab)=c("treatment","mediation","outcome")
              mergtab_sele=rbind(mergtab_sele,mergtab)
              mergtablist[[rowi]]=mergtab
}
save(mergtab_sele,mergtablist,file="start_mediation_list_all_transubset.RData")
# run for each omcis tri
restab=c()
revreslt=c()
for(trii in seq(nrow(omicschtab))){
              mergetab=mergtablist[[trii]]
              omicsvec=unlist(omicschtab[trii,])
              # overlapped samples
              samplist=sapply(listomics[omicsvec],rownames)
              subejects=intersect(Reduce(intersect,samplist),rownames(temptab_cl))
              inputmat=listomics[[omicsvec[1]]][subejects,,drop=FALSE]
              mediationmat=listomics[[omicsvec[2]]][subejects,,drop=FALSE]
              outcomemat=listomics[[omicsvec[3]]][subejects,,drop=FALSE]
              covarmat=temptab_cl[subejects,covar_cols]
              #
              tabstat_list<-foreach(rowi=seq(nrow(mergetab)),.packages=c("mediation"),.options.snow=opts)%dopar%{
                            set.seed(rowi)
                            locrecrd=mergetab[rowi,]
                            locdf=as.data.frame(cbind(data.frame(outcome=outcomemat[,locrecrd[,"outcome"]],mediation=mediationmat[,locrecrd[,"mediation"]],treatment=inputmat[,locrecrd[,"treatment"]]),covarmat))
                            locdf=na.omit(locdf)
                            if(dim(locdf)[1]<=10){
                                          return(c())
                            }
                            mediator_lm=lm(mediation ~ treatment + sex + bmi + ageY,data=locdf)
                            outcome_lm=lm(outcome ~ treatment + mediation + sex + bmi + ageY,data=locdf)
                            md_model=mediate(mediator_lm,outcome_lm,treat="treatment",mediator="mediation")
                            modelsumr=summary(md_model)
                            data.frame(treatomics=omicsvec[1],mediaomics=omicsvec[2],outomics=omicsvec[3],treatment=locrecrd[,"treatment"],mediation=locrecrd[,"mediation"],outcome=locrecrd[,"outcome"],mediationeffect=modelsumr$d0,mediationeffect_p=modelsumr$d0.p,directeffect=modelsumr$z0,directeffect_p=modelsumr$z0.p,proportion=modelsumr$n0,proportion_p=modelsumr$n0.p,step1coef=summary(mediator_lm)$coefficients["treatment","Estimate"],step2coef=summary(outcome_lm)$coefficients["mediation","Estimate"],step1coef_p=summary(mediator_lm)$coefficients["treatment","Pr(>|t|)"],step2coef_p=summary(outcome_lm)$coefficients["mediation","Pr(>|t|)"])
              }
              stattab=Reduce("rbind",tabstat_list)
              # reverse run
              tabstat_list_rev<-foreach(rowi=seq(nrow(mergetab)),.packages=c("mediation"),.options.snow=opts)%dopar%{
                            set.seed(rowi)
                            locrecrd=mergetab[rowi,]
                            locdf=as.data.frame(cbind(data.frame(outcome=outcomemat[,locrecrd[,"outcome"]],mediation=mediationmat[,locrecrd[,"mediation"]],treatment=inputmat[,locrecrd[,"treatment"]]),covarmat))
                            locdf=na.omit(locdf)
                            if(dim(locdf)[1]<=10){
                                          return(c())
                            }
                            mediator_lm=lm(outcome ~ treatment + sex + bmi + ageY,data=locdf)
                            outcome_lm=lm(mediation ~ treatment + outcome + sex + bmi + ageY,data=locdf)
                            md_model=mediate(mediator_lm,outcome_lm,treat="treatment",mediator="outcome")
                            modelsumr=summary(md_model)
                            data.frame(reatomics=omicsvec[1],mediaomics=omicsvec[2],outomics=omicsvec[3],treatment=locrecrd[,"treatment"],mediation=locrecrd[,"mediation"],outcome=locrecrd[,"outcome"],mediationeffect_rev=modelsumr$d0,mediationeffect_rev_p=modelsumr$d0.p)
              }
              stattab_rev=Reduce("rbind",tabstat_list_rev)
              # pvalue adjust and clean up
              outcomelist=unique(stattab[,"outcome"])
              stattab_padj=c()
              for(outcomele in outcomelist){
                            loctab=stattab[stattab[,"outcome"]==outcomele,]
                            stattab_rev_loc=stattab_rev[stattab_rev[,"outcome"]==outcomele,]
                            stattab_rev_loc$mediationeffect_rev_p_fdr=p.adjust(stattab_rev_loc[,"mediationeffect_rev_p"],method="fdr")
                            remtab=stattab_rev_loc[stattab_rev_loc$mediationeffect_rev_p_fdr<0.05,]
                            indvec=c()
                            for(rowi in seq(nrow(remtab))){
                                          locrecrd=remtab[rowi,]
                                          matchind=which(loctab[,"treatment"]==locrecrd[,"treatment"]&loctab[,"mediation"]==locrecrd[,"mediation"]&loctab[,"outcome"]==locrecrd[,"outcome"])
                                          indvec=c(indvec,matchind)
                            }
                            if(length(indvec)>0){
                                          loctab=loctab[-indvec,]
                            }
                            loctab$padj=p.adjust(loctab[,"mediationeffect_p"],method="fdr")
                            stattab_padj=rbind(stattab_padj,loctab)
                            revreslt=rbind(revreslt,stattab_rev_loc)
              }
              restab=rbind(restab,stattab_padj)
}
save(restab,revreslt,file="mediation_res_transubset.RData")
