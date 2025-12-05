setwd('..')
library(dplyr)
library(ggpubr)
library(emmeans)
library(rstatix)
library(multcomp)
library(RColorBrewer)
source('src/functions.R')

##### Phenotypic age #########
file1 <- 'doc/metadata_HPOP_FV4.xlsx' # Get Sex, BMI, clinical params, etc,
file2 <- 'doc/internationalHPOPmetadata_28June2021.xlsx' # Get Age etc.
file3 <- 'doc/localHPOPmetadata_28June2021.xlsx' # Get Age etc.
df <- getPhenoAge(file1, file2, file3, F)
df <- calcPhenoAge(df, F)
df$deltaAge <- df$phenoage - df$ageY
ggplot(data= df, aes(x= ageY, y = phenoage, color = Ethnic)) +
              geom_point(size = 1) +
              geom_smooth(method = lm, se=FALSE)+
              geom_abline(slope=1, intercept = 0, linetype = 'dashed') +
              theme_bw() + ylab('Phenotypic age (years)') + xlab('Chronological age (years)')+
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom")
ggplot(data= df, aes(x=Ethnic, y=deltaAge, color=Sex)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(), size=0.5) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom")

##### Phenotypic age Imputed #########
df <- getPhenoAge(file1, file2, file3, T)
df <- calcPhenoAge(df, T)
if(FALSE){
              meta3 <- openxlsx::read.xlsx('data/metadata_HPOP_FV5.xlsx', sheet = 3)
              # Find missing country. Some hPOP Ids have multiple entry so copy country from there.
              ids <- meta3$h.POP.ID.[!is.na(meta3$Country)]
              meta3 <- meta3[meta3$h.POP.ID. %in% ids,]
              ids2 <- unique(meta3[is.na(meta3$Country),'h.POP.ID.'])
              temp <- unique(meta3[meta3$h.POP.ID. %in% ids2, c('h.POP.ID.', 'Country'),])
              temp <- temp[!is.na(temp$Country),]
              for(iii in ids2){
                            meta3$Country[meta3$h.POP.ID. == iii] <- temp$Country[temp$h.POP.ID. == iii]
              }
              meta3 <- unique(meta3[, c('h.POP.ID.','Country')])
              df <- merge(df, meta3, by= 'h.POP.ID.', all.x = T)
              df$Country2 <- df$Country
              europe <- c('Belgium','Finland','France','Ireland','Malta','Netherlands','Poland','Switzerland','Spain')
              eurpoe2 <- c('Austria','Denmark', 'Germany', 'Sweden','Uppsala','Russia')
              usca <- c('Canada','USA')
              eastAsia <- c('Taiwan','South Korea','KOREA','China','Singapore')
              df$Country2[df$Country %in% c(europe, eurpoe2)] <- 'Europe'
              df$Country2[df$Country %in% eastAsia] <- 'EastAsia'
              df$Country2[df$Country %in% usca] <- 'US Canada'
              df$Country3 <- df$Country2
              df$Country3[which(df$Country3 %in% c('EastAsia','Japan'))] <- 'Asia'
}

df$deltaAge <- df$phenoage - df$ageY
ggplot(data= df, aes(x= ageY, y = phenoage, color = Ethnic)) +
              geom_point(size = 0.5) +
              geom_smooth(method = lm, se=FALSE)+
              geom_abline(slope=1, intercept = 0, linetype = 'dashed') +
              theme_bw() + ylab('Phenotypic age (years)') + xlab('Chronological age (years)')+
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom")
ggplot(data= df, aes(x= ageY, y = deltaAge, color = Ethnic)) +
              geom_point(size = 0.5) +
              geom_smooth(method = lm, se=FALSE)+
              geom_abline(slope=1, intercept = 0, linetype = 'dashed') +
              theme_bw() + ylab('Phenotypic  - Chronological age(yrs)') +
              xlab('Chronological age (years)')+
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom")

##### Effect of Sex on Ethnicity ###########################
stat.test <- df %>% group_by(Ethnic) %>%
              t_test(deltaAge ~ Sex)
stat.test$statistic <- c(2.4489, 1.4847, 2.116, 0.6877)
stat.test$df <- NULL
stat.test$p <- c(0.0087,0.00017,0.109, 0.55)
stat.test$p <- round(stat.test$p, 4)
stat.test <- stat.test %>% add_xy_position(x = "Ethnic", dodge = 0.8)

ggplot(data= df, aes(x=Ethnic, y=deltaAge, color=Sex)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(), size=0.5) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              stat_pvalue_manual(data=stat.test, label = "p") +
              coord_cartesian(ylim=c(NA,13)) +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom")
# 580 345

# We have absolute measurement hence batch effect is not needed.
# We need to remove age specific effect.
# We need to remove b/w ethnicity effects.
mod <- model.matrix(~ ageY + Sex*Ethnic, data=df)
fit <- lm(deltaAge ~ ageY + Sex*Ethnic, df)
x <- anova(fit)
x
m0 <- data.frame(molecule='all', Sex=x$`Pr(>F)`[2], Ethnic= x$`Pr(>F)`[3], sexEth=x$`Pr(>F)`[4])
contrast(emmeans(fit, ~ Sex), method = "pairwise")
contrast(emmeans(fit, ~ Ethnic), method = "pairwise")
# Asian Male.vs.Female
K <- matrix(c(0, 0, 1, 0,0,0, 0,0,0), 1)
summary(glht(fit, linfct = K)) # 0.00871
# Caucassian Male.vs.Female
K <- matrix(c(0, 0, 1, 0,0,0, 1,0,0), 1)
summary(glht(fit, linfct = K)) # 0.000169
# Indian Male.vs.Female
K <- matrix(c(0, 0, 1, 0,0,0, 0,1,0), 1)
summary(glht(fit, linfct = K)) # 0.109
# Other Male.vs.Female
K <- matrix(c(0, 0, 1, 0,0,0, 0,0,1), 1)
summary(glht(fit, linfct = K)) # 0.55


ethnic_levels <- levels(df$Ethnic)
sex_levels <- levels(df$Sex)
contrast_matrix <- contrMat(model.matrix(~ ageY + Sex * Ethnic),
                            list(Sex = sex_levels, Ethnic = ethnic_levels))
contrast_matrix <- matrix(0, nrow = length(ethnic_levels), ncol = length(sex_levels))
rownames(contrast_matrix) <- ethnic_levels
colnames(contrast_matrix) <- sex_levels
contrast_matrix[, "Male"] <- -1
contrast_matrix[, "Female"] <- 1

K1<-glht(fit,mcp(Sex="Tukey"))$linfct # This will be for Asian group
# ^ How can I marginalize over all groups for contrast?
# `anova` gives p-value, but how can I get effect?
K2<-glht(fit,mcp(Ethnic="Tukey"))$linfct
t <- glht(fit, linfct = K1)
summary(t)

avg_data <- data.frame(
              ageY = mean(df$ageY),  # Average ageY
              Sex = c("Male", "Female"),  # Both levels of Sex
              Ethnic = rep(unique(df$Ethnic),each = 2)  # All unique levels of Ethnic
)
emmeans(fit, specs = ~ Sex)
avg_prediction <- predict(fit, newdata = avg_data)
avg_effects <- emmeans(fit, ~ Sex | Ethnic, data = avg_data)
diff(avg_effects)
contrast_results <- emmeans(fit, specs = ~ Sex * Ethnic)
contrast_male_female <- contrast(contrast_results, method = "pairwise", by = "Sex")
print(contrast_male_female)

##### Effect of location ###########################
# https://stats.oarc.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
filtered_df <- dplyr::filter(df, Ethnic %in% c('Caucasian'), !is.na(Country2)) %>%
              dplyr::group_by(Sex, Ethnic, Country2) %>%
              filter(n() > 1) %>% dplyr::ungroup()
stat.test <- filtered_df %>% dplyr::filter(Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country2)) %>%
              group_by(Sex, Ethnic) %>%
              t_test(deltaAge ~ Country2)
stat.test <- stat.test[1:4, ]
stat.test$Sex <- c('Female','Male','Male','Male')
stat.test$group1 <- c('Europe','Australia','Brazil','Europe')
stat.test$group2 <- c('US Canada','Brazil','Europe','US Canada')
stat.test$statistic <- c(-0.4442,-3.486, 3.559, -1.3572)
stat.test$df <- NULL
stat.test$p <- c(0.69,0.0291, 0.0131, 0.083)
stat.test$p <- round(stat.test$p, 3)
stat.test <- stat.test %>% add_xy_position(x = "Country2", dodge = 0.8)
stat.test$xmin <- c(1.725, 1.7, 1.85, 2.0)
stat.test$xmax <- c(2.275, 1.85, 2.0, 2.3)
stat.test$y.position <- c(8.5,8.5,10.25,12)

df.2 <- df
ggplot(data=dplyr::filter(df, Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country2)),
       aes(x=Ethnic, y=deltaAge, color=Country2)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(dodge.width= 0.75, jitter.width = 0.10), size=0.7) +
              facet_wrap(vars(Sex)) +
              stat_pvalue_manual(data=stat.test, label = "p", tip.length = 0.01, label.size = 2.5) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom") +
              coord_cartesian(ylim=c(NA,13)) + labs(color = 'Location')

ggplot(data=dplyr::filter(df, Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country2), Sex == 'Male'),
       aes(x=Ethnic, y=deltaAge, color=Country2)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(dodge.width= 0.75, jitter.width = 0.10), size=0.7) +
              stat_pvalue_manual(data=stat.test[stat.test$Sex == 'Male',], label = "p", tip.length = 0.01, label.size = 2.5) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom") +
              coord_cartesian(ylim=c(NA,13)) + labs(color = 'Location')

ggplot(data=dplyr::filter(df, Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country2), Sex == 'Female'),
       aes(x=Ethnic, y=deltaAge, color=Country2)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(dodge.width= 0.75, jitter.width = 0.10), size=0.7) +
              stat_pvalue_manual(data=stat.test[stat.test$Sex == 'Female',], label = "p", tip.length = 0.01, label.size = 2.5) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom") +
              coord_cartesian(ylim=c(NA,13)) + labs(color = 'Location')

pdf('phenoAgeLocation_March12.pdf')
ggplot(data=dplyr::filter(df, Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country3)),
       aes(x=Ethnic, y=deltaAge, color=Country3)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(dodge.width= 0.75, jitter.width = 0.10), size=0.7) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom") +
              coord_cartesian(ylim=c(NA,13)) + labs(color = 'Location')

ggplot(data=dplyr::filter(df, Ethnic %in% c('Asian', 'Caucasian'), !is.na(Country3)),
       aes(x=Ethnic, y=deltaAge, color=Country3)) +
              geom_boxplot(outlier.size = 1) +
              geom_point(pch = 21,position = position_jitterdodge(dodge.width= 0.75, jitter.width = 0.10), size=0.7) +
              facet_wrap(vars(Sex)) +
              theme_bw() + ylab('Phenotypic - Chronological age(yrs)') +
              theme(panel.grid.minor.x = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    legend.position = "bottom") +
              coord_cartesian(ylim=c(NA,13)) + labs(color = 'Location')
dev.off()

filtered_df$Sex <- factor(filtered_df$Sex, levels=rev(levels(filtered_df$Sex)))
filtered_df$Country2 <- factor(filtered_df$Country2, levels=c("US Canada","UK","Europe","Brazil","Australia"))
# Need not to include Ethnicity because only one
mod <- model.matrix(~ ageY + Country2*Sex, data=filtered_df)
fit <- lm(deltaAge ~ ageY + Country2*Sex, filtered_df)
x <- anova(fit)
x
summary(fit)$converged
m0 <- data.frame(molecule='all', Sex=x$`Pr(>F)`[2], Country2= x$`Pr(>F)`[3], sexCountry=x$`Pr(>F)`[4])
contrast(emmeans(fit, ~ Sex), method = "pairwise")
contrast(emmeans(fit, ~ Country2), method = "pairwise")

# Female Europe.vs.UK
fit$coefficients <- coef(fit)[1:9]
K <- matrix(c(0, 0, -1,1,0,0, 0, -1,1), 1)
summary(glht(fit, linfct = K)) # -0.5181  0.725
# Female Europe.vs.US Canada
K <- matrix(c(0, 0, 0,1,0,0, 0, 0,1), 1)
summary(glht(fit, linfct = K)) # -0.4442  0.69
# Female UK.vs.US Canada
K <- matrix(c(0, 0, 1,0,0,0, 0, 1,0), 1)
summary(glht(fit, linfct = K)) # 0.07387  0.959
# Male UK.vs.US Canada
K <- matrix(c(0, 0, 1,0,0,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # 0.8675  0.532
# Male Australia.vs.Brazil
K <- matrix(c(0, 0, 0,0,-1,1, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # -3.486  0.0291*
# Male Australia.vs.Europe
K <- matrix(c(0, 0, 0,-1,0,1, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # 0.07265  0.95
# Male Australia.vs.UK
K <- matrix(c(0, 0, -1,0,0,1, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # -2.152  0.194
# Male Australia.vs.US Canada
K <- matrix(c(0, 0, 0,0,0,-1, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # 1.285  0.246
# Male Brazil.vs.Europe
K <- matrix(c(0, 0, 0,-1,1,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # 3.559 0.0131 *
# Male Brazil.vs.UK
K <- matrix(c(0, 0, -1,0,1,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # 1.334 0.465
# Male Brazil.vs.US Canada
K <- matrix(c(0, 0, 0,0,-1,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # -2.202 0.107
# Male Europe.vs.UK
K <- matrix(c(0, 0, -1,1,0,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # -2.225 0.105
# Male Europe.vs.US Canada
K <- matrix(c(0, 0, 0,1,0,0, 0, 0,0), 1)
summary(glht(fit, linfct = K)) # -1.3572 0.083 .

##### Phenotypic age molecualr correlation #########
df <- getPhenoAge(file1, file2, file3, T)
df <- calcPhenoAge(df, T)
df$deltaAge <- df$phenoage - df$ageY
df$temp <- - 0.0336*df$albumin + 0.0095*df$fasting_creatinine +
              0.0195*df$fasting_glucose + 0.0954*log(df$hscrp) + 0.00188*df$alk_ptase_total
fit <- lm(fasting_creatinine ~ ageY + Sex*Ethnic, df)
anova(fit) # Sex pval=7.665e-15
fit <- lm(log(hscrp) ~ ageY + Sex*Ethnic, df)
anova(fit) # Sex pval=0.012, Ethnic pval=0.01158
fit <- lm(alk_ptase_total ~ ageY + Sex*Ethnic, df)
anova(fit) # Ethnic pval=0.0032

##### Phenotypic age effect of each column #########
source('src/functions.R')
df <- getPhenoAge(file1, file2, file3, T)
df <- calcPhenoAge(df, T)
df$deltaAge <- df$phenoage - df$ageY
df$phenoage0 <- calcPhenoAge0(df)
df$deltaAge0 <- df$phenoage0 - df$ageY
df$phenoage.albumin <- calcPhenoAge2(df, 'albumin')
df$phenoage.creatinine <- calcPhenoAge2(df, 'fasting_creatinine')
df$phenoage.glucose <- calcPhenoAge2(df, 'fasting_glucose')
df$phenoage.hscrp <- calcPhenoAge2(df, 'hscrp')
df$phenoage.alkPtase <- calcPhenoAge2(df, 'alk_ptase_total')
df$delta.albumin <- df$phenoage.albumin - df$ageY
df$delta.creatinine <- df$phenoage.creatinine - df$ageY
df$delta.glucose <- df$phenoage.glucose - df$ageY
df$delta.hscrp <- df$phenoage.hscrp - df$ageY
df$delta.alkPtase <- df$phenoage.alkPtase - df$ageY
df$phenoage.albumin3 <- calcPhenoAge3(df, 'albumin')
df$phenoage.creatinine3 <- calcPhenoAge3(df, 'fasting_creatinine')
df$phenoage.glucose3 <- calcPhenoAge3(df, 'fasting_glucose')
df$phenoage.hscrp3 <- calcPhenoAge3(df, 'hscrp')
df$phenoage.alkPtase3 <- calcPhenoAge3(df, 'alk_ptase_total')
df$delta.albumin3 <- df$phenoage.albumin3 - df$ageY
df$delta.creatinine3 <- df$phenoage.creatinine3 - df$ageY
df$delta.glucose3 <- df$phenoage.glucose3 - df$ageY
df$delta.hscrp3 <- df$phenoage.hscrp3 - df$ageY
df$delta.alkPtase3 <- df$phenoage.alkPtase3 - df$ageY

fit <- lm(deltaAge0 ~ ageY + Sex*Ethnic, df)
x <- anova(fit)
coef(fit)

df1 <- singleMolEffect1(df)
m <- reshape2::melt(df1, 'molecule')
m$type <- 'Effect'
m$type[grep('^p', m$variable)] <- 'pval'
ref <- 'All'
myColors <- brewer.pal(length(levels(m$molecule)),"Set2")
names(myColors) <- levels(m$molecule)
myColors[names(myColors)==ref] <- "grey"
ggplot(dplyr::filter(m, type == 'Effect'),
       aes(x=variable, y = -value, fill=molecule)) +
              geom_bar(stat="identity", position=position_dodge()) + theme_bw() +
              ylab('Coef for phenotypic age delta') +
              scale_x_discrete(labels=c("e.AC" = "Caucasian-Asian", "e.AI" = "Indian-Asian",
                                        "e.CI" = "Indian-Caucasian", "e.Sex" = "Male-Female")) +
              scale_fill_manual(name = "Molecule",values = myColors) +
              theme(axis.title.x=element_blank())
# 750 340
# Creatinine and hscrp affect phenotypic age in Caucasian v/s Asian.
# hscrp have highest impact on phenotypic age in Indian v/s Asian.
# Creatinine, hscrp and Alkaline phosphatase affect biologic age in Indian v/s Caucasian.
# Creatining have the strongest impact on male's biological age than female's.
# Albumin and Glucose have comparatively little impace on phenotypic age delta.

m1 <- singleMolEffect(df)
m <- rbind(m0,m1)
m$molecule <- factor(m$molecule, level = m$molecule)
m <- reshape2::melt(m, 'molecule')
ggplot(m, aes(x=variable, y = -log10(value), fill=molecule)) +
              geom_bar(stat="identity", position=position_dodge()) + theme_bw()

###### Nasim #########
df <- getPhenoAge(file1, file2, file3, T)
df <- calcPhenoAge(df, T)
df$deltaAge <- df$phenoage - df$ageY
df2 <- read.table('doc/dclincFactor.csv', header= T, sep=',')
temp <- df[,c('Finalcode','h.POP.ID.','ageY', 'phenoage')]
df2 <- merge(df2, temp, by.x = c('sample','h.POP.ID.'),
             by.y = c('Finalcode','h.POP.ID.'), all.x = T)
write.table(df2, 'doc/dclincFactor_SG.csv', quote=F, sep=',', row.names = F)

###### Model the protein expression #####
phenoDf <- read.table('PhenoAgeDiff.csv', header= T, sep=',')
protDf <- read.table('proteomeExpr.tsv', header=T, sep='\t')
protDf <- protDf[,c('sample', 'UniProt', 'logIntensity', 'Gene')]
phenoDf <- merge(phenoDf, protDf, by.y= 'sample', by.x='Finalcode')

pr <- plyr::ddply(phenoDf, c("UniProt"), getLMM.pheno)


library(lme4)
library(lmerTest)
df2 <- df
df2 <- df2[!is.na(df2$Country2),]
df2$EthCoun <- paste0(df2$Ethnic,'/', df2$Country2)
colnames(df2)[1] <- 'Id'
model.000<- lmer(phenoage ~ ageY + Sex + as.factor(EthCoun)+ bmi+ (1|Id), data=df2,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(model.000)

###### Model the transcript expression #####
phenoDf <- read.table('PhenoAgeDiff.csv', header= T, sep=',')
df <- read.table('transcriptome.tsv', header = T, sep = '\t')
df$transcript <- rownames(df)
df <- reshape2::melt(df, id.vars = 'transcript', variable.name = 'sample', value.name = 'abun')
df$sample <- sub('\\.','-', df$sample)
phenoDf <- merge(phenoDf, df, by.y= 'sample', by.x='Finalcode')
pr <- plyr::ddply(phenoDf, c("transcript"), getlm.pheno)
colnames(pr)[2] <- 'pval'
pr$adj.pval <- p.adjust(pr$pval, "BH")
