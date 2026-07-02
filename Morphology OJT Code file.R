install.packages("readxl")
library(readxl)
##import the data 
MorphologyData=read_excel("C:\\Users\\Mahi2\\Downloads\\Morphology data 12-02-2026.xlsx")
str(MorphologyData)#to check the structure of data ,data is proper imported or not 
head(MorphologyData)#use this command to confirm headers and values are properly place or not 
dim(MorphologyData)
summary(MorphologyData)
MorphologyData$genotype=as.factor(MorphologyData$genotype)
class(MorphologyData$genotype)
MorphologyData$Replication=as.factor(MorphologyData$Replication)
class(MorphologyData$Replication)
##because these are categorical variable 

##check missing values
sum(is.na(MorphologyData)) #no missing values is here
#understanding the dataset 
summary(MorphologyData)
hist(MorphologyData$PH) ##data approximately follow norml distribution
#there is no any abnormality in PH dataset 
boxplot(MorphologyData$PH)
qqnorm(MorphologyData$PH)
qqline(MorphologyData$PH)#normality assumption are satisfied

### RCBD Anova
PH=MorphologyData$PH
genotype=MorphologyData$genotype
Replication=MorphologyData$Replication
model<-aov(PH~genotype+Replication,data=MorphologyData)
##assumption of anova
#1 normality of residual
res<-residuals(model)
shapiro.test(res)
qqnorm(res)
hist(res)
#2homosgeneity of variance(equal variance)
library(car)
leveneTest(PH~genotype,data=MorphologyData)
#3 residual vs fitted plot
plot(model,which = 1)
summary(model)

###genetic variability 
#anova table extract 
anova_tab<-summary(model)[[1]]
anova_tab

##MSG OR MSE 
#MSG=gennotype ki total variation 
#MSE=Environmental/error variation 
MSG=anova_tab["genotype","Mean Sq"]
MSG
MSE=anova_tab["Residuals","Mean Sq"]
MSE

##now calculatinng the of PH
mean_PH<-mean(MorphologyData$PH)
mean_PH

##define the Replication
r<-3
r

##calculating the GV genetic variation 
GV<-(MSG-MSE)/r
GV #there is huge variance present in the genotype but we can not say it is high or low
    #that why we calculate the GCV
##calculating the phenotypic variance 
PV=GV+MSE
PV
#PV>GV because in PV there is environmenta effect is also add

## to say  these variation is low or  high we calculate GCV /PCV
##GCV
GCV<-(sqrt(GV)/mean_PH)*100
GCV ## 15.54 it mean plant hight exhibited genetic variabitily 

##PCV
PCV<-(sqrt(PV)/mean_PH)*100
PCV

##PCV>GCV GCV=15.65% PCV=17.67% 
# 17.67-15.65=2.02%
#plant height is mainly control by genetic factor , and environmental effects are comparetively low

## heritability
H2<-(GV/PV)*100
H2 #high hertability

##GA 
k<-2.06  #selection intensity constant 

GA<-k*sqrt(PV)*(H2/100)
GA

#GAM
GAM<-(GA/mean_PH)*100
GAM

results


###correlation###

##phenotypic correlation
data_num<-MorphologyData[,3:15]
phenotypic_cor<-cor((data_num))
phenotypic_cor
round(phenotypic_cor,3)

#genotypic correlation
library(dplyr)

geno_mean <- MorphologyData %>%
  group_by(genotype) %>%
  summarise(across(PH:`Y/P`, mean))
head(geno_mean)


geno_data <- geno_mean[, -1]

head(geno_data)
genotypic_cor <- cor(geno_data)

round(genotypic_cor, 3)

#visualization
library(corrplot)
#1. Phenotypic Correlation Heatmap
corrplot(phenotypic_cor,
         method="color",
         type="lower",
         addCoef.col="black")
#####

corrplot(phenotypic_cor,
         method = "color",
         type = "lower",
         tl.cex = 0.8)     # name size control
#2. genotypic Correlation Heatmap
corrplot(genotypic_cor,
         method="color",
         type="lower",
         addCoef.col="black")
corrplot(genotypic_cor,
         method = "color",
         type = "lower",
         tl.cex = 0.8)  

###path analysis 
##phenotyypic path anaysis
R <- phenotypic_cor[1:12,1:12]

R
dim(R)
r_y <- phenotypic_cor[1:12,13]

r_y
direct_effect <- solve(R) %*% r_y

round(direct_effect,4)
#residual effect

R2 <- sum(direct_effect * r_y)

Residual_effect <- sqrt(1 - R2)

R2
Residual_effect

##genotyypic path analysis
Rg <- genotypic_cor[1:12,1:12]

rg_y <- genotypic_cor[1:12,13]

direct_geno <- solve(Rg) %*% rg_y

round(direct_geno,4)
R2g <- sum(direct_geno * rg_y)

Residual_g <- sqrt(1 - R2g)

R2g
Residual_g


##visual representation 

##bar plot(optional)
direct_df <- data.frame(
  Trait = rownames(direct_geno),
  DirectEffect = as.vector(direct_geno)
)

barplot(direct_df$DirectEffect,
        names.arg = direct_df$Trait,
        las = 2,
        main = "Genotypic Direct Effects on Yield",
        ylab = "Direct Effect")
abline(h=0)

######
# Direct effects ko dataframe me convert karo

direct_df <- data.frame(
  Trait = rownames(direct_geno),
  DirectEffect = as.vector(direct_geno)
)

# Bar Plot

barplot(
  direct_df$DirectEffect,
  names.arg = direct_df$Trait,
  las = 2,
  ylab = "Direct Effect",
  main = "Genotypic Path Coefficients",
  col = "lightblue"
)

abline(h = 0, lty = 2)

###PCA ##
pca <- prcomp(data_num, scale. = TRUE)

summary(pca)
eigenvalues <- pca$sdev^2

round(eigenvalues,3)
round(pca$rotation[,1:4],3)

install.packages("factoextra")

library(factoextra)
fviz_eig(pca)
fviz_pca_biplot(pca,repel = TRUE)

##applying the regression on pca 
scores<-as.data.frame(pca$x)
#adding th eyeild variable 
scores$Yeild<-MorphologyData$`Y/P`
#fitting the regression 
pc1<-scores$PC1
pc2<-scores$PC2
pc3<-scores$PC3
pc4<-scores$PC4
model<-lm(Yeild~pc1+pc2+pc3+pc4,data=scores)
summary(model)#regressio model is highly significant h 



##CLUSTER ANALYSIS 
cluster_data <- MorphologyData %>%
  group_by(genotype) %>%
  summarise(across(PH:`Y/P`, mean))
data_scaled <- scale(cluster_data[,-1])
d <- dist(data_scaled,
          method="euclidean")
ward_cluster <- hclust(d,
                       method="ward.D2")
ward_cluster 
clusters <- cutree(ward_cluster,
                   k=4)
table(clusters)
subset(cluster_data,
       clusters==4)
aggregate(cluster_data[,2:14],
          by=list(Cluster=clusters),
          mean)
subset(cluster_data,
       clusters==3)

##inter cluster 
# Numeric traits from genotype means
data_num <- cluster_data[,2:14]

# Standardization
data_scaled <- scale(data_num)

# Covariance Matrix
S <- cov(data_scaled)

# Number of genotypes
n <- nrow(data_scaled)

# Mahalanobis D² Matrix
D2 <- matrix(0,n,n)

for(i in 1:n){
  for(j in 1:n){
    
    D2[i,j] <- mahalanobis(
      data_scaled[i,],
      data_scaled[j,],
      S
    )
    
  }
}

# Cluster numbers from Ward clustering
cluster_data$Cluster <- clusters

# Inter-cluster Distance Matrix
K <- length(unique(clusters))

interD <- matrix(0,K,K)

for(i in 1:K){
  
  for(j in 1:K){
    
    idx1 <- which(clusters==i)
    idx2 <- which(clusters==j)
    
    interD[i,j] <- mean(
      D2[idx1,idx2]
    )
    
  }
}

round(interD,2)
library(pheatmap)

rownames(interD) <- paste0("C",1:K)
colnames(interD) <- paste0("C",1:K)

pheatmap(interD,
         display_numbers = TRUE)

