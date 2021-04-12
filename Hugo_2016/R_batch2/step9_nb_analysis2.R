
library(openxlsx)
library(ggplot2)
library(gridExtra)

#--------------------------------------------------------------------
# Step 9: Neoantigen Burden Analysis 
 # 1. compare nb of our anlaysis versus those reported in the paper
 # 2. associatoin nb with response
 # 3. Plots 
 # 4. Mulitnomial model 
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# load mb/nb results and clinical information
#--------------------------------------------------------------------

sample_nb = read.table("output/neoAg_burden.txt", sep="\t", 
                       header=TRUE, as.is=TRUE)
dim(sample_nb)
colnames(sample_nb)[1] ="matchID"
sample_nb[1:2,]


sample_mb = read.table(file = "hugo_patient_mb_info.txt", 
                       header = TRUE, sep = " ")
dim(sample_mb)
# remove Vand_Pt27_2 
sample_mb = sample_mb[which(sample_mb$sample!="Vand_Pt27_2"),]
#dim(sample_mb)
colnames(sample_mb)[1] ="matchID"
sample_mb[1:2,]

clinic = read.delim("patient_info_with_mutations.txt", sep="\t")
dim(clinic)
clinic[1:2,1:26]

#--------------------------------------------------------------------
# only consider the patients with mb information
# add mb data into clinic data
#--------------------------------------------------------------------

table(tolower(clinic$matchID) %in% tolower(sample_mb$matchID))

patient = intersect(clinic$matchID, sample_mb$matchID)
length(patient)
clinic[which(! clinic$matchID %in% patient),1:26]

clinic = clinic[match(patient, clinic$matchID),]
dim(clinic)
head(clinic[,1:26])

mat= match(patient, sample_mb$matchID)
clinic[["mb.calc"]] = log10(sample_mb$n_mutations_neoAg[mat] + 1)
head(clinic[,c(1:26, 3421)])

#--------------------------------------------------------------------
# merge clinic data and nb data
#--------------------------------------------------------------------

sample_nb[1:2,]
for(i in 2:7){
  sample_nb[,i] = log10(sample_nb[,i] + 1)
}
sample_nb[1:2,]

clinic = merge(clinic[,c(1:26,3421)], sample_nb, by="matchID", all.x = TRUE)
dim(clinic)
clinic[1:2,]


#--------------------------------------------------------------------
# set nb to be 0 if mb is 0
#--------------------------------------------------------------------

table(is.na(clinic$mhci))
table(is.na(clinic$mb.calc))

summary(clinic$mb.calc[is.na(clinic$mhci)])
table(clinic$mb.calc == 0)

#--------------------------------------------------------------------
# compare the mb/nb of our anlaysis vs. prevous work
#--------------------------------------------------------------------

#clinic$Mutation.Load = log10(clinic$Mutation.Load + 1)
#clinic$`Neo-antigen.Load` = log10(clinic$`Neo-antigen.Load` + 1)
#clinic$`Neo-peptide.Load` = log10(clinic$`Neo-peptide.Load` + 1)

#pdf("figures/compare_mb_nb_hugo_vs_us.pdf", width=9, height=6)
#par(mfrow=c(2,3), mar=c(5,4,1,1), bty="n")
#plot(clinic$Mutation.Load, clinic$`Neo-antigen.Load`, 
#     xlab="hugo et al. MB", ylab="hugo et al. NeoAg")
#plot(clinic$Mutation.Load, clinic$`Neo-peptide.Load`, 
#     xlab="hugo et al. MB", ylab="hugo et al. Neo-pep")
#plot(clinic$`Neo-antigen.Load`, clinic$`Neo-peptide.Load`,
#     xlab="hugo et al. NeoAg", ylab="hugo et al. Neo-pep")

#plot(clinic$Mutation.Load, clinic$mb.pre, 
#     xlab="hugo et al. MB", ylab="our MB")
#abline(0, 1, col="brown")
#plot(clinic$`Neo-antigen.Load`, clinic$mhci,
#     xlab="hugo et al. NeoAg", ylab="our Neo-pep")
#abline(0, 1, col="brown")
#plot(clinic$`Neo-peptide.Load`, clinic$mhci,
#     xlab="hugo et al. Neo-pep", ylab="our Neo-pep")
#abline(0, 1, col="brown")
#dev.off()

#--------------------------------------------------------------------
# compare with response data
#--------------------------------------------------------------------
table(clinic$irRECIST)

#change 
clinic$irRECIST = as.character(clinic$irRECIST)

clinic$Response2 = clinic$irRECIST
clinic$Response2[which(clinic$irRECIST == "Progressive Disease")] = "PDSD"
clinic$Response2[which(clinic$irRECIST == "Complete Response")] = "PRCR"
clinic$Response2[which(clinic$irRECIST == "Partial Response")] = "PRCR"


table(clinic$irRECIST, clinic$Response2)

anova(lm(clinic$mb ~ clinic$Response2))

anova(lm(clinic$mb.calc ~ clinic$Response2))

anova(lm(clinic$mhci ~ clinic$Response2))
anova(lm(clinic$mhci.lt.ref ~ clinic$Response2))
anova(lm(clinic$mhci.lt.ref.weighted ~ clinic$Response2))

anova(lm(clinic$mhciipan ~ clinic$Response2))
anova(lm(clinic$mhciipan.lt.ref ~ clinic$Response2))
anova(lm(clinic$mhciipan.lt.ref.weighted ~ clinic$Response2))


x = clinic$mhci.lt.ref.weighted
y = clinic$mhciipan.lt.ref.weighted
z = x + y
anova(lm(z ~ clinic$Response2))

#--------------------------------------------------------------------
# boxplot
#--------------------------------------------------------------------
clinic$Response2[which(clinic$irRECIST == "Progressive Disease")] = "PDSD"
clinic$Response2[which(clinic$irRECIST == "Complete Response")] = "PRCR"
clinic$Response2[which(clinic$irRECIST == "Partial Response")] = "PRCR"

clinic$sum.nb = clinic$mhci.lt.ref.weighted + clinic$mhciipan.lt.ref.weighted
clinic$sum.nb.raw = clinic$mhci + clinic$mhciipan
table(clinic$Response2, clinic$irRECIST)

ggplot(clinic, aes(x=Response2, y=sum.nb)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "Sum of HLA-I and HLA-II") 

clinic1 = clinic[which(clinic$irRECIST!="Partial Response"),]
table(clinic1$Response2)

ggplot(clinic, aes(x=irRECIST, y=sum.nb)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "Sum of HLA-I and HLA-II")

ggplot(clinic, aes(x=irRECIST, y=sum.nb.raw)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "Sum of HLA-I and HLA-II")
#--------------------------------------------------------------------
# check the way to combine MHC-I and MHC-II
#--------------------------------------------------------------------

g0 = ggplot(clinic, aes(x=mhci, y=mhciipan, 
                        shape=Response2, color=Response2)) + geom_point()

g1 = ggplot(clinic, aes(x=mhci.lt.ref, y=mhciipan.lt.ref, 
                        shape=Response2, color=Response2)) + geom_point()

g2 = ggplot(clinic, aes(x=mhci.lt.ref.weighted, y=mhciipan.lt.ref.weighted, 
                   shape=Response2, color=Response2)) + geom_point()

ymin = 1
ymax = 2
x = clinic$mhci.lt.ref.weighted
y = clinic$mhciipan.lt.ref.weighted
length(x)
length(y)

w2kp = which(y > ymin & y < ymax)
length(w2kp)

load("../../Riaz2017/final_anno_ds/riaz_model1.rda")
#lm1 = lm(y ~ -1 + x, subset=w2kp)
beta0 = 0
beta1 = lm1$coef[1]
beta1 

g2 = g2 + geom_segment(aes(x = (ymin-beta0)/beta1, y = ymin, 
                           xend = (ymax-beta0)/beta1, yend = ymax), 
                       col="black")

xnew = x[which(y > ymin & y < ymax)]
yhat = predict(lm1, data.frame(x=xnew))

resid1 = abs(y[which(y > ymin & y < ymax)] - yhat)

df1 = data.frame(resid=resid1, resposne=clinic$Response2[which(y > ymin & y < ymax)])
g3 = ggplot(df1, aes(x=resposne, y=resid)) + 
      geom_boxplot() + 
      geom_jitter(position=position_jitter(0.2), size=1) + 
      labs(title="", x="Response", y = "|Residuals|") 

l = mget(c("g0", "g1", "g2"))

ggsave(paste("figures/scatter_riazLine_nb_vs_response_ymax", toString(ymax), ".pdf", sep =""),
       width=4.5, height=9, 
       marrangeGrob(grobs = l, nrow=3, ncol=1))
ggsave(paste("figures/boxplot_riazLine_nb_res_vs_response_ymax", toString(ymax), ".pdf", sep =""), 
       g3, width=3, height=3)

#--------------------------------------------------------------------
# save updated clinical information
#--------------------------------------------------------------------

write.table(clinic, file = "output/clinic_info.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

sessionInfo()
q(save="no")

