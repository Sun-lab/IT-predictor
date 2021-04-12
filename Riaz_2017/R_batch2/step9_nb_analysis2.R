
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(quantreg)

theme_set(theme_classic())

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
sample_nb[1:2,]

sample_mb = read.table(file = "riaz_patient_mb_info.txt", 
                       header = TRUE, sep = " ")
dim(sample_mb)
sample_mb[1:2,]

sample_nb$PreOn = sub(".*_", "", sample_nb$sample)
sample_nb$Patient = sub("_.*", "", sample_nb$sample)
dim(sample_nb)
sample_nb[1:2,]

sample_mb$PreOn = sub(".*_", "", sample_mb$sample)
sample_mb$Patient = sub("_.*", "", sample_mb$sample)
dim(sample_mb)
sample_mb[1:2,]

clinic = read.xlsx("_supp/mmc2.xlsx", startRow=2)
dim(clinic)
clinic[1:2,]

#--------------------------------------------------------------------
# only consider the patients with mb information
# add mb data into clinic data
#--------------------------------------------------------------------

table(clinic$Patient %in% sample_mb$Patient)
table(tolower(clinic$Patient) %in% tolower(sample_mb$Patient))

patient = intersect(clinic$Patient, sample_mb$Patient)
length(patient)
clinic[which(! clinic$Patient %in% patient),]

clinic = clinic[match(patient, clinic$Patient),]
dim(clinic)
head(clinic)

sample_mb.pre = sample_mb[which(sample_mb$PreOn == "pre"),]
sample_mb.on  = sample_mb[which(sample_mb$PreOn == "on"),]

dim(sample_mb.pre)
dim(sample_mb.on)

mat.pre = match(patient, sample_mb.pre$Patient)
clinic[["mb.pre"]] = log10(sample_mb.pre$n_mutations_neoAg[mat.pre] + 1)

mat.on = match(patient, sample_mb.on$Patient)
clinic[["mb.on"]]  = log10(sample_mb.on$n_mutations_neoAg[mat.on] + 1)

#--------------------------------------------------------------------
# merge clinic data and nb data
#--------------------------------------------------------------------

sample_nb[1:2,]
sample_nb.pre = sample_nb[which(sample_nb$PreOn == "pre"),-c(1,8)]
sample_nb.on  = sample_nb[which(sample_nb$PreOn == "on"),-c(1,8)]

dim(sample_nb.pre)
dim(sample_nb.on)

names(sample_nb.on)[-7] = paste0(names(sample_nb.on)[-7], ".on")

sample_nb.pre[1:2,]
sample_nb.on[1:2,]

sample_nb.pre[,-7] = log10(sample_nb.pre[,-7] + 1)
sample_nb.on[,-7]  = log10(sample_nb.on[,-7] + 1)

clinic = merge(clinic, sample_nb.pre, by="Patient", all.x = TRUE)
dim(clinic)
clinic[1:2,]

clinic = merge(clinic, sample_nb.on, by="Patient", all.x = TRUE)
dim(clinic)
clinic[1:2,]

#--------------------------------------------------------------------
# set nb to be 0 if mb is 0
#--------------------------------------------------------------------

table(is.na(clinic$mhci))
table(is.na(clinic$mb.pre))

summary(clinic$mb.pre[is.na(clinic$mhci)])
table(clinic$mb.pre == 0)
table(clinic$mb.on == 0)

colSums(is.na(clinic))
grep("^mhc", names(clinic))
grep(".on$", names(clinic))

clinic[which(clinic$mb.pre == 0),15:20] = 0
clinic[which(clinic$mb.on == 0),21:26]  = 0

table(clinic$mb.pre == 0, clinic$mb.on == 0)

table(clinic$mb.pre == clinic$mb.on)
clinic[which(clinic$mb.pre == clinic$mb.on),]

#--------------------------------------------------------------------
# compare the mb/nb of our anlaysis vs. prevous work
#--------------------------------------------------------------------

clinic$Mutation.Load = log10(clinic$Mutation.Load + 1)
clinic$`Neo-antigen.Load` = log10(clinic$`Neo-antigen.Load` + 1)
clinic$`Neo-peptide.Load` = log10(clinic$`Neo-peptide.Load` + 1)

pdf("figures/compare_mb_nb_Riaz_vs_us.pdf", width=9, height=6)
par(mfrow=c(2,3), mar=c(5,4,1,1), bty="n")
plot(clinic$Mutation.Load, clinic$`Neo-antigen.Load`, 
     xlab="Riaz et al. MB", ylab="Riaz et al. NeoAg")
plot(clinic$Mutation.Load, clinic$`Neo-peptide.Load`, 
     xlab="Riaz et al. MB", ylab="Riaz et al. Neo-pep")
plot(clinic$`Neo-antigen.Load`, clinic$`Neo-peptide.Load`,
     xlab="Riaz et al. NeoAg", ylab="Riaz et al. Neo-pep")

plot(clinic$Mutation.Load, clinic$mb.pre, 
     xlab="Riaz et al. MB", ylab="our MB")
abline(0, 1, col="brown")
plot(clinic$`Neo-antigen.Load`, clinic$mhci,
     xlab="Riaz et al. NeoAg", ylab="our Neo-pep")
abline(0, 1, col="brown")
plot(clinic$`Neo-peptide.Load`, clinic$mhci,
     xlab="Riaz et al. Neo-pep", ylab="our Neo-pep")
abline(0, 1, col="brown")
dev.off()

#--------------------------------------------------------------------
# compare with response data
#--------------------------------------------------------------------

clinic$Response_orig = clinic$Response
clinic$Response3 = clinic$Response
clinic$Response3[which(clinic$Response == "NE")] = "PD"
clinic$Response3[which(clinic$Response == "CR")] = "PRCR"
clinic$Response3[which(clinic$Response == "PR")] = "PRCR"

clinic$Response = clinic$Response3
clinic$Response[which(clinic$Response3 == "SD")] = "PDSD"
clinic$Response[which(clinic$Response3 == "PD")] = "PDSD"

table(clinic$Response, clinic$Response3)
table(clinic$Response, clinic$Response)

clinic$Cohort[which(clinic$Cohort == "NIV3-NAIVE")] = "Ipi naive"
clinic$Cohort[which(clinic$Cohort == "NIV3-PROG")] = "Ipi prog"

anova(lm(clinic$Mutation.Load ~ clinic$Response))
anova(lm(clinic$`Neo-antigen.Load` ~ clinic$Response))
anova(lm(clinic$`Neo-peptide.Load` ~ clinic$Response))

anova(lm(clinic$mb.pre ~ clinic$Response))

anova(lm(clinic$mhci ~ clinic$Response))
anova(lm(clinic$mhci.lt.ref ~ clinic$Response))
anova(lm(clinic$mhci.lt.ref.weighted ~ clinic$Response))

anova(lm(clinic$mhciipan ~ clinic$Response))
anova(lm(clinic$mhciipan.lt.ref ~ clinic$Response))
anova(lm(clinic$mhciipan.lt.ref.weighted ~ clinic$Response))


x = clinic$mhci.lt.ref.weighted
y = clinic$mhciipan.lt.ref.weighted
z = x + y
anova(lm(z ~ clinic$Response))

clinic$mut_diff = clinic$mb.pre - clinic$mb.on
anova(lm(clinic$mut_diff ~ clinic$Response))

#--------------------------------------------------------------------
# boxplot
#--------------------------------------------------------------------

p0 = ggplot(clinic, aes(x=Response, y=mb.pre)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "log10(MB)") 

p1 = ggplot(clinic, aes(x=Response, y=mhci.lt.ref.weighted)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "log10(MHC-I NB weighted)") 

p2 = ggplot(clinic, aes(x=Response, y=mhciipan.lt.ref.weighted)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "log10(MHC-II NB weighted)") 

p3 = ggplot(clinic, aes(x=Response, y=mb.pre - mb.on)) + 
  geom_boxplot() + 
  geom_jitter(position=position_jitter(0.2), size=1) + 
  labs(title="", x="Response", y = "log10(MB) pre - on") 

pdf("figures/boxplot_nb_vs_response.pdf", width=6, height=6)
ggarrange(p0, p1, p2, p3, ncol = 2, nrow = 2)
dev.off()

#--------------------------------------------------------------------
# check the way to combine MHC-I and MHC-II
#--------------------------------------------------------------------

g0 = ggplot(clinic, aes(x=mhci, y=mhci.lt.ref, color=Response)) + 
  geom_point()

g1 = ggplot(clinic, aes(x=mhci, y=mhci.lt.ref.weighted, color=Response)) + 
  geom_point()

g2 = ggplot(clinic, aes(x=mhciipan, y=mhciipan.lt.ref, color=Response)) + 
  geom_point()

g3 = ggplot(clinic, aes(x=mhciipan, y=mhciipan.lt.ref.weighted, 
                        color=Response)) + 
  geom_point()

pdf("figures/scatter_nb_vs_response.pdf", width=8, height=6)
ggarrange(g0, g1, g2, g3, ncol = 2, nrow = 2)
dev.off()

#----------------------------------------------------------------
# fit linear regression 
#----------------------------------------------------------------

x = clinic$mhci
y = clinic$mhciipan

g2 = ggplot(clinic, aes(x=mhci, y=mhciipan, color=Response)) + 
  geom_point(size=1.2) + 
  labs(x="HLA-I neoantigen burden", y = "HLA-II neoantigen burden") 

w2kp = which(x > 2 & x < 3.5)

lm1  = lm(y ~ -1 + x, subset=w2kp)
summary(lm1)

beta0 = 0
beta1 = lm1$coef[1]
beta1

g2_lin = g2 + geom_segment(aes(x = 2, y = beta0 + beta1*2, 
                           xend = 3.5, yend = beta0 + beta1*3.5), 
                       col="black")

#----------------------------------------------------------------
# fit median quantile regression 
#----------------------------------------------------------------

med.fit = rq(y ~ -1 + x, tau = .5, subset = w2kp)
summary(med.fit)

med.beta0 = 0
med.beta1 = med.fit$coef[1]
med.beta1

g2_med = g2 + geom_segment(aes(x = 2, y = med.beta0 + med.beta1*2, 
                               xend = 3.5, yend = med.beta0 + med.beta1*3.5), 
                           col="black")

#----------------------------------------------------------------
# Prediction
#----------------------------------------------------------------

xnew = x[w2kp]

#Linear Model 
yhat = predict(lm1, data.frame(x=xnew))
resid1 = abs(y[w2kp] - yhat)

df1 = data.frame(Absolute_Residual=resid1, Response=clinic$Response[w2kp], 
                 Cohort=clinic$Cohort[w2kp])

g3_lin  = ggboxplot(df1, x = "Response", y = "Absolute_Residual", 
                    color = "Response", add = "jitter",
                    outlier.shape = NA, legend = "none") + 
  stat_compare_means(label = "p.format")

g3_lin2 = ggboxplot(df1, x = "Cohort", y = "Absolute_Residual",
               color = "Response", add = "jitter", 
               outlier.shape = NA, legend = "none") + 
  stat_compare_means(aes(group = Response), label = "p.format")

# Median quantile regression 
yhat_med = predict(med.fit, data.frame(x=xnew))
resid1_med = abs(y[w2kp] - yhat_med)

df1_med = data.frame(Absolute_Residual=resid1_med, 
                     Response=clinic$Response[w2kp], 
                     Cohort=clinic$Cohort[w2kp])

g3_med  = ggboxplot(df1_med, x = "Response", y = "Absolute_Residual", 
                    color = "Response", add = "jitter", 
                    outlier.shape = NA, legend = "none") + 
  stat_compare_means(label = "p.format")

g3_med2 = ggboxplot(df1_med, x = "Cohort", y = "Absolute_Residual",
                    color = "Response", add = "jitter", 
                    outlier.shape = NA, legend = "none") + 
  stat_compare_means(aes(group = Response), label = "p.format")

pdf("figures/boxplot_nb_res_vs_response_linear.pdf", width=9, height=3)
ggarrange(g2_lin, g3_lin, g3_lin2, ncol = 3, nrow = 1, widths=c(4,2,3), 
          labels = c("A", "B", "C"))
dev.off()

pdf("figures/boxplot_nb_res_vs_response_medquant.pdf", width=9, height=3)
ggarrange(g2_med, g3_med, g3_med2, ncol = 3, nrow = 1, widths=c(4,2,3),
          labels = c("A", "B", "C"))
dev.off()

#--------------------------------------------------------------------
# save updated clinical information
#--------------------------------------------------------------------

write.table(clinic, file = "output/clinic_info.txt", append = FALSE, 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

sessionInfo()
q(save="no")

