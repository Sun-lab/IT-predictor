#--------------------------------------------------------------------
# Step 6: Neoantigen Burden Analysis 
 # 1. Import Hugo and Riaz prediction ds
 # 3. Calculate nb 
 # 4. Attach mb and drug response by matchID
 # 5. Plots 
 # 6. Mulitnomial model 
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# 
#--------------------------------------------------------------------
riaz_hlai_sum_min_pred= load("Riaz2017_data/final_anno_ds/for_netmhci/pred/riaz_hlai_sum_min_pred.RData")

#--------------------------------------------------------------------
# Merge mb and response with nbdx
#--------------------------------------------------------------------
riaz_hlai_sum_min_pred_response = merge(riaz_hlai_sum_min_pred, 
                                        riaz.patient.mut1[,c("matchID","irRECIST3","mb")], by = "matchID")
hlaii_sum_min_pred_response = merge(hlaii_sum_min_pred, patient.mut1[,c(1,3:4)], by = "matchID")

#--------------------------------------------------------------------
# Step 3. Calculate nb
#--------------------------------------------------------------------
riaz_hlai_sum_min_pred_response$nb1.ind = riaz_hlai_sum_min_pred_response$min_pred<500
hlaii_sum_min_pred_response$nb2.ind = hlaii_sum_min_pred_response$min_pred <500

save(riaz_hlai_sum_min_pred_response, file = "Riaz2017_data/final_anno_ds/for_netmhci/pred/riaz_hlai_sum_min_pred_response.RData")
save(hlaii_sum_min_pred_response, file = "for_netmhciipan/pred/hlaii_sum_min_pred_response.RData")

#----MHC-I
library(dplyr)
nb1.ind.sum=aggregate(riaz_hlai_sum_min_pred_response$nb1.ind,
                      by=list(riaz_hlai_sum_min_pred_response$matchID), FUN=sum)
names(nb1.ind.sum)= c("matchID", "nburden1")


#----MHC-II
nb2.ind.sum=aggregate(hlaii_sum_min_pred_response$nb2.ind,
                      by=list(hlaii_sum_min_pred_response$matchID), FUN=sum)
names(nb2.ind.sum)= c("matchID", "nburden2")



#--------------------------------------------------------------------
# Step 4. Attach mb and drug response by matchID
#--------------------------------------------------------------------
# --- HLA-I
nb1.ind.sum.response = merge(nb1.ind.sum, riaz.patient.mut1[,c("matchID", "mb", 
                                                               "irRECIST3")], 
                             by = "matchID")
write.table(nb1.ind.sum.response, file = "Riaz2017_data/final_anno_ds/for_netmhci/pred/nb_hlai_9pep.txt", 
            quote = FALSE, col.names = TRUE, row.names=FALSE)

# --- HLA-II
nb2.ind.sum.response = merge(nb2.ind.sum, patient.mut1, by = "matchID")
write.table(nb2.ind.sum.response, file = "for_netmhciipan/pred/nb_hlaii_15pep.txt", 
            quote = FALSE, col.names = TRUE, row.names=FALSE)
#--------------------------------------------------------------------
# Step 5. Plots 
#--------------------------------------------------------------------
library(ggplot2)
#------------------------------ MB vs. Drug response 
ggplot(nb1.ind.sum.response, aes(x=irRECIST3, y=log(mb+1,10))) + 
  geom_jitter(position=position_jitter(0.2)) + 
  stat_summary(fun.y=median, geom="point", shape=18,size=3, color="red") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue") + 
  ggtitle("Figure 0. log10(mb +1) by Drug response")

#------------------------------ MB vs. NB 
# HLA-I 
set.seed(92)
ggplot(nb1.ind.sum.response, aes(x=log(mb+1,10), y=log(nburden1+1,10))) + 
  geom_point() + 
  ggtitle("Figure 1.1. log10(mb+1) vs. log10(Neoantigen Burden Defn 1 +1) ")


# HLA-II
ggplot(nb2.ind.sum.response, aes(x=log(mb+1,10), y=log(nburden2+1,10))) + 
  geom_jitter(position=position_jitter(0.2)) + 
  ggtitle("Figure 1.2. log10(mb+1) vs. log10(Neoantigen Burden Defn 2 +1) ")

#------------------------------ NB vs. Drug response
# HLA-I 
ggplot(nb1.ind.sum.response, aes(x=irRECIST3, y=log(nburden1+1,10))) + 
  geom_jitter(position=position_jitter(0.2)) + 
  stat_summary(fun.y=median, geom="point", shape=18,size=3, color="red") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue") + 
  ggtitle("Figure 2.1. log10(Neoantigen Burden Defn 1 +1) by Drug response")

# HLA-II 
ggplot(nb2.ind.sum.response, aes(x=irRECIST, y=log(nburden2+1,10))) + 
  geom_jitter(position=position_jitter(0.2)) + 
  stat_summary(fun.y=median, geom="point", shape=18,size=3, color="red") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue") + 
  ggtitle("Figure 2.2. log10(Neoantigen Burden Defn 2 +1) by Drug response")



#--------------------------------------------------------------------
# Step 6. Multinomial Model 
#--------------------------------------------------------------------
# ----------------------- Multinomial model w/o log(mb)
# -------------HLA-I
library(nnet)
nb1.ind.sum.response$irRECIST2 <- relevel(nb1.ind.sum.response$irRECIST, ref = "Progressive Disease")
model1 <- multinom(irRECIST2 ~ log10(nburden1+1), data = nb1.ind.sum.response)

#calculate z-score and p-value (Wald Test) of regression coefficients
z <- summary(model1)$coefficients/summary(model1)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

#output model summary
model1.sum= cbind( summary(model1)$coefficients, summary(model1)$standard.errors, p)
colnames(model1.sum)=c("Intercept", "log(burden)", "Int_SE", "log(burden)_SE", 
                       "Int_p", "log(burden)_p")

model1.sum


#------------- HLA-II
nb2.ind.sum.response$irRECIST2 <- relevel(nb2.ind.sum.response$irRECIST, ref = "Progressive Disease")
model2 <- multinom(irRECIST2 ~ log10(nburden2+1), data = nb2.ind.sum.response)

#calculate z-score and p-value (Wald Test) of regression coefficients
z <- summary(model2)$coefficients/summary(model2)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

#output model summary
model2.sum= cbind( summary(model2)$coefficients, summary(model2)$standard.errors, p)
colnames(model2.sum)=c("Intercept", "log(burden2)", "Int_SE", "log(burden2)_SE", 
                       "Int_p", "log(burden2)_p")

model2.sum


#--------------------------Multinomial model + log(mb)
# -------------HLA-I
model1.1<- multinom(irRECIST2 ~ log10(nburden1+1)  + log10(mb+1), data = nb1.ind.sum.response)

#calcualte p-values
z <- summary(model1.1)$coefficients/summary(model1.1)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

##summary 
#summary(model5)
model1.1.sum= cbind( summary(model1.1)$coefficients, summary(model1.1)$standard.errors, p)
colnames(model1.1.sum)=c("Intercept", "log(burden)", "log(mb)", "IntSE", "log(burden)SE",
                       "log(mb)SE", "Int_p", "log(burden)_p", "log(mb)_p")

model1.1.sum

# -------------HLA-II
model2.1<- multinom(irRECIST2 ~ log10(nburden2+1)  + log10(mb+1), data = nb2.ind.sum.response)

#calcualte p-values
z <- summary(model2.1)$coefficients/summary(model2.1)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

##summary 
#summary(model5)
model2.1.sum= cbind( summary(model2.1)$coefficients, summary(model2.1)$standard.errors, p)
colnames(model2.1.sum)=c("Intercept", "log(burden)", "log(mb)", "IntSE", "log(burden)SE",
                         "log(mb)SE", "Int_p", "log(burden)_p", "log(mb)_p")

model2.1.sum

#--------------------------------------------------------------------
# Step 7. ANOVA - overall p-value
#--------------------------------------------------------------------
#Model with just mb
patient.mut1$irRECIST2 <- relevel(patient.mut1$irRECIST, ref = "Progressive Disease")
model0 = multinom(irRECIST2 ~  log10(mb+1), data = patient.mut1)
summary(model0)

# ANOVA w/ nb1
anova(model0, model1.1)

#ANOVA w/ nb2 
anova(model0, model2.1)

#--------------------------------------------------------------------
# Step 8. Linear Model 
#--------------------------------------------------------------------
### Linear Model 
linearMod1=  lm(log(nburden1+1,10)~  factor(irRECIST2) + log(mb+1,10), data=nb1.ind.sum.response)  
summary(linearMod1)

linearMod2=  lm(log(nburden2+1,10)~  factor(irRECIST2) + log(mb+1,10), data=nb2.ind.sum.response)  
summary(linearMod2)

#--------------------------------------------------------------------
# Step 9. Combine CR and PR == Drug response
#--------------------------------------------------------------------
# combine CR and PR
nb1.ind.sum.response$irRECIST3 = ifelse(nb1.ind.sum.response$irRECIST == "Progressive Disease", 
                                        "Progressive Disease", "Drug Response")
model1a = multinom(irRECIST3 ~ log10(nburden1+1), data = nb1.ind.sum.response)
#calculate z-score and p-value (Wald Test) of regression coefficients
z <- summary(model1a)$coefficients/summary(model1a)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

#output model summary
model1a.sum= cbind( summary(model1a)$coefficients, summary(model1a)$standard.errors, p)
colnames(model1a.sum)=c("ParamEst", "SE", "p")

model1a.sum

#--- plot 
ggplot(nb1.ind.sum.response, aes(x=irRECIST3, y=log(nburden1+1,10))) + 
  geom_jitter(position=position_jitter(0.2)) + 
  stat_summary(fun.y=median, geom="point", shape=18,size=3, color="red") + 
  stat_summary(fun.y=mean, geom="point", shape=18, size=3, color="blue") + 
  ggtitle("Figure 2.1. log10(Neoantigen Burden Defn 1 +1) by Drug response")

# --- multinomial with mb
model1.1a<- multinom(irRECIST3 ~ log10(nburden1+1)  + log10(mb+1), data = nb1.ind.sum.response)

#calcualte p-values
z <- summary(model1.1a)$coefficients/summary(model1.1a)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

##summary 
#summary(model5)
model1.1a.sum= cbind( summary(model1.1a)$coefficients, summary(model1.1a)$standard.errors, p)
colnames(model1.1a.sum)=c("ParamEst", "SE", "p")

model1.1a.sum

# ANOVA w/ nb1
patient.mut1$irRECIST3 <- ifelse(patient.mut1$irRECIST == "Progressive Disease", "Progressive Disease", "Drug Response")
model0a = multinom(irRECIST3 ~  log10(mb+1), data = patient.mut1)
anova(model0a, model1.1a)
