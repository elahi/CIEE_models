source("Bdiv_model_multistress.R")
source("change_metrics.R")
#make sure that you have these packages
library(vegan)
library(picante)
library(cluster)
library(RColorBrewer)
library(FD)
library(ggplot2)

#this is the model
Stress = 1
system.time(model.results<-Bdiv_model(disp=0.01,Int_type="Comp",nSpecies=120,Stress = Stress,Stress2 = 5))
#disp is the proportion of biomass that disperses in each time step - recommend using values between 0 and 0.005
#for Int_type choose the type of interactions: 
#"NoInt" - no interactions; "Comp" - competition; "Mixed" -competition, mutualism, and parasitsm; "Trophic" - tri trophic food web

Abundance_ts<-model.results$Abundance_ts
Traits<-model.results$traits
Stressor<-model.results$Stressor
Stress.df<-model.results$Stress.df

sub_sample<-seq(1000,1000+Stress*1000,by=1000)
Abund_sub<-Abundance_ts[,,sub_sample]
Stressor_sub<-Stressor[sub_sample]

Stress.pre.post<-rbind(data.frame(Stressor1=0,Stressor2=rep(0,400)),Stress.df)

RDA<-rda(decostand(rbind(Abund_sub[,,1],Abund_sub[,,2]),method="hellinger")~.,Stress.pre.post)

Stress.pch<-Stress.pre.post$Stressor2[401:800]+1
Stress.pch[Stress.pch==6]<-2

pdf("Model fingerprint.pdf",8,8)
par(pty='s',las=1)
plot(RDA, type='n')
points(RDA,display="sites",col=rep(c(1,2),each=400), pch=Stress.pch)
text(RDA,display="bp", lwd=2)
legend("bottomleft",bty='n',c("Pre stress", "Post stress","Regional stress","Local and regional stress"),pch=c(1,1,1,2),col=c(1,2,2,2))
dev.off()

RsquareAdj(RDA)

anova.cca(RDA,permutations = 999)
anova.cca(RDA,by="axis",permutations = 999)

RDA_constrained1<-rda(decostand(rbind(Abund_sub[,,1],Abund_sub[,,2]),method="hellinger"),Stress.pre.post$Stressor2,Stress.pre.post$Stressor1)
RsquareAdj(RDA_constrained1)
RDA_constrained2<-rda(decostand(rbind(Abund_sub[,,1],Abund_sub[,,2]),method="hellinger"),Stress.pre.post$Stressor1,Stress.pre.post$Stressor2)
RsquareAdj(RDA_constrained2)

par(mfrow=c(3,1),mar=c(3,3,3,3))
plot(RDA, type='n',main="Full RDA")
points(RDA,display="sites",col=rep(c(1,2),each=400), pch=20)
text(RDA,display="bp", lwd=2)

plot(RDA_constrained1, type='n',main="Constrained by Stress 1")
points(RDA_constrained1,display="sites",col=rep(c(1,2),each=400), pch=20)
text(RDA_constrained1,display="bp", lwd=2)

plot(RDA_constrained2, type='n', main="Constrained by Stress 2")
points(RDA_constrained2,display="sites",col=rep(c(1,2),each=400), pch=20)
text(RDA_constrained2,display="bp", lwd=2)
par(mfrow=c(1,1))

Var_part<-varpart(decostand(rbind(Abund_sub[,,1],Abund_sub[,,2]),method="hellinger"),Stress.pre.post$Stressor1,Stress.pre.post$Stressor2)
plot(Var_part)
