#Definir o diretorio e carregar os pacotes
library("geomorph")
library("phytools")
library("geiger")
library("ggpubr")
library("viridis")
library("caper")
library("psych")

#MORPHOGEO LATERAL####
#carregar o arquivo TPS, definindo o speciesID
landmarks1 <- readland.tps(file = "input_lateralna2.tps", specID = "imageID")
landmarks <- estimate.missing(landmarks1, method = "TPS")

#Definir os agrupamentos (por spp, ou genero, ou dieta, etc)
Groups <- read.csv("input_list_lateral.csv", header = T, sep = ";")
Spp<-as.factor(Groups[,1])#list of specific groupings

#GPA
GPA <- gpagen(landmarks)
plotAllSpecimens(GPA$coords)
lands2d <- two.d.array(GPA$coords)
Csize2d <- log(GPA$Csize)
GPA3 <- arrayspecs(lands2d, 26, 2)
plotOutliers(GPA3, inspect.outliers = TRUE)
write.csv(GPA3, file ="output_GPA_lateral.csv")

#interlmkdist(A, lmks)
lmks <- matrix(c(1,22,3,4,5,8), ncol=2, byrow=TRUE, 
               dimnames = list(c("total", "rostraltube","postincisive"),c("start", "end")))
A <- landmarks
lineardists <- interlmkdist(A, lmks)
lineardists
linear <- aggregate(lineardists, by = list(Spp), FUN = mean) 
write.csv(lineardists, file = "output_linear_specimens_lateral.csv")
write.csv(linear, file = "output_linear_species_lateral.csv")

#size
CS_species <- aggregate(GPA$Csize, by = list(Spp), FUN = mean) 
write.csv(CS_species, file="output_cs_species_lateral.csv")
#logsize
logCS_species <- aggregate(Csize2d, by = list(Spp), FUN = mean) 
write.csv(logCS_species, file = "output_logcs_species_lateral.csv")
write.csv(Csize2d, file = "output_logCS_specimens_lateral.csv")

#Pegando o valor midio de landmarks por especie
Sppmean <-aggregate(lands2d, by = list(Spp), FUN=mean)[,2:53] 
Sppmean3d <-arrayspecs(Sppmean, 26, 2, sep=".")
write.csv(Sppmean, file="output_sppmeanlandmarks_lateral.csv")

#PCA por especime
pca<-gm.prcomp(lands2d, tol=0)
pcs<-pca$x
summary(pca)
Sppsd<-aggregate(pcs, by = list(Spp), FUN=sd)[,2:53]
Sppsd
write.csv(pcs, file = "output_pcs_specimens_lateral.csv")

#PCA por especie
pca <- gm.prcomp(Sppmean3d, tol=0)
summary(pca)
pcs<-pca$x
plot(pca)
pcs
write.csv(pcs, file = "output_pcs_species_lateral.csv")

#plotting mean shape
ref <- mshape(GPA3)
plot(ref)

#outline
lands<-readland.tps(file = "input_outline_lateral.TPS")
out<-warpRefOutline(file = "input_outline_lateral.txt", lands[,,1], ref)####

plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "points")

dev.off()
plotRefToTarget(ref,pca$shapes$shapes.comp1$min,outline=out$outline)
title("pc1 min")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max,outline=out$outline)
title("pc1 max")

plotRefToTarget(ref,pca$shapes$shapes.comp2$min,outline=out$outline)
title("pc2 min")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max,outline=out$outline)
title("pc2 max")


#MORPHOGEO DORSAL####
#loading TPS file, defining speciesID
landmarks1 <- readland.tps(file = "input_dorsalna.tps", specID = "imageID", negNA = TRUE)
landmarks <- estimate.missing(landmarks1, method = "TPS")
options(scipen=999)
#defining groups (spp)
Groups <- read.csv("input_list_dorsal.csv", header = T, sep = ";")
Spp<-as.factor(Groups[,1])#list of specific groupings
ind<-as.factor(Groups[,2])#list of individuals for symmetric component extraction

#GPA
GPA <- gpagen(landmarks)
plotAllSpecimens(GPA$coords)
lands2d <- two.d.array(GPA$coords)
GPA3 <- arrayspecs(lands2d, 21, 2)
plotOutliers(GPA3, inspect.outliers = TRUE)

#symmetric component
lm<-read.csv("input_lmpairs_dorsal.csv", header=TRUE, sep = ",")#landmarks pairs on each side of bilateral symmetry
lmm<-as.matrix(lm)

gdf <- geomorph.data.frame(shape = GPA3,
                           ind = ind)
sym <- bilat.symmetry(A = shape, ind = ind,
                      object.sym = TRUE,
                      land.pairs = lmm, data = gdf, RRPP = TRUE, iter = 149)
summary(sym)
sym<-sym$symm.shape
write.csv(sym, file="output_sym_dorsal.csv")

sym2 <- two.d.array(sym)
write.csv(sym2, file ="output_sym2_dorsal.csv")

#size
CS_specimens <- GPA$Csize
write.csv(CS_specimens, file="output_cs_specimens_dorsal.csv")
CS_species <- aggregate(GPA$Csize, by = list(Spp), FUN = mean) 
write.csv(CS_species, file="output_cs_species_dorsal.csv")

#logsize
Csize2d <- log(GPA$Csize)
write.csv(Csize2d, file="output_logcs_specimens_dorsal.csv")
logCS_species <- aggregate(Csize2d, by = list(Spp), FUN = mean) 
write.csv(logCS_species, file = "output_logcs_species_dorsal.csv")

#taking the average value of landmarks by species
Sppmean <-aggregate(sym2, by = list(Spp), FUN=mean)[,2:43] 
Sppmean3d <-arrayspecs(Sppmean, 21, 2, sep=".")
write.csv(Sppmean, file="output_sppmeanlandmarks_dorsal.csv")

#PCA por especime
pca<-gm.prcomp(sym2, tol=0)
summary(pca)
pcs<-pca$x
write.csv(pcs, file = "output_pcs_specimens_dorsal.csv")

#PCA (species)
pca <- gm.prcomp(Sppmean3d, tol=0)
summary(pca)
pcs<-pca$x
write.csv(pcs, file="output_pcs_species_dorsal.csv")

#plotting mean shape
ref <- mshape(Sppmean3d)
findMeanSpec(GPA$coords)
#loading outlines for plotting shape changes on each component
lands<-readland.tps(file = "input_outline_dorsal.TPS")
out<-warpRefOutline(file = "input_outline_dorsal.csv", lands[,,1], ref)####

#plotting shape changes by points
plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "points")

#plotting shape changes by vectors
plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "vector")

##plotting shape changes by outlines
dev.off()
plotRefToTarget(ref,pca$shapes$shapes.comp1$min,outline=out$outline)
title("pc1 min")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max,outline=out$outline)
title("pc1 max")

plotRefToTarget(ref,pca$shapes$shapes.comp2$min,outline=out$outline)
title("pc2 min")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max,outline=out$outline)
title("pc2 max")

#MORPHOGEO VENTRAL####
#loading TPS file, defining speciesID
landmarks1 <- readland.tps(file = "input_ventralna.tps", specID = "imageID", negNA = TRUE)
landmarks <- estimate.missing(landmarks1, method = "TPS")
options(scipen=999)
#defining groups (spp)
Groups <- read.csv("input_list_ventral.csv", header = T, sep = ";")
Spp<-as.factor(Groups[,1])#list of specific groupings
ind<-as.factor(Groups[,2])#list of individuals for symmetric component extraction

#GPA
GPA <- gpagen(landmarks)
plotAllSpecimens(GPA$coords)
lands2d <- two.d.array(GPA$coords)
GPA3 <- arrayspecs(lands2d, 46, 2)
plotOutliers(GPA3, inspect.outliers = TRUE)

#symmetric component
lm<-read.csv("input_lmpairs_ventral.csv", header=TRUE, sep = ",")#landmarks pairs on each side of bilateral symmetry
lmm<-as.matrix(lm)

gdf <- geomorph.data.frame(shape = GPA3,
                           ind = ind)
sym <- bilat.symmetry(A = shape, ind = ind,
                      object.sym = TRUE,
                      land.pairs = lmm, data = gdf, RRPP = TRUE, iter = 149)
summary(sym)
sym<-sym$symm.shape
write.csv(sym, file="output_sym_ventral.csv")

sym2 <- two.d.array(sym)
write.csv(sym2, file ="output_sym2_ventral.csv")

#size
CS_specimens <- GPA$Csize
write.csv(CS_specimens, file="output_cs_specimens_ventral.csv")
CS_species <- aggregate(GPA$Csize, by = list(Spp), FUN = mean) 
write.csv(CS_species, file="output_cs_species_ventral.csv")

#logsize
Csize2d <- log(GPA$Csize)
write.csv(Csize2d, file="output_logcs_specimens_ventral.csv")
logCS_species <- aggregate(Csize2d, by = list(Spp), FUN = mean) 
write.csv(logCS_species, file = "output_logcs_species_ventral.csv")

#taking the average value of landmarks by species
Sppmean <-aggregate(sym2, by = list(Spp), FUN=mean)[,2:93] 
Sppmean3d <-arrayspecs(Sppmean, 46, 2, sep=".")
write.csv(Sppmean, file="output_sppmeanlandmarks_ventral.csv")

#PCA por especime
pca<-gm.prcomp(sym2, tol=0)
summary(pca)
pcs<-pca$x
write.csv(pcs, file = "output_pcs_specimens_ventral.csv")

#PCA (species)
pca <- gm.prcomp(Sppmean3d, tol=0)
summary(pca)
pcs<-pca$x
write.csv(pcs, file="output_pcs_species_ventral.csv")

#plotting mean shape
ref <- mshape(Sppmean3d)
findMeanSpec(GPA$coords)
#loading outlines for plotting shape changes on each component
lands<-readland.tps(file = "input_outline_ventral.TPS")
out<-warpRefOutline(file = "input_outline_ventral.csv", lands[,,1], ref)####

#plotting shape changes by points
plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "points")

#plotting shape changes by vectors
plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "vector")

##plotting shape changes by outlines
plotRefToTarget(ref,pca$shapes$shapes.comp1$min,outline=out$outline)
title("pc1 min")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max,outline=out$outline)
title("pc1 max")

plotRefToTarget(ref,pca$shapes$shapes.comp2$min,outline=out$outline)
title("pc2 min")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max,outline=out$outline)
title("pc2 max")


#MORPHOGEO MANDIBLE####
#carregar o arquivo TPS, definindo o speciesID
landmarks1 <- readland.tps(file = "input_mandible2.tps", specID = "imageID")
landmarks <- estimate.missing(landmarks1, method = "TPS")

#Definir os agrupamentos (por spp, ou genero, ou dieta, etc)
Groups <- read.csv("input_list_mandible.csv", header = T, sep = ";")
Spp<-as.factor(Groups[,1])#list of specific groupings

#GPA
GPA <- gpagen(landmarks)
plotAllSpecimens(GPA$coords)
lands2d <- two.d.array(GPA$coords)
Csize2d <- log(GPA$Csize)
GPA3 <- arrayspecs(lands2d, 16, 2)
plotOutliers(GPA3, inspect.outliers = TRUE)
write.csv(GPA3, file ="output_GPA_mandible.csv")

#size
CS_specimens <- GPA$Csize
write.csv(CS_specimens, file="output_cs_specimens_mandible.csv")
CS_species <- aggregate(GPA$Csize, by = list(Spp), FUN = mean) 
write.csv(CS_species, file="output_cs_species_mandible.csv")

#logsize
Csize2d <- log(GPA$Csize)
write.csv(Csize2d, file="output_logcs_specimens_mandible.csv")
logCS_species <- aggregate(Csize2d, by = list(Spp), FUN = mean) 
write.csv(logCS_species, file = "output_logcs_species_mandible.csv")


#Pegando o valor medio de landmarks por especie
Sppmean <-aggregate(lands2d, by = list(Spp), FUN=mean)[,2:33] 
Sppmean3d <-arrayspecs(Sppmean, 16, 2, sep=".")
write.csv(Sppmean, file="output_sppmeanlandmarks_mandible.csv")

#PCA por especime
pca<-gm.prcomp(lands2d, tol=0)
pcs<-pca$x
summary(pca)
Sppsd<-aggregate(pcs, by = list(Spp), FUN=sd)[,2:33]
Sppsd
write.csv(pcs, file = "output_pcs_specimens_mandible.csv")

#PCA por especie
pca <- gm.prcomp(Sppmean3d, tol=0)
summary(pca)
pcs<-pca$x
pcs
write.csv(pcs, file = "output_pcs_species_mandible.csv")

#plotting mean shape
ref <- mshape(GPA3)
plot(ref)

#outline
lands<-readland.tps(file = "input_outline_mandible.TPS")
out<-warpRefOutline(file = "input_outline_mandible.txt", lands[,,1], ref)####

plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "points")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "points")

plotRefToTarget(ref,pca$shapes$shapes.comp1$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$min, method = "vector")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max, method = "vector")

plotRefToTarget(ref,pca$shapes$shapes.comp1$min,outline=out$outline)
title("pc1 min")
plotRefToTarget(ref,pca$shapes$shapes.comp1$max,outline=out$outline)
title("pc1 max")

plotRefToTarget(ref,pca$shapes$shapes.comp2$min,outline=out$outline)
title("pc2 min")
plotRefToTarget(ref,pca$shapes$shapes.comp2$max,outline=out$outline)
title("pc2 max")

#PHYLOMORPHOSPACE####
#phylomorpho
#DORSAL

data_plot <- read.csv("pcs_species_dorsal.csv", header = TRUE, row.names = 1, sep = ";")
tree <- read.nexus("new_oxy.nexus")
plot(tree)
names<-name.check(tree,data_plot)
names
tree1<-drop.tip(tree, tip = names$tree_not_data)
name.check(tree1,data_plot)
data <- data_plot[match(tree1$tip.label, rownames(data_plot)), ]
tree2 <- force.ultrametric(tree1)
is.ultrametric(tree2)
write.nexus(tree2, file = "new_oxy_dropped.nexus")

#invertendo os pcs (lembrar que foram invertidos quando for colocar os outlines)
data_plot[,1] <- -data_plot[,1]

#par(xpd=NA) e par(xpd = FALSE) é para deixar os nomes das espécies ultrapassarem o limite do plot
par(xpd = NA)
phylomorphospace(tree2, data_plot[,c(1, 2)],xlab="PC1 (53.7%)",ylab="PC2 (21.2%)", label = "horizontal")
par(xpd = FALSE)
title(main="PC1xPC2_dorsal",
      font.main=3)

#VENTRAL
data_plot <- read.csv("pcs_species_ventral.csv", header = TRUE, row.names = 1, sep = ";")
data_plot
data_plot[,1] <- -data_plot[,1]
names<-name.check(tree2,data_plot)
names
par(xpd = NA)
phylomorphospace(tree2, data_plot[,c(1, 2)],xlab="PC1 (53.7%)",ylab="PC2 (17.8%)", label = "horizontal")
par(xpd = FALSE)
title(main="PC1xPC2_ventral",
      font.main=3)

#LATERAL
data_plot <- read.csv("pcs_species_lateral.csv", header = TRUE, row.names = 1, sep = ";")
#invertendo os pcs (lembrar que foram invertidos quando for colocar os outlines)
data_plot[,1] <- -data_plot[,1]

par(xpd = NA)
phylomorphospace(tree2, data_plot[,c(1, 2)],xlab="PC1 (36.5%)",ylab="PC2 (20.0%)", label = "horizontal")
par(xpd = FALSE)
title(main="PC1xPC2_lateral",
      font.main=3)

#MANDIBLE
data_plot <- read.csv("pcs_species_mandible.csv", header = TRUE, row.names = 1, sep = ";")
names<-name.check(tree2,data_plot)
names
tree3<-drop.tip(tree2, tip = names$tree_not_data)
name.check(tree3,data_plot)

data_plot[,1] <- -data_plot[,1]
data_plot[,2] <- -data_plot[,2]

par(xpd = NA)
phylomorphospace(tree3, data_plot[,c(1, 2)],xlab="PC1 (38.6%)",ylab="PC2 (22.5%)", label = "horizontal")
par(xpd = FALSE)
title(main="PC1xPC2_mandible",
      font.main=3)

#MAPPING####

#testing correlation between centroid sizes before mapping
dat <- read.csv("output_cs_new.csv", header = T, sep = ";")

corr.test(dat[,c("dorsal","lateral","ventral","mandible")],
          method = "spearman")

ICC(dat[,c("dorsal", "lateral", "ventral", "mandible")])

#highly correlated, so we are going to use an average of all four centroid sizes
dat$CS_mean <- rowMeans(dat[,c("dorsal","lateral","ventral","mandible")], na.rm = TRUE)
dat$log_CS <- log(dat$CS_mean)

write.csv(dat, "output_cs_meanlog.csv", row.names = FALSE)

#mapping
data <- read.csv("linear_measurements2.csv", header = TRUE, row.names = 1, sep = ";")
names<-name.check(tree,data)
plot(tree)
names
tree1<-drop.tip(tree, tip = names$tree_not_data)
name.check(tree2,data)
data <- data[match(tree1$tip.label, rownames(data)), ]
tree2 <- force.ultrametric(tree1)
is.ultrametric(tree2)
write.nexus(tree2, file = "new_oxy_dropped.nexus")

#variables
# vetores nomeados
rt <- setNames(data$rostraltube_sc, rownames(data))
pi <- setNames(data$postincisive_sc, rownames(data))
size <- setNames(data$log_CS, rownames(data))

#paleta
vir_pal <- viridis(100)

# rostraltube_sc
obj_rt <- contMap(tree2, rt, plot = FALSE)
obj_rt <- setMap(obj_rt, vir_pal)   # <- JEITO CERTO
plot(obj_rt, lwd = 4.5, fsize = 0.8, legend = TRUE,
     main = "rostraltube_sc (viridis)")
title("Rostral tube (scaled)")

# postincisive_sc
obj_pi <- contMap(tree2, pi, plot = FALSE)# cria o objeto sem plotar
obj_pi <- setMap(obj_pi, vir_pal)# FORÇA a paleta no objeto
plot(obj_pi, lwd = 4.5, fsize = 0.8, legend = TRUE,
     main = "postincisive_sc (viridis)")
title("Post incisive nasal length (scaled)")

# size
obj_size <- contMap(tree2, size, plot = FALSE)# cria o objeto sem plotar
obj_size <- setMap(obj_size, vir_pal)# FORÇA a paleta no objeto
plot(obj_size, lwd = 4.5, fsize = 0.8, legend = TRUE,
     main = "postincisive_sc (viridis)")
title("Size (logCS)")

#PGLS####

options(scipen = 999)
dat <- read.csv("linear_measurements2.csv", header = TRUE, sep = ";")
## -------------------------------
## Comparative data
## -------------------------------

cd <- comparative.data(data=dat, phy=tree2, names.col=sp, vcv=TRUE, vcv.dim=3, warn.dropped = TRUE)

## -------------------------------
## PGLS
## -------------------------------

m_rt <- pgls(
  rt ~ pi,
  data = cd,
  lambda = "ML",
  param.CI = 0.95
)

summary(m_rt)

#alometry
m_allo1 <- pgls(
  rt ~ size,
  data = cd,
  lambda = "ML",
  param.CI = 0.95
)

summary(m_allo1)

m_allo2 <- pgls(
  pi ~ size,
  data = cd,
  lambda = "ML",
  param.CI = 0.95
)

summary(m_allo2)

## -------------------------------
## Outliers filogenéticos
## -------------------------------
res <- residuals(m_rt, phylo = TRUE)
res_std <- res / sqrt(var(res))[1]
names(res_std) <- rownames(m_rt$residuals)

outliers <- names(res_std)[abs(res_std) > 3]
outliers

## -------------------------------
## Diagnósticos e plots
## -------------------------------
dev.off()
plot.pgls(m_rt)

## Visualização simples (parcial, não filogenética)

plot(dat$rostraltube_sc, dat$postincisive_sc,
     xlab = "Rostral length (mean per species)",
     ylab = "Rostral tube (mean per species)",
     pch = 19)

