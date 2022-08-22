#### All code - Final 

#~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
### 2.3(.2) LiDAR data pre-processing:----
library(lidR)
library(rTLS)
library(leafR)

BGs <- c(1:5)
PCs <- c(1:45)

## BUGGs
# Height normalise:
mycsf <- csf(TRUE, 1, 1, time_step = 1)
for(i in 1:length(BGs)){ 
  PC <- readLAS(readLAS(paste0("BUGG",i,".las")))
  PC <- classify_ground(PC, mycsf)     
  PC <- normalize_height(PC, algorithm = knnidw(k = 6L, p = 2))
  writeLAS(PC, paste0("BUGG",i,"_normalised.las"))
}

# Remove ground:
for(i in 1:length(BGs)){
  Point_count<-readLAS(paste0("BUGG",i,"_normalised.las"))
  Point_count <- Point_count[which(Point_count@data$Classification != 2)]
  writeLAS(Point_count, paste0("BUGG",i,"_noGround.las"))
}

# Clip:
for(i in 1:length(BGs)){
  Point_count<-readLAS(paste0("BUGG",i,"_noGround.las"))
  Point_count <- clip_circle(Point_count, radius = 20, xcenter = 0, ycenter = 0)
  writeLAS(Point_count, paste0("BUGG",i,"_clipped.las"))
}

## PCs
# Height normalise:
for(i in 1:length(PCs)){ 
  PC <- readLAS(readLAS(paste0("BUGG",i,".las")))
  PC <- classify_ground(PC, mycsf)     
  PC <- normalize_height(PC, algorithm = knnidw(k = 6L, p = 2))
  writeLAS(PC, paste0("PC",i,"_normalised.las"))
}

# Remove ground:
for(i in 1:length(PCs)){
  Point_count<-readLAS(paste0("PC",i,"_normalised.las"))
  Point_count <- Point_count[which(Point_count@data$Classification != 2)]
  writeLAS(Point_count, paste0("PC",i,"_noGround.las"))
}

# Clip:
for(i in 1:length(BGs)){
  Point_count<-readLAS(paste0("PC",i,"_noGround.las"))
  Point_count <- clip_circle(Point_count, radius = 20, xcenter = 0, ycenter = 0)
  writeLAS(Point_count, paste0("PC",i,"_clipped.las"))
}
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###### 2.4. Habitat complexity analyses ------------------------------------------
### 2.4.1 Canopy height: 
## BUGGs
mean_canopy_BUGG <- canopy_c.o.v_BUGG <- vector()
max_canopy_BUGG <- c()
for(i in 1:5){
  PC <- readLAS(paste0("BUGG",i,"_clipped.las"))
  PC <- PC@data[,1:3]
  voxels <- voxels(PC, edge_length = c(1,1,150))
  for(d in 1:length(voxels$voxels)){
    if(voxels$voxels$N[d] >= 10){
      subset <- dplyr::filter(voxels$cloud, 
                              (voxels$cloud$X > (voxels$voxels$X[d] - 0.5)) & 
                                (voxels$cloud$X < (voxels$voxels$X[d] + 0.5)) & 
                                (voxels$cloud$Y > (voxels$voxels$Y[d] - 0.5)) &
                                (voxels$cloud$Y < (voxels$voxels$Y[d] + 0.5)))
      maxz <- max(subset$Z)
      max_canopy_BUGG <- c(max_canopy_BUGG, maxz) 
    }}
  mean_canopy_BUGG <- c(mean_canopy_BUGG, mean(max_canopy_BUGG))
  cov <- (sd(max_canopy_BUGG)/mean(max_canopy_BUGG))
  canopy_c.o.v_BUGG <- c(canopy_c.o.v_BUGG, cov)
}

BUGG.veg.metrics <- cbind(mean_canopy_BUGG, canopy_c.o.v_BUGG)

## PCs
mean_canopy.PC <- canopy_c.o.v.PC <- vector()
max_canopy.PC <- vector()
for(i in 1:length(PCs)){
  PC <- readLAS(paste0("LiDAR/PCs heightNorm/PC",i,"_noGround.las"))
  PC <- PC@data[,1:3]
  voxels <- voxels(PC, edge_length = c(1,1,150))
  for(d in 1:length(voxels$voxels)){
    if(voxels$voxels$N[d] >= 10){
      subset <- dplyr::filter(voxels$cloud, 
                              (voxels$cloud$X > (voxels$voxels$X[d] - 0.5)) & 
                                (voxels$cloud$X < (voxels$voxels$X[d] + 0.5)) & 
                                (voxels$cloud$Y > (voxels$voxels$Y[d] - 0.5)) &
                                (voxels$cloud$Y < (voxels$voxels$Y[d] + 0.5)))
      maxz <- max(subset$Z)
      max_canopy.PC <- c(max_canopy.PC, maxz) 
    }}
  mean_canopy.PC <- c(mean_canopy.PC, mean(max_canopy.PC))
  cov <- (sd(max_canopy.PC)/mean(max_canopy.PC))
  canopy_c.o.v.PC <- c(canopy_c.o.v.PC, cov)
}

PC.veg.metrics <- cbind(mean_canopy.PC, canopy_c.o.v.PC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2.4.2 LAD-based metrics: 
## BUGGs:
cc_handheld_BUGG <- All_HH_LAI_BUGG <-vector()
for(i in 1:5){
  VOXELS_LAD <- lad.voxels(paste0("LiDAR/PCs heightNorm/BUGG",i,"_clipped.las"))
  lad_profile <- lad.profile(VOXELS_LAD)
  plot_LAI <- lai(lad_profile)
  plot.mean <- mean(lad_profile$lad)
  plot.sd <- sd(lad_profile$lad)
  plot.cov <- plot.sd/plot.mean
  cc_handheld_BUGG[i] <- lai(lad_profile, min = 3, max = 150)
  understory.lai <- lai(lad_profile, min = 0, max = 5)
  mid.lai<-lai(lad_profile,min=5,max=10)
  canopy.lai<-lai(lad_profile,min=10,max=20)
  
  under <- dplyr::filter(lad_profile, height < 5)
  u.mean <- mean(under$lad)
  u.sd <- sd(under$lad)
  u.cov <- u.sd/u.mean
  mid <- dplyr::filter(lad_profile, (height > 5) & (height < 10))
  m.mean <- mean(mid$lad)
  m.sd <- sd(mid$lad)
  m.cov <- m.sd/m.mean
  canopy <- dplyr::filter(lad_profile, (height > 10) & (height < 20))
  c.mean <- mean(canopy$lad)
  c.sd <- sd(canopy$lad)
  c.cov <- c.sd/c.mean
  
  All_HH_LAI_BUGG<-rbind(All_HH_LAI_BUGG, c(plot_LAI,understory.lai,mid.lai,canopy.lai,u.cov,m.cov,c.cov)) 
}

colnames(All_HH_LAI_BUGG) <- c("plai","ulai","mlai","clai","uladcv",
                          "mladcv","cladcv")

BUGG.veg.metrics <- cbind(All_HH_LAI_BUGG, BUGG.veg.metrics)

## PCs
cc_handheld <- All_HH_LAI <-vector()
for(i in 1:length(PCs)){
  VOXELS_LAD <- leafR::lad.voxels(paste0("LiDAR/PCs heightNorm/PC",i,"_noGround.las"))
  lad_profile <- leafR::lad.profile(VOXELS_LAD)
  plot_LAI <- leafR::lai(lad_profile)
  plot.mean <- mean(lad_profile$lad)
  plot.sd <- sd(lad_profile$lad)
  plot.cov <- plot.sd/plot.mean
  cc_handheld[i] <- leafR::lai(lad_profile, min = 3, max = 150)
  understory.lai <- leafR::lai(lad_profile, min = 0, max = 5)
  mid.lai<-leafR::lai(lad_profile,min=5,max=10)
  canopy.lai<-leafR::lai(lad_profile,min=10,max=20)
  
  under <- dplyr::filter(lad_profile, height < 5)
  u.mean <- mean(under$lad)
  u.sd <- sd(under$lad)
  u.cov <- u.sd/u.mean
  mid <- dplyr::filter(lad_profile, (height > 5) & (height < 10))
  m.mean <- mean(mid$lad)
  m.sd <- sd(mid$lad)
  m.cov <- m.sd/m.mean
  canopy <- dplyr::filter(lad_profile, (height > 10) & (height < 20))
  c.mean <- mean(canopy$lad)
  c.sd <- sd(canopy$lad)
  c.cov <- c.sd/c.mean
  
  
  All_HH_LAI<-rbind(All_HH_LAI, c(plot_LAI,understory.lai,mid.lai,canopy.lai,u.cov,m.cov,c.cov)) 
}

PC.veg.metrics <- cbind(All_HH_LAI, PC.veg.metrics)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2.4.3 Fractal dimension 
# BUGGs:
## fractal geometry per strata block (understory, mid, canopy, high canopy):
layer_fractals_BUGG <- vector()
for(i in 1:5){
  PC <- readLAS(paste0("BUGG",i,"_clipped.las"))
  PC <- PC@data[,1:3]
  
  under <- dplyr::filter(PC, (PC[,3] > 0) & (PC[,3] <= 5))
  ufractals <- rTLS::voxels_counting(under, min_size = 0.1)
  udimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = ufractals)
  under_fractal <- udimension$coefficients[2]
  mid <- dplyr::filter(PC, (PC[,3] > 5) & (PC[,3] <= 10))
  mfractals <- rTLS::voxels_counting(mid, min_size = 0.1)
  mdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = mfractals)
  mid_fractal <- mdimension$coefficients[2]
  canopy <- dplyr::filter(PC, (PC[,3] > 10) & (PC[,3] <= 20))
  cfractals <- rTLS::voxels_counting(canopy, min_size = 0.1)
  cdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = cfractals)
  canopy_fractal <- cdimension$coefficients[2]
  if(max(PC[,3] <= 20)){high_fractal <- 0}
  else {
    high <- dplyr::filter(PC, (PC[,3] > 20) & (PC[,3] < 100))
    hfractals <- rTLS::voxels_counting(high, min_size = 0.1)
    hdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = hfractals, na.rm=T)
    high_fractal <- hdimension$coefficients[2]
  }
  layer_fractals_BUGG <- rbind(layer_fractals_BUGG, c(under_fractal, mid_fractal, canopy_fractal)) 
  print(i)
}

# Plot-wide:
plot_fractals_BUGG <- vector()
for(i in 1:5){
  PC <- readLAS(paste0("BUGG",i,"_clipped.las"))
  PC <- PC@data[,1:3]
  fractals <- rTLS::voxels_counting(PC, min_size = 0.1)
  dimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = fractals)
  plot_fractal <- dimension$coefficients[2]
  plot_fractals_BUGG <- c(plot_fractals_BUGG, plot_fractal)
  print(i)
}

BUGG.veg.metrics <- cbind(plot_fractals_BUGG, layer_fractals_BUGG, BUGG.veg.metrics)
# PCs:
## fractal geometry per strata block (understory, mid, canopy, high canopy):
layer_fractals <- vector()
for(i in 1:length(PCs)){
  PC <- readLAS(paste0("PC",i,"_noGround.las"))
  PC <- PC@data[,1:3]
  
  under <- dplyr::filter(PC, (PC[,3] > 0) & (PC[,3] <= 5))
  ufractals <- rTLS::voxels_counting(under, min_size = 0.1)
  udimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = ufractals)
  under_fractal <- udimension$coefficients[2]
  mid <- dplyr::filter(PC, (PC[,3] > 5) & (PC[,3] <= 10))
  mfractals <- rTLS::voxels_counting(mid, min_size = 0.1)
  mdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = mfractals)
  mid_fractal <- mdimension$coefficients[2]
  canopy <- dplyr::filter(PC, (PC[,3] > 10) & (PC[,3] <= 20))
  cfractals <- rTLS::voxels_counting(canopy, min_size = 0.1)
  cdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = cfractals)
  canopy_fractal <- cdimension$coefficients[2]
  if(max(PC[,3] <= 20)){high_fractal <- 0}
  else {
    high <- dplyr::filter(PC, (PC[,3] > 20) & (PC[,3] < 100))
    hfractals <- rTLS::voxels_counting(high, min_size = 0.1)
    hdimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = hfractals, na.rm=T)
    high_fractal <- hdimension$coefficients[2]
  }
  layer_fractals <- rbind(layer_fractals, c(under_fractal, mid_fractal, canopy_fractal, high_fractal)) 
  print(i)
}

layer_fractals <- layer_fractals[,1:3]
  
## Plot-wide fractal geometry:
plot_fractals <- vector()
for(i in 1:length(PCs)){
  PC <- readLAS(paste0("PC",i,"_noGround.las"))
  PC <- PC@data[,1:3]
  fractals <- rTLS::voxels_counting(PC, min_size = 0.1)
  dimension <- lm(log10(N_voxels) ~ log10(1/(Edge.X*Edge.Y*Edge.Z)), data = fractals)
  plot_fractal <- dimension$coefficients[2]
  plot_fractals <- c(plot_fractals, plot_fractal)
  print(i)
}

PC.veg.metrics <- cbind(layer_fractals, plot_fractals, PC.veg.metrics)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2.4.5. Principal components analysis 
pca.veg <- prcomp(all.veg[,-14],
                       center = TRUE,
                       scale. = TRUE)

# Biplot:
veg.biplot <- ggbiplot(pca.veg,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = all.veg$block,
                  ellipse = TRUE,
                  circle = FALSE,
                  ellipse.prob = 0.80)
veg.biplot <- veg.biplot + geom_text(label = obs.labels, nudge_x = -0.2, nudge_y = 0.15)
veg.biplot <- veg.biplot + scale_color_discrete(name = '')
veg.biplot <- veg.biplot + theme(legend.direction = 'vertical', legend.position = 'right')

print(veg.biplot)

## ANOVA & Tukey HSD:

coords.pca.veg <- as.data.frame(factoextra::get_pca_ind(pca.veg)$coord[,1:2])
coords.pca.veg <- cbind(coords.pca.veg, all.veg$block)

veg.aov1 <- aov(PC1 ~ block, data = coords.pca.veg)
TukeyHSD(veg.aov1)

veg.aov2 <- aov(PC2 ~ block, data = coords.pca.veg)
TukeyHSD(veg.aov2)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##### 2.5. Avian community data collection -------------------------------------
## ((2.5.1)) Detection filtering

library(tidyverse)
library(dplyr)
library(fs)
library(lubridate)
library(stringr)
library(janitor)

# Read in all BUGG detections from folder:
BUGGs <- fs::dir_ls("BUGG_Detections")
BUGG_reads <- list()
for(i in seq_along(BUGGs)){
  BUGG_reads[[i]] <- read.csv(file = BUGGs[[i]])
}
BUGG_reads <- setNames(BUGG_reads, BUGGs)

for(i in seq_along(BUGG_reads)){
BUGG_reads[[i]]$Date <- as.Date(BUGG_reads[[i]]$Date) 
}


BUGG_reads.bsq <- BUGG_reads[[1]]
length(unique(BUGG_reads.7.bsq$Date)) # check no. of days

BUGG_reads.hamb <- BUGG_reads[[2]]
BUGG_reads.hamb <- BUGG_reads.hamb %>% filter(Date < ymd(20220322))
length(unique(BUGG_reads.hamb$Date))

BUGG_reads.irle <- BUGG_reads[[3]]
BUGG_reads.irle <- BUGG_reads.irle %>% filter(Date < ymd(20220315))
length(unique(BUGG_reads.irle$Date))

BUGG_reads.irli <- BUGG_reads[[4]]
BUGG_reads.irli <- BUGG_reads.irli %>% filter(Date < ymd(20220314))
length(unique(BUGG_reads.irli$Date))

BUGG_reads.res <- BUGG_reads[[5]]
BUGG_reads.res <- BUGG_reads.res %>% filter(Date < ymd(20220316))
length(unique(BUGG_reads.res$Date))

BUGG_reads <- list(BUGG_reads.bsq, BUGG_reads.hamb, BUGG_reads.irle, BUGG_reads.irli, BUGG_reads.res)


# Function below takes a list of detections as outputted by the loop above 
# (each list element is a separate device/location), as well as a threshold 
# (no. of days a species must be present to be included in final output), 
# filters by this threshold + a confidence threshold (>=0.7),
# and outputs a list, the elements of which contain a vector of scientific names 
# of species present. 

make.species.list <- function(BUGG.readz, thresh){
  final.species <- list()
  daysall <- list()
  for(i in seq_along(BUGG.readz)){   
    BUGG.readz[[i]]$Date <- as.Date(BUGG.readz[[i]]$Date) 
    unq.days <- unique(BUGG.readz[[i]]$Date) 
    days <- list()
    for(j in seq_along(unq.days)){
      day <- BUGG.readz[[i]] %>% filter(Date == unq.days[j])  
      days[[j]] <- day                              
    }
    daysall[[i]] <- days 
  }
  
  daycounts.all <- list()
  for(i in seq_along(daysall)){
    unq.species <- unique(BUGG.readz[[i]]$Scientific.name)
    daycount <- rep(0, times = length(unq.species))
    for(j in seq_along(unq.species)){
      for(k in seq_along(daysall[[i]])){
        if((unq.species[j] %in% daysall[[i]][[k]]$Scientific.name) == T){
          daycount[j] <- daycount[j] + 1
        }
      }
    }
    daycounts.all[[i]] <- daycount
  }
  
  
  final.species.alltime <- list()
  for(i in seq_along(daycounts.all)){
    unq.species <- unique(BUGG.readz[[i]]$Scientific.name)
    all.species <- vector()
    for(j in seq_along(unq.species)){
      if(daycounts.all[[i]][[j]] >= thresh){
        all.species <- c(all.species, unq.species[j])
      }
    }
    final.species.alltime[[i]] <- all.species
  }
  return(final.species.alltime)
}

species.list <- make.species.list(BUGG_reads, 3)

species.bsq <- species.list[[1]]
species.hamb <- species.list[[2]]
species.irle <- species.list[[3]]
species.irli <- species.list[[4]]
species.res <- species.list[[5]]

# Make occurrence matrix:
spec.vec <- c(species.bsq, species.hamb, species.irle, species.irli, species.res)
unique.species <- unique(spec.vec)

species.sites <- list(species.bsq, species.hamb, species.irle, species.irli, species.res)

species.x.sites <- list(rep(0,111),rep(0,111),rep(0,111),rep(0,111),rep(0,111))
for(i in seq_along(species.sites)){
  for(j in seq_along(unique.species)){
    if(unique.species[j] %in% species.sites[[i]] == T){
      species.x.sites[[i]][[j]] <- 1
    }}}


##### 2.6. Functional diversity --------------------------------------------------
## 2.6.1. Trait data

# Load AVONET datasets:
AVONET1 <- read.csv('AVONET/AVONET1_BirdLife.csv')
AVONET2 <- read.csv('AVONET/AVONET2_eBird.csv')
colnames(AVONET2) <- colnames(AVONET1)

# Function to extract traits:
AVONET_extract_new2 <- function(newspecies){
  oldspecies <- data.frame(AVONET1[10215,])
  for(i in seq_along(newspecies)){
    print(i)
    if((newspecies[i] %in% AVONET1[,1]) == TRUE){
      m <- match(newspecies[i], AVONET1[,1])
      oldspecies <- rbind(oldspecies, AVONET1[m, ])
    } else {
      p <- match(newspecies[i], AVONET2[,1])
      oldspecies <- rbind(oldspecies, AVONET2[p, ])
    }
  }
  oldspecies <- oldspecies[-1,]
  return(oldspecies)
}

alltraits <- AVONET_extract_new2(unique.species)

# Function to extract PC axes from PCA (removal of multicollinearity) as new traits from biometric traits:
transform.biometry <- function(trait.matrix){
  trophic <- log(trait.matrix[,1:4])
  locom <- log(trait.matrix[,c(5:8,10)])
  
  PCA1.trophic <- prcomp(trophic, center = TRUE, scale. = TRUE)
  coords.trophic <- factoextra::get_pca_ind(PCA1.trophic)
  trophic.shape <- coords.trophic$coord[,2]
  trophic.size <- coords.trophic$coord[,1]
  
  PCA1.locom <- prcomp(locom, center = TRUE, scale. = TRUE)
  coords.locom <- factoextra::get_pca_ind(PCA1.locom)
  locom.shape <- coords.locom$coord[,2]
  locom.size <- coords.locom$coord[,1]
  
  pcsize <- cbind(trophic.size, locom.size)
  PCA2.size <- prcomp(pcsize, center = TRUE, scale. = TRUE)
  coords.size <- factoextra::get_pca_ind(PCA2.size)
  size <- coords.size$coord[,1]
  
  new.variables <- cbind(size, trophic.shape, locom.shape)
  rownames(new.variables) <- trait.matrix[,1]
  return(new.variables)
}

# Format trait matrix for FD calculation:

alltraits.short <- alltraits[,c(10:20,26:31)]
alltime.PCA <- transform.biometry(alltraits.short)

foraging.strata <- read.csv("Other data/foraging_strata_alltime.csv")

all.traits.FD <- cbind(alltime.PCA, alltime.trait.matrix.short[,12:17], foraging.strata[,2])
rownames(alltraits.FD) <- alltraits$Species1
alltraits.FD[,4] <- as.factor(alltraits.FD[,4])
alltraits.FD[,5] <- as.factor(alltraits.FD[,5])
alltraits.FD[,6] <- as.factor(alltraits.FD[,6])
alltraits.FD[,7] <- as.factor(alltraits.FD[,7])
alltraits.FD[,8] <- as.factor(alltraits.FD[,8])
alltraits.FD[,9] <- as.factor(alltraits.FD[,9])
alltraits.FD[,10] <- as.factor(alltraits.FD[,10])

colnames(all.traits.alltime)[8] <- "Niche.Trophic"
colnames(all.traits.alltime)[10] <- "Foraging.Preference"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 2.6.3. alpha-diversity indices
#### Follow mFD workflow:
library(mFD)

# Functional matrix: species x traits 
functional.matrix <- all.traits.FD
rownames(functional.matrix) <- gsub(" ", ".", rownames(functional.matrix)) 

# Data frame of trait types
trait.type <- as.data.frame(cbind(colnames(functional.matrix), c("Q", "Q", "Q", "N", "O", "O", "N", "N", "N","N")))
colnames(trait.type) <- c("trait_name", "trait_type")

# Species traits summary
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = trait.type,   
  sp_tr      = functional.matrix, 
  stop_if_NA = TRUE)

# Assemblages * species 
assemblages_species_sample <- species.x.sites
colnames(assemblages_species_sample) <- gsub(" ", ".", colnames(assemblages_species_sample)) 
row.names(assemblages_species_sample) <- c("hamburgo", "bosque", "irlandae", "irlandai", "restoration")

# Summary of the assemblages * species dataframe
asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w = as.matrix(assemblages_species_sample))

# Species richness per assemblage
asb_sp_summ$"asb_sp_richn"

# Functional distance between species
sp_dist <- mFD::funct.dist(
  sp_tr         = functional.matrix,
  tr_cat        = trait.type,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

# Compute multidimensional functional spaces and their quality
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

# Quality metrics of spaces
round(fspaces_quality$"quality_fspaces", 3)  

# Illustrating the quality of the multidimensional functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d", "pcoa_7d", "pcoa_8d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# Test correlation between functional axes and traits
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"
tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = functional.matrix, 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")], 
  plot           = TRUE)

# Print traits with significant effect:
tr_faxes$"tr_faxes_stat"[which(tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

tr_faxes$"tr_faxes_plot"

# Plot functional space
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  faxes           = NULL,
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot$patchwork

# Compute functional diversity indices and plot them
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")],
  asb_sp_w         = as.matrix(assemblages_species_sample),
  ind_vect         = c("feve", "fric", "fdis"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values

# Plot functional diversity indices
plots_alpha1 <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("hamburgo"),
  ind_nm                   = c("feve", "fric", "fdis"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#7cad02"),
  color_vert               = c(pool = "grey50", asb1 = "#7cad02"),
  fill_sp                  = c(pool = NA, asb1 = "#7cad02"),
  fill_vert                = c(pool = NA, asb1 = "#7cad02"),
  color_ch                 = c(pool = NA, asb1 = "#7cad02"),
  fill_ch                  = c(pool = "white", asb1 = "#7cad02"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
aplot1 <- plots_alpha1$fric$PC1_PC2

plots_alpha2 <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("irlandai"),
  ind_nm                   = c("feve", "fric", "fdis"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#05bfc4"),
  color_vert               = c(pool = "grey50", asb1 = "#05bfc4"),
  fill_sp                  = c(pool = NA, asb1 = "#05bfc4"),
  fill_vert                = c(pool = NA, asb1 = "#05bfc4"),
  color_ch                 = c(pool = NA, asb1 = "#05bfc4"),
  fill_ch                  = c(pool = "white", asb1 = "#05bfc4"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
aplot2 <- plots_alpha2$fric$PC1_PC2

plots_alpha3 <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("irlandae"),
  ind_nm                   = c("feve", "fric", "fdis"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#c47d05"),
  color_vert               = c(pool = "grey50", asb1 = "#c47d05"),
  fill_sp                  = c(pool = NA, asb1 = "#c47d05"),
  fill_vert                = c(pool = NA, asb1 = "#c47d05"),
  color_ch                 = c(pool = NA, asb1 = "#c47d05"),
  fill_ch                  = c(pool = "white", asb1 = "#c47d05"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
aplot3 <- plots_alpha3$fric$PC1_PC2

plots_alpha4 <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("restoration"),
  ind_nm                   = c("feve", "fric", "fdis"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#c77cff"),
  color_vert               = c(pool = "grey50", asb1 = "#c77cff"),
  fill_sp                  = c(pool = NA, asb1 = "#c77cff"),
  fill_vert                = c(pool = NA, asb1 = "#c77cff"),
  color_ch                 = c(pool = NA, asb1 = "#c77cff"),
  fill_ch                  = c(pool = "white", asb1 = "#c77cff"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
aplot4 <- plots_alpha4$fric$PC1_PC2

plots_alpha5 <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("bosque"),
  ind_nm                   = c("feve", "fric", "fdis"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#f8766d"),
  color_vert               = c(pool = "grey50", asb1 = "#f8766d"),
  fill_sp                  = c(pool = NA, asb1 = "#f8766d"),
  fill_vert                = c(pool = NA, asb1 = "#f8766d"),
  color_ch                 = c(pool = NA, asb1 = "#f8766d"),
  fill_ch                  = c(pool = "white", asb1 = "#f8766d"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22),
  shape_centroid_fdiv      = c(asb1 = 22),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 
aplot5 <- plots_alpha5$fric$PC1_PC2


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2.6.2. beta-diversity indices
## Calculate functional beta diversity:
beta.dissim <- beta.fd.multidim(sp_faxes_coord = sp_faxes_coord[,1:5], asb_sp_occ = assemblages_species_sample, 
                                check_input = T, details_returned = T, betapart_para = T)
bdiv.matrix <- beta.dissim$pairasb_fbd_indices$jac_diss

# Plot:
# Make correlation plot:
bdiv.matrix.mod <- as.data.frame(as.matrix(bdiv.matrix))
colnames(bdiv.matrix.mod) <- rownames(bdiv.matrix.mod) <- c("e.hamburgo","a.bosque","d.irlandae","c.irlandai","b.restoration")

bdiv.matrix.mod <- bdiv.matrix.mod %>%  as_tibble(rownames = "A") %>% 
  pivot_longer(-A, names_to = "B", values_to = "dissim")

# Heatplot: 
tick.labels <- c("Bosque","Restoration","Irlanda 1","Irlanda 2","Hamburgo")
ramp <- c(rep("#08306B", times = 11), rev(RColorBrewer::brewer.pal(9, "Blues")))

heatplot <- bdiv.matrix.mod %>% ggplot(aes(x=A, y=B, fill= dissim )) +
  geom_tile() +
  geom_text(aes(label = round(dissim, digits = 3)), color = "black", size = 4) +
  scale_fill_gradientn(colors = ramp, 
                       guide = guide_colourbar(
                         barheight = 15,
                         title = expression(paste("",beta, " dissim.")), 
                         title.position = "top", title.hjust = -3)) +
  xlab(label = NULL) +
  ylab(label = NULL) +
  ggtitle(label = "a") +
  scale_x_discrete(labels= tick.labels) +
  scale_y_discrete(labels= tick.labels) +
  theme(plot.title = element_text(size = 15),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13),
        legend.title = element_text(size=13),
        legend.text = element_text(size=13))


##### 2.7. Habitat-trait associations --------------------------------------------
### 2.7.1. Response of F-beta to habitat complexity

## Regression with euclidean distances in PC space:
# Get coords from vegetation pca:
coords.pca.veg <- factoextra::get_pca_ind(pca.veg)$coord[,1:2] 

# Get average coords (centroid) for each block:
coords.veg.irlandai <- as.data.frame(coords.pca.veg[c(1,5,6,11,12,13,15,16,17,18,30,47), ])
coords.veg.irlandae <- as.data.frame(coords.pca.veg[c(2,3,4,7,8,9,10,14,27,28,29,48), ])
coords.veg.hamburgo <- as.data.frame(coords.pca.veg[30:45,])
coords.veg.bosque <- as.data.frame(coords.pca.veg[23:26,])
coords.veg.restoration <- as.data.frame(coords.pca.veg[19:22,])

PCA.mean.irli <- c(mean(coords.veg.irlandai$Dim.1), mean(coords.veg.irlandai$Dim.2))
PCA.mean.irle <- c(mean(coords.veg.irlandae$Dim.1), mean(coords.veg.irlandae$Dim.2))
PCA.mean.hamb <- c(mean(coords.veg.hamburgo$Dim.1), mean(coords.veg.hamburgo$Dim.2))
PCA.mean.bsq <- c(mean(coords.veg.bosque$Dim.1), mean(coords.veg.bosque$Dim.2))
PCA.mean.res <- c(mean(coords.veg.restoration$Dim.1), mean(coords.veg.restoration$Dim.2))

# Combine into dataframe
veg.pc.means <- rbind(PCA.mean.hamb, PCA.mean.bsq, PCA.mean.irle, PCA.mean.irli, PCA.mean.res)
colnames(veg.pc.means) <- c("PC1","PC2")
rownames(veg.pc.means) <- c("hamburgo","bosque","irlandae","irlandai","restoration")

# Funct. for euclidean distance:
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# Get combinations of blocks to match dissimilarity matrix:
combn(nrow(veg.pc.means),2)
veg.euclid <- combn(nrow(veg.pc.means), 2, function(x) euclidean(veg.pc.means[x[1],1:2], veg.pc.means[x[2],1:2]))
combos <- c("h.b","h.ie","h.ii","h.r","b.ie","b.ii","b.r","ie.ii","ie.r","ii.r")

# Make dataframe of functional dissimilarity values and corresponding PC distances:
regr <- cbind(veg.euclid, as.vector(bdiv.matrix))
rownames(regr) <- combos
colnames(regr) <- c("veg.euclid", "bdiv")
regr <- as.data.frame(regr)

# Construct regression:
habitat.community.model1 <- lm(bdiv ~ veg.euclid, regr)
summary(habitat.community.model1)
# Plot:
regr.plot1 <- regr %>% ggplot(aes(x = veg.euclid, y = bdiv)) +
  geom_point() +  
  stat_smooth(method = "lm", color = "#2171B5") +
  xlab(label = "Habitat structural dissimilarity") +
  ylab(label = expression(paste("Functional ",beta, " diversity (dissimilarity)"))) +
  ggtitle(label = "b") +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size=13),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10)) 

## Regression with individual PC metrics:

# Get differences for each veg metric individually: 
plai.diff <- combn(length(veg.means.all$plai), 2, function(x) (veg.means.all$plai[x[1]] - veg.means.all$plai[x[2]]))
ulai.diff <- combn(length(veg.means.all$ulai), 2, function(x) (veg.means.all$ulai[x[1]] - veg.means.all$ulai[x[2]]))
mlai.diff <- combn(length(veg.means.all$mlai), 2, function(x) (veg.means.all$mlai[x[1]] - veg.means.all$mlai[x[2]]))
clai.diff <- combn(length(veg.means.all$clai), 2, function(x) (veg.means.all$clai[x[1]] - veg.means.all$clai[x[2]]))
uladcv.diff <- combn(length(veg.means.all$uladcv), 2, function(x) (veg.means.all$uladcv[x[1]] - veg.means.all$uladcv[x[2]]))
mladcv.diff <- combn(length(veg.means.all$mladcv), 2, function(x) (veg.means.all$mladcv[x[1]] - veg.means.all$mladcv[x[2]]))
cladcv.diff <- combn(length(veg.means.all$cladcv), 2, function(x) (veg.means.all$cladcv[x[1]] - veg.means.all$cladcv[x[2]]))
ufrd.diff <- combn(length(veg.means.all$ufrd), 2, function(x) (veg.means.all$ufrd[x[1]] - veg.means.all$ufrd[x[2]]))
mfrd.diff <- combn(length(veg.means.all$mfrd), 2, function(x) (veg.means.all$mfrd[x[1]] - veg.means.all$mfrd[x[2]]))
cfrd.diff <- combn(length(veg.means.all$cfrd), 2, function(x) (veg.means.all$cfrd[x[1]] - veg.means.all$cfrd[x[2]]))
pfrd.diff <- combn(length(veg.means.all$pfrd), 2, function(x) (veg.means.all$pfrd[x[1]] - veg.means.all$pfrd[x[2]]))
mch.diff <- combn(length(veg.means.all$mch), 2, function(x) (veg.means.all$mch[x[1]] - veg.means.all$mch[x[2]]))
chcv.diff <- combn(length(veg.means.all$chcv), 2, function(x) (veg.means.all$chcv[x[1]] - veg.means.all$chcv[x[2]]))

model.data <- cbind(as.vector(bdiv.matrix), plai.diff, ulai.diff, mlai.diff, clai.diff, uladcv.diff, 
                    mladcv.diff, cladcv.diff, ufrd.diff, mfrd.diff, cfrd.diff, pfrd.diff, mch.diff, chcv.diff)
colnames(model.data)[1] <- "bdiv"
model.data <- as.data.frame(model.data)

habitat.community.model2 <- lm(bdiv ~ plai.diff, model.data)
habitat.community.model3 <- lm(bdiv ~ ulai.diff, model.data)
habitat.community.model4 <- lm(bdiv ~ mlai.diff, model.data)
habitat.community.model5 <- lm(bdiv ~ clai.diff, model.data)
habitat.community.model6 <- lm(bdiv ~ uladcv.diff, model.data)
habitat.community.model7 <- lm(bdiv ~ mladcv.diff, model.data)
habitat.community.model8 <- lm(bdiv ~ cladcv.diff, model.data)
habitat.community.model9 <- lm(bdiv ~ ufrd.diff, model.data)
habitat.community.model10 <- lm(bdiv ~ mfrd.diff, model.data)
habitat.community.model11 <- lm(bdiv ~ cfrd.diff, model.data)
habitat.community.model12 <- lm(bdiv ~ pfrd.diff, model.data)
habitat.community.model13 <- lm(bdiv ~ mch.diff, model.data)
habitat.community.model14 <- lm(bdiv ~ chcv.diff, model.data)

# Again (test!!!!) but with indiv. BUGG values:

plai.diff <- combn(length(BUGG.veg.metrics$plai), 2, function(x) (BUGG.veg.metrics$plai[x[1]] - BUGG.veg.metrics$plai[x[2]]))
ulai.diff <- combn(length(BUGG.veg.metrics$ulai), 2, function(x) (BUGG.veg.metrics$ulai[x[1]] - BUGG.veg.metrics$ulai[x[2]]))
mlai.diff <- combn(length(BUGG.veg.metrics$mlai), 2, function(x) (BUGG.veg.metrics$mlai[x[1]] - BUGG.veg.metrics$mlai[x[2]]))
clai.diff <- combn(length(BUGG.veg.metrics$clai), 2, function(x) (BUGG.veg.metrics$clai[x[1]] - BUGG.veg.metrics$clai[x[2]]))
uladcv.diff <- combn(length(BUGG.veg.metrics$uladcv), 2, function(x) (BUGG.veg.metrics$uladcv[x[1]] - BUGG.veg.metrics$uladcv[x[2]]))
mladcv.diff <- combn(length(BUGG.veg.metrics$mladcv), 2, function(x) (BUGG.veg.metrics$mladcv[x[1]] - BUGG.veg.metrics$mladcv[x[2]]))
cladcv.diff <- combn(length(BUGG.veg.metrics$cladcv), 2, function(x) (BUGG.veg.metrics$cladcv[x[1]] - BUGG.veg.metrics$cladcv[x[2]]))
ufrd.diff <- combn(length(BUGG.veg.metrics$ufrd), 2, function(x) (BUGG.veg.metrics$ufrd[x[1]] - BUGG.veg.metrics$ufrd[x[2]]))
mfrd.diff <- combn(length(BUGG.veg.metrics$mfrd), 2, function(x) (BUGG.veg.metrics$mfrd[x[1]] - BUGG.veg.metrics$mfrd[x[2]]))
cfrd.diff <- combn(length(BUGG.veg.metrics$cfrd), 2, function(x) (BUGG.veg.metrics$cfrd[x[1]] - BUGG.veg.metrics$cfrd[x[2]]))
pfrd.diff <- combn(length(BUGG.veg.metrics$pfrd), 2, function(x) (BUGG.veg.metrics$pfrd[x[1]] - BUGG.veg.metrics$pfrd[x[2]]))
mch.diff <- combn(length(BUGG.veg.metrics$mch), 2, function(x) (BUGG.veg.metrics$mch[x[1]] - BUGG.veg.metrics$mch[x[2]]))
chcv.diff <- combn(length(BUGG.veg.metrics$chcv), 2, function(x) (BUGG.veg.metrics$chcv[x[1]] - BUGG.veg.metrics$chcv[x[2]]))

model.data <- cbind(as.vector(bdiv.matrix), plai.diff, ulai.diff, mlai.diff, clai.diff, uladcv.diff, 
                    mladcv.diff, cladcv.diff, ufrd.diff, mfrd.diff, cfrd.diff, pfrd.diff, mch.diff, chcv.diff)
colnames(model.data)[1] <- "bdiv"
model.data <- as.data.frame(model.data)

habitat.community.model2 <- lm(bdiv ~ plai.diff, model.data)
habitat.community.model3 <- lm(bdiv ~ ulai.diff, model.data)
habitat.community.model4 <- lm(bdiv ~ mlai.diff, model.data)
habitat.community.model5 <- lm(bdiv ~ clai.diff, model.data)
habitat.community.model6 <- lm(bdiv ~ uladcv.diff, model.data)
habitat.community.model7 <- lm(bdiv ~ mladcv.diff, model.data)
habitat.community.model8 <- lm(bdiv ~ cladcv.diff, model.data)
habitat.community.model9 <- lm(bdiv ~ ufrd.diff, model.data)
habitat.community.model10 <- lm(bdiv ~ mfrd.diff, model.data)
habitat.community.model11 <- lm(bdiv ~ cfrd.diff, model.data)
habitat.community.model12 <- lm(bdiv ~ pfrd.diff, model.data)
habitat.community.model13 <- lm(bdiv ~ mch.diff, model.data)
habitat.community.model14 <- lm(bdiv ~ chcv.diff, model.data)



#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## 2.7.2. Binary trait-habitat associations (Fourth-corner analysis)
library(ade4)

# Run 4th-corner:
four.comb.BUGG1  <- ade4::fourthcorner(veg.means.all, species.x.sites,
                                       all.traits.FD, modeltype = 6, 
                                       p.adjust.method.G = "none", 
                                       p.adjust.method.D = "none", 
                                       nrepet = 3999)
# Plot:
plot(four.comb.BUGG6, alpha = 0.05, stat = "D2")

# With adjusted p-values:
four.comb.BUGG2  <- ade4::fourthcorner(veg.means.all, species.x.sites,
                                       all.traits.FD, modeltype = 6, 
                                       p.adjust.method.G = "fdr", 
                                       p.adjust.method.D = "fdr", 
                                       nrepet = 3999)
plot(four.comb.BUGG2, alpha = 0.05, stat = "D2")

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 



