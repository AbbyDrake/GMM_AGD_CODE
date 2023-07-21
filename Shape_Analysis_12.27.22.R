#load libraries
library(geometry)
library(geomorph)
library(RRPP)

setwd("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/")
dir()

#Load Datafile
Data_File<-
  read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/All_Data_12.27.22.txt", header=T)

#get coordinates from file
Symmdata<-as.matrix(Data_File[1:length(Data_File$Id),14:dim(Data_File)[2]]) 
head(Symmdata)

# Log centroid size
Landmark.lnCS<-Data_File$Log_Centroid_Size  

#PCA
Landmark.PCA<-prcomp(Symmdata,center=T)
summary(Landmark.PCA)
# Importance of components:
#   PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10    PC11    PC12    PC13
# Standard deviation     0.09713 0.04724 0.03290 0.02854 0.02254 0.01712 0.01603 0.01512 0.01385 0.01309 0.01262 0.01144 0.01121
# Proportion of Variance 0.54700 0.12938 0.06277 0.04721 0.02945 0.01699 0.01490 0.01326 0.01112 0.00994 0.00924 0.00759 0.00729
# Cumulative Proportion  0.54700 0.67638 0.73916 0.78637 0.81582 0.83281 0.84771 0.86097 0.87209 0.88203 0.89127 0.89885 0.90614

##Correlation of PCs with log centroid size
cor.test(Landmark.PCA$x[,1],Landmark.lnCS,method="pearson") #r = 0.6977262          p-value < 2.2e-16 
cor.test(Landmark.PCA$x[,2],Landmark.lnCS,method="pearson") #r = 0.3134604    p-value < 2.2e-16 

##Convex Hull Volume of PCs 1-3; Procrustes Variance; Procrustes Distances##

#DOMESTIC CATS
Cat_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_4"] == "DOM_CAT"), ]
Cat_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Cat_Subset_Landmark_mat[1:length(Cat_Subset_Landmark_mat$Id),14:dim(Cat_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
DOMCAT_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
DOMCAT_vol$vol
# 0.001551999

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.01003003

#Find maximum Procrustes distance  
Proc_dist_Cats<-as.matrix(dist(Shape_data))
Proc_dist_Cats
max(Proc_dist_Cats) # 0.3248903; Persian (Schwager_Hildegard_Blue_Ocean_Henkel): Siamese 451299
#write.csv(Proc_dist_Cats, "Proc_dist_Shape_Cats.csv") 

#DOMESTIC DOGS
Dog_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_4"] == "DOM_DOG"), ]
Dog_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Dog_Subset_Landmark_mat[1:length(Dog_Subset_Landmark_mat$Id),14:dim(Dog_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
DOMDOG_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
DOMDOG_vol$vol
# 0.00713864

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.01307608

#Find maximum Procrustes distance 
Proc_dist_Dogs<-as.matrix(dist(Shape_data))
Proc_dist_Dogs
max(Proc_dist_Dogs) # 0.4512234; Pekingese (409): Borzoi (318)
#write.csv(Proc_dist_Dogs, "Proc_dist_Shape_Dogs.csv") 

#CARNIVORA
Carnivora_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_4"] == "Carnivora"), ]
Carnivora_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Carnivora_Subset_Landmark_mat[1:length(Carnivora_Subset_Landmark_mat$Id),14:dim(Carnivora_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
Carnivora_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
Carnivora_vol$vol
# 0.007566942

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.0153711

#Find maximum Procrustes distance
Proc_dist_Carnivora<-as.matrix(dist(Shape_data))
Proc_dist_Carnivora
max(Proc_dist_Carnivora) # 0.3950533; Odobenus rosmarus: Eupleres goudoti
#write.csv(Proc_dist_Carnivora, "Proc_dist_Shape_Carnivora.csv") 

#CANIDAE
Canidae_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_3"] == "Canidae"), ]
Canidae_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Canidae_Subset_Landmark_mat[1:length(Canidae_Subset_Landmark_mat$Id),14:dim(Canidae_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
Canidae_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
Canidae_vol$vol
# 0.001166492

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.003254366

#Find maximum Procrustes distance  
Proc_dist_Canidae<-as.matrix(dist(Shape_data))
Proc_dist_Canidae
max(Proc_dist_Canidae) # 0.2145899; Speothos venaticus: Canis lupus (17289)
#write.csv(Proc_dist_Canidae, "Proc_dist_Shape_Canidae.csv") 

#FELIDAE
Felidae_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_3"] == "Felidae"), ]
Felidae_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Felidae_Subset_Landmark_mat[1:length(Felidae_Subset_Landmark_mat$Id),14:dim(Felidae_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
Felidae_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
Felidae_vol$vol
# 0.001764006

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.007011887

#Find maximum Procrustes distance  
Proc_dist_Felidae<-as.matrix(dist(Shape_data))
Proc_dist_Felidae
max(Proc_dist_Felidae) # 0.2334182; Herpailurus yagouaroundi (ZMB21295): Acinonyx	jubatus (27897)
#write.csv(Proc_dist_Felidae, "Proc_dist_Shape_Felidae.csv") 

#WILDCATS
Wildcat_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_2"] == "Wildcat"), ]
Wildcat_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Wildcat_Subset_Landmark_mat[1:length(Wildcat_Subset_Landmark_mat$Id),14:dim(Wildcat_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
Wildcat_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
Wildcat_vol$vol
# 0.0001932972

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.002450796

#Find maximum Procrustes distance  
Proc_dist_Wildcats<-as.matrix(dist(Shape_data))
Proc_dist_Wildcats
max(Proc_dist_Wildcats) # 0.131812
#write.csv(Proc_dist_Wildcats, "Proc_dist_Shape_Wildcats.csv") 

#WOLVES
Wolf_Subset_Landmark_mat <- Data_File[which(Data_File[,"Group_2"] == "Wolf"), ]
Wolf_Subset_Landmark_mat

#get coordinates from file
Shape_data<-as.matrix(Wolf_Subset_Landmark_mat[1:length(Wolf_Subset_Landmark_mat$Id),14:dim(Wolf_Subset_Landmark_mat)[2]]) 
head(Shape_data)

#PCA
Shape.PCA<-prcomp(Shape_data,center=T)

#Volume of Convex Hull of PCs 1-3
Wolf_vol<-convhulln(
  Shape.PCA$x[,1:3],
  options = NULL,
  output.options = "FA",
  return.non.triangulated.facets = FALSE
)
Wolf_vol$vol
# 0.0002058303

#Procrustes Variance
Proc.Var<-morphol.disparity(Shape_data~1)
Proc.Var
# 0.001685166

#Find maximum Procrustes distance within Wolves 
Proc_dist_Wolves<-as.matrix(dist(Shape_data))
Proc_dist_Wolves
max(Proc_dist_Wolves) # 0.1213065
#write.csv(Proc_dist_Wolves, "Proc_dist_Shape_Wolves.csv") 

######### Creating warped skull shapes for PCs 1 and 2
setwd("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/")
dir()

#Load Datafile
Data_File<-
  read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/All_Data_12.27.22.txt", header=T)

#get coordinates from file
Symmdata<-as.matrix(Data_File[1:length(Data_File$Id),14:dim(Data_File)[2]]) 
head(Symmdata)

#Load the shape file and coordinate data
InitialPLY<-read.ply("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/Average_Skull.ply", ShowSpecimen = FALSE, addNormals = FALSE)
InitialPLY.coor<-read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/47_Landmarks_Average_Skull.txt", header=FALSE)
InitialPLY.coor
InitialPLY.coor<-as.matrix(InitialPLY.coor)

InitialPLY.coor.array<-matrix(unlist(InitialPLY.coor), nrow = 3, ncol = 47)
InitialPLY.coor.array<-t(InitialPLY.coor.array)
InitialPLY.coor.array

#Calculate the mean shape using mshape
Symmdata_Mean<- apply(Symmdata,2,mean)
Symmdata_Mean
Symmdata_Mean.coor<-as.matrix(Symmdata_Mean)
class(Symmdata_Mean)

Symmdata_Mean.coor.array<-matrix(unlist(Symmdata_Mean.coor), nrow=3, ncol=47)
Symmdata_Mean.coor.array<-t(Symmdata_Mean.coor.array)
Symmdata_Mean.coor.array
write.table(Symmdata_Mean.coor.array, file = "Symmdata_Mean.coor.array.txt")

#Warp intial ply to mean shape
AvgShape<-warpRefMesh(InitialPLY, InitialPLY.coor.array, Symmdata_Mean.coor.array)
writePLY("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/avgskullshape.ply")

#PC1_Min
PC1minPLY.coor<-read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC1neg0.15.txt", header=FALSE)
PC1minPLY.coor
PC1minPLY.coor<-as.matrix(PC1minPLY.coor)

PC1min<-warpRefMesh(AvgShape, Symmdata_Mean.coor.array, PC1minPLY.coor, color = "azure3")

writePLY("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC1_min_shape.ply")

#PC1_Max
PC1maxPLY.coor<-read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC1pos0.27.txt", header=FALSE)
PC1maxPLY.coor
PC1maxPLY.coor<-as.matrix(PC1maxPLY.coor)

PC1max<-warpRefMesh(AvgShape, Symmdata_Mean.coor.array, PC1maxPLY.coor, color = "azure3")

writePLY("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC1_max_shape.ply")

#PC2_Min
PC2minPLY.coor<-read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC2neg0.15.txt", header=FALSE)
PC2minPLY.coor
PC2minPLY.coor<-as.matrix(PC2minPLY.coor)

PC2Min<-warpRefMesh(AvgShape, Symmdata_Mean.coor.array, PC2minPLY.coor, color = "azure3")

writePLY("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC2_min_shape.ply")

#PC2_Max
PC2maxPLY.coor<-read.table("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC2pos0.17.txt", header=FALSE)
PC2maxPLY.coor
PC2maxPLY.coor<-as.matrix(PC2maxPLY.coor)

PC2max<-warpRefMesh(AvgShape, Symmdata_Mean.coor.array, PC2maxPLY.coor, color = "azure3")

writePLY("C:/Users/abbyg/Desktop/Research/Cat_Skulls/December_Analysis/Shape_Models/PC2_max_shape.ply")
