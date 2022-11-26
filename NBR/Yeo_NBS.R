#set working directory
setwd("~/Shania_graph")
#load packages
packages <- c("NBR","lattice","R.matlab","tseries","lattice","readODS","parallel","nlme","readxl")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}
lapply(packages, library, character.only = TRUE)

#load the demographic file
Data <- read_ods("NBS_demog.ods", col_names = TRUE)
#loading mat file with functional connectivity values
Dat <- readMat("yeo_nbs/YeoNBS.mat")

#converting list into matrix#check the structure of the list and take the one with the matrix
new <- Dat[['Mat']]
A <- Dat$Mat


# Set max-absolute value in order to set a color range centered in zero.
avg_mx <- apply(new, 1:2, mean, na.rm=TRUE)
flim <- max(abs(avg_mx)[is.finite(avg_mx)])
levelplot(avg_mx, main = "Average", ylab = "ROI", xlab = "ROI",
          at = seq(-flim, flim, length.out = 100))

#load data with names of brain regions
brain_labs <- read_excel("Yeo_roi.xlsx", col_names = TRUE)
brain_labs <- colnames(brain_labs[1:17])

#check the no. of permutations required for the current analysis
null_ses_str <- NBR_final$nudist[,2]
obs_ses_str <- NBR_final$fwe$Age[,4]
nperm <- length(null_ses_str)
cumpval <- cumsum(null_ses_str >= obs_ses_str)/(1:nperm)
##Plot p value stability
plot(cumpval, type="l", xlim = c(0,1000), ylim = c(0,1), las = 1, xlab = "Permutation index", ylab = "p-value", main = "Cumulative p-value for Session strength")
abline(h=0.05, col="red", lty=2)
mepval <- 2*sqrt(cumpval*(1-cumpval)/1:nperm)
lines(cumpval+mepval, col = "chartreuse4")
lines(cumpval-mepval, col = "chartreuse4")
#Always check whether the number of observations in net and idata is same
set.seed(18900217)
before <- Sys.time()
 NBR_final <- nbr_lme(
  net = new,
  nnodes = 17,
  idata = Data,
  thrP = 0.05,
  thrT = NULL,
  nperm = 1000,
  nudist = T,
  mod = "~ Age*Group+Gender+FD+MRI_Scanner+Med",
  rdm = "~ 1+Wave|ID",
  cores = detectCores(),
  na.action = na.exclude
)
after <- Sys.time()
show(after-before)
show(NBR_final)

#extracting fwe value of age/group/interaction
B <- show(nbr_result$fwe$Age)
length(nbr_result)
# Plot significant component
edge_mat <- array(0, dim(avg_mx))
edge_mat[nbr_result$components$Age [,2:3]] <- 1
levelplot(edge_mat, col.regions = rev(heat.colors(100)),
          main = "Component", ylab = "ROI", xlab = "ROI")
NBR_fwe_age <- data.frame(nbr_result$fwe$Age)
NBR_fwe_Group <- data.frame(nbr_result$fwe$Group)
NBR_fwe_AgeGroup <- data.frame(nbr_result$fwe$`Age:Group`)
#Components
NBR_component_age <- data.frame(nbr_result$components$Age)
NBR_component_Group <- data.frame(nbr_result$components$Group)
NBR_component_AgeGroup <- data.frame(nbr_result$components$`Age:Group`)



if(!require("corrplot")) install.packages("corrplot"); library(corrplot)
if(!require("circlize")) install.packages("circlize"); library(circlize)
if(!require("scales")) install.packages("scales"); library(scales)
if(!require("plot.matrix")) install.packages("plot.matrix"); library(plot.matrix)
edge_mat <- array(0, dim(avg_mx))
edge_mat[NBR_final$components$`Group` [,2:3]] <- NBR_final$components$`Group`[, "strn"]
##For age
a <- NBR_final$components$'Age'
colnames(edge_mat) <- row.names(edge_mat) <- brain_labs
a<- a[-10,]
edge_mat[a[,2:3]] <- a[,"strn"]
####corrplot upper matrix

corrplot(edge_mat,is.corr = FALSE, col.lim = c(-2,2), tl.col = "black", method = "square", type = 'upper', diag = TRUE, col = colorRampPalette(c("blue", "ivory","red"))(100))
colours()
###chord plot
col.pal = c(visual1 = "purple", visual2 = "red", SAL = "lightpink", FPN3 = "maroon", DMN4 = "lightpink3", limbic1 = "darkolivegreen", limbic2 = "mediumseagreen", FPN2 = "orange", SMN2 = "lightgreen", FPN1 = "grey", DAN2 = "forestgreen", DMN2 = "navy", DMN1 ="blue")
col_fun = colorRamp2(c(-0.2,0,0.2), c("blue", "white", "red"))
chordDiagram(edge_mat, col = col_fun, grid.col = col.pal)
chordDiagram(edge_mat, col = col_fun, grid.col = col.pal, annotationTrack = "grid")


