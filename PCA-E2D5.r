library(flowWorkspace)
library(ggplot2)
library(reshape2)
library(rgl)

# Original samples for experiment 2 donor 5 (E2D5)
fcs_file <- "E2D5PI.fcs"
flow_frame <- read.FCS(fcs_file)
markers <- colnames(flow_frame)
markers <- markers[!grepl("Tet", markers)]
markers <- markers[!grepl("blank", markers)]
markers <- markers[!grepl("DNA", markers)]
markers <- markers[!grepl("CD3", markers)]
markers <- markers[!grepl("LiveDead", markers)]
markers <- markers[!grepl("MuCD45", markers)]
markers <- markers[!grepl("Ba138", markers)]
markers <- setdiff(markers, c("Time", "Cell_length"))

# Filters out markers irrelevant to the analysis of Newell et al. (2012)
x_original <- exprs(flow_frame[, markers])
colnames(x_original) <- sapply(strsplit(markers, "\\("), head, n = 1)

# The CD8+ T-cell subsets were provided as separate FCS files rather than via a
# workspace with the requisite gates.
fcs_subsets <- c("E2D5Naive.fcs", "E2D5SLEC.fcs", "E2D5CM.fcs", "E2D5EM.fcs")
fs_subset <- read.flowSet(files = fcs_subsets)

# Filters out markers irrelevant to the analysis of Newell et al. (2012)
fs_subset <- fs_subset[, markers]
colnames(fs_subset) <- sapply(strsplit(markers, "\\("), head, n = 1)

df_subset <- rbind(
    data.frame(exprs(fs_subset[[1]]), Population = "Naive"),
    data.frame(exprs(fs_subset[[2]]), Population = "SLEC"),
    data.frame(exprs(fs_subset[[3]]), Population = "CM"),
    data.frame(exprs(fs_subset[[4]]), Population = "EM")
)
# ggplot(df_subset, aes(x = CCR7, y = CD45RA, color = Population)) + geom_point()

# Below, we replicate the so-called  3D-PCA analysis from Newell et al. (2012).
pca_out <- prcomp(x_original, scale = TRUE)
pca_rotation <- pca_out$rotation

# Extracts PCA variance
pca_prop <- summary(pca_out)$importance[2, ]
pca_cumprop <- summary(pca_out)$importance[3, ]

# Replicates Figure 3B
barplot(pca_prop[1:10], ylim = c(0, 0.8))
lines(pca_cumprop[1:10], type = "b")

# Replicates Figure 3C-3E
pca_3D_out <- pca_rotation[, 1:3]
pca_3D_out <- melt(pca_3D_out)
colnames(pca_3D_out) <- c("Marker", "PC", "Loading")

p <- ggplot(pca_3D_out, aes(x = Marker, weight = Loading)) + geom_bar()
p <- p + ylab("Principal Component Loading") + facet_wrap(~ PC, ncol = 1)
p + theme(axis.text.x = element_text(angle = 90))

# Scales subsetted cells using those from the original data
x_subset <- scale(x = as.matrix(subset(df_subset, select = -Population)),
                  center = pca_out$center, scale = pca_out$scale)

# Rotates the cells using the PCA obtained from the original data
x_subset <- x_subset %*% pca_rotation

# Plots the first 3 PCs of the original data using 'rgl'
# Then, overlays the subsetted cells
# The colors mostly follow Newell et al. (2012).
# I changed yellow to orange, which is easier to see.
# Also, the original observations are small and in gray rather than black;
# otherwise, they overwhelm the plot.
plot_colors <- df_subset$Population
levels(plot_colors) <- list(green = "Naive", red = "SLEC", orange = "CM", blue = "EM")


pca_original <- pca_out$x[, 1:3]
pca_subsets <- x_subset[, 1:3]


save(pca_original, pca_subsets, plot_colors, file = "cache/PCA-CyTOF.RData")


plot3d(pca_original, alpha = 0.15)
plot3d(pca_subsets, col = plot_colors, add = TRUE)
