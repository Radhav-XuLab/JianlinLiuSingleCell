#Fishplot for tumor evolution
library(timescape)
timepoint <- c("Normal", "Mutant Mammary Gland", "Primary Tumor")
clone_id <- c(1, 2, 3)
clonal_prev <- c(1, 1, 1)
clonal_prev <- as.data.frame(cbind(timepoint, clone_id, clonal_prev))
source <- c(1, 2)
target <- c(2, 3)
tree_edges <- as.data.frame(cbind(source, target))
pert_name <- c("Metastasis")
prev_tp <- c("Liver")
perturbations <- as.data.frame(cbind(pert_name, prev_tp))
timescape(clonal_prev = clonal_prev, tree_edges = tree_edges,
          perturbations = perturbations, yaxis_title = "Clonal Frequency",
          height=300)

#Fishplot for #153 tumor evolution using Arhgef11 and Plekha5
library(timescape)
#Clone1 of PT cells with Arhgef11 mutation only
Clone1_P <- c("153PT1", "153PT6", "153PT9", "153PT9", "153PT11", "153PT14", "153PT15", "153PT17")
#Clone2 of PT cells with Arhgef11 and Plekha5 mutation
Clone2_P <- c("153PT2", "153PT3", "153PT4", "153PT7", "153PT8", "153PT12", "153PT13", "153PT18")
#Clone1 of LMT cells with Arhgef11 and Plekha5 mutation
Clone1_L <- c("153LMT1", "153LMT2", "153LMT5", "153LMT6", "153LMT7", "153LMT9", "153LMT10",  "153LMT12", "153LMT13", "153LMT15", "153LMT18", "153LMT19", "153LMT20")
#Clone2 of LMT cells with Arhgef11  mutation only
Clone2_L <- c("153LMT11", "153LMT14")
timepoint <- c("Mutant Mammary Gland", "Primary Tumor", "Primary Tumor", "Metastatic Tumor", "Metastatic Tumor")
clone_id <- c(1, 2, 3, 2, 3)
clonal_prev <- c(1,
                 length(Clone1_P)/(length(Clone1_P) + length(Clone2_P)),
                 length(Clone2_P)/(length(Clone1_P) + length(Clone2_P)),
                 length(Clone2_L)/(length(Clone1_L) + length(Clone2_L)),
                 length(Clone1_L)/(length(Clone1_L) + length(Clone2_L))
)
clonal_prev <- as.data.frame(cbind(timepoint, clone_id, clonal_prev))
source <- c(1, 2, 2)
target <- c(2, 3, 4)
tree_edges <- as.data.frame(cbind(source, target))
pert_name <- c("Metastasis")
prev_tp <- c("Liver")
perturbations <- as.data.frame(cbind(pert_name, prev_tp))
timescape(clonal_prev = clonal_prev, tree_edges = tree_edges,
          perturbations = perturbations, yaxis_title = "Clonal Frequency",
          height=300)