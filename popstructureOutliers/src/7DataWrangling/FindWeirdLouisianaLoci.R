#####################
# Spurious Outliers #
#####################

data <- read.table("data/large_outputs/AllOutlier_WildForAssocEnviAssoc_MergedData_Lotterhos.txt",
                   sep = " ",
                   header = TRUE, 
                   stringsAsFactors = FALSE)

# louisiana.wild - column indices that correspond to the wild LA populations
louisiana.wild <- which(wild$Pop.ID == "CL" | wild$Pop.ID == "SL")

# louisiana.wild_spurious.outliers - all weird loci
# This is where all other individuals are 0's where LA wild individuals are 
#   1 or 2. I couldn't get this to work by using c(1:12) and I am not sure why.
#   I didn't have the patience to figure it out so I found a workaround...
louisiana.wild_spurious.outliers <- which(rowSums(wild$G[, -louisiana.wild]) == 0 & 
                                          wild$G[, louisiana.wild[1]] != 0 & 
                                          wild$G[, louisiana.wild[2]] != 0 & 
                                          wild$G[, louisiana.wild[3]] != 0 & 
                                          wild$G[, louisiana.wild[4]] != 0 & 
                                          wild$G[, louisiana.wild[5]] != 0 & 
                                          wild$G[, louisiana.wild[6]] != 0 & 
                                          wild$G[, louisiana.wild[7]] != 0 & 
                                          wild$G[, louisiana.wild[8]] != 0 & 
                                          wild$G[, louisiana.wild[9]] != 0 & 
                                          wild$G[, louisiana.wild[10]] != 0 & 
                                          wild$G[, louisiana.wild[11]] != 0 & 
                                          wild$G[, louisiana.wild[12]] != 0)
