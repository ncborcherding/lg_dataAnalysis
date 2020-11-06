#Chains for MAIT and INKT cells
scoreMAIT <- function(membership, species = NULL) {
    comp <- list(mouse = list(v = "TRAV1", j = "TRAJ33", length = 12), 
                human = list(v = "TRAV1-2", j = c("TRAJ33", "TRAJ20", "TRAJ12"), length = 12))
    score <- data.frame("cdr3" = membership$cdr3, score = 0)
    for (i in seq_len(nrow(membership))) {
        v <- membership$v[i]
        j <- membership$j[i]
        length <- nchar(membership$cdr3)[i]
        if(comp[[species]]$v == v &  j %in% comp[[species]]$j & length %in% comp[[species]]$length){
            score$score[i] <- 1
        }
        else {
            next()
        }
    }
    return(score)
}

scoreINKT <- function(membership, species = NULL) {
    comp <- list(mouse = list(v = "TRAV11", j = "TRAJ18", length = 15), 
                 human = list(v = "TRAV10", j = c("TRAJ18", "TRBV25"), length = c(14,15,16)))
    score <- data.frame("cdr3" = membership$cdr3, score = 0)
    for (i in seq_len(nrow(membership))) {
        v <- membership$v[i]
        j <- membership$j[i]
        length <- nchar(membership$cdr3)[i]
        if(comp[[species]]$v == v &  j %in% comp[[species]]$j & length %in% comp[[species]]$length){
            score$score[i] <- 1
        }
        else {
            next()
        }
    }
    return(score)
}

# Getting positions from matrix of global and local convergence.

# out is the matrix
# integer.start is the low end of values to scan
# integer.stop is the high end of values to scan
# order is to change column order in a mirrored matrix to eliminate redundancy
getPostitions <- function(out, integer.start, integer.stop, order = FALSE) {
    form <- paste0("[", integer.start, "-", integer.stop, "]")
    positions <- which(matrix(grepl(form, out), ncol=ncol(out)), arr.ind = TRUE)
    if (order == TRUE) {
        for (y in seq_len(nrow(positions))) {
            a <- positions[y,1]
            b <- positions[y,2]
            if (a < b) {
                tmp <- a
                positions[y,1] <- b
                positions[y,2] <- tmp
            }
        }
        positions <- unique(positions)
    }
    return(positions)
}

# Ensure the convergence analysis is only examining similarity between CDR3 sequences derived from
# the same v-gene

# positions is the final information of AA relations
# TCR is the information gathered for the initial analysis
checkVgenes <- function(positions, TCR) {
    count <- NULL
    for (i in seq_len(nrow(positions))) {
        a <- positions[i,"To"]
        b <- positions[i,"From"]
        d <- TCR[which(TCR[,"Var1"] == a), "v"]
        e <- TCR[which(TCR[,"Var1"] == b), "v"]
        if (d[1] == e[1]) {
            next()
        } else {
            count <- c(count,i)
        } 
    }
    if (length(count) > 0) {
        positions <- positions[-count,]
    }
    return(positions)
}

# Chain is either TCRA or TCRB
chainCheck <- function(chain)  {
    if (chain == "TCRB") {
        cdr3 <- "cdr3_aa2"
        vgenes <- "TCR2"
    } else if (chain == "TCRA") {
        cdr3 <- "cdr3_aa1"
        vgenes <- "TCR1"
    }
    return(list(cdr3, vgenes))
}

#T his function will calculate the edit distance between all unique aa-based
# CDR3 sequences, saving those that are 0-1 away from one another and then
# eliminating the interactions that do not share the same Vgene. Returning
# a relational data frame without redundancy.

# chain specifies which chain to pull the information from
# edit.distance is the cut.off for levenshtein distance to be taken 
# into account for clustering
globalConvergence <- function(tmp, edit.distance = NULL) {
    
    TCR <- getTCR(tmp)
    TCR_dist <- as.matrix(stringdistmatrix(TCR$Var1, method = "lv"))
    TCR_dist[TCR_dist >= 1+edit.distance] <- NA
    for (i in seq_len(ncol(TCR_dist))) {
        TCR_dist[i,i] <- NA
    }
    positions_LV <- getPostitions(TCR_dist, 0,1, order = T)
    if (nrow(positions_LV) == 0) {
        return(positions_LV)
    }
    positions_LV <- data.frame(positions_LV, "Type" = "LV")
    colnames(positions_LV)[1:2] <- c("To", "From")
    
    positions_LV$To <- as.character(TCR$Var1)[positions_LV$To]
    positions_LV$From <- as.character(TCR$Var1)[positions_LV$From]
    return(positions_LV)
}

getTCR <- function(tmp) {
    chains <- list(one = c("TCR1","cdr3_aa1"), two = c("TCR2","cdr3_aa2"))
    TCR <- NULL
    for (i in 1:2) {
        sub <- as.data.frame(na.omit(table(tmp[,chains[[i]]][2])))
        gene <- unique(tmp[,c(chains[[i]][2], chains[[i]][1])])
        colnames(gene) <- c("cdr3_aa", "TCR")

        gene$v <- str_split(gene$TCR, "[.]", simplify = T)[,1]
        gene$j <- str_split(gene$TCR, "[.]", simplify = T)[,2]
        sub <- merge(sub, gene, by.x = "Var1", by.y = "cdr3_aa")
        TCR <- rbind.data.frame(TCR, sub)
    }
    return(TCR)
}

localConvergence <- function(tmp, motif.length = 3) {
    load("./data/AA_combinations.rda")
    TCR <- getTCR(tmp)
    positions <- NULL
    vgenes <- unique(TCR[,4])
    for (j in seq_along(vgenes)) {
        subset <- TCR[TCR$v %in% vgenes[j],]
        out <- matrix(nrow = nrow(subset), ncol =length(AA_combinations), 0)
        colnames(out) <- AA_combinations
        rownames(out) <- subset$Var1
    
        for (i in seq_len(nrow(subset))) {
            tmp <- as.character(subset$Var1[i])
            index <- seq_len(nchar(tmp))
            index <- index[-c(1,2,3, nchar(tmp)-1, nchar(tmp)-2, nchar(tmp)-3)]
            for (j in index) {
                string <- substr(tmp, j+3, j+(motif.length-1)+3)
                if (string %in% AA_combinations & nchar(string) == motif.length) {
                    position <- which(colnames(out) == string)
                    out[i,position] <- out[i,position] + 1
                } else {
                    next()
                }
            }
        }

    out[out == 0] <- NA
    positions_motif <- getPostitions(out, 0,9)
    motifs <- unique(positions_motif[,2])
    df.edge <- NULL
    for (k in seq_along(motifs)) {
        correspond <- positions_motif[positions_motif[,2] == motifs[k],]
        if (length(correspond) == 2) {
            next()
        }
        out <- t(combn(c(correspond[,1], correspond[,1]), 2))
        out <- unique(out)
        out <- out[out[,1] != out[,2],]
        To <- as.character(subset$Var1)[out[,1]]
        From <- as.character(subset$Var1)[out[,2]]
        out <- data.frame(To,From, spec.Motif = AA_combinations[motifs[k]])
        df.edge <- rbind(df.edge, out)
        }
    positions <- rbind(positions, df.edge)
    }
    positions <- data.frame(positions, "Type" = "Motif")
    positions<- unique(positions)
    return(positions)
}

# for bootstrapping the global convergence
sampleControl_gc <- function(i) {
    con <- controls[sample(nrow(controls), nrow(TCR)/2),]
    con2 <- rbind.data.frame(pbmc[pbmc[,c("cdr3_aa1")] %in% con$Var1,],
                            pbmc[pbmc[,c("cdr3_aa2")] %in% con$Var1,])
    gc_con <- globalConvergence(con2, edit.distance = edit.distance)
    gc_con <- checkVgenes(gc_con, controls)
    y <- nrow(gc_con)
}

# for bootstrapping the local convergence
sampleControl_lc <- function(i) {
    con <- controls[sample(nrow(controls), nrow(TCR)/2),]
    con2 <- rbind.data.frame(pbmc[pbmc[,c("cdr3_aa1")] %in% con$Var1,],
                             pbmc[pbmc[,c("cdr3_aa2")] %in% con$Var1,])
    lc_con <- localConvergence(con2)
    y <- unlist(table(lc_con$spec.Motif))
    return(y)
}

# reorganization of the combined list from scRepertoire
parsingContigList <- function(combined, group = NULL) {
    if (is.null(group)) {
        group <- "sample"
    }
    tmp <- bind_rows(combined)
    group.list <- split(tmp, f = tmp[,group])
    return(group.list)
}



#Function to calculate local and global convergence over contig information
#processed by scRepertoire. Will call T cell clusters by local and global
#convergence

#combined is the object from combineTCR()
#chain is either TCRB or TCRA
#group is how to calculate convergence, if null, will calculate across each sample/list element
#motif.length The window to evaluate amino acide motifs
#num.cores how many cores to use during the bootstrapping process
#boot.straps the number of bootstraps to calculate
#edit.distance is the cut.off for levenshtein distance to be taken 
# into account for clustering
#fc.motif the fold-change cut off for the motif/local convergence analysis
#p.value.motif the one-side p-value cut-off for the motif/local convergence analysis
calculateConvergence <- function(combined, group = NULL, 
                                 motif.length = 3, num.cores =2, boot.straps = 1000, 
                                 edit.distance = 1, p.value.motif = 0.001, fc.motif = 5,
                                 score.cluster = NULL) {
    tmp.list <- parsingContigList(combined, group = group)
    load("./data/pbmcControls.rda")
    controls <- getTCR(pbmc)
    new.list <- list()
    for (x in seq_along(tmp.list)) {
        le <- x
        
        tmp <- tmp.list[[x]]
        TCR <- getTCR(tmp)
        message(paste("Calculating Global Convergence in:", names(tmp.list)[x]))
        positions_LV <- globalConvergence(tmp, edit.distance = edit.distance)
        TCR_dist <- as.matrix(stringdistmatrix(TCR$Var1, method = "lv"))
        positions_LV <- checkVgenes(positions_LV, TCR)
        message(paste("Calculating Local Convergence in:", names(tmp.list)[x]))
        positions_motif <- localConvergence(tmp, motif.length=motif.length)
        
        message(paste("Bootstrapping Random Global Convergence in:", names(tmp.list)[x]))
        bootstrap_gc <- pbmclapply(1:boot.straps, sampleControl_gc, mc.cores = num.cores)
        
        message(paste("Bootstrapping Random Local Convergence in:", names(tmp.list)[x]))
        bootstrap_lc <- pbmclapply(1:boot.straps, sampleControl_lc, mc.cores = num.cores)
        
        message(paste("Filtering Bootstrapped Local Convergence in:", names(tmp.list)[x]))
        #Make single data frame of motifs in bootstrap values
        x <- lapply(bootstrap_lc, function(x){data.frame(x)})
        df <- suppressWarnings(Reduce(function(x, y) merge(x=x, y=y, by="Var1", all.x=T, all.y=T), x))
        df[is.na(df)] <- 0
        
        ##for each motif, calculate p-value and fold-change
        table <- as.data.frame(table(positions_motif$spec.Motif))
        
        for (i in seq_len(nrow(table))) {
            vector <- df[table[,"Var1"][i], 2:ncol(df)]
            number <- table[i,"Freq"]
            table$p[i] <- length(vector[vector >= number])/length(vector)
            table$fc[i] <- number/median(as.numeric(vector))
        }
        #Filter out motifs based on p-value and fold-change
        table <- table[table$p <= p.value.motif & table$fc >= fc.motif,]
        positions_motif <- positions_motif[positions_motif$spec.Motif %in% table$Var1,]
        edges <- rbind.data.frame(positions_LV, positions_motif[,c(1,2,4)])
        vertices <- unique(as.data.frame(TCR$Var1))
        
        g = graph_from_data_frame(edges, directed = FALSE, vertices = vertices)
        
        clusters <- components(g, mode = "strong")
        membership <- data.frame(cdr3 = names(clusters$membership), TCRcluster = clusters$membership)
        membership$TCRcluster <- paste0(names(tmp.list)[le], "_", membership$TCRcluster)
        membership <- merge(membership, TCR, by.x = "cdr3", by.y = "Var1")
        
        motifs <- positions_motif %>% 
            group_by(To) %>% 
            summarise(.groups = "keep", filtered.motifs = paste(unique(spec.Motif), collapse = ','))
        names.motif <- unique(motifs$filtered.motifs) 
        pLC <- NULL
        for (i in seq_along(names.motif)) {
            mot <- unlist(str_split(names.motif[i], ","))
            sum <- sum(table[table$Var1 %in% mot,]$p)
            pLC <- c(pLC, sum)
        }
        names(pLC) <- names.motif
        pLC <- as.data.frame(pLC)
        membership <- merge(membership, motifs, by.x = "cdr3", by.y = "To")
        membership <- merge(membership, pLC, by.x = "filtered.motifs", by.y = "row.names")
        iNKT <- scoreINKT(membership, species)
        MAIT <- scoreMAIT(memership, species)
        new.list[[le]] <- list(contigs = tmp, clusters = membership,  edit.distances = TCR_dist, 
                                iNKT = iNKT, MAIT = MAIT)
    }
    names(new.list) <- names(tmp.list)
    if (score.cluster == TRUE) {
        message(paste("Calculating Cluster Probabilities"))
        new.list <- scoreCluster(new.list)
    }
    return(new.list)
}

vgeneDiversity <- function(sub) {
    genes <- seq_len(length(unique(sub[,"vgene"])))
    names(genes) <- unique(sub[,"vgene"])
    vgenes <- genes[sub[,"vgene"]] 
    vgene.diversity <- diversity(vgenes, index = "simpson")
    return(vgene.diversity)
}
vgeneDiversity.perm <- function(i) {
    cont.perm <- tmp[sample(nrow(tmp), nrow(sub)),]
    diversity <- vgeneDiversity(cont.perm)
    return(diversity)
}

length.perm <- function(i) {
    cont.perm <- tmp[sample(nrow(tmp), nrow(sub)),]
    diversity <- diversity(cont.perm[,"length"], index = "simpson")
    return(diversity)
}

freq.perm <- function(i) {
    cont.perm <- tmp[sample(nrow(tmp), nrow(sub)),]
    max.freq <- length(unique(names(table(cont.perm$cdr3_aa2))))/nrow(sub)
    return(max.freq)
}

probabilityFun <- function(raw, perm) {
    prob <- length(perm[perm >= raw])/ length(perm)
    return(prob)
    
}


scoreCluster <- function(convergence) {
    out <- NULL
    for (i in seq_along(convergence)) {
        tmp <- convergence[[i]][[1]]
        cluster <- unique(tmp$TCRcluster)
        for (j in seq_along(cluster)) {
            sub <- tmp[tmp$TCRcluster == cluster[j],]
            
            #Probability of randomly length diversity
            #somewhat biased because of the selection of vgenes
            vgene.diversity <- vgeneDiversity(sub)
            vgene.diversity.perm <- unlist(lapply(1:1000, vgeneDiversity.perm))
            pVgene <- probabilityFun(vgene.diversity, vgene.diversity.perm)
            
            #Probability of randomly length diversity
            length.diversity <- diversity(sub[,"length"], index = "simpson")
            length.diversity.perm <- unlist(lapply(1:1000, vgeneDiversity.perm))
            pLength <- probabilityFun(length.diversity, length.diversity.perm)
            
            #Probability of randomly sampled frequences
            freq.size <-  length(unique(names(table(sub[,cdr3]))))/nrow(sub)
            freq.size.perm <- unlist(lapply(1:1000, freq.perm))
            pFreq <- probabilityFun(freq.size, freq.size.perm)
            
            #mean motif p-values for cluster
            pLC <- mean(sub[,ncol(sub)])
            
            #probability of global similarity
            pGC <- probabilityFun(nrow(positions_LV), unlist(bootstrap_gc))
            
            if(any(list(pVgene, pLength, pFreq, pLC, pGC) == 0)) {
                vars <- which(list(pVgene, pLength, pFreq, pLC, pGC) == 0)
                for (k in vars) {
                    p <- list(pVgene, pLength, pFreq, pLC, pGC)[[k]] 
                    varia <- c("pVgene", "pLength", "pFreq", "pLC", "pGC")[k] 
                    p <- 0.001
                    assign(varia, p)
                }
            }
            prob <- (pLC * pGC * pVgene *pLength * pFreq) / 
                ((pLC * pGC * pVgene *pLength * pFreq) + 
                     ((1-pLC) * (1-pGC) * (1-pVgene) * (1-pLength) * (1-pFreq)))
        }
        sub$Cluster.Enrichment <- prob
        out <- rbind(out, sub)    
    }
    convergence[[i]][[1]] <- out
    
}
