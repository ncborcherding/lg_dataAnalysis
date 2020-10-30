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
        d <- TCR[which(TCR[,"Var1"] == a), "vgene"]
        e <- TCR[which(TCR[,"Var1"] == b), "vgene"]
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

# Chain specifies which chain to pull the information from
globalConvergence <- function(tmp, chain = "TCRB") {
    
    TCR <- getTCR(tmp, chain)
    TCR_dist <- as.matrix(stringdistmatrix(TCR$Var1, method = "lv"))
    TCR_dist[TCR_dist >= 2] <- NA
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

getTCR <- function(tmp, chain) {
    cdr3 <- chainCheck(chain)[[1]]
    vgene <- chainCheck(chain)[[2]]
    TCR <- as.data.frame(na.omit(table(tmp[,cdr3])))
    TCR$length <- nchar(as.character(TCR$Var1))
    vgene <- unique(tmp[,c(cdr3, vgene)])
    colnames(vgene) <- c("cdr3", "TCR")
    vgene$TCR <- str_split(vgene$TCR, "[.]", simplify = T)[,1]
    vgene <- unique(vgene)
    TCR <- merge(TCR, vgene, by.x = "Var1", by.y = "cdr3")
    colnames(TCR)[ncol(TCR)] <- "vgene"
    return(TCR)
}

localConvergence <- function(tmp, chain = "TCRB", motif.length = 3) {
    load("./data/AA_combinations.rda")
    TCR <- getTCR(tmp,chain)
    positions <- NULL
    vgenes <- unique(TCR[,4])
    for (j in seq_along(vgenes)) {
        subset <- TCR[TCR$vgene %in% vgenes[j],]
        out <- matrix(nrow = nrow(subset), ncol =length(AA_combinations), 0)
        colnames(out) <- AA_combinations
        rownames(out) <- subset$Var1
    
        for (i in seq_len(nrow(subset))) {
            tmp <- as.character(subset$Var1[i])
            for (j in seq_len(nchar(tmp))) {
                string <- substr(tmp, j, j+(motif.length-1))
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
    con <- controls[sample(nrow(controls), nrow(TCR)),]
    cdr3 <- chainCheck(chain)[[1]]
    con2 <- pbmc[pbmc[,cdr3] %in% con$Var1,]
    gc_con <- globalConvergence(con2, chain = "TCRB")
    gc_con <- checkVgenes(gc_con, con)
    y <- nrow(gc_con)
}

# for bootstrapping the local convergence
sampleControl_lc <- function(i) {
    con <- controls[sample(nrow(controls), nrow(TCR)),]
    cdr3 <- chainCheck(chain)[[1]]
    con2 <- pbmc[pbmc[,cdr3] %in% con$Var1,]
    lc_con <- localConvergence(con2, chain = "TCRB")
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