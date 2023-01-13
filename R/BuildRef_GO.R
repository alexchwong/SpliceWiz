.process_ontology <- function(
        reference_path, genome_type, verbose = TRUE
) {
    hasPackage <- 
        .check_package_installed("DBI", "1.0.0", "silent") && 
        .check_package_installed("GO.db", "3.12.0", "silent")

    if(!hasPackage & genome_type %in% c("hg38", "hg19", "mm10", "mm9")) {
        .log(paste("Packages DBI and GO.db are required",
            "for gene ontology annotations, skipping..."
        ), "message")
        return()
    }
    if(genome_type %in% c("hg38", "hg19")) {
        species <- "Homo sapiens"
    } else if(genome_type %in% c("mm10", "mm9")) {
        species <- "Mus musculus"
    } else {
        if(verbose) 
            .log("Gene ontology not prepared for this reference", "message")
        return()
    }
    ontDT <- .get_geneGO(species, verbose)
    fst::write.fst(ontDT, file.path(reference_path, "fst", "Ontology.fst"))
}

.ora_internal <- function(
    collate_path,
    genes, universe,
    ontologyType = "BP",
    pAdjustMethod="BH",
    ...
) {
    .check_package_installed("fgsea", "1.0.0", "silent")

    reference_path <- file.path(collate_path, "Reference")
    ontFile <- file.path(reference_path, "fst/Ontology.fst")
    if(!file.exists(ontFile)) .log(paste(
        "No gene ontology resource found in", reference_path))

    ont <- as.data.table(fst::read.fst(ontFile))
    ontUse <- ont[get("ontology") == ontologyType]
    if(nrow(ontUse) == 0) .log(paste(
        ontologyType, "not found as a gene ontology category"
    ))

    pathways <- split(ontUse$ensembl_id, ontUse$go_id)
    
    foraRes <- as.data.table(fgsea::fora(
        pathways, 
        genes=genes, 
        universe=universe,
        ...
    ))

    foraHeader <- data.table(
        go_id = foraRes$pathway,
        go_term = ontUse$go_term[match(
            foraRes$pathway,
            ontUse$go_id
        )]
    )
    
    final <- cbind(foraHeader, foraRes[, -1])
    
    # Fold enrichment
    final$expected <- ceiling(length(unique(genes)) * final$size / 
        length(unique(universe)))
    final$foldEnrichment <- final$overlap / (final$expected + 0.001)
    final$foldEnrichment <- round(final$foldEnrichment, 2)
    final$padj <- p.adjust(final$pval, pAdjustMethod)
    return(final)
}

.subset_EventNames_by_GO <- function(
    EventNames,
    go_id,
    collate_path
) {
    reference_path <- file.path(collate_path, "Reference")
    ontFile <- file.path(reference_path, "fst/Ontology.fst")
    if(!file.exists(ontFile)) .log(paste(
        "No gene ontology resource found in", reference_path))
    
    mapperFile <- file.path(collate_path, "rowEvent.mapGenes.fst")
    if(!file.exists(mapperFile)) .log(paste(
        "Required file", mapperFile, "not found."
    ))
    
    allEvents <- read.fst(mapperFile)
    allEvents <- allEvents[allEvents$EventName %in% EventNames]
    
    ont <- as.data.table(fst::read.fst(ontFile), 
        columns = c("go_id", "ensembl_id"))
    
    GO_gene_id <- ont$ensembl_id[ont$go_id == go_id]
    
    allEvents <- allEvents[
        allEvents$gene_id %in% GO_gene_id |
        allEvents$gene_id %in% GO_gene_id_b,]
    
    if(nrow(allEvents) == 0) return(NULL)
    return(allEvents$EventName)
}

.extract_gene_ids_for_GO <- function(
    enrichedEventNames,
    universeEventNames = NULL,
    collate_path
) {
    mapperFile <- file.path(collate_path, "rowEvent.mapGenes.fst")
    if(!file.exists(mapperFile)) .log(paste(
        "Required file", mapperFile, "not found."
    ))
    
    allEvents <- read.fst(mapperFile)
    
    uniqueEventnames <- unique(c(enrichedEventNames, universeEventNames))
    if(!all(uniqueEventnames %in% allEvents$EventName)) {
        nonMatches <- uniqueEventnames[!(uniqueEventnames %in% 
            allEvents$EventName)]
        .log(paste("One or more EventName(s) not found in reference!",
            "Culprit examples:", head(nonMatches)
        ))
    }
    
    genes <- unique(c(
        allEvents$gene_id[match(
            enrichedEventNames, allEvents$EventName)],
        allEvents$gene_id_b[match(
            enrichedEventNames, allEvents$EventName)]
    ))
    
    if(!is.null(universeEventNames)) {
        universe <- unique(c(
            allEvents$gene_id[match(
                universeEventNames, allEvents$EventName)],
            allEvents$gene_id_b[match(
                universeEventNames, allEvents$EventName)]
        ))    
    } else {
        GeneFile <- file.path(collate_path, "Reference/fst/Genes.fst")
        if(!file.exists(GeneFile))
            .log(paste("Gene reference", GeneFile, "not found"))

        Genes <- read.fst(GeneFile, columns = "gene_id")
        universe <- Genes$gene_id
    }
    
    # Remove novelGene ID's (they don't exist in ontologies)
    genes <- genes[!grepl("novelGene", genes)]
    universe <- universe[!grepl("novelGene", universe)]
    
    final <- list(
        genes = genes,
        universe = universe
    )
    return(final)
}

goASE <- function(
    enrichedEventNames,
    universeEventNames = NULL,
    collate_path,
    ontologyType = c("BP", "MF", "CC"),
    pAdjustMethod = c("BH", "holm", "hochberg", "hommel", 
        "bonferroni", "BY", "fdr", "none"),
    ...
) {
    ontologyType <- match.arg(ontologyType)
    if(ontologyType == "") 
        .log("ontologyType must be one of `BP`, `MF` or `CC`")
    pAdjustMethod <- match.arg(pAdjustMethod)
    if(pAdjustMethod == "") 
        .log(paste("pAdjustMethod must a valid option.",
            "See details section in ?p.adjust"))

    geneIds <- .extract_gene_ids_for_GO(
        enrichedEventNames,
        universeEventNames,
        collate_path
    )
    
    res <- .ora_internal(collate_path, geneIds$genes, geneIds$universe, 
        ontologyType, pAdjustMethod, ...)
    
    return(res)
}

################################################################################

# Global functions for gene ontology

.fetch_orgDB <- function(
    species = c("Homo sapiens", "Mus musculus"),
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    species <- match.arg(species)
    if(species == "") 
        .log("Species for orgDB must be Homo sapiens or Mus musculus")

    ah_orgList <- subset(ah, ah$rdataclass == "OrgDb")
    ah_orgDb <- subset(ah_orgList, ah_orgList$species == species)

    cache_loc <- AnnotationHub::cache(ah_orgDb[1])

    return(cache_loc)
}

.fetch_GOterms <- function() {
    .check_package_installed("GO.db", "3.12.0")
    xx <- as.list(GO.db::GOTERM)
    go_DT <- data.table::rbindlist(
        lapply(xx, function(x) data.frame(
            go_id = x@GOID, go_term = x@Term
        ))
    )
    return(go_DT)
}

.get_geneGO <- function(
    species = c("Homo sapiens", "Mus musculus"),
    verbose = TRUE,
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    .check_package_installed("DBI", "1.0.0")

    if(verbose) .log("Retrieving gene GO-term pairings", "message")
    cache_loc <- .fetch_orgDB(species, localHub, ah)
    dbcon <- DBI::dbConnect(
        DBI::dbDriver("SQLite"), 
        dbname = cache_loc
    )
    genes_DT <- as.data.table(DBI::dbGetQuery(dbcon, paste(
      "SELECT *",
      "FROM go",
      "LEFT JOIN ensembl",
      "ON go._id = ensembl._id"
    )))
    DBI::dbDisconnect(dbcon)

    if(verbose) .log("Retrieving GO terms from GO.db", "message")
    GO_DT <- .fetch_GOterms()

    final_DT <- genes_DT[GO_DT, on = "go_id"]
    final_DT <- final_DT[complete.cases(final_DT)]
    final_DT[, c("_id", "_id.1") := list(NULL, NULL)]

    return(final_DT)
}