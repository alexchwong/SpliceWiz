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
    reference_path,
    genes, universe,
    ontologyType = "BP",
    ...
) {
    .check_package_installed("fgsea", "1.0.0", "silent")

    ontFile <- file.path(reference_path, "fst/Ontology.fst")
    if(!file.exists()) .log(paste(
        "No gene ontology resource found in", reference_path))

    ont <- as.data.table(fst::read.fst(ontFile))
    ontUse <- ont[get("ontology") == ontologyType]
    if(nrows(ontUse) == 0) .log(paste(
        ontologyType, "not found as a gene ontology category"
    ))

    pathways <- split(ont$ensembl_id, ont$go_id)
    
    foraRes <- as.data.table(fora(
        pathways, 
        genes=genes, 
        universe=universe,
        ...
    ))

    foraHeader <- data.table(
        go_id = foraRes$pathways,
        go_term <- ont$go_term[match(
            foraRes$pathways,
            ont$go_id
        )]
    )
    
    final <- cbind(foraHeader, foraRes[, -1])
    return(final)
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