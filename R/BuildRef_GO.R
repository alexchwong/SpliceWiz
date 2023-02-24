.process_ontology <- function(
        reference_path, species, verbose = TRUE
) {
    hasPackage <- 
        .check_package_installed("DBI", "1.0.0", "silent") && 
        .check_package_installed("GO.db", "3.12.0", "silent")

    if(!hasPackage) {
        .log(paste("Packages DBI and GO.db are required",
            "for gene ontology annotations, skipping..."
        ), "message")
        return()
    }
    if(!is_valid(species)) {
        if(verbose) 
            .log("Gene ontology not prepared for this reference", "message")
        return()
    } else {
        ontDT <- .get_geneGO(species, verbose)
    }
    
    if(!is.null(ontDT))
        fst::write.fst(ontDT, file.path(reference_path, "fst", "Ontology.fst"))
}

.getOntologySpecies <- function(genome_type) {
    if(genome_type %in% c("mm10", "mm9")) {
        return("Mus musculus")
    } else if(genome_type %in% c("hg38", "hg19")) {
        return("Homo sapiens")
    } else {
        return("")
    }
}

.ora_internal <- function(
    se, genes, universe, ontologyType = "BP", pAdjustMethod="BH", ...
) {
    .check_package_installed("fgsea", "1.0.0", "silent")

    ont <- ref(se)$ontology
    if(is.null(ont)) 
        .log("No gene ontology reference found for this NxtSE object")
    
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

#' @describeIn Gene-ontology-methods Returns a list of ASEs enriched in a given
#' gene ontology category
#' @export
subset_EventNames_by_GO <- function(
    EventNames, go_id, se
) {
    ont <- ref(se)$ontology
    if(is.null(ont)) 
        .log("No gene ontology reference found for this NxtSE object")
    
    allEvents <- as.data.frame(rowData(se))
    allEvents <- allEvents[, c("EventName", "gene_id", "gene_id_b")]
    allEvents <- allEvents[allEvents$EventName %in% EventNames,]
    
    GO_gene_id <- ont$ensembl_id[ont$go_id == go_id]
    
    allEvents <- allEvents[
        allEvents$gene_id %in% GO_gene_id |
        allEvents$gene_id_b %in% GO_gene_id,]
    
    if(nrow(allEvents) == 0) return(NULL)
    return(allEvents$EventName)
}

#' @describeIn Gene-ontology-methods Produces a list containing enriched and
#'   universe `gene_id`s of given enriched and background ASE `EventName`s
#' @export
extract_gene_ids_for_GO <- function(
    enrichedEventNames,
    universeEventNames = NULL,
    se
) {
    ont <- ref(se)$ontology
    if(is.null(ont)) 
        .log("No gene ontology reference found for this NxtSE object")
    
    allEvents <- as.data.frame(rowData(se))
    allEvents <- allEvents[, c("EventName", "gene_id", "gene_id_b")]
    
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
        Genes <- ref(se)$geneList
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


#' Gene ontology (over-representation) analysis using enriched genes of 
#' top alternative splicing events
#'
#' Genes containing differential alternative splicing events (ASEs) may be
#' enriched in key functional pathways. This can be identified using a simple
#' over-representation analysis. Biologists can identify key pathways of
#' interest in order to focus on studying ASEs belonging to genes of functional
#' interest.
#'
#' @param EventNames,go_id In `subset_EventNames_by_GO()`, a vector of ASE 
#'   `EventName`s to subset against the given `go_id`.
#' @param enrichedEventNames A vector of `EventName`s. This is typically one
#'   or more `EventName`s of differential ASEs
#' @param universeEventNames A vector of `EventName`s, typically the 
#'   `EventName`s of all ASEs that were tested. If left as `NULL`, all genes
#'   are considered background genes.
#' @param se The NxtSE object containing the GO reference and the `EventName`s
#' @param ontologyType One of either `"BP"` - biological pathways, `"MF"` -
#'   molecular function, or `"CC"` - cellular component.
#' @param pAdjustMethod The method for p-value adjustment for FDR.
#'   See `?p.adjust`
#' @param ... Additional arguments to be passed to `fgsea::fora()`
#' @return 
#' For `goASE()`, A data table containing the following:
#'    * go_id: Gene ontology ID
#'    * go_term: Gene ontology term
#'    * pval: Raw p values
#'    * padj: Adjusted p values
#'    * overlap: Number of enriched genes (of enriched ASEs)
#'    * size: Number of background genes (of background ASEs)
#'    * overlapGenes: A list of `gene_id`'s from genes of enriched ASEs
#'    * expected: The number of overlap genes expected by random 
#' 
#' For `extract_gene_ids_for_GO()`, A list containing the following:
#'    * genes: A vector of enriched `gene_id`s
#'    * universe: A vector of background `gene_id`s
#'
#' For `subset_EventNames_by_GO()`, A vector of all ASE `EventName`s belonging
#'   to the given gene ontology `go_id`
#'
#' @examples
#' \dontrun{
#' # Below example code of how to use output of differential ASEs for GO analysis
#' # It will not work because the reference must be either human / mouse, or a
#' # valid `ontologySpecies` given to buildRef()
#' 
#' se <- SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' 
#' colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' 
#' require("limma")
#' res_limma <- ASE_limma(se, "treatment", "A", "B")
#'
#' # Perform gene ontology analysis of the first 10 differential ASEs
#' 
#' go_df <- goASE(
#'   enrichedEventNames = res_limma$EventName[1:10],
#'   universeEventNames = res_limma$EventName,
#'   se = se
#' )
#'
#' # Return a list of all ASEs belonging to the top enriched category
#'
#' GOsubset_EventName <- subset_EventNames_by_GO(
#'   EventNames = res_limma$EventName,
#'   go_id = go_df$go_id[1],
#'   se = se
#' )
#' 
#' # Return a list of all ASEs belonging to the top enriched category.
#' # - typically used if one wishes to export `gene_id` for use in other gene
#' #   ontology tools
#'
#' gene_id_list <- extract_gene_ids_for_GO(
#'   enrichedEventNames = res_limma$EventName[1:10],
#'   universeEventNames = res_limma$EventName,
#'   se = se
#' )
#'
#' }
#' @name Gene-ontology-methods
#' @md
NULL

#' @describeIn Gene-ontology-methods Performs over-representation gene ontology
#' analysis using a given list of enriched / background ASEs
#' @export
goASE <- function(
    enrichedEventNames,
    universeEventNames = NULL,
    se,
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

    geneIds <- extract_gene_ids_for_GO(
        enrichedEventNames,
        universeEventNames,
        se
    )
    
    res <- .ora_internal(se, geneIds$genes, geneIds$universe, 
        ontologyType, pAdjustMethod, ...)
    
    return(res)
}

################################################################################

# Global functions for gene ontology

.check_GO_species <- function(
    species = "",
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    ah_orgList <- subset(ah, ah$rdataclass == "OrgDb")
    ah_orgListEns <- query(ah_orgList, "Ensembl")
    
    supportedSpecies <- unique(ah_orgListEns$species)
    if(!(species %in% supportedSpecies)) {
        .log(paste(
            species, 
            "not supported in AnnotationHub. Supported species:",
            supportedSpecies
        ), message)
        return("")
    } else {
        return(species)
    }
}

.fetch_orgDB <- function(
    species = "",
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    ah_orgList <- subset(ah, ah$rdataclass == "OrgDb")
    ah_orgListEns <- query(ah_orgList, "Ensembl")
    
    supportedSpecies <- unique(ah_orgListEns$species)
    if(!(species %in% supportedSpecies))
        .log(paste(
            species, 
            "not supported in AnnotationHub. Supported species:",
            supportedSpecies
        ))
        
    ah_orgDb <- subset(ah_orgListEns, ah_orgList$species == species)
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
    species <- .check_GO_species(species, localHub, ah)
    if(species == "") return(NULL)

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