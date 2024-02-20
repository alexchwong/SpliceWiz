#' Gene ontology (over-representation) analysis using enriched genes of 
#' top alternative splicing events
#'
#' Genes containing differential alternative splicing events (ASEs) may be
#' enriched in key functional pathways. This can be identified using a simple
#' over-representation analysis. Biologists can identify key pathways of
#' interest in order to focus on studying ASEs belonging to genes of functional
#' interest.
#'
#' @details
#' Users can perform GO analysis using either the GO annotation compiled via
#' building the SpliceWiz reference using `buildRef()` , or via a 
#' custom-supplied gene ontology annotation. This is done by
#' supplying their own GO annotations as an argument to `ontologyRef`. This
#' should be coerceable to a `data.frame` containing the following columns:
#'
#' * `gene_id` Gene ID's matching that used by the SpliceWiz reference
#' * `go_id` Gene ontology ID terms, of the form `GO:XXXXXXX`
#'
#' @param enrichedGenes A vector of `gene_id` representing the list of enriched
#'   genes. To generate a list of valid `gene_id`, see [viewGenes]
#' @param universeGenes (default `NULL`) A vector of `gene_id` representing the 
#'  list of background genes.
#' @param ontologyRef A valid gene ontology reference. This can be generated
#'   either using `viewGO(reference_path)` or `ref(se)$ontology`. This field
#'   is required for `goGenes()` and optional for `goASE()`. See details.
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
#' @param res_go For `plotGO`, the gene ontology results data object returned
#'   by the `goASE()` function.
#' @param plot_x,plot_size,plot_color What parameters should be plotted on
#'   the x-axis, bubble-size, or bubble-color? Valid options are
#'   `c("log10FDR", "foldEnrichment", "nGenes"). Defaults are
#'   `"log10FDR", "nGenes", "foldEnrichment"` for x-axis, bubble size/color,
#'   respectively
#' @param filter_n_terms (default `20`) How many top terms to plot.
#' @param filter_padj,filter_pvalue Whether given GO results should be filtered
#'   by adjusted p value (FDR) or nominal p value, respectively, prior to plot
#' @param trim_go_term (default `50`) For long GO terms, description will be
#'   trimmed by first n characters, where `trim_go_term = n`
#' @param ... Additional arguments to be passed to `fgsea::fora()`
#' @return 
#' For `goASE()` and `goGenes()`, a data table containing the following:
#'    * go_id: Gene ontology ID
#'    * go_term: Gene ontology term
#'    * pval: Raw p values
#'    * padj: Adjusted p values
#'    * overlap: Number of enriched genes (of enriched ASEs)
#'    * size: Number of background genes (of background ASEs)
#'    * overlapGenes: A list of `gene_id`'s from genes of enriched ASEs
#'    * expected: The number of overlap genes expected by random 
#' 
#' For `extract_gene_ids_for_GO()`, a list containing the following:
#'    * genes: A vector of enriched `gene_id`s
#'    * universe: A vector of background `gene_id`s
#'
#' For `subset_EventNames_by_GO()`, a vector of all ASE `EventName`s belonging
#'   to the given gene ontology `go_id`
#'
#' @examples
#' # Generate example reference with `Homo sapiens` gene ontology
#'
#' ref_path <- file.path(tempdir(), "Reference_withGO")
#' buildRef(
#'     reference_path = ref_path,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf(),
#'     ontologySpecies = "Homo sapiens"
#' )
#'
#' # Perform GO analysis using first 1000 genes
#' ontology <- viewGO(ref_path)
#' allGenes <- sort(unique(ontology$gene_id))
#'
#' exampleGeneID <- allGenes[1:1000]
#' exampleBkgdID <- allGenes
#' 
#' go_df <- goGenes(
#'     enrichedGenes = exampleGeneID, 
#'     universeGenes = exampleBkgdID, 
#'     ontologyRef = ontology
#' )
#'
#' # Plots the top 12 GO terms
#' 
#' plotGO(go_df, filter_n_terms = 12)
#'
#' # Below example code of how to use output of differential ASEs for GO analysis
#' # It will not work with the example dataset because the reference must be 
#' # either human / mouse, or a  valid `ontologySpecies` given to buildRef()
#' # We hope the example code is simple enough to understand for users to adapt
#' # to their own workflows.
#'
#' \dontrun{
#' 
#' se <- SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' 
#' colData(se)$treatment <- rep(c("A", "B"), each = 3)
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
#' @seealso
#' [Build-Reference-methods] on how to generate gene ontology annotations\cr\cr
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
    ontologyRef = NULL,
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

    ont <- NULL
    if(is.null(ontologyRef)) {
        ont <- ref(se)$ontology
        if(is.null(ont)) 
            .log("No gene ontology reference found for this NxtSE object")
        ont <- .validate_GO_ref(ont)
    } else {
        ont <- .validate_GO_ref(ontologyRef)
    }   

    res <- .ora_internal(ont, geneIds$genes, geneIds$universe, 
        ontologyType, pAdjustMethod, ...)
    
    return(res)
}

#' @describeIn Gene-ontology-methods Performs GO analysis given the set of
#'   enriched and (optionally) the background (universe) genes.
#' @export
goGenes <- function(
    enrichedGenes, universeGenes = NULL, ontologyRef, 
    ontologyType = c("BP", "MF", "CC"),
    pAdjustMethod = c("BH", "holm", "hochberg", "hommel", 
        "bonferroni", "BY", "fdr", "none"),
    ...
) {
    ont <- .validate_GO_ref(ontologyRef)
    ontologyType <- match.arg(ontologyType)
    if(ontologyType == "") 
        .log("ontologyType must be one of `BP`, `MF` or `CC`")
    pAdjustMethod <- match.arg(pAdjustMethod)
    if(pAdjustMethod == "") 
        .log(paste("pAdjustMethod must a valid option.",
            "See details section in ?p.adjust"))

    .ora_internal(
        ont, enrichedGenes, universeGenes, 
        ontologyType, pAdjustMethod, ...
    )
}

#' @describeIn Gene-ontology-methods Produces a list containing enriched and
#'   universe `gene_id`s of given enriched and background ASE `EventName`s
#' @export
extract_gene_ids_for_GO <- function(
    enrichedEventNames,
    universeEventNames = NULL,
    se
) {
    # ont <- ref(se)$ontology
    # if(is.null(ont)) 
        # .log("No gene ontology reference found for this NxtSE object")
    # ont <- .validate_GO_ref(ont)
    
    allEvents <- as.data.frame(rowData(se))
    allEvents <- allEvents[, c("EventName", "gene_id", "gene_id_b")]
    
    uniqueEventnames <- unique(c(enrichedEventNames, universeEventNames))
    if(!all(uniqueEventnames %in% allEvents$EventName)) {
        nonMatches <- uniqueEventnames[!(uniqueEventnames %in% 
            allEvents$EventName)]
        .log(paste("One or more EventName(s) not found in reference!",
            "Culprit examples:", paste(head(nonMatches), collapse = "; ")
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

#' @describeIn Gene-ontology-methods Returns a list of ASEs enriched in a given
#' gene ontology category
#' @export
subset_EventNames_by_GO <- function(
    EventNames, go_id, se
) {
    ont <- ref(se)$ontology
    if(is.null(ont)) 
        .log("No gene ontology reference found for this NxtSE object")
    ont <- .validate_GO_ref(ont)
    
    allEvents <- as.data.frame(rowData(se))
    allEvents <- allEvents[, c("EventName", "gene_id", "gene_id_b")]
    allEvents <- allEvents[allEvents$EventName %in% EventNames,]
    
    GO_gene_id <- ont$gene_id[ont$go_id == go_id]
    
    allEvents <- allEvents[
        allEvents$gene_id %in% GO_gene_id |
        allEvents$gene_id_b %in% GO_gene_id,]
    
    if(nrow(allEvents) == 0) return(NULL)
    return(allEvents$EventName)
}

#' @describeIn Gene-ontology-methods Produces a lollipop plot based on the
#' given gene ontology results object
#' @export
plotGO <- function(
    res_go = NULL,
    plot_x = c("log10FDR", "foldEnrichment", "nGenes"),
    plot_size = c("nGenes", "foldEnrichment", "log10FDR"),
    plot_color = c("foldEnrichment", "nGenes", "log10FDR"),
    filter_n_terms = 20,
    filter_padj = 1, # don't filter by default
    filter_pvalue = 1, # don't filter by default
    trim_go_term = 50
) {
    if(is.null(res_go)) return(NULL)
    
    res <- .format_GO_result(res_go, trim_go_term)
    return(
        .generate_ggplot_GO(
            res, plot_x, plot_size, plot_color,
            filter_n_terms, filter_padj, filter_pvalue, trim_go_term
        )
    )
}

################################################################################

### INTERNALS - buildRef - ###

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
    
    genes <- viewGenes(reference_path)
    
    if(!is_valid(species)) {
        if(verbose) 
            .log("Gene ontology not prepared for this reference", "message")
        return()
    } else {
        ontDT <- .get_geneGO(species, genes$gene_id, genes$gene_name, verbose)
    }
    
    if(!("gene_id" %in% colnames(ontDT))) {
        ontDT[, c("gene_id") := genes$gene_id[match(
            get("gene_name"), genes$gene_name
        )]]
        ontDT <- ontDT[!is.na("gene_id")]
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

.build_GO_table <- function(ont) {
    if(!is(ont, "data.table")) ont <- as.data.table(ont)
    # Reduce to 2-column ont (+/- evidence)
    if("evidence" %in% names(ont)) {
        genes_DT <- ont[, c("gene_id", "go_id", "evidence"), with = FALSE]    
    } else {
        genes_DT <- ont[, c("gene_id", "go_id"), with = FALSE]    
    }
    genes_DT <- unique(genes_DT)

    GO_DT <- .fetch_GOterms()
    genes_DT[, c("go_term", "ontology") := list(
        GO_DT$go_term[match(get("go_id"), GO_DT$go_id)],
        GO_DT$Ontology[match(get("go_id"), GO_DT$go_id)]
    )]

    genes_DT <- genes_DT[complete.cases(genes_DT)]
    return(genes_DT)
}

.validate_GO_ref <- function(ont) {
    ont_cols <- colnames(ont)
    if("ensembl_id" %in% ont_cols) {
        # support legacy SpliceWiz references
        colnames(ont)["ensembl_id" %in% ont_cols] <- "gene_id"
    }
    if(!all(c("go_id", "ontology", "go_term", "gene_id") %in% ont_cols)) {
        # See if minimal GO is satisfied; repair if required
        if( all(c("gene_id", "go_id") %in% ont_cols) ) {
            ont <- .build_GO_table(ont)        
        } else {
            .log("Given ontology reference is not valid")        
        }
    }
    return(ont)
}

### INTERNALS - fgsea::fora() wrapper - ###

.ora_internal <- function(
    ont, genes, universe, ontologyType = "BP", pAdjustMethod="BH", ...
) {
    .check_package_installed("fgsea", "1.0.0", "silent")

    if(is.null(ont)) 
        .log("Invalid gene ontology reference: `ont`")
    
    # should be done upstream
    # ont <- .validate_GO_ref(ont)
    
    if(!is(ont, "data.table")) {
        ont <- as.data.table(ont)
    }
    
    ontUse <- ont[get("ontology") == ontologyType]
    if(nrow(ontUse) == 0) .log(paste(
        ontologyType, "not found as a gene ontology category"
    ))
    
    ontUse <- unique(ontUse, by = c("gene_id", "go_id"))
    
    pathways <- split(ontUse$gene_id, ontUse$go_id)
    
    genes <- .gencode_correct_id_batch(genes)
    universe <- .gencode_correct_id_batch(universe)
    
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

################################################################################

# Global functions for gene ontology

# Check if given species is available
.check_GO_species <- function(
    species = "",
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    supportedSpecies <- getAvailableGO(localHub, ah)
    
    if(!(species %in% supportedSpecies)) {
        .log(paste(
            species, 
            "not supported in AnnotationHub."
        ), message)
        return("")
    } else {
        return(species)
    }
}

.check_cached_resource <- function(cache_loc, ah_id, ah) {
    dbcon <- DBI::dbConnect(
        DBI::dbDriver("SQLite"), 
        dbname = cache_loc
    )
    is_resource_legit <- tryCatch(
        {
            genes_DT <- as.data.table(DBI::dbGetQuery(dbcon, paste(
              "SELECT *",
              "FROM go",
              "LEFT JOIN ensembl",
              "ON go._id = ensembl._id"
            )))
            TRUE
        }, error = function(e) {
            FALSE
        }
    )
    DBI::dbDisconnect(dbcon)
    if(!is_resource_legit) {
        removeResources(ah, ah_id)
        cache_loc <- AnnotationHub::cache(ah_id)
    }
    return(cache_loc)
}

# Retrieves cache file name of orgDB resource of given species
.fetch_orgDB <- function(
    species = "",
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    supportedSpecies <- getAvailableGO(localHub, ah)
    
    if(!(species %in% supportedSpecies))
        .log(paste(
            species, 
            "not supported in AnnotationHub. Supported species:",
            supportedSpecies
        ))

    ah_orgList <- subset(ah, ah$rdataclass == "OrgDb")
    ah_orgListEns <- query(ah_orgList, "Ensembl")
    ah_orgDb <- subset(ah_orgListEns, ah_orgListEns$species == species)
    cache_loc <- AnnotationHub::cache(ah_orgDb[1])
    cache_loc <- .check_cached_resource(cache_loc, names(ah_orgDb[1]), ah)
    return(cache_loc)
}

.fetch_orgDB_cache <- function(ah, species) {
    ah_orgList <- subset(ah, ah$rdataclass == "OrgDb")
    ah_orgDb <- ah_orgList[ah_orgList$species == species]
    cache_loc <- AnnotationHub::cache(ah_orgDb[1])
    
    return(cache_loc)
}

# Get GO.db's gene ontology terms
.fetch_GOterms <- function() {
    .check_package_installed("GO.db", "3.12.0")
    # xx <- as.list(GO.db::GOTERM)
    # go_DT <- data.table::rbindlist(
        # lapply(xx, function(x) data.frame(
            # go_id = x@GOID, go_term = x@Term
        # ))
    # )
    xx <- AnnotationDbi::toTable(GO.db::GOTERM)
    go_DT <- data.table(
        go_id = xx$go_id,
        go_term = xx$Term,
        Ontology = xx$Ontology
    )
    return(go_DT)
}

# Compile gene ontology annotations
.get_geneGO <- function(
    species = c("Homo sapiens", "Mus musculus"),
    gene_ids, gene_names,
    verbose = TRUE,
    localHub = FALSE, ah = AnnotationHub(localHub = localHub)
) {
    .check_package_installed("DBI", "1.0.0")
    species <- .check_GO_species(species, localHub, ah)
    if(species == "") return(NULL)

    if(verbose) .log("Retrieving gene GO-term pairings", "message")
    cache_loc <- .fetch_orgDB_cache(ah, species)
    dbcon <- DBI::dbConnect(
        DBI::dbDriver("SQLite"), 
        dbname = cache_loc
    )
    
    coverage <- 0
    genes_DT <- data.table()
    tableNames <- DBI::dbListTables(dbcon)
    
    # Test if 'ensembl' is a table
    if("ensembl" %in% tableNames) {
        genes_DT <- as.data.table(DBI::dbGetQuery(dbcon, paste(
          "SELECT *",
          "FROM go",
          "LEFT JOIN ensembl",
          "ON go._id = ensembl._id"
        )))
        coverage <- length(intersect(gene_ids, genes_DT$ensembl_id)) / 
            length(gene_ids)
    }

    if(coverage < 0.1) {
        genes_DT <- as.data.table(DBI::dbGetQuery(dbcon, paste(
          "SELECT *",
          "FROM go",
          "LEFT JOIN genes",
          "ON go._id = genes._id"
        )))
        coverage <- length(intersect(gene_ids, genes_DT$gene_id)) / 
            length(gene_ids)

        if(coverage < 0.1) {
            genes_DT <- as.data.table(DBI::dbGetQuery(dbcon, paste(
              "SELECT *",
              "FROM go",
              "LEFT JOIN gene_info",
              "ON go._id = gene_info._id"
            )))
            coverage <- length(intersect(gene_names, genes_DT$symbol)) / 
                length(gene_names)

            if(coverage < 0.1) {
                .log(paste(
                    "Gene ontology failed to match gene_id or gene_name entries,",
                    "skipping gene ontology..."
                ), "warning")
            } else {
                setnames(genes_DT, "symbol", "gene_name", skip_absent=TRUE)
                
                dup_symbols <- gene_names[duplicated(gene_names)]
                uniq_symbols <- gene_names[!(gene_names %in% dup_symbols)]
                
                genes_DT <- genes_DT[get("gene_name") %in% uniq_symbols]
                genes_DT[, c("gene_id") := gene_ids[
                    match(get("gene_name"), gene_names)
                ]]
            }
        }
    } else {
        setnames(genes_DT, "ensembl_id", "gene_id", skip_absent=TRUE)
    }
        
    DBI::dbDisconnect(dbcon)
    
    genes_DT <- .build_GO_table(copy(genes_DT))

    if(verbose) .log(paste(
        "Gene ontology coverage accounts for",
        round(coverage * 100, 2), "% of all genes in reference"
    ), "message")

    return(genes_DT)
}

### INTERNALS - plotting GO - ###

.generate_plot_GO_labels <- function(axis_term) {
    if(axis_term == "log10FDR") {
        return("Enrichment FDR (-log10)")
    } else if(axis_term == "nGenes") {
        return("Number of Enriched Genes")
    } else {
        return("Fold Enrichment")
    }
}

.format_GO_result <- function(res, trim_go_term = 50) {
    res$go_term <- substr(res$go_term, 1, trim_go_term + 1)
    res[nchar(get("go_term")) > 50, c("go_term") :=
        paste0(substr(get("go_term"), 1, trim_go_term-3), "...")]
    res$Term <- paste(res$go_term, res$go_id, sep = "~")
    res$Term <- factor(res$Term, res$Term, ordered = TRUE)
    
    res$FDR = res$padj
    res$log10FDR = -log10(res$FDR)
    res$nGenes = res$overlap

    return(res)
}

.generate_ggplot_GO <- function(
    res,
    plot_x = c("log10FDR", "foldEnrichment", "nGenes"),
    plot_size = c("nGenes", "foldEnrichment", "log10FDR"),
    plot_color = c("foldEnrichment", "nGenes", "log10FDR"),
    filter_n_terms = 20,
    filter_padj = 1, # don't filter by default
    filter_pvalue = 1, # don't filter by default
    trim_go_term = 50
) {
    if(nrow(res) == 0) return(NULL)
    
    plot_x <- match.arg(plot_x)
    plot_size <- match.arg(plot_size)
    plot_color <- match.arg(plot_color)

    res_use <- res[seq_len(filter_n_terms)]
    res_use <- res_use[get("pval") <= filter_pvalue]    
    res_use <- res_use[get("padj") <= filter_padj]

    res_use[, c("Information") := paste(
        get("Term"),
        paste0(plot_x, ": ", get(plot_x)),
        paste0(plot_size, ": ", get(plot_size)),
        paste0(plot_color, ": ", get(plot_color)),
        sep = "\n"
    )]

    setorderv(res_use, plot_x)

    p <- ggplot(res_use, aes(text = get("Information"))) + 
        geom_segment(data = res_use, mapping = aes(
            x = 0, xend = get(plot_x), 
            y = get("Term"), yend = get("Term"), 
            color = get(plot_color)
        )) +
        geom_point(data = res_use, mapping = aes(
            x = get(plot_x), y = get("Term"), 
            size = get(plot_size), color = get(plot_color)
        )) +
        scale_colour_gradient(low = "blue", high = "red") +
        scale_y_discrete(limits=rev) +
        labs(
            x = .generate_plot_GO_labels(plot_x),
            y = "Gene Ontology Term",
            color = .generate_plot_GO_labels(plot_color), 
            size = .generate_plot_GO_labels(plot_size)
        )
    
    return(p)
}

#

.gencode_correct_id_indiv <- function(geneId) {
    if(
        substr(geneId, 1, 3) == "ENS" &&
        length(tstrsplit(geneId, split = ".", fixed = TRUE)) == 2
    ) {
        return(tstrsplit(geneId, split = ".", fixed = TRUE)[[1]])
    } else {
        return(geneId)
    }
}

.gencode_correct_id_batch <- function(geneIds) {
    return(unname(vapply(geneIds, .gencode_correct_id_indiv, character(1))))
}