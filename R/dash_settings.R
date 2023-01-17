# 
setreactive_system <- function() {
    reactiveValues(
        threads_initialized = FALSE,
        n_threads = 1
    )
}

# Settings for "Make New SpliceWiz Reference"
setreactive_newref <- function() {
    reactiveValues(
        newref_path = "",
        newref_fasta = "",
        newref_gtf = "",
        newref_AH_fasta = "",
        newref_AH_gtf = "",
        newref_mappa = "",
        newref_NPA = "",
        newref_bl = "",
        
        ui_newrefAH_Species = "",
        ui_newrefAH_Version_Trans = "",
        ui_newrefAH_Trans = "",
        ui_newrefAH_Assembly = "",
        ui_newrefAH_Version_Genome = "",
        ui_newrefAH_Genome = "",
        
        ui_newref_genome_type = ""
    )
}

# Settings for "Make New SpliceWiz Experiment"
setreactive_expr <- function() {
    reactiveValues(
        ref_path = "",
        bam_path = "",
        sw_path = "",
        anno_file = "",
        collate_path = "",
        collate_path_prompt = "",

        df.files = c(),
        df.anno = c(),
        disallow_df_update = FALSE,

        df.files_savestate = c(),
        df.anno_savestate = c(),

        ref_settings = c(),
        expr_path = "",
        selected_rows = c(),
        df = c(),
        se = NULL
    )
}

# Settings for navQC
setreactive_QC <- function() {
    reactiveValues(
        QC = c()
    )
}

# Settings for navFilter
setreactive_filtered_SE <- function() {
    reactiveValues(
        filterSummary = NULL,
        filters = list(),
        filtered_SE = NULL
    )
}

# Settings for navAnalyse
setreactive_DE <- function() {
    reactiveValues(
        res = NULL,
        res_settings = list(),
        method = NULL,
        batchVar1 = NULL,
        batchVar2 = NULL,
        DE_Var = NULL,
        nom_DE = NULL,
        denom_DE = NULL,
        IRmode_DE = NULL,
        dof = 1,
        
        DT_DE_rows_all = NULL,
        DT_DE_rows_selected = NULL,
        command_selected = NULL,
        filters = list()
    )
}

setreactive_GO <- function() {
    reactiveValues(
        res = NULL,

        # gene_ids that can optionally be saved to file by user
        gene_ids = NULL,
        univ_ids = NULL,

        # Save GO results
        filteredVolc = NULL,
        
        resGO = NULL,

        final_plot = NULL
    )
}

# Settings for Diag and Volcano
setreactive_Diag <- function() {
    # NB same code as Volcano
    reactiveValues(
        useDE = NULL,

        plot_ini = FALSE,
        plotly_click = NULL,
        final_plot = NULL,
        ggplot = NULL,
        selected = NULL
    )
}

# Settings for Heatmap
setreactive_Heat <- function() {
    reactiveValues(
        useDE = NULL,
        
        ggplot = NULL,
        final_plot = NULL
    )
}

# Settings for Coverage plots
setreactive_Cov <- function() {
    reactiveValues(
        useDE = NULL,
        
        view_chr = "",
        view_start = "",
        view_end = "",
        data_start = 0,
        data_end = 0,

        view_strand = "*",

        event.ranges = NULL,

        plotly_relayout = NULL,
        plot_ini = FALSE,
        
        final_plot = NULL,
        
        trigger = NULL,
        plot_params = NULL
    )
}