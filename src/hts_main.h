#ifndef CODE_MAIN_HTS
#define CODE_MAIN_HTS

#include "includedefine.h"
#include "SpliceWiz.h"    // Full Rcpp functionality incl RcppProgress

#ifdef _OPENMP
#include <omp.h>
#endif

#include "covTools.h"         // For COV I/O
#include "FastaReader.h"
#include "GZTools.h"          // For gzip I/O
#include "synthReadGenerator.h"

#include "swEngine_hts.h"

#ifdef SPLICEWIZ
  int SpliceWizCore_hts(
    swEngine_hts Engine,
    std::string const &bam_file, 
    std::string const &s_output_txt, 
    std::string const &s_output_cov,
    bool const verbose,
    int const read_pool
  );

  int SpliceWizMain_hts(
      std::string bam_file, std::string reference_file, std::string output_file, 
      bool verbose = true, int n_threads = 1, int read_pool = 1000000
  );

  int SpliceWizMain_hts_naked(
    std::string bam_file, std::string reference_file, std::string output_file,
    bool verbose, int n_threads, int read_pool
  );

  int c_BAM2COV_hts(
    std::string bam_file, std::string output_file, 
    bool verbose = true, int n_threads = 1, int read_pool = 1000000
  );
  
  int c_BAM2COV_hts_naked(
    std::string bam_file, std::string output_file, 
    bool verbose = true, int n_threads = 1, int read_pool = 1000000
  );  
  
#endif

#endif