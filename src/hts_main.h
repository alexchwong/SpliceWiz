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
  int SpliceWizCore_htsMulti(
    swEngine_hts &Engine,
    std::vector<std::string> const &bam_file, 
    std::vector<std::string> const &s_output_txt, 
    std::vector<std::string> const &s_output_cov,
    bool const verbose,
    int const read_pool
  );

  int BAM2COVCore_hts(
    swEngine_hts &Engine,
    std::string const &bam_file, 
    std::string const &s_output_cov,
    bool const verbose,
    int const read_pool
  );

  int doStatsCore_hts(
    swEngine_hts &Engine,
    std::string const &bam_file, 
    std::string const &s_output_txt,
    bool const verbose,
    int const read_pool
  );

  int SpliceWizMain_hts(
    std::string reference_file, 
    StringVector bam_files, StringVector output_files,
    int max_threads = 1, bool verbose = true, int read_pool = 1000000
  );
  
  int SpliceWizMain_multi_hts(
      std::string reference_file, 
      StringVector bam_files, 
      StringVector output_files,
      int max_threads = 1, 
      bool verbose = true, 
      int read_pool = 1000000
  );

  int c_BAM2COV_hts(
    std::string bam_file, std::string output_file, 
    bool verbose = true, int n_threads = 1, int read_pool = 1000000
  );

  int c_doStats_hts(
    std::string bam_file, std::string output_file, 
    bool verbose = true, int n_threads = 1, int read_pool = 1000000
  );

  
#endif

#endif