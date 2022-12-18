#ifndef CODE_MAIN
#define CODE_MAIN

#include "includedefine.h"
#include "SpliceWiz.h"    // Full Rcpp functionality incl RcppProgress

#ifdef _OPENMP
#include <omp.h>
#endif

#include "BAM2blocks.h"       // For BB
#include "covTools.h"         // For COV I/O
#include "FastaReader.h"
#include "GZTools.h"          // For gzip I/O
#include "ReadBlockProcessor_CoverageBlocks.h"  // includes FragmentsMap and others
#include "ReadBlockProcessor_TandemJunctions.h"

#include "synthReadGenerator.h"

#include "swEngine.h"

bool checkFileExists(const std::string& name);

int Has_OpenMP();
int Test_OpenMP_For();

bool c_Check_Cov(std::string s_in);

// #############################################################################

#ifdef SPLICEWIZ
// Rcpp-only functions
  List c_RLE_From_Cov(
    std::string s_in, std::string seqname, int start, int end, int strand
  );

  StringVector c_Cov_Seqnames(std::string s_in);

  List c_RLEList_From_Cov(std::string s_in, int strand);

  List c_gunzip_DF(std::string s_in, StringVector s_header_begin);
  
#endif

int c_gunzip(std::string s_in, std::string s_out);

#ifdef SPLICEWIZ
  int SpliceWizMain(
      std::string bam_file, std::string reference_file, std::string output_file, 
      bool verbose = true, int n_threads = 1, bool multiRead = false
  );

  int SpliceWizMain_multi(
      std::string reference_file, StringVector bam_files, StringVector output_files,
      int max_threads = 1, bool verbose = true, bool multiRead = false
  );

  int c_GenerateMappabilityReads(
    std::string genome_file, std::string out_fa,
    int read_len, int read_stride, int error_pos
  );

  int c_GenerateMappabilityRegions(
    std::string bam_file, std::string output_file, 
    int threshold, int includeCov = 0, bool verbose = true,
    int n_threads = 1
  );

  int c_BAM2COV(
    std::string bam_file, std::string output_file, 
    bool verbose = true, int n_threads = 1, bool multiRead = false
  );

  int c_doStats(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, bool multiRead
  );

#else
  int SpliceWizMain(
      std::string bam_file, std::string reference_file, std::string s_output_txt,
      std::string s_output_cov, int n_threads = 1
  );

  int c_GenerateMappabilityReads(
    std::string genome_file, std::string out_fa,
    int read_len, int read_stride, int error_pos
  );

  int c_GenerateMappabilityRegions(
    std::string bam_file, std::string s_output_txt, 
    int threshold, int n_threads = 1, std::string s_output_cov = ""
  );	

  int c_BAM2COV(
    std::string bam_file, std::string output_file, 
    int n_threads = 1, bool multiRead = false
  );

  int c_doStats(
      std::string bam_file, std::string output_file, int n_threads, bool multiRead
  );

  int main(int argc, char * argv[]);
#endif

#endif