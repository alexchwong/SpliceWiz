/* hts_main.h Htslib-based SpliceWiz functions - for benchmarking

Copyright (C) 2021 Alex Chit Hei Wong
Copyright (C) 2016 William Ritchie
  - original: https://github.com/williamritchie/IRFinder/tree/IRFinder-1.3.1)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.  */

#ifdef WITH_HTSLIB

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
    bool const verbose, bool const skipCOV,
    int const read_pool = 1000000
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
    int max_threads = 1, bool verbose = true, 
    bool const skipCOV = false, int read_pool = 1000000
  );
  
  int SpliceWizMain_multi_hts(
      std::string reference_file, 
      StringVector bam_files, 
      StringVector output_files,
      int max_threads = 1, 
      bool verbose = true, bool const skipCOV = false, 
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

#endif