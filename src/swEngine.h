/* swEngine.h Core SpliceWiz processBAM engine

Copyright (C) 2021 Alex Chit Hei Wong

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

#ifndef CODE_SWENGINE
#define CODE_SWENGINE

#include "includedefine.h"
#include "SpliceWiz.h"    // Full Rcpp functionality incl RcppProgress

#ifdef _OPENMP
#include <omp.h>
#endif

#include <chrono>

#include "BAM2blocks.h"       // For BB
#include "covTools.h"         // For COV I/O
#include "FastaReader.h"
#include "GZTools.h"          // For gzip I/O
#include "ReadBlockProcessor_CoverageBlocks.h"  // includes FragmentsMap and others
#include "ReadBlockProcessor_TandemJunctions.h"

class swEngine {
  private:
    std::vector<std::string> ref_names; 
    std::vector<std::string> ref_alias;
    std::vector<uint32_t> ref_lengths;
    std::string CB_string;
    std::string SP_string;
    std::string ROI_string;
    std::string JC_string;
    std::string TJ_string;
    unsigned int n_threads_to_use;

    std::vector<CoverageBlocksIRFinder> oCB;
    std::vector<SpansPoint> oSP;
    std::vector<FragmentsInROI> oROI;
    std::vector<FragmentsInChr> oChr;
    std::vector<JunctionCount> oJC;
    std::vector<TandemJunctions> oTJ;
    std::vector<FragmentsMap> oFM;
    std::vector<BAM2blocks> BBchild;

    bool refLoaded;
    bool CBloaded;
    bool SPloaded;
    bool ROIloaded;
    bool Chrloaded;
    bool JCloaded;
    bool TJloaded;
    bool FMloaded;
    bool BAMLoaded;
  public:
    swEngine();
    void clear();

    int Set_Threads(int n_threads);
    bool checkFileExists(const std::string& name);
    int ReadChrAlias(std::istringstream &IN);
    int readReference(std::string &reference_file, bool const verbose = FALSE);
    
    int loadReference();
    int loadReference(
      bool loadCB,
      bool loadSP,
      bool loadROI,
      bool loadChr,
      bool loadJC,
      bool loadTJ,
      bool loadFM
    );
    int refreshReference();
    int associateBAM(
      pbam_in * _IN, 
      std::vector<string> chr_name,
      std::vector<uint32_t> chr_len
    );

    int SpliceWizMultiCore(
      std::vector<std::string> &bam_file, 
      std::vector<std::string> &s_output_txt, 
      std::vector<std::string> &s_output_cov,
      bool const verbose,
      bool const multithreadedRead = false
    );

    int BAM2COVcore(
      std::string const &bam_file, 
      std::string const &s_output_cov,
      bool const verbose,
      bool const multithreadedRead = false
    );

    int doStatsCore(
      std::string const &bam_file, 
      std::string const &s_output_txt,
      bool const verbose,
      bool const multithreadedRead = false
    );

    int MappabilityRegionsCore(
      std::string const &bam_file,
      std::string const &s_output_txt,
      std::string const &s_output_cov,
      int threshold,
      bool const includeCov,
      bool const verbose,
      bool const multithreadedRead
    );
};

#endif