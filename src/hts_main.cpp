/* hts_main.cpp Htslib benchmarking

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

#include <chrono>
#include "hts_main.h"

// [[Rcpp::export]]
int SpliceWizMain_hts(
    std::string bam_file, std::string reference_file, std::string output_file,
    bool verbose, int n_threads, int read_pool
) {
  
  std::string s_output_txt = output_file + ".txt.gz";
  std::string s_output_cov = output_file + ".cov";
  
  swEngine_hts Engine;
   
  std::string s_bam = bam_file;
  std::string s_ref = reference_file;  
  if(!Engine.checkFileExists(s_bam)) {
    cout << "File " << s_bam << " does not exist!\n";
    return(-1);
  } 
  if(!Engine.checkFileExists(s_ref)) {
    cout << "File " << s_ref << " does not exist!\n";
    return(-1);
  } 

  int use_threads = Engine.Set_Threads(n_threads);
  int ret = Engine.readReference(s_ref, verbose);
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid SpliceWiz reference\n";
    return(ret);
  }
  
  if(verbose) {
    cout << "Running SpliceWiz (htslib) " << s_bam;
    // if(Has_OpenMP() != 0) cout << " with OpenMP ";
    cout << "using " << use_threads << " threads"
      << "\n" << "Reference: " << s_ref << "\n"
      << "Output file: " << s_output_txt << "\t" << s_output_cov << "\n\n"
      << "Reading reference file\n";
  }

  // main:
  auto start = chrono::steady_clock::now();
  auto check = start;
  int ret2 = Engine.SpliceWizCore(
      s_bam, s_output_txt, s_output_cov,
      verbose, read_pool
  );
    
  if(ret2 == -2) {
    cout << "Process interrupted running SpliceWiz on " << s_bam << '\n';
    return(ret2);
  } else if(ret2 == -1) {
    cout << "Error encountered processing " << s_bam << "\n";
  } else {
    check = chrono::steady_clock::now();
    auto time_sec = chrono::duration_cast<chrono::seconds>(check - start).count();
    cout << s_bam << " processed (" << time_sec << " seconds)\n";
  }
  
  return(ret2);
}

// [[Rcpp::export]]
int c_BAM2COV_hts(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, int read_pool
){
  
  swEngine_hts Engine;
  int nthr = Engine.Set_Threads(n_threads);
  if(verbose) {
    cout << "Running BAM2COV (htslib) " << bam_file;
    cout << " using " << nthr << " threads\n";
  }
  
  int ret = Engine.BAM2COVcore(bam_file, output_file, verbose, read_pool);
  return(ret);
}

// [[Rcpp::export]]
int c_BAM2COV_hts_serial(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, int read_pool
){
  
  swEngine_hts Engine;
  int nthr = Engine.Set_Threads(n_threads);
  if(verbose) {
    cout << "Running BAM2COV (htslib) " << bam_file;
    cout << " using " << nthr << " threads\n";
  }
  
  int ret = Engine.BAM2COVcore_serial(bam_file, output_file, verbose, read_pool);
  return(ret);
}