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

// SpliceWiz core:
int SpliceWizCore_hts(
    swEngine_hts Engine,
    std::string const &bam_file, 
    std::string const &s_output_txt, 
    std::string const &s_output_cov,
    bool const verbose,
    int const read_pool
) {

  if(!Engine.checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
 
	if(verbose) cout << "Processing BAM file " << bam_file << "\n";
  
  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);

  hts_tpool *pool;
  const int queue_size = 0;
  if (Engine.n_threads_to_use > 1) {
      pool = hts_tpool_init(Engine.n_threads_to_use);
      bgzf_thread_pool(fp, pool, queue_size);
  }

  // Abort here if BAM corrupt
  if(header->n_targets <= 0){
    cout << bam_file << " - contains no chromosomes mapped\n";
    return(-1);
  }
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  for (int i = 0; i < header->n_targets; ++i) {
    s_chr_names.push_back(header->target_name[i]);
    u32_chr_lens.push_back(header->target_len[i]);
  }
  
  // Compile here a list of chromosomes; use BAM chromosomes for order
  // Add reference-only chromosomes at the end
  std::vector<std::string> bam_chr_name;
  std::vector<uint32_t> bam_chr_len;
  for(unsigned int i = 0; i < s_chr_names.size(); i++) {
    for(unsigned int j = 0; j < Engine.ref_alias.size(); j++) {
      if( 0==strncmp(
            Engine.ref_alias.at(j).c_str(), 
            s_chr_names.at(i).c_str(), 
            s_chr_names.at(i).size()
          ) && s_chr_names.at(i).size() == Engine.ref_alias.at(j).size()
      ) {
        bam_chr_name.push_back(Engine.ref_names.at(j));
        bam_chr_len.push_back(u32_chr_lens.at(i));
        break;
      }
    }
    if(i == bam_chr_name.size()) {
      bam_chr_name.push_back(s_chr_names.at(i));
      bam_chr_len.push_back(u32_chr_lens.at(i));
    }
  }
  // Now fill in reference chromosomes not in BAM:
  for(unsigned int i = 0; i < Engine.ref_names.size(); i++) {
    auto it = std::find(bam_chr_name.begin(), bam_chr_name.end(), Engine.ref_names.at(i));
    if(it == bam_chr_name.end()) {
      bam_chr_name.push_back(Engine.ref_names.at(i));
      bam_chr_len.push_back(Engine.ref_lengths.at(i));      
    }
  }
  
  std::vector<CoverageBlocksIRFinder*> oCB(Engine.n_threads_to_use);
  std::vector<SpansPoint*> oSP(Engine.n_threads_to_use);
  std::vector<FragmentsInROI*> oROI(Engine.n_threads_to_use);
  std::vector<FragmentsInChr*> oChr(Engine.n_threads_to_use);
  std::vector<JunctionCount*> oJC(Engine.n_threads_to_use);
  std::vector<TandemJunctions*> oTJ(Engine.n_threads_to_use);
  std::vector<FragmentsMap*> oFM(Engine.n_threads_to_use);
  std::vector<htsBAM2blocks*> BBchild(Engine.n_threads_to_use);

  // Multi-threaded results container initialization
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
  #endif
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    oCB.at(i) = new CoverageBlocksIRFinder(Engine.CB_string);
    oSP.at(i) = new SpansPoint(Engine.SP_string);
    oROI.at(i) = new FragmentsInROI(Engine.ROI_string);
    oChr.at(i) = new FragmentsInChr;
    oJC.at(i) = new JunctionCount(Engine.JC_string);
    oTJ.at(i) = new TandemJunctions(Engine.TJ_string);
    oFM.at(i) = new FragmentsMap;
    BBchild.at(i) = new htsBAM2blocks(bam_chr_name, bam_chr_len);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &(*oJC.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &(*oJC.at(i)), std::placeholders::_1) );

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&TandemJunctions::ChrMapUpdate, &(*oTJ.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&TandemJunctions::ProcessBlocks, &(*oTJ.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &(*oChr.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &(*oChr.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &(*oSP.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &(*oSP.at(i)), std::placeholders::_1) );
        
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &(*oROI.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &(*oROI.at(i)), std::placeholders::_1) );
    
    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &(*oCB.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &(*oCB.at(i)), std::placeholders::_1) );

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i)->initializeChrs();
  }
  
  // Read pool size
  unsigned int pool_cap = (unsigned int)read_pool;
  
  // Initialize bam1_t vector
  std::vector<bam1_t *> bpool;
  for(unsigned int i = 0; i < pool_cap; i++) {
    bpool.push_back(bam_init1());
  }
  
  // Pre-partition for n threads
  std::vector<int> pool_starts;
  std::vector<int> pool_ends;
  int est_tp_size = 1 + (pool_cap / Engine.n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    if(poolPos + est_tp_size > (int)pool_cap) {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(pool_cap - 1);
      poolPos = pool_cap;
    } else {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(poolPos + est_tp_size - 1);
      poolPos += est_tp_size;
    }
  }
  for(unsigned int i = pool_starts.size(); i < Engine.n_threads_to_use; i++) {
    pool_starts.push_back(-1);
    pool_ends.push_back(-1);
  }
  
  // BAM processing loop
  bool error_detected = false;
  off_t prevPos = 0; 
  off_t curPos = 0;
#ifdef SPLICEWIZ
  Progress p(GetFileSize(bam_file), verbose);
  while(!p.check_abort()) {
    curPos = htell(fp->fp);
    p.increment(curPos - prevPos);
    prevPos = curPos;
    
#else
  while(!p.check_abort()) {
#endif
    
    // Load n reads here and partition by thread
    unsigned int pool_size = 0;
    for(unsigned int i = 0; i < bpool.size(); i++) {
      int ret = bam_read1(fp, bpool.at(i));
      if(ret < 0) {
        break;
      } else {
        pool_size++;
      }
    }

    // End of file
    if(pool_size == 0) {
      break;
    }

    // Multi-threaded process reads
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      if(pool_starts.at(i) >= 0 && pool_starts.at(i) < (int)pool_size) {
        int true_end = pool_size;
        if(true_end > pool_ends.at(i)) {
          true_end = pool_ends.at(i) + 1; // [first, last)
        }

        int pa_ret = BBchild.at(i)->processAll(
          bpool, pool_starts.at(i), true_end
        );
        if(pa_ret == -1) {
          
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          error_detected = true;
        }
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < Engine.n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
    }
  }

  for(unsigned int i = 0; i < bpool.size(); i++) {
    bam_destroy1(bpool.at(i));
  }
  bam_hdr_destroy(header);
  bgzf_close(fp);

#ifdef SPLICEWIZ
  if(p.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      delete oJC.at(i);
      delete oTJ.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    if(error_detected) {
      return(-1);
    }
	// Process aborted; stop processBAM for all requests
    return(-2);
  }


  if(Engine.n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(Engine.n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < Engine.n_threads_to_use) {
          BBchild.at(i)->processSpares(*BBchild.at(i_new));
          BBchild.at(i)->processStats(*BBchild.at(i_new));
          delete BBchild.at(i_new);

          oJC.at(i)->Combine(*oJC.at(i_new));
          oTJ.at(i)->Combine(*oTJ.at(i_new));
          oChr.at(i)->Combine(*oChr.at(i_new));
          oSP.at(i)->Combine(*oSP.at(i_new));
          oROI.at(i)->Combine(*oROI.at(i_new));
          oCB.at(i)->Combine(*oCB.at(i_new));
          oFM.at(i)->Combine(*oFM.at(i_new));
          
          delete oJC.at(i_new);
          delete oTJ.at(i_new);
          delete oChr.at(i_new);
          delete oSP.at(i_new);
          delete oROI.at(i_new);
          delete oCB.at(i_new);
          delete oFM.at(i_new);          
        }
      }
    }
  }

  // Write Coverage Binary file:
  std::ofstream ofCOV;
  ofCOV.open(s_output_cov, std::ofstream::binary);
  covWriter outCOV;
  outCOV.SetOutputHandle(&ofCOV);
  oFM.at(0)->WriteBinary(&outCOV, verbose, Engine.n_threads_to_use);
  ofCOV.close();

// Write output to file:  
  if(verbose) cout << "Writing output file\n";

  std::ofstream out;                            
  out.open(s_output_txt, std::ios::binary);  // Open binary file
  GZWriter outGZ;                               
  outGZ.SetOutputHandle(&out); // GZ compression

  int outret = outGZ.writeline("BAM_report\tValue"); 
  if(outret != Z_OK) {
    cout << "Error writing gzip-compressed output file\n";
    out.close();
    delete oJC.at(0);
    delete oTJ.at(0);
    delete oChr.at(0);
    delete oSP.at(0);
    delete oROI.at(0);
    delete oCB.at(0);
    delete oFM.at(0);
    delete BBchild.at(0);
    return(-1);
  }

  // Output stuff here

// Write stats here:
  std::string myLine;
  BBchild.at(0)->WriteOutput(myLine);
  outGZ.writestring(myLine); outGZ.writeline("");

  int directionality = oJC.at(0)->Directional(myLine);
  outGZ.writeline("Directionality\tValue"); 
  outGZ.writestring(myLine); outGZ.writeline("");

  // Generate output but save this to strings:
  std::string myLine_ROI;
  std::string myLine_JC;
  std::string myLine_TJ;
  std::string myLine_SP;
  std::string myLine_Chr;
  std::string myLine_ND;
  std::string myLine_Dir;
  std::string myLine_QC;
  
  oROI.at(0)->WriteOutput(myLine_ROI, myLine_QC);
  oJC.at(0)->WriteOutput(myLine_JC, myLine_QC);
  oTJ.at(0)->WriteOutput(myLine_TJ, myLine_QC);
  oSP.at(0)->WriteOutput(myLine_SP, myLine_QC);
  oChr.at(0)->WriteOutput(myLine_Chr, myLine_QC);
  oCB.at(0)->WriteOutput(myLine_ND, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), Engine.n_threads_to_use);
  if (directionality != 0) {
    oCB.at(0)->WriteOutput(myLine_Dir, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), Engine.n_threads_to_use, directionality); // Directional.
	}

  outGZ.writeline("QC\tValue"); outGZ.writestring(myLine_QC); outGZ.writeline("");
	
  outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
  outGZ.writestring(myLine_ROI); outGZ.writeline("");
  
  outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
  outGZ.writestring(myLine_JC); outGZ.writeline("");

  outGZ.writeline("TJ_seqname\tstart1\tend1\tstart2\tend2\tstrand\ttotal\tpos\tneg");
  outGZ.writestring(myLine_TJ); outGZ.writeline("");
  
  outGZ.writeline("SP_seqname\tcoord\ttotal\tpos\tneg");
  outGZ.writestring(myLine_SP); outGZ.writeline("");
  
  outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
  outGZ.writestring(myLine_Chr); outGZ.writeline("");
  
  outGZ.writestring(myLine_ND); outGZ.writeline("");
  
  if (directionality != 0) {
    outGZ.writestring(myLine_Dir); outGZ.writeline("");
  }
  outGZ.flush(true);
  out.flush(); out.close();
  
  // destroy objects:

  delete oJC.at(0);
  delete oTJ.at(0);
  delete oChr.at(0);
  delete oSP.at(0);
  delete oROI.at(0);
  delete oCB.at(0);
  delete oFM.at(0);
  delete BBchild.at(0);

  return(0);
}

// SpliceWiz core:
int SpliceWizCore_htsMulti(
    swEngine_hts Engine,
    std::vector<std::string> const &bam_file, 
    std::vector<std::string> const &s_output_txt, 
    std::vector<std::string> const &s_output_cov,
    bool const verbose,
    int const read_pool
) {

  // check all BAM files exist
  for(unsigned int i = 0; i < bam_file.size(); i++) {
    if(!Engine.checkFileExists(bam_file.at(i))) {
      cout << "File " << bam_file.at(i) << " does not exist!\n";
      return(-1);
    }
  }

  // Initialize results container
  cout << "Allocating thread memory to " << Engine.n_threads_to_use 
    << " threads for SpliceWiz (htslib)...";
  auto start = chrono::steady_clock::now();
  auto check = start;
    Engine.loadReference();
  check = chrono::steady_clock::now();
  auto time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
  cout << "initialized (" << time_sec << " milliseconds)\n"; 

 
  // swResultsContainer res(n_threads_to_use);
  for(unsigned int z = 0; z < bam_file.size(); z++) {
    Engine.refreshReference();
    
    if(!Engine.checkFileExists(bam_file.at(z))) {
      cout << "File " << bam_file.at(z) << " does not exist!\n";
      continue;
    }
    if(verbose) cout << "Processing BAM file " << bam_file.at(z) << "\n";
    start = chrono::steady_clock::now();
  
    BGZF *fp = bgzf_open(bam_file.at(z).c_str(), "r");
    bam_hdr_t *header = bam_hdr_read(fp);

    hts_tpool *pool;
    const int queue_size = 0;
    if (Engine.n_threads_to_use > 1) {
        pool = hts_tpool_init(Engine.n_threads_to_use);
        bgzf_thread_pool(fp, pool, queue_size);
    }

    // Abort here if BAM corrupt
    if(header->n_targets <= 0){
      cout << bam_file.at(z) << " - contains no chromosomes mapped\n";
      return(-1);
    }
    std::vector<std::string> s_chr_names;
    std::vector<uint32_t> u32_chr_lens;
    for (int i = 0; i < header->n_targets; ++i) {
      s_chr_names.push_back(header->target_name[i]);
      u32_chr_lens.push_back(header->target_len[i]);
    }
  
    // Compile here a list of chromosomes; use BAM chromosomes for order
    // Add reference-only chromosomes at the end
    std::vector<std::string> bam_chr_name;
    std::vector<uint32_t> bam_chr_len;
    for(unsigned int i = 0; i < s_chr_names.size(); i++) {
      for(unsigned int j = 0; j < Engine.ref_alias.size(); j++) {
        if( 0==strncmp(
              Engine.ref_alias.at(j).c_str(), 
              s_chr_names.at(i).c_str(), 
              s_chr_names.at(i).size()
            ) && s_chr_names.at(i).size() == Engine.ref_alias.at(j).size()
        ) {
          bam_chr_name.push_back(Engine.ref_names.at(j));
          bam_chr_len.push_back(u32_chr_lens.at(i));
          break;
        }
      }
      if(i == bam_chr_name.size()) {
        bam_chr_name.push_back(s_chr_names.at(i));
        bam_chr_len.push_back(u32_chr_lens.at(i));
      }
    }
    // Now fill in reference chromosomes not in BAM:
    for(unsigned int i = 0; i < Engine.ref_names.size(); i++) {
      auto it = std::find(bam_chr_name.begin(), bam_chr_name.end(), Engine.ref_names.at(i));
      if(it == bam_chr_name.end()) {
        bam_chr_name.push_back(Engine.ref_names.at(i));
        bam_chr_len.push_back(Engine.ref_lengths.at(i));      
      }
    }

    Engine.associateBAM(bam_chr_name, bam_chr_len);
  
    // Read pool size
    unsigned int pool_cap = (unsigned int)read_pool;
    
    // Initialize bam1_t vector
    std::vector<bam1_t *> bpool;
    for(unsigned int i = 0; i < pool_cap; i++) {
      bpool.push_back(bam_init1());
    }
    
    // Pre-partition for n threads
    std::vector<int> pool_starts;
    std::vector<int> pool_ends;
    int est_tp_size = 1 + (pool_cap / Engine.n_threads_to_use);
    int poolPos = 0;
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      if(poolPos + est_tp_size > (int)pool_cap) {
        pool_starts.push_back(poolPos);
        pool_ends.push_back(pool_cap - 1);
        poolPos = pool_cap;
      } else {
        pool_starts.push_back(poolPos);
        pool_ends.push_back(poolPos + est_tp_size - 1);
        poolPos += est_tp_size;
      }
    }
    for(unsigned int i = pool_starts.size(); i < Engine.n_threads_to_use; i++) {
      pool_starts.push_back(-1);
      pool_ends.push_back(-1);
    }
  
    // BAM processing loop
    bool error_detected = false;
    off_t prevPos = 0; 
    off_t curPos = 0;
#ifdef SPLICEWIZ
    Progress p(GetFileSize(bam_file.at(z)), verbose);
    while(!p.check_abort()) {
      curPos = htell(fp->fp);
      p.increment(curPos - prevPos);
      prevPos = curPos;
    
#else
    while(!p.check_abort()) {
#endif
      
      // Load n reads here and partition by thread
      unsigned int pool_size = 0;
      for(unsigned int i = 0; i < bpool.size(); i++) {
        int ret = bam_read1(fp, bpool.at(i));
        if(ret < 0) {
          break;
        } else {
          pool_size++;
        }
      }

      // End of file
      if(pool_size == 0) {
        break;
      }

      // Multi-threaded process reads
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
      #endif
      for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
        if(pool_starts.at(i) >= 0 && pool_starts.at(i) < (int)pool_size) {
          int true_end = pool_size;
          if(true_end > pool_ends.at(i)) {
            true_end = pool_ends.at(i) + 1; // [first, last)
          }

          int pa_ret = Engine.BBchild.at(i)->processAll(
            bpool, pool_starts.at(i), true_end
          );
          if(pa_ret == -1) {
            
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            error_detected = true;
          }
        }
      }
      
      if(error_detected) break;

      // combine unpaired reads after each fillReads / processAlls
      for(unsigned int i = 1; i < Engine.n_threads_to_use; i++) {
        Engine.BBchild.at(0)->processSpares(*Engine.BBchild.at(i));
      }
    }

    for(unsigned int i = 0; i < bpool.size(); i++) {
      bam_destroy1(bpool.at(i));
    }
    bam_hdr_destroy(header);
    bgzf_close(fp);

  #ifdef SPLICEWIZ
    if(p.check_abort() || error_detected) {
      // interrupted:
  #else
    if(error_detected) {
  #endif
      // no need to delete results object, it will go out of scope
      if(error_detected) {
        // return(-1);
        cout << "Error encountered processing " << bam_file.at(z) << "\n";
        continue;
      }
    // Process aborted; stop processBAM for all requests
      cout << "Process interrupted running SpliceWiz on " 
        << bam_file.at(z) << '\n';
      return(-1);
    }

  if(Engine.n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(Engine.n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < Engine.n_threads_to_use) {
          Engine.BBchild.at(i)->processSpares(*Engine.BBchild.at(i_new));
          Engine.BBchild.at(i)->processStats(*Engine.BBchild.at(i_new));

          Engine.oJC.at(i)->Combine(*Engine.oJC.at(i_new));
          Engine.oTJ.at(i)->Combine(*Engine.oTJ.at(i_new));
          Engine.oChr.at(i)->Combine(*Engine.oChr.at(i_new));
          Engine.oSP.at(i)->Combine(*Engine.oSP.at(i_new));
          Engine.oROI.at(i)->Combine(*Engine.oROI.at(i_new));
          Engine.oCB.at(i)->Combine(*Engine.oCB.at(i_new));
          Engine.oFM.at(i)->Combine(*Engine.oFM.at(i_new));       
        }
      }
    }
  }

  // Write Coverage Binary file:
    std::ofstream ofCOV;
    ofCOV.open(s_output_cov.at(z), std::ofstream::binary);
    covWriter outCOV;
    outCOV.SetOutputHandle(&ofCOV);
    Engine.oFM.at(0)->WriteBinary(&outCOV, verbose, Engine.n_threads_to_use);
    ofCOV.close();

  // Write output to file:  
    if(verbose) cout << "Writing output file\n";

    std::ofstream out;                            
    out.open(s_output_txt.at(z), std::ios::binary);  // Open binary file
    GZWriter outGZ;                               
    outGZ.SetOutputHandle(&out); // GZ compression

    int outret = outGZ.writeline("BAM_report\tValue"); 
    if(outret != Z_OK) {
      cout << "Error writing gzip-compressed output file\n";
      out.close();
      continue;
      // return(-1);
    }

  // Output stuff here

  // Write stats here:
    std::string myLine;
    Engine.BBchild.at(0)->WriteOutput(myLine);
    outGZ.writestring(myLine); outGZ.writeline("");

    int directionality = Engine.oJC.at(0)->Directional(myLine);
    outGZ.writeline("Directionality\tValue"); 
    outGZ.writestring(myLine); outGZ.writeline("");

    // Generate output but save this to strings:
    std::string myLine_ROI;
    std::string myLine_JC;
    std::string myLine_TJ;
    std::string myLine_SP;
    std::string myLine_Chr;
    std::string myLine_ND;
    std::string myLine_Dir;
    std::string myLine_QC;
    
    Engine.oROI.at(0)->WriteOutput(myLine_ROI, myLine_QC);
    Engine.oJC.at(0)->WriteOutput(myLine_JC, myLine_QC);
    Engine.oTJ.at(0)->WriteOutput(myLine_TJ, myLine_QC);
    Engine.oSP.at(0)->WriteOutput(myLine_SP, myLine_QC);
    Engine.oChr.at(0)->WriteOutput(myLine_Chr, myLine_QC);
    Engine.oCB.at(0)->WriteOutput(myLine_ND, myLine_QC, 
      *Engine.oJC.at(0), *Engine.oSP.at(0), *Engine.oFM.at(0), 
      Engine.n_threads_to_use);
    if (directionality != 0) {
      Engine.oCB.at(0)->WriteOutput(myLine_Dir, myLine_QC, 
        *Engine.oJC.at(0), *Engine.oSP.at(0), *Engine.oFM.at(0), 
        Engine.n_threads_to_use, directionality); // Directional.
    }

    outGZ.writeline("QC\tValue"); outGZ.writestring(myLine_QC); outGZ.writeline("");
    
    outGZ.writeline("ROIname\ttotal_hits\tpositive_strand_hits\tnegative_strand_hits");
    outGZ.writestring(myLine_ROI); outGZ.writeline("");
    
    outGZ.writeline("JC_seqname\tstart\tend\tstrand\ttotal\tpos\tneg");
    outGZ.writestring(myLine_JC); outGZ.writeline("");

    outGZ.writeline("TJ_seqname\tstart1\tend1\tstart2\tend2\tstrand\ttotal\tpos\tneg");
    outGZ.writestring(myLine_TJ); outGZ.writeline("");
    
    outGZ.writeline("SP_seqname\tcoord\ttotal\tpos\tneg");
    outGZ.writestring(myLine_SP); outGZ.writeline("");
    
    outGZ.writeline("ChrCoverage_seqname\ttotal\tpos\tneg");
    outGZ.writestring(myLine_Chr); outGZ.writeline("");
    
    outGZ.writestring(myLine_ND); outGZ.writeline("");
    
    if (directionality != 0) {
      outGZ.writestring(myLine_Dir); outGZ.writeline("");
    }
    outGZ.flush(true);
    out.flush(); out.close();
  
    check = chrono::steady_clock::now();
    time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
    cout << bam_file.at(z) << " processed (" << time_sec << " milliseconds)\n";
  }

  return(0);
}

// BAM2COV htslib core:
int BAM2COVCore_hts(
    swEngine_hts Engine,
    std::string const &bam_file, 
    std::string const &s_output_cov,
    bool const verbose,
    int const read_pool
) {
  if(!Engine.checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 

  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);
  
  hts_tpool *pool;
  const int queue_size = 0;
  if (Engine.n_threads_to_use > 1) {
      pool = hts_tpool_init(Engine.n_threads_to_use);
      bgzf_thread_pool(fp, pool, queue_size);
  }
  
  // Abort here if BAM corrupt
  if(header->n_targets <= 0){
    cout << bam_file << " - contains no chromosomes mapped\n";
    return(-1);
  }
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  for (int i = 0; i < header->n_targets; ++i) {
    s_chr_names.push_back(header->target_name[i]);
    u32_chr_lens.push_back(header->target_len[i]);
  }

  std::vector<FragmentsMap*> oFM(Engine.n_threads_to_use);
  std::vector<htsBAM2blocks*> BBchild(Engine.n_threads_to_use);

  // Multi-threaded results container initialization
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
  #endif
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    oFM.at(i) = new FragmentsMap;
    BBchild.at(i) = new htsBAM2blocks(s_chr_names, u32_chr_lens);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i)->initializeChrs();
  }
  
  // Read pool size
  unsigned int pool_cap = (unsigned int)read_pool;
  
  // Initialize bam1_t vector
  std::vector<bam1_t *> bpool;
  for(unsigned int i = 0; i < pool_cap; i++) {
    bpool.push_back(bam_init1());
  }
  
  // Pre-partition for n threads
  std::vector<int> pool_starts;
  std::vector<int> pool_ends;
  int est_tp_size = 1 + (pool_cap / Engine.n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    if(poolPos + est_tp_size > (int)pool_cap) {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(pool_cap - 1);
      poolPos = pool_cap;
    } else {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(poolPos + est_tp_size - 1);
      poolPos += est_tp_size;
    }
  }
  for(unsigned int i = pool_starts.size(); i < Engine.n_threads_to_use; i++) {
    pool_starts.push_back(-1);
    pool_ends.push_back(-1);
  }
  
  // BAM processing loop
  bool error_detected = false;
  off_t prevPos = 0; 
  off_t curPos = 0;
#ifdef SPLICEWIZ
  Progress p(GetFileSize(bam_file), verbose);
  while(!p.check_abort()) {
    curPos = htell(fp->fp);
    p.increment(curPos - prevPos);
    prevPos = curPos;
    
#else
  while(!p.check_abort()) {
#endif
    
    // Load n reads here and partition by thread
    unsigned int pool_size = 0;
    for(unsigned int i = 0; i < bpool.size(); i++) {
      int ret = bam_read1(fp, bpool.at(i));
      if(ret < 0) {
        break;
      } else {
        pool_size++;
      }
    }

    // End of file
    if(pool_size == 0) {
      break;
    }

    // Multi-threaded process reads
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      if(pool_starts.at(i) >= 0 && pool_starts.at(i) < (int)pool_size) {
        int true_end = pool_size;
        if(true_end > pool_ends.at(i)) {
          true_end = pool_ends.at(i) + 1; // [first, last)
        }

        int pa_ret = BBchild.at(i)->processAll(
          bpool, pool_starts.at(i), true_end
        );
        if(pa_ret == -1) {
          
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          error_detected = true;
        }
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < Engine.n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
    }
  }

  for(unsigned int i = 0; i < bpool.size(); i++) {
    bam_destroy1(bpool.at(i));
  }
  bam_hdr_destroy(header);
  bgzf_close(fp);

#ifdef SPLICEWIZ
  if(p.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    if(error_detected) {
      return(-1);
    }
	// Process aborted; stop processBAM for all requests
    return(-2);
  }


  if(Engine.n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(Engine.n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < Engine.n_threads_to_use) {
          BBchild.at(i)->processSpares(*BBchild.at(i_new));
          BBchild.at(i)->processStats(*BBchild.at(i_new));
          delete BBchild.at(i_new);

          oFM.at(i)->Combine(*oFM.at(i_new));
          delete oFM.at(i_new);          
        }
      }
    }
  }

  // Write Coverage Binary file:
  std::ofstream ofCOV;
  ofCOV.open(s_output_cov, std::ofstream::binary);
  covWriter outCOV;
  outCOV.SetOutputHandle(&ofCOV);
  oFM.at(0)->WriteBinary(&outCOV, verbose, Engine.n_threads_to_use);
  ofCOV.close();

  delete oFM.at(0);
  delete BBchild.at(0);

  return(0);
}

// SpliceWiz core:
int doStatsCore_hts(
    swEngine_hts Engine,
    std::string const &bam_file, 
    std::string const &s_output_txt,
    bool const verbose,
    int const read_pool
) {
  if(!Engine.checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 

  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);
  
  hts_tpool *pool;
  const int queue_size = 0;
  if (Engine.n_threads_to_use > 1) {
      pool = hts_tpool_init(Engine.n_threads_to_use);
      bgzf_thread_pool(fp, pool, queue_size);
  }
  
  // Abort here if BAM corrupt
  if(header->n_targets <= 0){
    cout << bam_file << " - contains no chromosomes mapped\n";
    return(-1);
  }
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  for (int i = 0; i < header->n_targets; ++i) {
    s_chr_names.push_back(header->target_name[i]);
    u32_chr_lens.push_back(header->target_len[i]);
  }

  std::vector<htsBAM2blocks*> BBchild(Engine.n_threads_to_use);

  // Multi-threaded results container initialization
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
  #endif
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    BBchild.at(i) = new htsBAM2blocks(s_chr_names, u32_chr_lens);

    BBchild.at(i)->initializeChrs();
  }
  
  // Read pool size
  unsigned int pool_cap = (unsigned int)read_pool;
  
  // Initialize bam1_t vector
  std::vector<bam1_t *> bpool;
  for(unsigned int i = 0; i < pool_cap; i++) {
    bpool.push_back(bam_init1());
  }
  
  // Pre-partition for n threads
  std::vector<int> pool_starts;
  std::vector<int> pool_ends;
  int est_tp_size = 1 + (pool_cap / Engine.n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
    if(poolPos + est_tp_size > (int)pool_cap) {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(pool_cap - 1);
      poolPos = pool_cap;
    } else {
      pool_starts.push_back(poolPos);
      pool_ends.push_back(poolPos + est_tp_size - 1);
      poolPos += est_tp_size;
    }
  }
  for(unsigned int i = pool_starts.size(); i < Engine.n_threads_to_use; i++) {
    pool_starts.push_back(-1);
    pool_ends.push_back(-1);
  }
  
  // BAM processing loop
  bool error_detected = false;
  off_t prevPos = 0; 
  off_t curPos = 0;
#ifdef SPLICEWIZ
  Progress p(GetFileSize(bam_file), verbose);
  while(!p.check_abort()) {
    curPos = htell(fp->fp);
    p.increment(curPos - prevPos);
    prevPos = curPos;
    
#else
  while(!p.check_abort()) {
#endif
    
    // Load n reads here and partition by thread
    unsigned int pool_size = 0;
    for(unsigned int i = 0; i < bpool.size(); i++) {
      int ret = bam_read1(fp, bpool.at(i));
      if(ret < 0) {
        break;
      } else {
        pool_size++;
      }
    }

    // End of file
    if(pool_size == 0) {
      break;
    }

    // Multi-threaded process reads
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(Engine.n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      if(pool_starts.at(i) >= 0 && pool_starts.at(i) < (int)pool_size) {
        int true_end = pool_size;
        if(true_end > pool_ends.at(i)) {
          true_end = pool_ends.at(i) + 1; // [first, last)
        }

        int pa_ret = BBchild.at(i)->processAll(
          bpool, pool_starts.at(i), true_end
        );
        if(pa_ret == -1) {
          
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          error_detected = true;
        }
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < Engine.n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
    }
  }

  for(unsigned int i = 0; i < bpool.size(); i++) {
    bam_destroy1(bpool.at(i));
  }
  bam_hdr_destroy(header);
  bgzf_close(fp);

#ifdef SPLICEWIZ
  if(p.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < Engine.n_threads_to_use; i++) {
      delete BBchild.at(i);
    }
    if(error_detected) {
      return(-1);
    }
	// Process aborted; stop processBAM for all requests
    return(-2);
  }


  if(Engine.n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(Engine.n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < Engine.n_threads_to_use) {
          BBchild.at(i)->processSpares(*BBchild.at(i_new));
          BBchild.at(i)->processStats(*BBchild.at(i_new));
          delete BBchild.at(i_new);
        }
      }
    }
  }

// Write output to file:  
  if(verbose) cout << "Writing output file\n";

  std::ofstream out;                            
  out.open(s_output_txt, std::ios::binary);  // Open binary file
  GZWriter outGZ;                               
  outGZ.SetOutputHandle(&out); // GZ compression

  int outret = outGZ.writeline("BAM_report\tValue"); 
  if(outret != Z_OK) {
    cout << "Error writing gzip-compressed output file\n";
    out.close();
    delete BBchild.at(0);
    return(-1);
  }
  
// Write stats here:
  std::string myLine;
  BBchild.at(0)->WriteOutput(myLine);
  outGZ.writestring(myLine); outGZ.writeline("");

  outGZ.flush(true);
  out.flush(); out.close();
  
  delete BBchild.at(0);

  return(0);
}

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

  int ret2 = SpliceWizCore_hts(Engine,
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
    auto time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
    cout << s_bam << " processed (" << time_sec << " milliseconds)\n";
  }
  
  return(ret2);
}

// [[Rcpp::export]]
int SpliceWizMain_multi_hts(
    std::string reference_file, 
    StringVector bam_files, StringVector output_files,
    int max_threads, bool verbose, int read_pool
){
	if(bam_files.size() != output_files.size() || bam_files.size() < 1) {
		cout << "bam_files and output_files are of different sizes\n";
		return(1);	
	}
	std::vector< std::string > v_bam;
	std::vector< std::string > v_out_txt;
	std::vector< std::string > v_out_cov;
  for(int z = 0; z < bam_files.size(); z++) {
		v_bam.push_back(string(bam_files(z)));
		v_out_txt.push_back(string(output_files(z)) + ".txt.gz");
		v_out_cov.push_back(string(output_files(z)) + ".cov");
	}

  swEngine_hts Engine;
  Engine.Set_Threads(max_threads);

  std::string s_ref = reference_file;
  if(!Engine.checkFileExists(s_ref)) {
    cout << "File " << s_ref << " does not exist!\n";
    return(-1);
  } 

  if(verbose) cout << "Reading reference file\n";
  int ret = Engine.readReference(s_ref, verbose);
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid SpliceWiz reference\n";
    return(ret);
  }
  int ret2 = SpliceWizCore_htsMulti(
    Engine,
    v_bam, v_out_txt, v_out_cov,
    verbose, read_pool
  );

  return(ret2);
}


// [[Rcpp::export]]
int c_BAM2COV_hts(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, int read_pool
){
  swEngine_hts Engine;
  unsigned int n_threads_to_use = (unsigned int)Engine.Set_Threads(n_threads);
  std::string s_output_cov = output_file;
  std::string s_bam = bam_file;
  
  if(verbose) {
    cout << "Running BAM2COV (htslib) on " << s_bam;
    cout << " using " << n_threads_to_use << " threads\n";
  }
  
  // main:
  auto start = chrono::steady_clock::now();
  auto check = start;
  
  int ret2 = BAM2COVCore_hts(Engine,
      s_bam, s_output_cov,
      verbose, read_pool
  );

  if(ret2 == -2) {
    cout << "Process interrupted running BAM2COV on " << s_bam << '\n';
    return(ret2);
  } else if(ret2 == -1) {
    cout << "Error encountered processing " << s_bam << "\n";
  } else {
    check = chrono::steady_clock::now();
    auto time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
    cout << s_bam << " processed (" << time_sec << " milliseconds)\n";
  }
  
  return(ret2);
}

// [[Rcpp::export]]
int c_doStats_hts(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, int read_pool
){
  swEngine_hts Engine;
  unsigned int n_threads_to_use = (unsigned int)Engine.Set_Threads(n_threads);
  std::string s_output_txt = output_file;
  std::string s_bam = bam_file;
  
  if(verbose) {
    cout << "Running doStats (htslib) on " << s_bam;
    cout << " using " << n_threads_to_use << " threads\n";
  }
  
  // main:
  auto start = chrono::steady_clock::now();
  auto check = start;
  
  int ret2 = doStatsCore_hts(Engine,
      s_bam, s_output_txt,
      verbose, read_pool
  );

  if(ret2 == -2) {
    cout << "Process interrupted running BAM2COV on " << s_bam << '\n';
    return(ret2);
  } else if(ret2 == -1) {
    cout << "Error encountered processing " << s_bam << "\n";
  } else {
    check = chrono::steady_clock::now();
    auto time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
    cout << s_bam << " processed (" << time_sec << " milliseconds)\n";
  }
  
  return(ret2);
}