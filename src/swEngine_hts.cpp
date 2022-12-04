/* swEngine_hts.cpp SpliceWiz processBAM engine (htslib)

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

#include "swEngine_hts.h"

long GetFileSize(std::string filename){
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

swEngine_hts::swEngine_hts() {
  CB_string = "";
  SP_string = "";
  ROI_string = "";
  JC_string = "";
  TJ_string = "";
  n_threads_to_use = 1;
}

bool swEngine_hts::checkFileExists(const std::string& name) {
    std::ifstream f;
    f.open(name);
    if(f){
      // cout << "File " << name << " exists\n";
      return(true);
    }
    // cout << "File " << name << " doesn't exist\n";
    return(false);
}

int swEngine_hts::Set_Threads(int n_threads) {
#ifdef _OPENMP
  int use_threads = 1;
  if(n_threads > 1 && n_threads < omp_get_thread_limit()) {
    use_threads = n_threads;
	} else if(n_threads >= omp_get_thread_limit()) {
		use_threads = omp_get_thread_limit();
		if(use_threads < 1) {
			use_threads = 1;
		}
	}
	// omp_set_num_threads(use_threads);
  n_threads_to_use = use_threads;
  return(use_threads);
#else
  n_threads_to_use = 1;
	return(1);
#endif
}

int swEngine_hts::ReadChrAlias(std::istringstream &IN) {
  ref_names.clear();
  ref_alias.clear();
  ref_lengths.clear();
  
  std::string myLine;
  myLine.reserve(1000);
  std::string myChr;
  myChr.reserve(100);
  std::string myAlias;
  myAlias.reserve(100);
  std::string myLength;
  myLength.reserve(100);
  
  while(!IN.eof() && !IN.fail()) {
    getline(IN, myLine, '\n');
    if (IN.eof() || IN.fail()) {
      if (myLine.length() == 0) {
        // This line is empty - just a blank line at the end of the file.
        // Checking at this stage allows correct handling of files both with and without a trailing \n after the last record.
        break;
      }else{
        // Error line in input, ignore.
        break;
      }
    }
    std::istringstream lineStream;
    lineStream.str(myLine);
    getline(lineStream, myChr, '\t');
    getline(lineStream, myLength, '\t');
    getline(lineStream, myAlias, '\t');
    if(myChr.size() > 0) {
      ref_names.push_back(myChr);
      ref_lengths.push_back((uint32_t)stoul(myLength));
      ref_alias.push_back(myAlias);      
    }
  }
  // cout << "Debug:" << ref_names.size() << " chromosome aliases loaded\n";
  return(0);
}

// SpliceWiz reference reader
int swEngine_hts::readReference(std::string &reference_file, bool verbose) { 
  if(!checkFileExists(reference_file)) {
    cout << "File " << reference_file << " does not exist!\n";
    return(-1);
  }

  GZReader * gz_in = new GZReader;
  int ret = gz_in->LoadGZ(reference_file, true);   // streamed mode
  if(ret != 0) return(-1);
  
  // Allows reference blocks to be read in any order
  std::string headerCover ("ref-cover.bed");
  std::string headerSpans ("ref-read-continues.ref");
  std::string headerROI ("ref-ROI.bed");
  std::string headerSJ ("ref-sj.ref");
  std::string headerTJ ("ref-tj.ref");
  std::string headerChr ("ref-chrs.ref");
  std::string headerEOF ("EOF");
  
  bool doneCover = false;
  bool doneSpans = false;
  bool doneROI = false;
  bool doneSJ = false;
  bool doneTJ = false;
  bool doneChrs = false;
  
  std::string myLine;
  std::string myBuffer;
  
  getline(gz_in->iss, myLine, '#');    // discard anything before the first hash
  getline(gz_in->iss, myLine, '\n');   // Get block name
  
  // Check non-empty ref block name
  if(myLine.size() == 0) {
    cout << "Invalid SpliceWiz reference detected\n";
    return(-1);
  }

  while(myLine.find(headerEOF)==std::string::npos) {
    // getline(gz_in->iss, myBuffer, '#');  // this is the data block

    if(myLine.find(headerCover)!=std::string::npos && !doneCover) {
      getline(gz_in->iss, CB_string, '#');
      doneCover = true;
    } else if(myLine.find(headerSpans)!=std::string::npos && !doneSpans) {
      getline(gz_in->iss, SP_string, '#');
      doneSpans = true;
    } else if(myLine.find(headerROI)!=std::string::npos && !doneROI) {
      getline(gz_in->iss, ROI_string, '#');
      doneROI = true;
    } else if(myLine.find(headerSJ)!=std::string::npos && !doneSJ) {
      getline(gz_in->iss, JC_string, '#');
      doneSJ = true;
    } else if(myLine.find(headerTJ)!=std::string::npos && !doneTJ) {
      getline(gz_in->iss, TJ_string, '#');
      doneTJ = true;
    } else if(myLine.find(headerChr)!=std::string::npos && !doneChrs) {
      getline(gz_in->iss, myBuffer, '#');
      std::istringstream inChrAlias;
      inChrAlias.str(myBuffer);
      ReadChrAlias(inChrAlias);
      doneChrs = true;
    } else {
      cout << "Error: Invalid SpliceWiz reference block detected\n";
      return(-1);
    }
    // Get next data block name
    getline(gz_in->iss, myLine, '\n');
  }

  delete gz_in;
  
  if(!doneCover || !doneSpans || !doneROI || !doneSJ) {
    cout << "Error: Incomplete SpliceWiz reference detected\n";
    return(-1);
  } else if(!doneTJ) {
    cout << "Note: Tandem junction reference not detected. " <<
      "Rebuild reference using SpliceWiz v0.99.3 or above.\n";
  }
  return(0);
}

// SpliceWiz core:
int swEngine_hts::SpliceWizCore(
    std::string const &bam_file, 
    std::string const &s_output_txt, 
    std::string const &s_output_cov,
    bool const verbose,
    int const read_pool
) {

  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
 
	if(verbose) cout << "Processing BAM file " << bam_file << "\n";
  
  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);

  hts_tpool *pool;
  const int queue_size = 64;
  if (n_threads_to_use > 1) {
      pool = hts_tpool_init(n_threads_to_use);
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
    for(unsigned int j = 0; j < ref_alias.size(); j++) {
      if( 0==strncmp(
            ref_alias.at(j).c_str(), 
            s_chr_names.at(i).c_str(), 
            s_chr_names.at(i).size()
          ) && s_chr_names.at(i).size() == ref_alias.at(j).size()
      ) {
        bam_chr_name.push_back(ref_names.at(j));
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
  for(unsigned int i = 0; i < ref_names.size(); i++) {
    auto it = std::find(bam_chr_name.begin(), bam_chr_name.end(), ref_names.at(i));
    if(it == bam_chr_name.end()) {
      bam_chr_name.push_back(ref_names.at(i));
      bam_chr_len.push_back(ref_lengths.at(i));      
    }
  }
  
  std::vector<CoverageBlocksIRFinder*> oCB(n_threads_to_use);
  std::vector<SpansPoint*> oSP(n_threads_to_use);
  std::vector<FragmentsInROI*> oROI(n_threads_to_use);
  std::vector<FragmentsInChr*> oChr(n_threads_to_use);
  std::vector<JunctionCount*> oJC(n_threads_to_use);
  std::vector<TandemJunctions*> oTJ(n_threads_to_use);
  std::vector<FragmentsMap*> oFM(n_threads_to_use);
  std::vector<htsBAM2blocks*> BBchild(n_threads_to_use);

  // Multi-threaded results container initialization
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
  #endif
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oCB.at(i) = new CoverageBlocksIRFinder(CB_string);
    oSP.at(i) = new SpansPoint(SP_string);
    oROI.at(i) = new FragmentsInROI(ROI_string);
    oChr.at(i) = new FragmentsInChr;
    oJC.at(i) = new JunctionCount(JC_string);
    oTJ.at(i) = new TandemJunctions(TJ_string);
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
  int est_tp_size = 1 + (pool_cap / n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
  for(unsigned int i = pool_starts.size(); i < n_threads_to_use; i++) {
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
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
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
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
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


  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < n_threads_to_use) {
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
  oFM.at(0)->WriteBinary(&outCOV, verbose, n_threads_to_use);
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
  oCB.at(0)->WriteOutput(myLine_ND, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), n_threads_to_use);
  if (directionality != 0) {
    oCB.at(0)->WriteOutput(myLine_Dir, myLine_QC, *oJC.at(0), *oSP.at(0), *oFM.at(0), n_threads_to_use, directionality); // Directional.
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
int swEngine_hts::BAM2COVcore(
    std::string const &bam_file,
    std::string const &s_output_cov,
    bool const verbose,
    bool const read_pool
) {
  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
	if(verbose) cout << "BAM2COV: " << bam_file << "\n";
  
  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);
  
  hts_tpool *pool;
  const int queue_size = 64;
  if (n_threads_to_use > 1) {
      pool = hts_tpool_init(n_threads_to_use);
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

  std::vector<FragmentsMap*> oFM(n_threads_to_use);
  std::vector<htsBAM2blocks*> BBchild(n_threads_to_use);

  // Multi-threaded results container initialization
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
  #endif
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
  int est_tp_size = 1 + (pool_cap / n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
  for(unsigned int i = pool_starts.size(); i < n_threads_to_use; i++) {
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
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
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
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    if(error_detected) {
      return(-1);
    }
	// Process aborted; stop processBAM for all requests
    return(-2);
  }


  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine objects (multi-threaded):
    int n_rounds = ceil(log(n_threads_to_use) / log(2));
    for(int j = n_rounds; j > 0; j--) {
      int n_bases = (int)pow(2, j-1);
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_bases) schedule(static,1)
      #endif
      for(int i = 0; i < n_bases; i++) {
        int i_new = i + n_bases;
        if((unsigned int)i_new < n_threads_to_use) {
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
  oFM.at(0)->WriteBinary(&outCOV, verbose, n_threads_to_use);
  ofCOV.close();

  delete oFM.at(0);
  delete BBchild.at(0);

  return(0);
}

// do nothing except decompress BAM reads:
int swEngine_hts::doNothing(
    std::string const &bam_file,
    bool const verbose,
    bool const read_pool
) {
  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
	if(verbose) cout << "doNothing: " << bam_file << "\n";
  
  BGZF *fp = bgzf_open(bam_file.c_str(), "r");
  bam_hdr_t *header = bam_hdr_read(fp);
  
  hts_tpool *pool;
  const int queue_size = 64;
  if (n_threads_to_use > 1) {
      pool = hts_tpool_init(n_threads_to_use);
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
  int est_tp_size = 1 + (pool_cap / n_threads_to_use);
  int poolPos = 0;
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
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
  for(unsigned int i = pool_starts.size(); i < n_threads_to_use; i++) {
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
  }

  for(unsigned int i = 0; i < bpool.size(); i++) {
    bam_destroy1(bpool.at(i));
  }
  bam_hdr_destroy(header);
  bgzf_close(fp);

  return(0);
}