/* swEngine.cpp SpliceWiz processBAM engine

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

#include "swEngine.h"
#include <sys/stat.h>
#include <chrono>

swEngine::swEngine() {
  CB_string = "";
  SP_string = "";
  ROI_string = "";
  JC_string = "";
  TJ_string = "";
  n_threads_to_use = 1;
  
  refLoaded = false;
  BAMLoaded = false;

  CBloaded = false;
  SPloaded = false;
  ROIloaded = false;
  Chrloaded = false;
  JCloaded = false;
  TJloaded = false;
  FMloaded = false;
}

bool swEngine::checkFileExists(const std::string& name) {
    std::ifstream f;
    f.open(name);
    if(f){
      // cout << "File " << name << " exists\n";
      return(true);
    }
    // cout << "File " << name << " doesn't exist\n";
    return(false);
}

int swEngine::Set_Threads(int n_threads) {
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

int swEngine::ReadChrAlias(std::istringstream &IN) {
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
int swEngine::readReference(std::string &reference_file, bool verbose) { 
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

int swEngine::loadReference() {
  if(!refLoaded) {
    oCB.resize(n_threads_to_use);
    oSP.resize(n_threads_to_use);
    oROI.resize(n_threads_to_use);
    oChr.resize(n_threads_to_use);
    oJC.resize(n_threads_to_use);
    oTJ.resize(n_threads_to_use);
    oFM.resize(n_threads_to_use);
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      oCB.at(i).initialize(CB_string);
      oSP.at(i).initialize(SP_string);
      oROI.at(i).initialize(ROI_string);
      oJC.at(i).initialize(JC_string);
      oTJ.at(i).initialize(TJ_string);
    }
    
    refLoaded = true;
    CBloaded = true;
    SPloaded = true;
    ROIloaded = true;
    Chrloaded = true;
    JCloaded = true;
    TJloaded = true;
    FMloaded = true;
  }
  return(0);
}

int swEngine::loadReference(
  bool loadCB,
  bool loadSP,
  bool loadROI,
  bool loadChr,
  bool loadJC,
  bool loadTJ,
  bool loadFM
) {
  if(!refLoaded) {
    if(loadCB) oCB.resize(n_threads_to_use);
    if(loadSP) oSP.resize(n_threads_to_use);
    if(loadROI) oROI.resize(n_threads_to_use);
    if(loadChr) oChr.resize(n_threads_to_use);
    if(loadJC) oJC.resize(n_threads_to_use);
    if(loadTJ) oTJ.resize(n_threads_to_use);
    if(loadFM) oFM.resize(n_threads_to_use);
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      if(loadCB) oCB.at(i).initialize(CB_string);
      if(loadSP) oSP.at(i).initialize(SP_string);
      if(loadROI) oROI.at(i).initialize(ROI_string);
      if(loadJC) oJC.at(i).initialize(JC_string);
      if(loadTJ) oTJ.at(i).initialize(TJ_string);
    }
    
    refLoaded = true;
    CBloaded = true;
    SPloaded = true;
    ROIloaded = true;
    Chrloaded = true;
    JCloaded = true;
    TJloaded = true;
    FMloaded = true;
  }
  return(0);
}

int swEngine::refreshReference() {
  if(refLoaded) {
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      oCB.at(i).Reset();
      oSP.at(i).Reset();
      oROI.at(i).Reset();
      oJC.at(i).Reset();
      oTJ.at(i).Reset();
      oChr.at(i).Reset();
      oFM.at(i).Reset();
    }
  }

  if(BAMLoaded) {
    BBchild.clear();
    BAMLoaded = false;
  }
  
  return(0);
}

int swEngine::associateBAM(
  pbam_in * _IN, 
  std::vector<string> chr_name,
  std::vector<uint32_t> chr_len
) {
  BBchild.clear();
  BBchild.resize(n_threads_to_use);
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    BBchild.at(i).initialize(chr_name, chr_len);

    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &(oJC.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &(oJC.at(i)), std::placeholders::_1) );

    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&TandemJunctions::ChrMapUpdate, &(oTJ.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&TandemJunctions::ProcessBlocks, &(oTJ.at(i)), std::placeholders::_1) );
    
    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&FragmentsInChr::ChrMapUpdate, &(oChr.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&FragmentsInChr::ProcessBlocks, &(oChr.at(i)), std::placeholders::_1) );
    
    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&SpansPoint::ChrMapUpdate, &(oSP.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&SpansPoint::ProcessBlocks, &(oSP.at(i)), std::placeholders::_1) );
        
    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&FragmentsInROI::ChrMapUpdate, &(oROI.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&FragmentsInROI::ProcessBlocks, &(oROI.at(i)), std::placeholders::_1) );
    
    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&CoverageBlocks::ChrMapUpdate, &(oCB.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&CoverageBlocks::ProcessBlocks, &(oCB.at(i)), std::placeholders::_1) );

    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(oFM.at(i)), std::placeholders::_1) );
    
    BBchild.at(i).openFile(_IN);
  }
  return(0);
}

// SpliceWiz core:
int swEngine::SpliceWizMultiCore(
    std::vector<std::string> &bam_file, 
    std::vector<std::string> &s_output_txt, 
    std::vector<std::string> &s_output_cov,
    bool const verbose,
    bool const multithreadedRead
) {

  // check all BAM files exist
  for(unsigned int i = 0; i < bam_file.size(); i++) {
    if(!checkFileExists(bam_file.at(i))) {
      cout << "File " << bam_file.at(i) << " does not exist!\n";
      return(-1);
    }
  }
  
  // Initialize results container
  cout << "Allocating memory to " << n_threads_to_use
    << " threads for SpliceWiz (ompBAM)...";
  auto start = chrono::steady_clock::now();
  auto check = start;
    loadReference();
  check = chrono::steady_clock::now();
  auto time_sec = chrono::duration_cast<chrono::milliseconds>(check - start).count();
  cout << "initialized (" << time_sec << " milliseconds)\n"; 

  for(unsigned int z = 0; z < bam_file.size(); z++) {
    refreshReference();
    
    if(!checkFileExists(bam_file.at(z))) {
      cout << "File " << bam_file.at(z) << " does not exist!\n";
      continue;
    }
    if(verbose) cout << "Processing BAM file " << bam_file.at(z) << "\n";
    start = chrono::steady_clock::now();
    
    pbam_in inbam((size_t)5e8, (size_t)1e9, 5, multithreadedRead);
    inbam.openFile(bam_file.at(z), n_threads_to_use);
    
    // Abort here if BAM corrupt
    std::vector<std::string> s_chr_names;
    std::vector<uint32_t> u32_chr_lens;
    int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
    if(chrcount < 1) {
      cout << bam_file.at(z) << " - contains no chromosomes mapped\n";
      continue;
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

    associateBAM(&inbam, bam_chr_name, bam_chr_len);

    // BAM processing loop
    bool error_detected = false;
  #ifdef SPLICEWIZ
    // for compatibility with 32-bit file sizes
    uint64_t filesize = inbam.GetFileSize();
    int divFactor = 1;
    if(filesize > 4294967295) divFactor = 1000000;

    Progress pbar(filesize/divFactor, verbose);
    while(0 == inbam.fillReads() && !pbar.check_abort()) {
      pbar.increment(inbam.IncProgress()/divFactor);
      

  #else
    while(0 == inbam.fillReads()) {
  #endif
      
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
      #endif
      for(unsigned int i = 0; i < n_threads_to_use; i++) {
        int pa_ret = BBchild.at(i).processAll(i);
        if(pa_ret == -1) {
          
          #ifdef _OPENMP
          #pragma omp critical
          #endif
          error_detected = true;
        }
      }
      
      if(error_detected) break;

      // combine unpaired reads after each fillReads / processAlls
      for(unsigned int i = 1; i < n_threads_to_use; i++) {
        BBchild.at(0).processSpares(BBchild.at(i));
      }
    }

  #ifdef SPLICEWIZ
    if(pbar.check_abort() || error_detected) {
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
            BBchild.at(i).processSpares(BBchild.at(i_new));
            BBchild.at(i).processStats(BBchild.at(i_new));

            oJC.at(i).Combine(oJC.at(i_new));
            oTJ.at(i).Combine(oTJ.at(i_new));
            oChr.at(i).Combine(oChr.at(i_new));
            oSP.at(i).Combine(oSP.at(i_new));
            oROI.at(i).Combine(oROI.at(i_new));
            oCB.at(i).Combine(oCB.at(i_new));
            oFM.at(i).Combine(oFM.at(i_new));         
          }
        }
      }
    }

    // Write Coverage Binary file:
    std::ofstream ofCOV;
    ofCOV.open(s_output_cov.at(z), std::ofstream::binary);
    covWriter outCOV;
    outCOV.SetOutputHandle(&ofCOV);
    oFM.at(0).WriteBinary(&outCOV, verbose, n_threads_to_use);
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
    BBchild.at(0).WriteOutput(myLine);
    outGZ.writestring(myLine); outGZ.writeline("");

    int directionality = oJC.at(0).Directional(myLine);
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
    
    oROI.at(0).WriteOutput(myLine_ROI, myLine_QC);
    oJC.at(0).WriteOutput(myLine_JC, myLine_QC);
    oTJ.at(0).WriteOutput(myLine_TJ, myLine_QC);
    oSP.at(0).WriteOutput(myLine_SP, myLine_QC);
    oChr.at(0).WriteOutput(myLine_Chr, myLine_QC);
    oCB.at(0).WriteOutput(
      myLine_ND, myLine_QC, 
      oJC.at(0), oSP.at(0), 
      oFM.at(0), n_threads_to_use
    );
    if (directionality != 0) {
      oCB.at(0).WriteOutput(myLine_Dir, myLine_QC, 
        oJC.at(0), oSP.at(0), oFM.at(0), 
        n_threads_to_use, directionality
      ); // Directional.
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

// BAM2COV core:
int swEngine::doStatsCore(
    std::string const &bam_file,
    std::string const &s_output_txt,
    bool const verbose,
    bool const multithreadedRead
) {
  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
	if(verbose) cout << "doStats (ompBAM): " << bam_file << "\n";
  
  pbam_in inbam((size_t)5e8, (size_t)1e9, 5, multithreadedRead);
  inbam.openFile(bam_file, n_threads_to_use);
  
  // Abort here if BAM corrupt
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  if(chrcount < 1) {
    cout << bam_file << " - contains no mapped chromosomes\n";
    return(-1);
  }

  BBchild.clear();
  BBchild.resize(n_threads_to_use);
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    BBchild.at(i).initialize(s_chr_names, u32_chr_lens);
    BBchild.at(i).openFile(&inbam);
  }
  
  // BAM processing loop
  bool error_detected = false;
#ifdef SPLICEWIZ
  // for compatibility with 32-bit file sizes
  uint64_t filesize = inbam.GetFileSize();
  int divFactor = 1;
  if(filesize > 4294967295) divFactor = 1000000;

  Progress pbar(filesize/divFactor, verbose);
  while(0 == inbam.fillReads() && !pbar.check_abort()) {
    pbar.increment(inbam.IncProgress()/divFactor);
    
#else
  while(0 == inbam.fillReads()) {
#endif
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      int pa_ret = BBchild.at(i).processAll(i);
      if(pa_ret == -1) {
        
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        error_detected = true;
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0).processSpares(BBchild.at(i));
    }
  }

#ifdef SPLICEWIZ
  if(pbar.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      // delete BBchild.at(i);
    }
    if(error_detected) {
      return(-1);
    }
	// Process aborted; stop processBAM for all requests
    return(-2);
  }

  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0).processSpares(BBchild.at(i));
      BBchild.at(0).processStats(BBchild.at(i));
      // delete BBchild.at(i);
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
    // delete BBchild.at(0);
    return(-1);
  }
  
// Write stats here:
  std::string myLine;
  BBchild.at(0).WriteOutput(myLine);
  outGZ.writestring(myLine); outGZ.writeline("");

  outGZ.flush(true);
  out.flush(); out.close();
  
  // delete BBchild.at(0);

  return(0);
}

// SpliceWiz Mappability:
int swEngine::MappabilityRegionsCore(
    std::string const &bam_file,
    std::string const &s_output_txt,
    std::string const &s_output_cov,
    int threshold,
    bool const includeCov,
    bool const verbose,
    bool const multithreadedRead
) {
  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
	if(verbose) cout 
    << "Calculating Mappability Exclusions: " 
    << bam_file << "\n";
  
  pbam_in inbam((size_t)5e8, (size_t)1e9, 5, multithreadedRead);
  inbam.openFile(bam_file, n_threads_to_use);
  
  // Abort here if BAM corrupt
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  if(chrcount < 1) {
    cout << bam_file << " - contains no mapped chromosomes\n";
    return(-1);
  }

  // std::vector<FragmentsMap*> oFM;
  // std::vector<BAM2blocks*> BBchild;

  BBchild.clear();
  BBchild.resize(n_threads_to_use);
  oFM.resize(n_threads_to_use);
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    // oFM.push_back(FragmentsMap());
    BBchild.at(i).initialize(s_chr_names, u32_chr_lens);

    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i).openFile(&inbam);
  }
  
  // BAM processing loop
  bool error_detected = false;
#ifdef SPLICEWIZ
  // for compatibility with 32-bit file sizes
  uint64_t filesize = inbam.GetFileSize();
  int divFactor = 1;
  if(filesize > 4294967295) divFactor = 1000000;

  Progress pbar(filesize/divFactor, verbose);
  while(0 == inbam.fillReads() && !pbar.check_abort()) {
    pbar.increment(inbam.IncProgress()/divFactor);
    

#else
  while(0 == inbam.fillReads()) {
#endif
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      int pa_ret = BBchild.at(i).processAll(i);
      if(pa_ret == -1) {
        
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        error_detected = true;
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0).processSpares(BBchild.at(i));
    }
  }

#ifdef SPLICEWIZ
  if(pbar.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      // delete oFM.at(i);
      // delete BBchild.at(i);
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
          BBchild.at(i).processSpares(BBchild.at(i_new));
          BBchild.at(i).processStats(BBchild.at(i_new));
          // delete BBchild.at(i_new);

          oFM.at(i).Combine(oFM.at(i_new));
          // delete oFM.at(i_new);          
        }
      }
    }
  }

  if(includeCov) {
    // Write Coverage Binary file:
    std::ofstream ofCOV;
    ofCOV.open(s_output_cov, std::ofstream::binary);
    covWriter outCOV;
    outCOV.SetOutputHandle(&ofCOV);
    oFM.at(0).WriteBinary(&outCOV, verbose, n_threads_to_use);
    ofCOV.close();    
  }

  std::ofstream outFragsMap;
  outFragsMap.open(s_output_txt, std::ifstream::out);
	
  oFM.at(0).WriteOutput(&outFragsMap, threshold, verbose, n_threads_to_use);
  outFragsMap.flush(); outFragsMap.close();

  // delete oFM.at(0);
  // delete BBchild.at(0);

  return(0);
}

// SpliceWiz BAM2COV:
int swEngine::BAM2COVcore(
    std::string const &bam_file,
    std::string const &s_output_cov,
    bool const verbose,
    bool const multithreadedRead
) {
  if(!checkFileExists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
	if(verbose) cout 
    << "SpliceWiz BAM2COV: " 
    << bam_file << "\n";
  
  pbam_in inbam((size_t)5e8, (size_t)1e9, 5, multithreadedRead);
  inbam.openFile(bam_file, n_threads_to_use);
  
  // Abort here if BAM corrupt
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  if(chrcount < 1) {
    cout << bam_file << " - contains no mapped chromosomes\n";
    return(-1);
  }

  BBchild.clear();
  BBchild.resize(n_threads_to_use);
  oFM.resize(n_threads_to_use);
  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    // oFM.push_back(FragmentsMap());
    BBchild.at(i).initialize(s_chr_names, u32_chr_lens);

    BBchild.at(i).registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i).registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i).openFile(&inbam);
  }
  
  // BAM processing loop
  bool error_detected = false;
#ifdef SPLICEWIZ
  // for compatibility with 32-bit file sizes
  uint64_t filesize = inbam.GetFileSize();
  int divFactor = 1;
  if(filesize > 4294967295) divFactor = 1000000;

  Progress pbar(filesize/divFactor, verbose);
  while(0 == inbam.fillReads() && !pbar.check_abort()) {
    pbar.increment(inbam.IncProgress()/divFactor);
    

#else
  while(0 == inbam.fillReads()) {
#endif
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      int pa_ret = BBchild.at(i).processAll(i);
      if(pa_ret == -1) {
        
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        error_detected = true;
      }
    }
    
    if(error_detected) break;

    // combine unpaired reads after each fillReads / processAlls
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0).processSpares(BBchild.at(i));
    }
  }

#ifdef SPLICEWIZ
  if(pbar.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      // delete oFM.at(i);
      // delete BBchild.at(i);
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
          BBchild.at(i).processSpares(BBchild.at(i_new));
          BBchild.at(i).processStats(BBchild.at(i_new));

          oFM.at(i).Combine(oFM.at(i_new));         
        }
      }
    }
  }


  // Write Coverage Binary file:
  std::ofstream ofCOV;
  ofCOV.open(s_output_cov, std::ofstream::binary);
  covWriter outCOV;
  outCOV.SetOutputHandle(&ofCOV);
  oFM.at(0).WriteBinary(&outCOV, verbose, n_threads_to_use);
  ofCOV.close();    

  return(0);
}