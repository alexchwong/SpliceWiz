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

uint64_t GetFileSize(std::string filename){
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

int swEngine_hts::loadReference() {
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

int swEngine_hts::loadReference(
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

int swEngine_hts::refreshReference() {
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

int swEngine_hts::associateBAM(
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
    
    BBchild.at(i).initializeChrs();
  }
  return(0);
}

// SWcore and B2C core performed insde hts_main.cpp