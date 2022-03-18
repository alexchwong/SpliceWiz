/* main.cpp Exported functions

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

#include "main.h"

// [[Rcpp::export]]
int Has_OpenMP() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
#endif
}

int Set_Threads(int n_threads) {
#ifdef _OPENMP
  int use_threads = 1;
	if(n_threads > 0 && n_threads <= omp_get_thread_limit()) {
    use_threads = n_threads;
	} else {
		use_threads = omp_get_thread_limit();
		if(use_threads < 1) {
			use_threads = 1;
		}
	}
	omp_set_num_threads(use_threads);
  return(use_threads);
#else
	return(1);
#endif
}

// [[Rcpp::export]]
int Test_OpenMP_For() {
// Returns -1 if no OpenMP; 0 if OpenMP works; 1 if compiled with OpenMP and only 1 thread; and 2 if compiled with OpenMP but parallel loops don't work
  
#ifdef _OPENMP
  bool test_res = false;
  if(omp_get_thread_limit() == 1) {
    return 1;
  } else {
    bool std_false = (bool)omp_in_parallel();
    // Test 2 threads
    #pragma omp parallel for num_threads(2) schedule(static,1)
    for(unsigned int i = 0; i < 2; i++) {
      #pragma omp critical
      test_res = (bool)omp_in_parallel();   // returns true iff in parallel loop
    }
    if(test_res && !std_false) return 0;
  }
  return 2;
#else
  return -1;
#endif
}

inline bool see_if_file_exists(const std::string& name) {
    std::ifstream f;
    f.open(name);
    if(f){
      // cout << "File " << name << " exists\n";
      return(true);
    }
    // cout << "File " << name << " doesn't exist\n";
    return(false);
}

// [[Rcpp::export]]
bool c_Check_Cov(std::string s_in) {
	// Checks if given file is a valid COV file
	
  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);	

  covReader inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(false);
  }

  int ret = inCov.ReadHeader();
  if(ret == -1){
		inCov_stream.close();	
    return(false);
  }	
	
  inCov_stream.close();	
	return(true);
}

// ########################### MAPPABILITY INTERNAL FN #########################

std::string GenerateReadError(
    char * input_read, 
    const unsigned int read_len, 
    const unsigned int error_pos,
    const unsigned int direction, 
    const size_t error_seed
) {
  
  // Copy https://github.com/williamritchie/IRFinder/blob/master/bin/util/generateReadsError.pl

  char * new_read_inter = new char[read_len + 1];
  new_read_inter[read_len] = '\0';
  memcpy(&new_read_inter[0], input_read, read_len);  

  char error_nuc = '\0';  // set this as something to avoid warning at compile
  if(error_seed % 3 == 0) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'G'; break;
      case 'C':
        error_nuc = 'A'; break;
      case 'G':
        error_nuc = 'T'; break;
      case 'T':
        error_nuc = 'C'; break;
      case 'a':
        error_nuc = 'g'; break;
      case 'c':
        error_nuc = 'a'; break;
      case 'g':
        error_nuc = 't'; break;
      case 't':
        error_nuc = 'c'; break;
      default:
        error_nuc = 'N';
    }
  } else if(error_seed % 3 == 1) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'T'; break;
      case 'C':
        error_nuc = 'G'; break;
      case 'G':
        error_nuc = 'C'; break;
      case 'T':
        error_nuc = 'A'; break;
      case 'a':
        error_nuc = 't'; break;
      case 'c':
        error_nuc = 'g'; break;
      case 'g':
        error_nuc = 'c'; break;
      case 't':
        error_nuc = 'a'; break;
      default:
        error_nuc = 'N';
    }
  } else {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'C'; break;
      case 'C':
        error_nuc = 'T'; break;
      case 'G':
        error_nuc = 'A'; break;
      case 'T':
        error_nuc = 'G'; break;
      case 'a':
        error_nuc = 'c'; break;
      case 'c':
        error_nuc = 't'; break;
      case 'g':
        error_nuc = 'a'; break;
      case 't':
        error_nuc = 'g'; break;
      default:
        error_nuc = 'N';
    }
  }
  
  memcpy(&new_read_inter[error_pos - 1], &error_nuc, 1);
  
  char * new_read = new char[read_len + 1];
  new_read[read_len] = '\0';
  if(direction == 0) {
    memcpy(&new_read[0], new_read_inter, read_len);  
  } else {
    for(unsigned int i = 0; i < read_len; i++) {
      switch(new_read_inter[i]) {   
        case 'A':
          new_read[read_len - i - 1] = 'T'; break;
        case 'T':
          new_read[read_len - i - 1] = 'A'; break;
        case 'G':
          new_read[read_len - i - 1] = 'C'; break;
        case 'C':
          new_read[read_len - i - 1] = 'G'; break;
        case 'a':
          new_read[read_len - i - 1] = 't'; break;
        case 't':
          new_read[read_len - i - 1] = 'a'; break;
        case 'g':
          new_read[read_len - i - 1] = 'c'; break;
        case 'c':
          new_read[read_len - i - 1] = 'g'; break;
        default :
          new_read[read_len - i - 1] = 'N';
      }         
    }
  }

  string return_str = string(new_read);
  delete[] new_read;
  return(return_str);
}

// Replicate old PERL script; return true if N's constitute less than half of length
bool checkDNA(char * input_read, unsigned int read_len) {
  unsigned int numN = 0;
  for(unsigned int i = 0; i < read_len; i++) {
    if(input_read[i]!='A' && input_read[i]!='T' && input_read[i]!='G' && input_read[i]!='C' &&
      input_read[i]!='a' && input_read[i]!='t' && input_read[i]!='g' && input_read[i]!='c') {
      numN++;
    }
  }
  return(numN < read_len / 2);
}

// #############################################################################

#ifdef SPLICEWIZ
// Below are Rcpp-only functions

// [[Rcpp::export]]
List c_RLE_From_Cov(std::string s_in, std::string seqname, int start, int end, int strand) {
// Returns an RLE covering the region described above
// s_in: The coverage file
// strand: 0 = -, 1 = +, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  if(!see_if_file_exists(s_in)) {
    cout << "File " << s_in << " does not exist!\n";
    return(NULL_RLE);
  }
  
  if(start > end || start < 0){
    return(NULL_RLE);
  }

  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covReader inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(NULL_RLE);
  }
	
  int ret = inCov.ReadHeader();
  if(ret < 0){
		cout << s_in << " appears to not be valid COV file... exiting\n";
		inCov_stream.close();	
    return(NULL_RLE);
  }
  
  // Find corresponding seqname
  unsigned int ref_index = 0;
  std::vector<chr_entry> chrs;
  inCov.GetChrs(chrs);
  while(0 != seqname.compare(0, seqname.size(), chrs.at(ref_index).chr_name)) {
    if(ref_index == chrs.size()) break;
    ref_index++;
  }
  if(ref_index == chrs.size()) {
    inCov_stream.close();	
    return(NULL_RLE);
  }
  
  int eff_end = 0;
  if(end == 0) {
    eff_end = chrs.at(ref_index).chr_len;
  } else {
    eff_end = end;
  }
  
  std::vector<int> values;
  std::vector<unsigned int> lengths;
  // Push first value
  values.push_back(0);
  lengths.push_back((unsigned int)start);
  
  inCov.FetchRLE(seqname, (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);

  inCov_stream.close();
  // Push last value
  if((uint32_t)eff_end < (uint32_t)chrs.at(ref_index).chr_len) {
    values.push_back(0);
    lengths.push_back((uint32_t)chrs.at(ref_index).chr_len - eff_end);
  }
    
  List RLE = List::create(
    _["values"] = values,
    _["lengths"] = lengths 
  );

  return(RLE);
}

// [[Rcpp::export]]
StringVector c_Cov_Seqnames(
  std::string s_in
) {
  
  StringVector s_out;
  
  if(!see_if_file_exists(s_in)) {
    cout << "File " << s_in << " does not exist!\n";
    return(s_out);
  }

  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covReader inCov;
  inCov.SetInputHandle(&inCov_stream);


  if(inCov.fail()){
    cout << "File " << s_in << " reading failed!\n";
		inCov_stream.close();
    return(s_out);
  }
  
  int ret = inCov.ReadHeader();
  if(ret == -1){
		cout << s_in << " appears to not be valid COV file... exiting";
		inCov_stream.close();	
    return(s_out);
  }
  
  std::vector<chr_entry> chrs;
  inCov.GetChrs(chrs);

  for(unsigned int i = 0; i < chrs.size(); i++) {
    s_out.push_back(chrs.at(i).chr_name);
  }

  return(s_out);
}

// [[Rcpp::export]]
List c_RLEList_From_Cov(std::string s_in, int strand) {
  // Returns an RLEList
  // s_in: The coverage file
  // strand: 0 = -, 1 = +, 2 = *
  
  List NULL_RLE = List::create(
    _["values"] = 0,
    _["lengths"] = 0 
  );
  
  if(!see_if_file_exists(s_in)) {
    cout << "File " << s_in << " does not exist!\n";
    return(NULL_RLE);
  }

  List RLEList;
  
  std::ifstream inCov_stream;
  inCov_stream.open(s_in, std::ifstream::binary);
  
  covReader inCov;
  inCov.SetInputHandle(&inCov_stream);

  if(inCov.fail()){
		inCov_stream.close();
    return(NULL_RLE);
  }
  
  int ret = inCov.ReadHeader();
  if(ret == -1){
		cout << s_in << " appears to not be valid COV file... exiting";
		inCov_stream.close();	
    return(NULL_RLE);
  }
  
  std::vector<chr_entry> chrs;
  inCov.GetChrs(chrs);

  for (unsigned int i = 0; i < chrs.size(); i++) {
    uint32_t eff_end = chrs.at(i).chr_len;
    uint32_t start = 0;
    
    std::vector<int> values;
    std::vector<unsigned int> lengths;

    inCov.FetchRLE(chrs.at(i).chr_name, (uint32_t)start, (uint32_t)eff_end, strand, &values, &lengths);
    
    List RLE = List::create(
      _["values"] = values,
      _["lengths"] = lengths 
    );
    RLEList.push_back(RLE, chrs.at(i).chr_name);
  }

  inCov_stream.close();
  
  return(RLEList);
}

// [[Rcpp::export]]
List c_gunzip_DF(std::string s_in, StringVector s_header_begin) {
  List Final_final_list;

  if(!see_if_file_exists(s_in)) {
    cout << "File " << s_in << " does not exist!\n";
    return(Final_final_list);
  }
  
  GZReader gz_in;
  int ret = gz_in.LoadGZ(s_in, false, true);   // lazy mode
  if(ret != 0) return(Final_final_list);

  // Look for first line of data to return
  std::string myLine;
  std::string myEntry;
  unsigned int q = 0;
  char delim = '\n';
	bool check_line = true;
	
  for(int z = 0; z < s_header_begin.size(); z++) {
    std::string header = string(s_header_begin(z));
    std::vector<std::string> columns;
    while(!gz_in.eof()) {
      gz_in.getline(myLine, delim); q++;
      myLine.erase( std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end() ); // remove \r 
      
			if(check_line == true) {
				if(strncmp(myLine.c_str(), header.c_str(), header.size()) == 0) {
					// read columns
					std::istringstream column_iss;
					column_iss.str(myLine);
					
					while(!column_iss.eof() && !column_iss.fail()) {
						getline(column_iss, myEntry, '\t');
						columns.push_back(myEntry);
					}
					
					break;
				} else {
					check_line = false;
				}
			} else {
				// screen for empty lines - then reactivate check_line
				if (myLine.length() == 0) {
					check_line = true;
				}
			}

    }

    // use a map of string vectors
    std::map< std::string, std::vector<std::string> > column_data;
      
    while(!gz_in.eof()) {
      gz_in.getline(myLine, delim); q++;
      myLine.erase( std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end() ); // remove \r 
      if (myLine.length() == 0) {
        break;  // End at detection of empty line
      }
      std::istringstream entry_iss;
      entry_iss.str(myLine);
      unsigned int j = 0;
      while(!entry_iss.eof() && !entry_iss.fail() && j < columns.size()) {
        getline(entry_iss, myEntry, '\t');
        column_data[columns.at(j)].push_back(myEntry);
        j++;
      }
      if(j > columns.size()) {
        cout << "Detecting extra rows at line" << q << '\n';
        // ignore for now
      } else if(j != columns.size()) {
        cout << "Missing entries detected at line" << q << '\n';
        // attempt to repair by putting blank entries
        for(unsigned int k = j; k < columns.size(); k++) {
          column_data[columns.at(k)].push_back("");
        }
      }
    }
    List final_list;
    for(unsigned int i = 0; i < columns.size(); i++) {
      final_list.push_back(column_data[columns.at(i)], columns.at(i));
    }
    Final_final_list.push_back(final_list, header);
  }
	gz_in.closeGZ();
  return(Final_final_list);
}

#endif
// End Rcpp-only functions

// [[Rcpp::export]]
int c_gunzip(std::string s_in, std::string s_out) {
  
  if(!see_if_file_exists(s_in)) {
    cout << "File " << s_in << " does not exist!\n";
    return(-1);
  }

  GZReader gz_in;
  int ret = gz_in.LoadGZ(s_in, true);   // streamed mode
  if(ret != 0) return(ret);
	
  std::ofstream out;
  out.open(s_out, std::ofstream::binary);
  std::string myLine;
  
  while(!gz_in.iss.eof()) {
    getline(gz_in.iss, myLine, '\n');
    out << myLine << "\n";
  }
  out.flush(); out.close();
  return(0);
}

int ReadChrAlias(std::istringstream &IN,
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths
) {
  ref_names.clear();
  ref_alias.clear();
  
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

// SpliceWiz reference reader:
int readReference(std::string &reference_file, 
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths,
    CoverageBlocksIRFinder &CB_template, 
    SpansPoint &SP_template, 
    FragmentsInROI &ROI_template,
    JunctionCount &JC_template, 
    bool verbose
) { 
  (void)(verbose);

  if(!see_if_file_exists(reference_file)) {
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
  std::string headerChr ("ref-chrs.ref");
  std::string headerEOF ("EOF");
  
  bool doneCover = false;
  bool doneSpans = false;
  bool doneROI = false;
  bool doneSJ = false;
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
    getline(gz_in->iss, myBuffer, '#');  // this is the data block

    if(myLine.find(headerCover)!=std::string::npos && !doneCover) {
      std::istringstream inCoverageBlocks;
      inCoverageBlocks.str(myBuffer);
      CB_template.loadRef(inCoverageBlocks);
      doneCover = true;
    } else if(myLine.find(headerSpans)!=std::string::npos && !doneSpans) {
      SP_template.setSpanLength(5,4);
      std::istringstream inSpansPoint;
      inSpansPoint.str(myBuffer);
      SP_template.loadRef(inSpansPoint);
      doneSpans = true;
    } else if(myLine.find(headerROI)!=std::string::npos && !doneROI) {
      std::istringstream inFragmentsInROI;
      inFragmentsInROI.str(myBuffer);
      ROI_template.loadRef(inFragmentsInROI);
      doneROI = true;
    } else if(myLine.find(headerSJ)!=std::string::npos && !doneSJ) {
      std::istringstream inJuncCount;
      inJuncCount.str(myBuffer);
      JC_template.loadRef(inJuncCount);
      doneSJ = true;
    } else if(myLine.find(headerChr)!=std::string::npos && !doneChrs) {
      std::istringstream inChrAlias;
      inChrAlias.str(myBuffer);
      ReadChrAlias(inChrAlias, ref_names, ref_alias, ref_lengths);
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
  }
  return(0);
}

// SpliceWiz core:
int SpliceWizCore(std::string const &bam_file, 
    std::string const &s_output_txt, std::string const &s_output_cov,
    std::vector<std::string> &ref_names, 
    std::vector<std::string> &ref_alias,
    std::vector<uint32_t> &ref_lengths,
    CoverageBlocksIRFinder const &CB_template, 
    SpansPoint const &SP_template, 
    FragmentsInROI const &ROI_template,
    JunctionCount const &JC_template,
    bool const verbose,
    int n_threads
) {
  unsigned int n_threads_to_use = (unsigned int)n_threads;   // Should be sorted out in calling function
 
  if(!see_if_file_exists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 
 
  std::string myLine;
	if(verbose) cout << "Processing BAM file " << bam_file << "\n";
  
  pbam_in inbam((size_t)5e8, (size_t)1e9, 5);

  inbam.openFile(bam_file, n_threads_to_use);
  
  // Abort here if BAM corrupt
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  int chrcount = inbam.obtainChrs(s_chr_names, u32_chr_lens);
  if(chrcount < 1) {
    cout << bam_file << " - contains no chromosomes mapped\n";
    return(-1);
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
  
  std::vector<CoverageBlocksIRFinder*> oCB;
  std::vector<SpansPoint*> oSP;
  std::vector<FragmentsInROI*> oROI;
  std::vector<FragmentsInChr*> oChr;
  std::vector<JunctionCount*> oJC;
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oCB.push_back(new CoverageBlocksIRFinder(CB_template));
    oSP.push_back(new SpansPoint(SP_template));
    oROI.push_back(new FragmentsInROI(ROI_template));
    oChr.push_back(new FragmentsInChr);
    oJC.push_back(new JunctionCount(JC_template));
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks(bam_chr_name, bam_chr_len));

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&JunctionCount::ChrMapUpdate, &(*oJC.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&JunctionCount::ProcessBlocks, &(*oJC.at(i)), std::placeholders::_1) );
    
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

    BBchild.at(i)->openFile(&inbam);
  }
  
  // BAM processing loop
  bool error_detected = false;
#ifdef SPLICEWIZ
  Progress p(inbam.GetFileSize(), verbose);
  while(0 == inbam.fillReads() && !p.check_abort()) {
    p.increment(inbam.IncProgress());
    
#else
  while(0 == inbam.fillReads()) {
#endif
    
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      int pa_ret = BBchild.at(i)->processAll(i);
      if(pa_ret == -1) {
        
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        error_detected = true;
      }
    }
    
    if(error_detected) break;
  }

#ifdef SPLICEWIZ
  if(p.check_abort() || error_detected) {
    // interrupted:
#else
  if(error_detected) {
#endif
    
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oJC.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }


  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
      delete BBchild.at(i);
    }
  // Combine objects:
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      oJC.at(0)->Combine(*oJC.at(i));
      oChr.at(0)->Combine(*oChr.at(i));
      oSP.at(0)->Combine(*oSP.at(i));
      oROI.at(0)->Combine(*oROI.at(i));
      oCB.at(0)->Combine(*oCB.at(i));
      oFM.at(0)->Combine(*oFM.at(i));
      
      delete oJC.at(i);
      delete oChr.at(i);
      delete oSP.at(i);
      delete oROI.at(i);
      delete oCB.at(i);
      delete oFM.at(i);
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

// Write stats here:
  BBchild.at(0)->WriteOutput(myLine);
  // If first write failed, then output error and fail early
  int outret = outGZ.writeline("BAM_report\tValue"); 
  if(outret != Z_OK) {
    cout << "Error writing gzip-compressed output file\n";
    out.close();
    return(-1);
  }
  outGZ.writestring(myLine); outGZ.writeline("");

  int directionality = oJC.at(0)->Directional(myLine);
  outGZ.writeline("Directionality\tValue"); outGZ.writestring(myLine); outGZ.writeline("");

  // Generate output but save this to strings:
  std::string myLine_ROI;
  std::string myLine_JC;
  std::string myLine_SP;
  std::string myLine_Chr;
  std::string myLine_ND;
  std::string myLine_Dir;
  std::string myLine_QC;
  
  oROI.at(0)->WriteOutput(myLine_ROI, myLine_QC);
	oJC.at(0)->WriteOutput(myLine_JC, myLine_QC);
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
  delete oChr.at(0);
  delete oSP.at(0);
  delete oROI.at(0);
  delete oCB.at(0);
  delete oFM.at(0);
  delete BBchild.at(0);

  return(0);
}

#ifdef SPLICEWIZ
// [[Rcpp::export]]
int SpliceWizMain(
    std::string bam_file, std::string reference_file, std::string output_file,
    bool verbose, int n_threads
) {
  
  std::string s_output_txt = output_file + ".txt.gz";
  std::string s_output_cov = output_file + ".cov";
#else
int SpliceWizMain(
    std::string bam_file, std::string reference_file, std::string s_output_txt,
    std::string s_output_cov, int n_threads
){
	
  bool verbose = true;
#endif
  
  int use_threads = Set_Threads(n_threads);
  
  std::string s_bam = bam_file;
  std::string s_ref = reference_file;
  
  if(!see_if_file_exists(s_bam)) {
    cout << "File " << s_bam << " does not exist!\n";
    return(-1);
  } 
  if(!see_if_file_exists(s_ref)) {
    cout << "File " << s_ref << " does not exist!\n";
    return(-1);
  } 
  
  if(verbose) {
    cout << "Running SpliceWiz on " << s_bam;
    if(Has_OpenMP() != 0) cout << " with OpenMP ";
    cout << "using " << use_threads << " threads"
      << "\n" << "Reference: " << s_ref << "\n"
      << "Output file: " << s_output_txt << "\t" << s_output_cov << "\n\n"
      << "Reading reference file\n";
  }

  CoverageBlocksIRFinder * CB_template = new CoverageBlocksIRFinder;
  SpansPoint * SP_template = new SpansPoint;
  FragmentsInROI * ROI_template = new FragmentsInROI;
  JunctionCount * JC_template = new JunctionCount;
  
  std::vector<std::string> ref_names;
  std::vector<std::string> ref_alias;
  std::vector<uint32_t> ref_lengths;
  
  int ret = 0;
  
  ret = readReference(s_ref, ref_names, ref_alias, ref_lengths,
    *CB_template, *SP_template, *ROI_template, *JC_template, verbose
  );
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid NxtIRF-generated SpliceWiz reference\n";
    return(ret);
  }
  // main:
  ret = SpliceWizCore(s_bam, s_output_txt, s_output_cov,
    ref_names, ref_alias, ref_lengths,
    *CB_template, *SP_template, *ROI_template, *JC_template, verbose, use_threads);
    
  if(ret != 0) cout << "Process interrupted running SpliceWiz on " << s_bam << '\n';
  
  delete CB_template;
  delete SP_template;
  delete ROI_template;
  delete JC_template;
  return(ret);
}

#ifdef SPLICEWIZ
// [[Rcpp::export]]
int SpliceWizMain_multi(
    std::string reference_file, StringVector bam_files, StringVector output_files,
    int max_threads, bool verbose
){
	
	int use_threads = Set_Threads(max_threads);

	if(bam_files.size() != output_files.size() || bam_files.size() < 1) {
		cout << "bam_files and output_files are of different sizes\n";
		return(1);	
	}
	
	std::vector< std::string > v_bam;
	std::vector< std::string > v_out;
  for(int z = 0; z < bam_files.size(); z++) {
		v_bam.push_back(string(bam_files(z)));
		v_out.push_back(string(output_files(z)));
	}

  std::string s_ref = reference_file;
  cout << "Reading reference file\n";
  
  CoverageBlocksIRFinder * CB_template = new CoverageBlocksIRFinder;
  SpansPoint * SP_template = new SpansPoint;
  FragmentsInROI * ROI_template = new FragmentsInROI;
  JunctionCount * JC_template = new JunctionCount;
  
  std::vector<std::string> ref_names;
  std::vector<std::string> ref_alias;
  std::vector<uint32_t> ref_lengths;
  
  int ret = 0;
  
  ret = readReference(s_ref, ref_names, ref_alias, ref_lengths,
    *CB_template, *SP_template, *ROI_template, *JC_template, false);
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid NxtIRF-generated SpliceWiz reference\n";
    return(ret);
  }

	cout << "Running SpliceWiz with OpenMP using " << use_threads << " threads\n";

  for(unsigned int z = 0; z < v_bam.size(); z++) {
    std::string s_bam = v_bam.at(z);
		std::string s_output_txt = v_out.at(z) + ".txt.gz";
		std::string s_output_cov = v_out.at(z) + ".cov";
    
    int ret2 = SpliceWizCore(s_bam, s_output_txt, s_output_cov,
      ref_names, ref_alias, ref_lengths,
      *CB_template, *SP_template, *ROI_template, *JC_template, verbose, use_threads);
    if(ret2 != 0) {
      cout << "Process interrupted running SpliceWiz on " << s_bam << '\n';
      delete CB_template;
      delete SP_template;
      delete ROI_template;
      delete JC_template;
      return(ret);
    } else {
      cout << s_bam << " processed\n";
    }
	}

  delete CB_template;
  delete SP_template;
  delete ROI_template;
  delete JC_template;
  return(0);
}

#endif

// ############################ MAPPABILITY READS AND REGIONS ##################

// [[Rcpp::export]]
int c_GenerateMappabilityReads(
  std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos
) {

  if(!see_if_file_exists(genome_file)) {
    cout << "File " << genome_file << " does not exist!\n";
    return(-1);
  } 
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  
  std::ofstream outFA;
  
  // STDOUT output is only allowed in GALAXY MODE
  // Allows writing to standard output if filename is '-'
#ifndef SPLICEWIZ
  int is_stdout = 0;
  if(out_fa == "-") {
    is_stdout = 1;
  } else {
    outFA.open(out_fa, std::ios::binary);    
  }
#else
  outFA.open(out_fa, std::ios::binary);    
#endif
      
  unsigned int direction = 0;
  char * read = new char[read_len + 1];
  size_t num_reads = 0;
  
  string chr;
  string sequence;

  FastaReader inFA;
  inFA.SetInputHandle(&inGenome);
  inFA.Profile();

#ifdef SPLICEWIZ
  Progress p(inFA.total_size);
#endif

  while(!inGenome.eof() && !inGenome.fail()) {

    inFA.ReadSeq();
    sequence = inFA.sequence;
    chr = inFA.seqname;
    char * buffer = new char[sequence.length() + 1];
    std::strcpy (buffer, sequence.c_str());
#ifdef SPLICEWIZ
    size_t seq_progress = 0;
#endif
    for(
        unsigned int bufferPos = 1; 
        bufferPos < sequence.length() - read_len + 1; 
        bufferPos += read_stride
    ) {
      memcpy(read, &buffer[bufferPos - 1], read_len);
      num_reads += 1;
      if(checkDNA(read, read_len)) {       
        std::string write_seq = GenerateReadError(
          read, read_len, error_pos, direction, num_reads
        ) ;

#ifndef SPLICEWIZ
        if(is_stdout == 1) {
          cout << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" 
            << std::to_string(bufferPos)
            << '\n' << write_seq << '\n';  
        } else {
          outFA << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" 
            << std::to_string(bufferPos)
            << '\n' << write_seq << '\n'; 
        }
#else
        outFA << (direction == 0 ? ">RF!" : ">RR!") << chr << "!" 
          << std::to_string(bufferPos)
          << '\n' << write_seq << '\n'; 
#endif

        direction = (direction == 0 ? 1 : 0);
        
        // update progress bar
#ifdef SPLICEWIZ
        if(num_reads % 100000 == 0) {
          p.increment(bufferPos - seq_progress);
          seq_progress = bufferPos;
        }
#endif
      }
    }
#ifdef SPLICEWIZ
    p.increment(sequence.length() - seq_progress);
#endif
    delete[] buffer;
  }
  delete[] read;
  
  inGenome.close();

#ifndef SPLICEWIZ
  if(is_stdout == 0) {
    outFA.flush();
    outFA.close();
  }
#else  
  outFA.flush();
  outFA.close();
#endif
  
  cout << num_reads << " synthetic reads generated\n";
  return(0);
}


#ifdef SPLICEWIZ
// [[Rcpp::export]]
int c_GenerateMappabilityRegions(
    std::string bam_file, std::string output_file, 
    int threshold, int includeCov, bool verbose,
    int n_threads
){
  
  std::string s_output_txt = output_file + ".txt";
  std::string s_output_cov = output_file + ".cov";
#else
int c_GenerateMappabilityRegions(
    std::string bam_file, std::string s_output_txt, 
    int threshold, int n_threads, std::string s_output_cov
){	
	bool verbose = true;
#endif

  if(!see_if_file_exists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 

  int use_threads = Set_Threads(n_threads);
  unsigned int n_threads_to_use = (unsigned int)use_threads;
 
  std::string myLine;
	if(verbose) cout << "Calculating Mappability Exclusions from aligned synthetic reads in BAM file " << bam_file << "\n";

  pbam_in inbam((size_t)5e8, (size_t)1e9, 5);
  inbam.openFile(bam_file, n_threads_to_use);

  // Assign children:
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i)->openFile(&inbam);
  }
  
  // BAM processing loop
#ifdef SPLICEWIZ
  Progress p(inbam.GetFileSize(), verbose);
  while(0 == inbam.fillReads() && !p.check_abort()) {
    p.increment(inbam.IncProgress());
  
#else
  while(0 == inbam.fillReads()) {
#endif
  
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      BBchild.at(i)->processAll(i, true);
    }
  }

#ifdef SPLICEWIZ
  if(p.check_abort()) {
    // interrupted:
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }
#endif
  
  inbam.closeFile();

  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
      delete BBchild.at(i);
    }
  // Combine objects:
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      oFM.at(0)->Combine(*oFM.at(i));

      delete oFM.at(i);
    }
  }

#ifdef SPLICEWIZ
  if(includeCov == 1) {
#else
  if(!s_output_cov.empty()) {
#endif
   // Write Coverage Binary file:
    std::ofstream ofCOV;
    ofCOV.open(s_output_cov, std::ofstream::binary);  
    covWriter outCOV;
    outCOV.SetOutputHandle(&ofCOV);
    oFM.at(0)->WriteBinary(&outCOV, verbose, n_threads_to_use);
    ofCOV.close();
  }
  
  std::ofstream outFragsMap;
  outFragsMap.open(s_output_txt, std::ifstream::out);
	
  oFM.at(0)->WriteOutput(&outFragsMap, threshold, verbose);
  outFragsMap.flush(); outFragsMap.close();

  delete oFM.at(0);
  delete BBchild.at(0);

  return(0);
}


#ifdef SPLICEWIZ
// [[Rcpp::export]]
int c_BAM2COV(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads
){
  
#else
int c_BAM2COV(
    std::string bam_file, std::string output_file, int n_threads
){	
	bool verbose = true;
#endif

  std::string s_output_cov = output_file;

  if(!see_if_file_exists(bam_file)) {
    cout << "File " << bam_file << " does not exist!\n";
    return(-1);
  } 

  int use_threads = Set_Threads(n_threads);
  unsigned int n_threads_to_use = (unsigned int)use_threads;
 
  std::string myLine;
	if(verbose) cout << "Creating COV file from " << bam_file << "\n";

  pbam_in inbam((size_t)5e8, (size_t)1e9, 5);
  inbam.openFile(bam_file, n_threads_to_use);

  // Assign children:
  std::vector<FragmentsMap*> oFM;
  std::vector<BAM2blocks*> BBchild;

  for(unsigned int i = 0; i < n_threads_to_use; i++) {
    oFM.push_back(new FragmentsMap);
    BBchild.push_back(new BAM2blocks);

    BBchild.at(i)->registerCallbackChrMappingChange( std::bind(&FragmentsMap::ChrMapUpdate, &(*oFM.at(i)), std::placeholders::_1) );
    BBchild.at(i)->registerCallbackProcessBlocks( std::bind(&FragmentsMap::ProcessBlocks, &(*oFM.at(i)), std::placeholders::_1) );

    BBchild.at(i)->openFile(&inbam);
  }
  
  // BAM processing loop
#ifdef SPLICEWIZ
  Progress p(inbam.GetFileSize(), verbose);
  while(0 == inbam.fillReads() && !p.check_abort()) {
    p.increment(inbam.IncProgress());
  
#else
  while(0 == inbam.fillReads()) {
#endif
  
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
    #endif
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      BBchild.at(i)->processAll(i);
    }
  }

#ifdef SPLICEWIZ
  if(p.check_abort()) {
    // interrupted:
    for(unsigned int i = 0; i < n_threads_to_use; i++) {
      delete oFM.at(i);
      delete BBchild.at(i);
    }
    return(-1);
  }
#endif
  
  inbam.closeFile();

  if(n_threads_to_use > 1) {
    if(verbose) cout << "Compiling data from threads\n";
  // Combine BB's and process spares
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      BBchild.at(0)->processSpares(*BBchild.at(i));
      delete BBchild.at(i);
    }
  // Combine objects:
    for(unsigned int i = 1; i < n_threads_to_use; i++) {
      oFM.at(0)->Combine(*oFM.at(i));

      delete oFM.at(i);
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


// ################################## MAIN #####################################

#ifndef SPLICEWIZ

void print_usage(std::string exec) {
  cout << "Usage:\n\t"
    << exec << " about\n\t\tDisplays version and OpenMP status\n\t"
    << exec <<  " main (-t 4) in.bam SpliceWiz.ref.gz out.txt.gz out.cov\n\t\t"
    << "(runs SpliceWiz - optionally using 4 threads)\n\t"
    << exec <<  " bam2cov (-t 4) in.bam out.cov"
    << "(runs SpliceWiz's Bam to Cov utility - optionally using 4 threads)\n\t"
    << exec <<  " gen_map_reads genome.fa reads_out.fa 70 10\n\t\t"
    << "(where synthetic read length = 70, and read stride = 10)\n\t"
    << exec <<  " gen_map_regions (-t 4) aligned_reads.bam 4 map.bed {map.cov}\n\t\t"   
    << "(where threshold for low mappability = 4, - optionally using 4 threads\n";
}

// main
int main(int argc, char * argv[]) {
	// Command line usage:
    // spliceWiz main -t N sample.bam IRFinder.ref.gz Output.txt.gz 
    // spliceWiz bam2cov -t N sample.bam sample.cov
    // spliceWiz gen_map_reads genome.fa reads_to_map.fa 70 10
    // spliceWiz gen_map_regions mappedreads.bam mappability.bed
  int n_thr = 1;
  int ret = 0;
  if(argc < 2) {
    print_usage(argv[0]);
    return 0;
  } else if(std::string(argv[1]) == "gen_map_reads") {
      if(argc < 3){
        print_usage(argv[0]);
        exit(1);
      }
      std::string s_genome,s_output;
      int read_len,read_stride,read_error;

      s_genome = argv[2];
      s_output = argv[3];

      if(argc > 4) {
        read_len = atoi(argv[4]);
      } else {
        read_len = 70;
      }
      if(argc > 5) {
        read_stride = atoi(argv[5]);
      } else {
        read_stride = 10;
      }
      read_error = 1 + (read_len / 2);
      ret = c_GenerateMappabilityReads(
        s_genome, s_output, read_len, read_stride, read_error
      );
      if(ret == 0) exit(0);
      exit(1);
  } else if(std::string(argv[1]) == "gen_map_regions") {
      if(argc < 5){
        print_usage(argv[0]);
        exit(1);
      }
      
      std::string s_bam,s_output,s_cov;
      if(std::string(argv[2]) == "-t") {
        n_thr = atoi(argv[3]);
        s_bam = argv[4];
        int threshold = atoi(argv[5]);
        s_output = argv[6];
        if(argc > 7) {
          s_cov = argv[7];
          ret = c_GenerateMappabilityRegions(s_bam, s_output, threshold, n_thr, s_cov);
        } else {
          ret = c_GenerateMappabilityRegions(s_bam, s_output, threshold, n_thr);
        }
      } else {
        s_bam = argv[2];
        s_output = argv[3];
        int threshold = atoi(argv[4]);
        if(argc > 5) {
          s_cov = argv[5];
          ret = c_GenerateMappabilityRegions(s_bam, s_output, threshold, n_thr,s_cov);
        } else {
          ret = c_GenerateMappabilityRegions(s_bam, s_output, threshold, n_thr);
        }
      }
      exit(ret);;      
  } else if(std::string(argv[1]) == "about") {
      std::string version = "0.99.0";
      cout << "SpliceWiz version " << version << "\t";
      int ret = Test_OpenMP_For();
      if(ret == -1) {
        cout << "compiled without OpenMP\n";
      } else if(ret == 0) {
        cout << "compiled with OpenMP, threads avail = " << Has_OpenMP() << "\n";
      } else if(ret == 1) {
        cout << "compiled with OpenMP, threads avail = " << 1 << "\n";
      } else {
        cout << "compiled with OpenMP but parallel loops could not be initialised" << "\n";
      }
      exit(0);
  } else if(std::string(argv[1]) == "main") {
      if(argc < 6){
        print_usage(argv[0]);
        exit(1);
      }
      
      int n_thr = 1; std::string s_bam,s_ref,s_output_txt,s_output_cov;
      
      if(std::string(argv[2]) == "-t") {
        n_thr = atoi(argv[3]);
        s_bam = argv[4];
        s_ref = argv[5];
        s_output_txt = argv[6];		
        s_output_cov = argv[7];
      } else {
        s_bam = argv[2];
        s_ref = argv[3];
        s_output_txt = argv[4];		
        s_output_cov = argv[5];
      }
      ret = SpliceWizMain(s_bam, s_ref, s_output_txt, s_output_cov, n_thr);
      exit(ret);
  } else if(std::string(argv[1]) == "bam2cov") {
      if(argc < 4){
        print_usage(argv[0]);
        exit(1);
      }
      
      int n_thr = 1; std::string s_bam,s_output_cov;
      
      if(std::string(argv[2]) == "-t" && argc == 6) {
        n_thr = atoi(argv[3]);
        s_bam = argv[4];
        s_output_cov = argv[5];
      } else if(argc == 4){
        s_bam = argv[2];
        s_output_cov = argv[3];
      } else {
        print_usage(argv[0]);
        exit(1);
      }
      ret = c_BAM2COV(s_bam, s_output_cov, n_thr);
      exit(ret);
  } else {
    print_usage(argv[0]);
    exit(0);
  }
}

#endif

// Has_OpenMP