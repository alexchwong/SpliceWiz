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

#include <chrono>
#include "main.h"

// for malloc_trim (Linux only)
#ifdef __linux__
#include <malloc.h>
#endif

bool checkFileExists(const std::string& name) {
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
int Has_OpenMP() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 0;
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
  
  if(!checkFileExists(s_in)) {
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
  
  if(!checkFileExists(s_in)) {
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
  
  if(!checkFileExists(s_in)) {
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

  if(!checkFileExists(s_in)) {
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
  
  if(!checkFileExists(s_in)) {
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

#ifdef SPLICEWIZ
// [[Rcpp::export]]
int SpliceWizMain(
    std::string bam_file, std::string reference_file, std::string output_file,
    bool verbose, int n_threads, bool skipCOV, bool multiRead
) {
  
  std::string s_output_txt = output_file + ".txt.gz";
  std::string s_output_cov = output_file + ".cov";
#else
int SpliceWizMain(
    std::string bam_file, std::string reference_file, std::string s_output_txt,
    std::string s_output_cov, int n_threads
){
	
  bool verbose = true;
  bool multiRead = false;
  bool skipCOV = false;
#endif
   
  std::string s_bam = bam_file;
  std::string s_ref = reference_file;  

	std::vector< std::string > v_bam;
	std::vector< std::string > v_out_txt;
	std::vector< std::string > v_out_cov;

  v_bam.push_back(s_bam);
  v_out_txt.push_back(s_output_txt);
  v_out_cov.push_back(s_output_cov);

  swEngine Engine;
  Engine.Set_Threads(n_threads);

  if(verbose) cout << "Reading reference file\n";
  int ret = Engine.readReference(s_ref, verbose);
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid SpliceWiz reference\n";
    return(ret);
  }
  int ret2 = Engine.SpliceWizMultiCore(
    v_bam, v_out_txt, v_out_cov,
    verbose, skipCOV, multiRead
  );
  Engine.clear();
#ifdef __linux__
  malloc_trim(0);
#endif
  return(ret2);
}

#ifdef SPLICEWIZ
// [[Rcpp::export]]
int SpliceWizMain_multi(
    std::string reference_file, StringVector bam_files, StringVector output_files,
    int max_threads, bool verbose, bool skipCOV, bool multiRead
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

  std::string s_ref = reference_file;
  if(!checkFileExists(s_ref)) {
    cout << "File " << s_ref << " does not exist!\n";
    return(-1);
  } 

  swEngine Engine;
  Engine.Set_Threads(max_threads);

  if(verbose) cout << "Reading reference file\n";
  int ret = Engine.readReference(s_ref, verbose);
  if(ret != 0) {
    cout << "Reading Reference file failed. Check if SpliceWiz.ref.gz exists and is a valid SpliceWiz reference\n";
    return(ret);
  }
  int ret2 = Engine.SpliceWizMultiCore(
    v_bam, v_out_txt, v_out_cov,
    verbose, skipCOV, multiRead
  );
  Engine.clear();
#ifdef __linux__
  // cout << "malloc trim() = " << malloc_trim(0) << "\n";
  malloc_trim(0);
#endif

  return(ret2);
}
#endif

// ############################ MAPPABILITY READS AND REGIONS ##################

// [[Rcpp::export]]
int c_GenerateMappabilityReads(
  std::string genome_file, std::string out_fa,
	int read_len, int read_stride, int error_pos
) {

  if(!checkFileExists(genome_file)) {
    cout << "File " << genome_file << " does not exist!\n";
    return(-1);
  } 
  
  std::ifstream inGenome;
  inGenome.open(genome_file, std::ifstream::in);
  
  std::ofstream outFA;
  synthReadGenerator synth((unsigned int)read_len, error_pos);
  
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
      if(synth.checkDNA(read)) {       
        std::string write_seq = synth.GenerateReadError(
          read, direction, num_reads
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
  bool doCov = (includeCov == 1);

#else
int c_GenerateMappabilityRegions(
    std::string bam_file, std::string s_output_txt, 
    int threshold, int n_threads, std::string s_output_cov
){	
	bool verbose = true;
  bool doCov = !s_output_cov.empty();

#endif

  swEngine Engine;
  Engine.Set_Threads(n_threads);
  
  int ret = Engine.MappabilityRegionsCore(
    bam_file, s_output_txt, s_output_cov,
    threshold, doCov,
    verbose, false
  );

  Engine.clear();
#ifdef __linux__
  malloc_trim(0);
#endif

  return(ret);
}


#ifdef SPLICEWIZ
// [[Rcpp::export]]
int c_BAM2COV(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, bool multiRead
){
  
#else
int c_BAM2COV(
    std::string bam_file, std::string output_file, int n_threads, bool multiRead
){	
	bool verbose = true;
#endif

  swEngine Engine;
  int nthr = Engine.Set_Threads(n_threads);
  
  std::string s_bam = bam_file;
  
  if(verbose) {
    cout << "Running BAM2COV (ompBAM) " << s_bam;
    cout << " using " << nthr << " threads\n";
  }

  auto start = chrono::steady_clock::now();
  auto check = start;
  int ret2 = Engine.BAM2COVcore(s_bam, output_file, verbose, multiRead);

  Engine.clear();
#ifdef __linux__
  malloc_trim(0);
#endif

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

#ifdef SPLICEWIZ
// [[Rcpp::export]]
int c_doStats(
    std::string bam_file, std::string output_file, 
    bool verbose, int n_threads, bool multiRead
){
  
#else
int c_doStats(
    std::string bam_file, std::string output_file, int n_threads, bool multiRead
){	
	bool verbose = true;
#endif

  swEngine Engine;
  int nthr = Engine.Set_Threads(n_threads);
  
  std::string s_bam = bam_file;
  
  if(verbose) {
    cout << "Running doStats (ompBAM) " << s_bam;
    cout << " using " << nthr << " threads\n";
  }

  auto start = chrono::steady_clock::now();
  auto check = start;
  int ret2 = Engine.doStatsCore(s_bam, output_file, verbose, multiRead);
  
  Engine.clear();
#ifdef __linux__
  malloc_trim(0);
#endif

  if(ret2 == -2) {
    cout << "Process interrupted running doStats on " << s_bam << '\n';
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
    // spliceWiz main -t N sample.bam SpliceWiz.ref.gz Output.txt.gz 
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