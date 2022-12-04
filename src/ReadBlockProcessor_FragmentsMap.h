/* ReadBlockProcessor_FragmentsMap.h Reads fragment coverage

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

#ifndef CODE_READBLOCKPROCESSOR_FRAGMENTSMAP
#define CODE_READBLOCKPROCESSOR_FRAGMENTSMAP

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ReadBlockProcessor.h"
#include "FragmentBlocks.h"

#include "covTools.h"

#include "SpliceWiz.h"

class FragmentsMap : public ReadBlockProcessor {
  // Counts mappability.
private:
  // 0 = -, 1 = +, 2 = both
  std::vector< std::vector< std::pair<unsigned int, int> > > chrName_vec_final[3];
  std::vector< std::vector< std::pair<unsigned int, int> > > chrName_vec_new[3];
  std::vector< std::vector< std::pair<unsigned int, int> > > temp_chrName_vec_new[3];

  uint32_t frag_count = 0;
	int sort_and_collapse_temp();

	bool final_is_sorted = false;
  
  vector<chr_entry> chrs;
public:
	void Combine(FragmentsMap &child);
	
  int sort_and_collapse_final(bool verbose, unsigned int n_threads_to_use);

  void ProcessBlocks(const FragmentBlocks &blocks);
  void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
  int WriteOutput(
    std::ostream *os, int threshold = 4, 
    bool verbose = false, unsigned int n_threads_to_use = 1
  ) ;
  int WriteBinary(
    covWriter *os, bool verbose = false, 
    unsigned int n_threads_to_use = 1
  ) ;
  
  void updateCoverageHist(std::map<unsigned int,unsigned int> &hist, unsigned int start, unsigned int end, unsigned int dir, const unsigned int &refID, bool debug = false) const;
};

#endif