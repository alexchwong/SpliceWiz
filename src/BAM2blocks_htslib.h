/* htsBAM2blocks_htslib.h Convert reads / fragments to FragmentBlocks (htslib)

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

#ifndef CODE_BAM2BLOCKS_HTS
#define CODE_BAM2BLOCKS_HTS

#include "includedefine.h"
#include "SpliceWiz.h"

#include "htslib/hfile.h"
#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/thread_pool.h"

#include "FragmentBlocks.h"

/* Little Endian .. for big endian each group of 4 bytes needs to be reversed before individual members are accessed. */


class htsBAM2blocks {
    FragmentBlocks oBlocks;

    std::vector< std::function<void(const std::vector<chr_entry> &)> > callbacksChrMappingChange;
    std::vector< std::function<void(const FragmentBlocks &)> > callbacksProcessBlocks;

    void cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len);

    unsigned int processPair(bam1_t * read1, bam1_t * read2);
    unsigned int processSingle(bam1_t * read1, bool mappability_mode = false);

    // Statistics.
    unsigned long cReadsProcessed;
    unsigned long long totalNucleotides;
    unsigned long cShortPairs;
    unsigned long cIntersectPairs;
    unsigned long cLongPairs;
    unsigned long cSingleReads;
    unsigned long cPairedReads;
    unsigned long cErrorReads;
    unsigned long cSkippedReads;
    unsigned long cChimericReads;

    bam1_t * reads[2];

    std::vector<chr_entry> chrs;

    std::map< std::string, bam1_t* > * spare_reads;
    bam1_t * SupplyRead(std::string& read_name);    
    // int realizeSpareReads();
    int read_name(bam1_t * b, std::string & dest);

// Disable copy construction / assignment (doing so triggers compile errors)
    htsBAM2blocks(const htsBAM2blocks &t);
    htsBAM2blocks & operator = (const htsBAM2blocks &t);
  public:
  	htsBAM2blocks();
  	htsBAM2blocks(
      std::vector<std::string> & ref_names, 
      std::vector<uint32_t> & ref_lengths
    );  // Define BB with defined references
    ~htsBAM2blocks();
    
  	unsigned int initializeChrs();

  	int processAll(
      std::vector<bam1_t*> & bpool, int starts, int ends,
      bool mappability_mode = false
    );

  	int processSpares(htsBAM2blocks& other);
  	int processStats(htsBAM2blocks& other);

  	int WriteOutput(std::string& output);

    void registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback );
    void registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback );
};


#endif
