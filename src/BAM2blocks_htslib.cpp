/* BAM2blocks_htslib.cpp Convert reads / fragments to FragmentBlocks

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

// WARNING: code is little endian only!

#include "BAM2blocks_htslib.h"
#include <chrono>

htsBAM2blocks::htsBAM2blocks() {
  oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

  cReadsProcessed = 0;
  totalNucleotides = 0;
  cShortPairs = 0;
  cIntersectPairs = 0;
  cLongPairs = 0;
  cSingleReads = 0;
  cPairedReads = 0;
  cErrorReads = 0;
  cSkippedReads = 0;
  cChimericReads = 0;
  
  spare_reads = new std::map< std::string, bam1_t* >;
  
  // Initialize vector of bam1_t pointers (2)
  reads[0] = NULL;
  reads[1] = NULL;
}

htsBAM2blocks::htsBAM2blocks(
    std::vector<std::string> & ref_names, 
    std::vector<uint32_t> & ref_lengths
) {
  oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

  cReadsProcessed = 0;
  totalNucleotides = 0;
  cShortPairs = 0;
  cIntersectPairs = 0;
  cLongPairs = 0;
  cSingleReads = 0;
  cPairedReads = 0;
  cErrorReads = 0;
  cSkippedReads = 0;
  cChimericReads = 0;
  
  if(ref_names.size() > 0) {
    for(unsigned int i = 0; i < ref_names.size(); i++) {
      chrs.push_back(chr_entry(i, ref_names.at(i), (int32_t)ref_lengths.at(i)));
    }
  }
  
  spare_reads = new std::map< std::string, bam1_t* >;

  // Initialize vector of bam1_t pointers (2)
  reads[0] = NULL;
  reads[1] = NULL;
}

htsBAM2blocks::~htsBAM2blocks() {
  for(auto it = spare_reads->begin(); it != spare_reads->end(); it++) {
    bam_destroy1(it->second);
  }
  delete spare_reads;
}

unsigned int htsBAM2blocks::initializeChrs() {

  for (auto & callback : callbacksChrMappingChange ) {
    callback(chrs);
  }
  return(0);
}

void htsBAM2blocks::cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len) {
  bool inBlock = true;
  int relpos = 0;
  int curblock = 0;
  starts.resize(1);  // Is this expensive or not -- does this call destroy on further items, or is it a single op, adjusting the end? If expensive we can revert to earlier behaviour where we keep track of how many blocks, just overwriting relevant parts of the vector.
  lens.resize(1);
  starts[curblock] = 0;
  lens[curblock] = 0;

  for (; n_cigar_op > 0; n_cigar_op--) {
    if (inBlock) {
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          // increment len of last block
          lens[curblock] += (*cigar >> 4);
          relpos += (*cigar >> 4);
          break;
        case 3:
          curblock++;
          relpos += (*cigar >> 4);
          // extend arrays. 
          starts.push_back(relpos);
          lens.push_back(0);
          inBlock = false;
          break;
      }
    }else{
      switch (*cigar & 15) {
        case 0: case 2: case 7: case 8:
          lens[curblock] = (*cigar >> 4);
          relpos += (*cigar >> 4);
          inBlock = true;
          break;
        case 3:
          // push start of next further out
          relpos += (*cigar >> 4);
          starts[curblock] = relpos;
          break;
        }
    }
    cigar++;
  }
  ret_genome_len = relpos;
}


//OK - translated - doesn't call the callbacks yet though.
unsigned int htsBAM2blocks::processPair(bam1_t * read1, bam1_t * read2) {
  // R1 is to the left of R2 (or equal starts).
  int r1_genome_len;
  //int r1_blocks;
  int r2_genome_len;

  bam1_t * r1 = read1;
  bam1_t * r2 = read2;
  
  // string debugstate;


  if (r1->core.flag & 0x40) {
    //this is first of pair.
    if (r1->core.flag & 0x10) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }else{
    if (r1->core.flag & 0x20) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }

  cigar2block(bam_get_cigar(r1), r1->core.n_cigar, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  cigar2block(bam_get_cigar(r2), r2->core.n_cigar, oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);

  bool merge_reads = false;
  bool swap_reads = false;
  bool goodPair = true;
  
  if (r1->core.pos + r1_genome_len < r2->core.pos) {
    cLongPairs++;
    //reads do not intersect
    oBlocks.readCount = 2;
    // debugstate.append( "-Long-");
  }else if (r1->core.pos + r1_genome_len >= r2->core.pos + r2_genome_len){
    if(r1->core.pos == r2->core.pos && r1_genome_len > r2_genome_len) {
      cIntersectPairs++;
      swap_reads = true;
      merge_reads = true;
    } else {
      cShortPairs++;
    }    
    // Read 2 is a short read & read 1 fully contains it (or perhaps just a trimmed read with two exactly complementary reads remaining).
    oBlocks.readCount = 1;
    // debugstate.append( "-Short-");    
  }else{
    // debugstate.append( "-Intersect-");
    cIntersectPairs++;
    
    oBlocks.readCount = 1;
    merge_reads = true;
  }
  
  if(swap_reads) {
    // recalculate blocks based on read1 <-> read2
    r2 = read1;
    r1 = read2;
    
    if (r1->core.flag & 0x40) {
      //this is first of pair.
      if (r1->core.flag & 0x10) {
        oBlocks.direction = 0;
      }else{
        oBlocks.direction = 1;
      }
    }else{
      if (r1->core.flag & 0x20) {
        oBlocks.direction = 0;
      }else{
        oBlocks.direction = 1;
      }
    }

    cigar2block(bam_get_cigar(r1), r1->core.n_cigar, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
    cigar2block(bam_get_cigar(r2), r2->core.n_cigar, oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);
  }
  
  if(merge_reads) {
    // We have two reads that intersect - construct just one complete fragment.

// Guaranteed assumptions:
//   Read 1 starts to the left of Read 2.
//   Read 2 end extends beyond the end of Read 1 end.
    int r1pos = r1->core.pos;
    int r2pos = r2->core.pos;
    for (unsigned int i = 0; i < oBlocks.rStarts[0].size(); i++) {
        if (r1pos + oBlocks.rStarts[0][i] + oBlocks.rLens[0][i] >= r2pos) {
          if (r1pos + oBlocks.rStarts[0][i] <= r2pos) {
            oBlocks.rLens[0][i] = r2pos - r1pos - oBlocks.rStarts[0][i] + oBlocks.rLens[1][0];
            //r1_blocks = i + r2_blocks;
            oBlocks.rStarts[0].resize(i + oBlocks.rStarts[1].size());
            oBlocks.rLens[0].resize(i + oBlocks.rStarts[1].size());
            // Maybe this can be optimised by using push_back below instead of running resize.
            for (unsigned int j = 1; j < oBlocks.rStarts[1].size(); j++) {
              i++;
              oBlocks.rLens[0][i] = oBlocks.rLens[1][j];
              oBlocks.rStarts[0][i] = oBlocks.rStarts[1][j] + r2pos - r1pos;
            }
            r1_genome_len = r2pos - r1pos + r2_genome_len;
            break;
          }else{
            //cerr << "Fault with this synthetic read, outputting each of the overlapping reads as singles: " << read1->read_name << endl;
            // This error is not worth reporting. The current version of STAR outputs a good number of these, concordance would be nice, but it is better to get at least one read illustrating the splice junction.
            goodPair = false;
            oBlocks.readCount = 2;
          }
        }
    }

    if (!goodPair) {
      oBlocks.readCount = 2;
    }
  }
  oBlocks.chr_id = r1->core.tid;
  oBlocks.readStart[0] = r1->core.pos;
  oBlocks.readEnd[0] = r1->core.pos + r1_genome_len;
  oBlocks.readName.resize(r1->core.l_qname - 1);
  oBlocks.readName.replace(0, r1->core.l_qname - 1, bam_get_qname(r1), r1->core.l_qname - 1); // is this memory/speed efficient?

  unsigned int totalBlockLen = 0;
  for (auto blockLen: oBlocks.rLens[0]) {
    totalBlockLen += blockLen;
  }
  if (oBlocks.readCount > 1) {
    oBlocks.readStart[1] = r2->core.pos;
    oBlocks.readEnd[1] = r2->core.pos + r2_genome_len;
    for (auto blockLen: oBlocks.rLens[1]) {
      totalBlockLen += blockLen;
    }
  }
  //DEBUG:
  // oBlocks.readName.append(debugstate);
  // oBlocks.readName.append(to_string(oBlocks.readCount));
// TODO - restructure -- we could instead do the manipulation from 2 reads-> 1 synthetic in a non-const callback.
//        not required until that future flexibility is needed if part of the framework is repurposed.
  for (auto & callback : callbacksProcessBlocks ) {
    callback(oBlocks);
  }
  return totalBlockLen;
}


unsigned int htsBAM2blocks::processSingle(bam1_t * read1, bool mappability_mode) {
  int r1_genome_len;

  // string debugstate;

  if (read1->core.flag & 0x10) {
    oBlocks.direction = 0;
  }else{
    oBlocks.direction = 1;
  }

  cigar2block(bam_get_cigar(read1), read1->core.n_cigar, oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  oBlocks.readCount = 1;

  oBlocks.chr_id = read1->core.tid;
  oBlocks.readStart[0] = read1->core.pos;
  oBlocks.readEnd[0] = read1->core.pos + r1_genome_len;
  oBlocks.readName.resize(read1->core.l_qname - 1);
  oBlocks.readName.replace(0, read1->core.l_qname - 1, bam_get_qname(read1), read1->core.l_qname - 1); // is this memory/speed efficient?
  
  // Below block only run from Mappability - only process reads if they are
  // mapped to the exact position from which synthetic reads were
  // generated in the genome
  if(mappability_mode) {
    // Return if not a perfect 70M match
    if(read1->core.n_cigar != 1) return(0);
    if( (*(bam_get_cigar(read1)) & 15) != 0) return(0);
    
    std::istringstream iss;
    std::string subline;
    iss.str(oBlocks.readName);
    std::getline(iss, subline, '!');  // ignore strand
    std::getline(iss, subline, '!');  // this is chr name
    if(0 != strncmp(
        subline.c_str(), chrs.at(oBlocks.chr_id).chr_name.c_str(), 
        subline.size()
      ) || 
      chrs.at(oBlocks.chr_id).chr_name.size() != subline.size()) 
    {
        return 0;
    }
    std::getline(iss, subline, '!');  // this is chr pos
    if(stoul(subline) != oBlocks.readStart[0] + 1) {
      return 0;
    }
  }
  
  //DEBUG:
  // oBlocks.readName.append(debugstate);
  // oBlocks.readName.append(to_string(oBlocks.readCount));
  //cout << "process pair - callbacks" << endl;  
  for (auto & callback : callbacksProcessBlocks ) {
    callback(oBlocks);
  }
  unsigned int totalBlockLen = 0;
  for (auto blockLen: oBlocks.rLens[0]) {
    totalBlockLen += blockLen;
  }
  return totalBlockLen;
}

// Prints statistics to file
int htsBAM2blocks::WriteOutput(std::string& output) {
  std::ostringstream oss;
  cErrorReads = spare_reads->size();
  oss << "Total reads processed\t" << cReadsProcessed << '\n';
  oss << "Total nucleotides\t" << totalNucleotides << '\n';
  oss << "Total singles processed\t" << cSingleReads << '\n';
  oss << "Total pairs processed\t" << cShortPairs+cIntersectPairs+cLongPairs << '\n';
  oss << "Short pairs\t" << cShortPairs << '\n';
  oss << "Intersect pairs\t" << cIntersectPairs << '\n';
  oss << "Long pairs\t" << cLongPairs << '\n';
  oss << "Skipped reads\t" << cSkippedReads << '\n';
  oss << "Chimeric reads\t" << cChimericReads << '\n';
  oss << "Error / Unpaired reads\t" << cErrorReads << '\n';
  output = oss.str();
  return(0);
}

// Returns a spare read; deletes same from BB's storage
bam1_t *  htsBAM2blocks::SupplyRead(std::string& read_name) {
  // Supplies the pointer to the last spare read
  // When called, transfer ownership of read to the parent BB
  if(spare_reads->size() == 0) return(NULL);
  auto it = spare_reads->begin();
  read_name = it->first;
  bam1_t * read = it->second;
  spare_reads->erase(it);
  cErrorReads -= 1;
  return(read);
}

// Summates statistics from child BB's
int htsBAM2blocks::processStats(htsBAM2blocks& other) {
															   
  cReadsProcessed += other.cReadsProcessed;
  totalNucleotides += other.totalNucleotides;
    
  cShortPairs += other.cShortPairs;
  cIntersectPairs += other.cIntersectPairs;
  cLongPairs += other.cLongPairs;
  cSingleReads += other.cSingleReads;
  cPairedReads += other.cPairedReads;
  cErrorReads += other.cErrorReads;
  cSkippedReads += other.cSkippedReads;
  cChimericReads += other.cChimericReads;
  
  other.cReadsProcessed = 0;
  other.totalNucleotides = 0;
    
  other.cShortPairs = 0;
  other.cIntersectPairs = 0;
  other.cLongPairs = 0;
  other.cSingleReads = 0;
  other.cPairedReads = 0;
  other.cErrorReads = 0;
  other.cSkippedReads = 0;
  other.cChimericReads = 0;
  return(0);
}

// Tries to match reads between other BB and self
int htsBAM2blocks::processSpares(htsBAM2blocks& other) {
  // Combines two BB's, and processes any matching paired reads
  bam1_t * spare_read;
  std::string read_name_s;
  while(1) {
    spare_read = other.SupplyRead(read_name_s);
        
    if(!spare_read) {
      break;
    }
    // cout << "Unmatched read: " 
      // << bam_get_qname(spare_read);
    
    auto it_read = spare_reads->find(read_name_s);
    if(it_read != spare_reads->end()){
      // cout << " matched\n";
      cPairedReads ++;
      if (spare_read->core.tid != it_read->second->core.tid) {
        cChimericReads += 1;
      } else {
        if (spare_read->core.pos <= it_read->second->core.pos) {    
          totalNucleotides += processPair(spare_read, it_read->second);
        } else{              
          totalNucleotides += processPair(it_read->second, spare_read);
        }
        cReadsProcessed+=2;
      }
      bam_destroy1(it_read->second);
      spare_reads->erase(read_name_s);
      bam_destroy1(spare_read);
      cErrorReads-=1;
    } else {
      // cout << " not matched\n";
      spare_reads->insert({read_name_s, spare_read});
    }
  }
  
  return(0);
}

// Saves bam1_t spare reads to dedicated buffer space
/*
int htsBAM2blocks::realizeSpareReads() {
  for (auto it = spare_reads->begin(); it != spare_reads->end(); it++) {
    if(!it->second->isReal()) {
      it->second->realize();
    }
  }
  return(0);
}
*/

// Extracts read name from bam1_t
int htsBAM2blocks::read_name(bam1_t * b, std::string & dest) {
  dest.clear();
  char *tmp = bam_get_qname(b);
  dest.assign(tmp);
  return(b->core.l_qname);
}

// Main function
int htsBAM2blocks::processAll(
  std::vector<bam1_t*> & bpool, 
  int starts, int ends,
  bool mappability_mode
) {
  // Reads from pbam_in until finished; do not create output
  unsigned int idx = 0;
  unsigned int nucs_proc = 0;
  
  // Use map pointer spare_reads:
  std::map< std::string, bam1_t* > * new_spare_reads;  

  bam1_t * b = NULL;
  std::string read_name_s;
  
  auto start = chrono::steady_clock::now();
  auto check = start;
  
  // while(1) {
  // cout << "htsBB profiling bpool start " << starts << ", end " << ends << '\n';

  if(starts > -1) {
    for(int i = starts; i < ends; i++) {
      check = chrono::steady_clock::now();
      if(chrono::duration_cast<chrono::seconds>(check - start).count() > 60) {
        cout << "Error: read processing appears very sluggish in " 
          << "one or more threads"
          << ". Suggest sort the BAM file by read name and try again\n"
          << "  e.g. use `samtools collate` or `sambamba sort -n`.\n"
          << "Alternatively, try to run NxtIRF/IRFinder using `n_threads = 1`\n";
        // realizeSpareReads();
        return(-1);
      }

      reads[idx] = bpool.at(i);

      if (reads[idx]->core.flag & 0x904) {
        // If is an unmapped / secondary / supplementary alignment -- discard/overwrite
        cSkippedReads ++;
      }else if (! (reads[idx]->core.flag & 0x1)) {
        // If is a single read -- process it as a single -- then discard/overwrite
        cSingleReads ++;
        nucs_proc = processSingle(reads[idx], mappability_mode);
        totalNucleotides += nucs_proc;
        if(nucs_proc > 0) cReadsProcessed++;
      }else{
        if(idx == 0 && spare_reads->size() == 0) {
          // If BAM is sorted by read name, then we don't need read size, simply use old system
          idx++;
        } else if(
            idx == 1 && spare_reads->size() == 0 && 
            reads[0]->core.l_qname == reads[1]->core.l_qname &&
            (0 == strncmp(
              bam_get_qname(reads[0]), 
              bam_get_qname(reads[1]), 
              reads[1]->core.l_qname)
            )) 
        {
          cPairedReads ++;
          if (reads[0]->core.pos <= reads[1]->core.pos) {
            totalNucleotides += processPair(reads[0], reads[1]);
          } else {
            totalNucleotides += processPair(reads[1], reads[0]);
          }
          cReadsProcessed+=2;
          idx = 0;
        } else {
          // Likely a coordinate sorted BAM file:
          for(unsigned int k = 0; k <= idx; k++) {
            read_name_s = bam_get_qname(reads[k]);
            auto it_read = spare_reads->find(read_name_s);
            
            if(it_read != spare_reads->end()){
              // Process matched read
              cPairedReads ++;
              if (reads[k]->core.tid != it_read->second->core.tid) {
                cChimericReads += 1;
              } else {
                if (reads[k]->core.pos <= it_read->second->core.pos) {    
                  totalNucleotides += processPair(reads[k], it_read->second);
                }else{           
                  totalNucleotides += processPair(it_read->second, reads[k]);
                }
                cReadsProcessed+=2;
                
                bam_destroy1(it_read->second);
                spare_reads->erase(read_name_s);
              }
            } else {
              // Bank unmatched read
              b = bam_init1();
              spare_reads->insert({read_name_s, bam_copy1(b, reads[k])});
            }
          }
          idx = 0;
        }
      }

      if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
        // Clean map by swapping for a new one
        new_spare_reads = new std::map< std::string, bam1_t* >;
        new_spare_reads->insert(spare_reads->begin(), spare_reads->end());
        spare_reads->swap(*new_spare_reads);
        delete new_spare_reads;
      }

    }
  }


  if(idx == 1 && spare_reads->size() == 0) {
    read_name_s = bam_get_qname(reads[0]);
    b = bam_init1();
    spare_reads->insert({read_name_s, bam_copy1(b, reads[0])});
  }
  cErrorReads = spare_reads->size();

  return(0);
}

void htsBAM2blocks::registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback ) {
  callbacksChrMappingChange.push_back(callback);
}

void htsBAM2blocks::registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback ) {  
  callbacksProcessBlocks.push_back(callback);
}


// Benchmark functions

// Main function
int htsBAM2blocks::processNothing(
  std::vector<bam1_t*> & bpool, 
  int starts, int ends,
  bool mappability_mode
) {
  // Reads from pbam_in until finished; do not create output
  unsigned int idx = 0;
  unsigned int nucs_proc = 0;
  
  // Use map pointer spare_reads:
  std::map< std::string, bam1_t* > * new_spare_reads;  

  bam1_t * b = NULL;
  std::string read_name_s;
  
  auto start = chrono::steady_clock::now();
  auto check = start;

  if(starts > -1) {
    for(int i = starts; i < ends; i++) {
      check = chrono::steady_clock::now();
      if(chrono::duration_cast<chrono::seconds>(check - start).count() > 60) {
        cout << "Error: read processing appears very sluggish in " 
          << "one or more threads"
          << ". Suggest sort the BAM file by read name and try again\n"
          << "  e.g. use `samtools collate` or `sambamba sort -n`.\n"
          << "Alternatively, try to run NxtIRF/IRFinder using `n_threads = 1`\n";
        // realizeSpareReads();
        return(-1);
      }

      // reads[idx] = bpool.at(i);
    }
  }

  return(0);
}

// Main function
int htsBAM2blocks::processTask1(
  std::vector<bam1_t*> & bpool, 
  int starts, int ends,
  bool mappability_mode
) {
  // Reads from pbam_in until finished; do not create output
  unsigned int idx = 0;
  unsigned int nucs_proc = 0;
  
  // Use map pointer spare_reads:
  std::map< std::string, bam1_t* > * new_spare_reads;  

  bam1_t * b = NULL;
  std::string read_name_s;
  
  auto start = chrono::steady_clock::now();
  auto check = start;

  if(starts > -1) {
    for(int i = starts; i < ends; i++) {
      check = chrono::steady_clock::now();
      if(chrono::duration_cast<chrono::seconds>(check - start).count() > 60) {
        cout << "Error: read processing appears very sluggish in " 
          << "one or more threads"
          << ". Suggest sort the BAM file by read name and try again\n"
          << "  e.g. use `samtools collate` or `sambamba sort -n`.\n"
          << "Alternatively, try to run NxtIRF/IRFinder using `n_threads = 1`\n";
        // realizeSpareReads();
        return(-1);
      }

      reads[idx] = bpool.at(i);

      if (reads[idx]->core.flag & 0x904) {
        // If is an unmapped / secondary / supplementary alignment -- discard/overwrite
        cSkippedReads ++;
      }else if (! (reads[idx]->core.flag & 0x1)) {
        // If is a single read -- process it as a single -- then discard/overwrite
        cSingleReads ++;
        // nucs_proc = processSingle(reads[idx], mappability_mode);
        totalNucleotides += nucs_proc;
        if(nucs_proc > 0) cReadsProcessed++;
      }else{
        if(idx == 0 && spare_reads->size() == 0) {
          // If BAM is sorted by read name, then we don't need read size, simply use old system
          idx++;
        } else if(
            idx == 1 && spare_reads->size() == 0 && 
            reads[0]->core.l_qname == reads[1]->core.l_qname &&
            (0 == strncmp(
              bam_get_qname(reads[0]), 
              bam_get_qname(reads[1]), 
              reads[1]->core.l_qname)
            )) 
        {
          cPairedReads ++;
          if (reads[0]->core.pos <= reads[1]->core.pos) {
            // totalNucleotides += processPair(reads[0], reads[1]);
          } else {
            // totalNucleotides += processPair(reads[1], reads[0]);
          }
          cReadsProcessed+=2;
          idx = 0;
        } else { /*
          // Likely a coordinate sorted BAM file:
          for(unsigned int k = 0; k <= idx; k++) {
            read_name_s = bam_get_qname(reads[k]);
            auto it_read = spare_reads->find(read_name_s);
            
            if(it_read != spare_reads->end()){
              // Process matched read
              cPairedReads ++;
              if (reads[k]->core.tid != it_read->second->core.tid) {
                cChimericReads += 1;
              } else {
                if (reads[k]->core.pos <= it_read->second->core.pos) {    
                  // totalNucleotides += processPair(reads[k], it_read->second);
                }else{           
                  totalNucleotides += processPair(it_read->second, reads[k]);
                }
                cReadsProcessed+=2;
                
                bam_destroy1(it_read->second);
                spare_reads->erase(read_name_s);
              }
            } else {
              // Bank unmatched read
              b = bam_init1();
              spare_reads->insert({read_name_s, bam_copy1(b, reads[k])});
            }
          }
          idx = 0;
        */}
      }
/*
      if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
        // Clean map by swapping for a new one
        new_spare_reads = new std::map< std::string, bam1_t* >;
        new_spare_reads->insert(spare_reads->begin(), spare_reads->end());
        spare_reads->swap(*new_spare_reads);
        delete new_spare_reads;
      }
*/
    }
  }

  return(0);
}