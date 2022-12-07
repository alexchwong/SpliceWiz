/* BAM2blocks.cpp Convert reads / fragments to FragmentBlocks

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

#include "BAM2blocks.h"
#include <chrono>

BAM2blocks::BAM2blocks() {
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
  
  spare_reads = new std::map< std::string, pbam1_t* >;
}

BAM2blocks::BAM2blocks(BAM2blocks&& rhs) {
  oBlocks = FragmentBlocks(); //Right syntax to call the default constructor on an object variable, declared but not initialised?

  cReadsProcessed = rhs.cReadsProcessed;
  totalNucleotides = rhs.totalNucleotides;
  cShortPairs = rhs.cShortPairs;
  cIntersectPairs = rhs.cIntersectPairs;
  cLongPairs = rhs.cLongPairs;
  cSingleReads = rhs.cSingleReads;
  cPairedReads = rhs.cPairedReads;
  cErrorReads = rhs.cErrorReads;
  cSkippedReads = rhs.cSkippedReads;
  cChimericReads = rhs.cChimericReads;
  spare_reads = rhs.spare_reads;
  
  rhs.cReadsProcessed = 0;
  rhs.totalNucleotides = 0;
  rhs.cShortPairs = 0;
  rhs.cIntersectPairs = 0;
  rhs.cLongPairs = 0;
  rhs.cSingleReads = 0;
  rhs.cPairedReads = 0;
  rhs.cErrorReads = 0;
  rhs.cSkippedReads = 0;
  rhs.cChimericReads = 0;
  rhs.spare_reads = nullptr;
}

BAM2blocks::BAM2blocks(
  std::vector<std::string> & ref_names, 
  std::vector<uint32_t> & ref_lengths
) {
  initialize(ref_names, ref_lengths);
}

void BAM2blocks::initialize(
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
  
  // spare_reads = new std::map< std::string, pbam1_t* >;
}

BAM2blocks::~BAM2blocks() {
  for(auto it = spare_reads->begin(); it != spare_reads->end(); it++) {
    delete it->second;
  }
  delete spare_reads;
}

unsigned int BAM2blocks::openFile(pbam_in * _IN) {
  // Pass pbam_in object to BB child object
  // Not thread safe!
  
  IN = _IN;
  
  // This block only run on default / empty BB objects
  std::vector<std::string> s_chr_names;
  std::vector<uint32_t> u32_chr_lens;
  IN->obtainChrs(s_chr_names, u32_chr_lens);
  if(chrs.size() == 0) {
    for(unsigned int i = 0; i < s_chr_names.size(); i++) {
      chrs.push_back(chr_entry(i, s_chr_names.at(i), (int32_t)u32_chr_lens.at(i)));
    }
  }

  for (auto & callback : callbacksChrMappingChange ) {
    callback(chrs);
  }
  return(0);
}

void BAM2blocks::cigar2block(uint32_t * cigar, uint16_t n_cigar_op, std::vector<int> &starts, std::vector<int> &lens, int &ret_genome_len) {
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
unsigned int BAM2blocks::processPair(pbam1_t * read1, pbam1_t * read2) {
  // R1 is to the left of R2 (or equal starts).
  int r1_genome_len;
  //int r1_blocks;
  int r2_genome_len;

  pbam1_t * r1 = read1;
  pbam1_t * r2 = read2;
  
  // string debugstate;


  if (r1->flag() & 0x40) {
    //this is first of pair.
    if (r1->flag() & 0x10) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }else{
    if (r1->flag() & 0x20) {
      oBlocks.direction = 0;
    }else{
      oBlocks.direction = 1;
    }
  }

  cigar2block(r1->cigar(), r1->n_cigar_op(), oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  cigar2block(r2->cigar(), r2->n_cigar_op(), oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);

  bool merge_reads = false;
  bool swap_reads = false;
  bool goodPair = true;
  
  if (r1->pos() + r1_genome_len < r2->pos()) {
    cLongPairs++;
    //reads do not intersect
    oBlocks.readCount = 2;
    // debugstate.append( "-Long-");
  }else if (r1->pos() + r1_genome_len >= r2->pos() + r2_genome_len){
    if(r1->pos() == r2->pos() && r1_genome_len > r2_genome_len) {
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
    
    if (r1->flag() & 0x40) {
      //this is first of pair.
      if (r1->flag() & 0x10) {
        oBlocks.direction = 0;
      }else{
        oBlocks.direction = 1;
      }
    }else{
      if (r1->flag() & 0x20) {
        oBlocks.direction = 0;
      }else{
        oBlocks.direction = 1;
      }
    }

    cigar2block(r1->cigar(), r1->n_cigar_op(), oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
    cigar2block(r2->cigar(), r2->n_cigar_op(), oBlocks.rStarts[1], oBlocks.rLens[1], r2_genome_len);
  }
  
  if(merge_reads) {
    // We have two reads that intersect - construct just one complete fragment.

// Guaranteed assumptions:
//   Read 1 starts to the left of Read 2.
//   Read 2 end extends beyond the end of Read 1 end.
    int r1pos = r1->pos();
    int r2pos = r2->pos();
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
  oBlocks.chr_id = r1->refID();
  oBlocks.readStart[0] = r1->pos();
  oBlocks.readEnd[0] = r1->pos() + r1_genome_len;
  oBlocks.readName.resize(r1->l_read_name() - 1);
  oBlocks.readName.replace(0, r1->l_read_name() - 1, r1->read_name(), r1->l_read_name() - 1); // is this memory/speed efficient?

  unsigned int totalBlockLen = 0;
  for (auto blockLen: oBlocks.rLens[0]) {
    totalBlockLen += blockLen;
  }
  if (oBlocks.readCount > 1) {
    oBlocks.readStart[1] = r2->pos();
    oBlocks.readEnd[1] = r2->pos() + r2_genome_len;
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


unsigned int BAM2blocks::processSingle(pbam1_t * read1, bool mappability_mode) {
  int r1_genome_len;

  // string debugstate;

  if (read1->flag() & 0x10) {
    oBlocks.direction = 0;
  }else{
    oBlocks.direction = 1;
  }

  cigar2block(read1->cigar(), read1->n_cigar_op(), oBlocks.rStarts[0], oBlocks.rLens[0], r1_genome_len);
  oBlocks.readCount = 1;

  oBlocks.chr_id = read1->refID();
  oBlocks.readStart[0] = read1->pos();
  oBlocks.readEnd[0] = read1->pos() + r1_genome_len;
  oBlocks.readName.resize(read1->l_read_name() - 1);
  oBlocks.readName.replace(0, read1->l_read_name() - 1, read1->read_name(), read1->l_read_name() - 1); // is this memory/speed efficient?
  
  // Below block only run from Mappability - only process reads if they are
  // mapped to the exact position from which synthetic reads were
  // generated in the genome
  if(mappability_mode) {
    // Return if not a perfect 70M match
    if(read1->n_cigar_op() != 1) return(0);
    if( (*(read1->cigar()) & 15) != 0) return(0);
    
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
int BAM2blocks::WriteOutput(std::string& output) {
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
pbam1_t *  BAM2blocks::SupplyRead(std::string& read_name) {
  // Supplies the pointer to the last spare read
  // When called, transfer ownership of read to the parent BB
  if(spare_reads->size() == 0) return(NULL);
  auto it = spare_reads->begin();
  read_name = it->first;
  pbam1_t * read = it->second;
  spare_reads->erase(it);
  cErrorReads -= 1;
  return(read);
}

// Summates statistics from child BB's
int BAM2blocks::processStats(BAM2blocks& other) {
															   
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
int BAM2blocks::processSpares(BAM2blocks& other) {
  // Combines two BB's, and processes any matching paired reads
  pbam1_t * spare_read;
  std::string read_name;
  while(1) {
    spare_read = other.SupplyRead(read_name);
    
    if(!spare_read) {
      break;
    }
    
    auto it_read = spare_reads->find(read_name);
    if(it_read != spare_reads->end()){
      cPairedReads ++;
      if (spare_read->refID() != it_read->second->refID()) {
        cChimericReads += 1;
      } else {
        if (spare_read->pos() <= it_read->second->pos()) {    
          totalNucleotides += processPair(&(*spare_read), &(*(it_read->second)));
        } else{              
          totalNucleotides += processPair(&(*(it_read->second)), &(*spare_read));
        }
        cReadsProcessed+=2;
      }
      delete (it_read->second);
      spare_reads->erase(read_name);
      delete spare_read;
      cErrorReads-=1;
    } else {
      spare_reads->insert({read_name, spare_read});
    }
  }
  
  return(0);
}

// Saves pbam1_t spare reads to dedicated buffer space
int BAM2blocks::realizeSpareReads() {
  for (auto it = spare_reads->begin(); it != spare_reads->end(); it++) {
    if(!it->second->isReal()) {
      it->second->realize();
    }
  }
  return(0);
}

// Main function
int BAM2blocks::processAll(unsigned int thread_number, bool mappability_mode) {
  // Reads from pbam_in until finished; do not create output
  unsigned int idx = 0;
  unsigned int nucs_proc = 0;

  bool any_reads_processed = false;
  
  // Use map pointer spare_reads:
  std::map< std::string, pbam1_t* > * new_spare_reads;  
  pbam1_t read;
  std::string read_name_s;
  pbam1_t * store_read;
  
  auto start = chrono::steady_clock::now();
  auto check = start;
  while(1) {
    check = chrono::steady_clock::now();
    if(chrono::duration_cast<chrono::seconds>(check - start).count() > 60) {
      cout << "Error: read processing appears very sluggish in thread " << thread_number
        << ". Suggest sort the BAM file by read name and try again\n"
        << "  e.g. use `samtools collate` or `sambamba sort -n`.\n"
        << "Alternatively, try to run NxtIRF/IRFinder using `n_threads = 1`\n";
      realizeSpareReads();
      return(-1);
    }
    
    read = IN->supplyRead(thread_number);
    if(!read.validate()) {
      if(idx == 1 && spare_reads->size() == 0) {
        reads[0].read_name(read_name_s);
        store_read = new pbam1_t;
        *(store_read) = reads[0];
        spare_reads->insert({read_name_s, store_read});
      }
      cErrorReads = spare_reads->size();
      realizeSpareReads();
      if(!any_reads_processed) return(1);
      return(0);   // This will happen if read fails - i.e. end of loaded buffer
    } else {
      any_reads_processed = true;
    }
    reads[idx] = read;

    if (reads[idx].flag() & 0x904) {
      // If is an unmapped / secondary / supplementary alignment -- discard/overwrite
      cSkippedReads ++;
    }else if (! (reads[idx].flag() & 0x1)) {
      // If is a single read -- process it as a single -- then discard/overwrite
      cSingleReads ++;
      nucs_proc = processSingle(&reads[idx], mappability_mode);
      totalNucleotides += nucs_proc;
      if(nucs_proc > 0) cReadsProcessed++;
    }else{
      if(idx == 0 && spare_reads->size() == 0) {
        // If BAM is sorted by read name, then we don't need read size, simply use old system
        idx++;
      } else if(idx == 1 && spare_reads->size() == 0 && 
          reads[0].l_read_name() == reads[1].l_read_name() &&
          (0 == strncmp(reads[0].read_name(), reads[1].read_name(), reads[1].l_read_name()))) {
        cPairedReads ++;
        if (reads[0].pos() <= reads[1].pos()) {
          totalNucleotides += processPair(&reads[0], &reads[1]);
        } else {
          totalNucleotides += processPair(&reads[1], &reads[0]);
        }
        cReadsProcessed+=2;
        idx = 0;
      } else {
        // Likely a coordinate sorted BAM file:
        for(unsigned int k = 0; k <= idx; k++) {
          reads[k].read_name(read_name_s);
          auto it_read = spare_reads->find(read_name_s);
          
          if(it_read != spare_reads->end()){
            // Process matched read
            cPairedReads ++;
            if (reads[k].refID() != it_read->second->refID()) {
              cChimericReads += 1;
            } else {
              if (reads[k].pos() <= it_read->second->pos()) {    
                totalNucleotides += processPair(&reads[k], &(*(it_read->second)));
              }else{           
                totalNucleotides += processPair(&(*(it_read->second)), &reads[k]);
              }
              cReadsProcessed+=2;
              delete (it_read->second);
              spare_reads->erase(read_name_s);
            }
          } else {
            // Bank unmatched read
            store_read = new pbam1_t;
            *(store_read) = reads[k];
            spare_reads->insert({read_name_s, store_read});
          }
        }
        idx = 0;
      }
    }

    if ( (cPairedReads + cSingleReads) % 1000000 == 0 ) {
      // Clean map by swapping for a new one
      new_spare_reads = new std::map< std::string, pbam1_t* >;
      new_spare_reads->insert(spare_reads->begin(), spare_reads->end());
      spare_reads->swap(*new_spare_reads);
      delete new_spare_reads;
    }
  }

  return(0);
}

void BAM2blocks::registerCallbackChrMappingChange( std::function<void(const std::vector<chr_entry> &)> callback ) {
  callbacksChrMappingChange.push_back(callback);
}

void BAM2blocks::registerCallbackProcessBlocks( std::function<void(const FragmentBlocks &)> callback ) {  
  callbacksProcessBlocks.push_back(callback);
}
