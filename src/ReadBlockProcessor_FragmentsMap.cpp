/* ReadBlockProcessor_FragmentsMap.cpp Reads fragment coverage

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

#include "ReadBlockProcessor_FragmentsMap.h"

void FragmentsMap::Reset() {
  chrs.resize(0);
  for(unsigned int i = 0; i < 3; i++) {
    chrName_vec_final[i].resize(0);
    chrName_vec_new[i].resize(0);
    temp_chrName_vec_new[i].resize(0);
  }
  frag_count = 0;
  final_is_sorted = false;
}

void FragmentsMap::ChrMapUpdate(const std::vector<chr_entry> &chrmap) {
  std::vector< std::pair<unsigned int, int> > empty_vector;
  empty_vector.push_back(std::make_pair (0,0));
  for(unsigned int j = 0; j < 3; j++) {   
    chrName_vec_final[j].resize(0);
    chrName_vec_new[j].resize(0);
    temp_chrName_vec_new[j].resize(0);
    for (unsigned int i = 0; i < chrmap.size(); i++) {
      chrName_vec_final[j].push_back(empty_vector);
      chrName_vec_new[j].push_back(empty_vector);
      temp_chrName_vec_new[j].push_back(empty_vector);
    }
  }

  for (unsigned int i = 0; i < chrmap.size(); i++) {
    chrs.push_back(chrmap.at(i));
  }
}

void FragmentsMap::ProcessBlocks(const FragmentBlocks &blocks) {
  for (int index = 0; index < blocks.readCount; index ++) {
    //Walk each block within each read.
    for (unsigned int j = 0; j < blocks.rLens[index].size(); j++) {
      // Stranded 
      (temp_chrName_vec_new[blocks.direction].at(blocks.chr_id)).push_back(std::make_pair( blocks.readStart[index] + blocks.rStarts[index][j], 1));
      (temp_chrName_vec_new[blocks.direction].at(blocks.chr_id)).push_back(std::make_pair( blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j], -1));

      // Unstranded 
      (temp_chrName_vec_new[2].at(blocks.chr_id)).push_back(std::make_pair( blocks.readStart[index] + blocks.rStarts[index][j], 1));
      (temp_chrName_vec_new[2].at(blocks.chr_id)).push_back(std::make_pair( blocks.readStart[index] + blocks.rStarts[index][j] + blocks.rLens[index][j], -1));
    }
  }
  frag_count += 1;
  if(frag_count % 1000000 == 0) {
    sort_and_collapse_temp();
  }
}

// Temporarily sorts the nested vector to reduce memory use; occurs every 1M reads
int FragmentsMap::sort_and_collapse_temp() {
  // Sort temp vectors and append to final:
  for(unsigned int j = 0; j < 3; j++) {
    unsigned int refID = 0;
    for (auto itChr=temp_chrName_vec_new[j].begin(); itChr!=temp_chrName_vec_new[j].end(); itChr++) {
      // sort
      if(itChr->size() > 0) {
        std::sort(
          itChr->begin(),
          itChr->end()
        );

        unsigned int loci = 0;
        int accum = 0;
        for(auto it_pos = itChr->begin(); it_pos != itChr->end(); it_pos++) {
          if(it_pos->first != loci) {
            if(accum != 0) chrName_vec_new[j].at(refID).push_back( std::make_pair(loci, accum) );
            loci = it_pos->first;
            accum = it_pos->second;
          } else {
            accum += it_pos->second;
          }
        }
        // final push
        chrName_vec_new[j].at(refID).push_back( std::make_pair(loci, accum) );

        // Clear temporary vector by swap trick
        // empty swap vector
        std::vector< std::pair<unsigned int, int> > empty_swap_vector;
        itChr->swap(empty_swap_vector);
      }      
      refID++;
    }
  }
  return(0);
}

// Final sorting. Also converts from loci/diff to loci/depth
int FragmentsMap::sort_and_collapse_final(
    bool verbose, unsigned int n_threads_to_use
) {
  if(!final_is_sorted) {
    sort_and_collapse_temp();
    if(verbose)  cout << "Performing final sort of fragment maps\n";

#ifdef _OPENMP
    #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
#endif
    for(unsigned int k = 0; k < 3 * chrs.size(); k++) {
      unsigned int j = k / chrs.size();
      unsigned int i = k - (j * chrs.size());
      
        auto itChr = &chrName_vec_new[j].at(i);
        auto itDest = &chrName_vec_final[j].at(i);
        itDest->resize(0);
        // sort
        std::sort(
          itChr->begin(),
          itChr->end()
        );
        
        // Progressors
        unsigned int   loci = 0;       // Current genomic coordinate
        unsigned int   old_loci = 0;       // Current genomic coordinate
        int           depth = 0;       // Current depth of cursor
        int           old_depth = 0;  // Previous depth of cursor
        
        for(auto it_pos = itChr->begin(); it_pos != itChr->end(); it_pos++) {
          if(it_pos->first != loci) {
            if(depth != old_depth) {  
              itDest->push_back( std::make_pair(old_loci, old_depth) );
              old_depth = depth;
              old_loci = loci;
            }
            loci = it_pos->first;
          }

          depth += it_pos->second;

          if(it_pos->first == 0) {
            old_depth = depth;  // ensure never trigger write when first time 
          }       
        }
        itDest->push_back( std::make_pair(old_loci, old_depth) );
        if(depth != old_depth) {
          itDest->push_back( std::make_pair(loci, depth) );
        }
        itChr->clear();
        
    }
    final_is_sorted = true;
  }
  return(0);
}

void FragmentsMap::Combine(FragmentsMap &child) {
  sort_and_collapse_temp();
  child.sort_and_collapse_temp();
  if(!final_is_sorted && !child.final_is_sorted) {
    for(unsigned int j = 0; j < 3; j++) {
      for(unsigned int i = 0; i < chrs.size(); i++) {
        chrName_vec_new[j].at(i).insert(chrName_vec_new[j].at(i).end(),
          child.chrName_vec_new[j].at(i).begin(), child.chrName_vec_new[j].at(i).end());
      }
    }
  } else if(final_is_sorted && child.final_is_sorted) {
    for(unsigned int j = 0; j < 3; j++) {
      for(unsigned int i = 0; i < chrs.size(); i++) {
        chrName_vec_final[j].at(i).insert(chrName_vec_final[j].at(i).end(),
          child.chrName_vec_final[j].at(i).begin(), child.chrName_vec_final[j].at(i).end());
      }
    }
    final_is_sorted = false;  // will request resort but keep incremental = false
  }
}

// updateCoverageHist from completed FragmentMap - directional:
void FragmentsMap::updateCoverageHist(std::map<unsigned int,unsigned int> &hist, unsigned int start, unsigned int end, unsigned int dir, const unsigned int &refID, bool debug) const {
  (void)(debug);
  
  if(refID >= chrName_vec_final[dir].size()) {
    hist.insert({0,0});
    return;
  }
  
  auto it_chr = &chrName_vec_final[dir].at(refID);
  auto it_pos = upper_bound(
      it_chr->begin(), 
      it_chr->end(), 
      make_pair(start, 0), 
      []( std::pair<unsigned int, int> const& a, std::pair<unsigned int, int> const& b ) { 
        return a.first < b.first; 
      }
  );
  
  if(it_pos == it_chr->end()) {
    // No coverage data
    hist[0] += end - start;
    return;
  }
  while(it_pos->first > start && it_pos != it_chr->begin()) {
    it_pos--; // shouldn't matter as the first vector pair should be at coord zero
  }
  int depth = it_pos->second;
  unsigned int cursor = start;
  while(cursor < end) {
    while(it_pos->first <= cursor && it_pos != it_chr->end()) {
      it_pos++;
    }
    if(it_pos == it_chr->end()) {
      hist[(unsigned int)depth] += end - cursor;
      break;
    }
    hist[(unsigned int)depth] += min(it_pos->first, end) - cursor;
    cursor = it_pos->first;
    depth = it_pos->second;
  }
}

int FragmentsMap::WriteBinary(
    covWriter *os, bool verbose, unsigned int n_threads_to_use
) {
  if(!final_is_sorted) {
    // Perform this separately as this is now multi-threaded
    sort_and_collapse_final(verbose, n_threads_to_use);   
  }
  if(verbose)  cout << "Writing COV file\n";

  os->InitializeCOV(chrs);

#ifdef SPLICEWIZ
  Progress p(3 * chrs.size(), verbose);
#endif
  for(unsigned int j = 0; j < 3; j++) {
    for(unsigned int i = 0; i < chrs.size(); i++) {
      unsigned int refID = chrs[i].refID;
      
      std::vector< std::pair<unsigned int, int> > * itDest;
      itDest = &chrName_vec_final[j].at(refID);
      
      os->WriteFragmentsMap(itDest, i, j, n_threads_to_use);
#ifdef SPLICEWIZ
      p.increment(1);
#endif
    }
  }
  
  os->WriteToFile();
  return(0);
}

int FragmentsMap::WriteOutput(std::ostream *os, 
    int threshold, bool verbose, unsigned int n_threads_to_use)  {

  // This is called on mappability

  std::vector<std::string> sort_chr_names;
  std::vector<int32_t> sort_chr_lens;
  for (auto chr = chrs.begin(); chr != chrs.end(); chr++) {
    sort_chr_names.push_back(chr->chr_name);
    sort_chr_lens.push_back(chr->chr_len);
  }

  unsigned int refID = 0;
  if(!final_is_sorted) {
    sort_and_collapse_final(verbose, n_threads_to_use);
  }
  if(verbose)  cout << "Writing Mappability Exclusions\n";
#ifdef SPLICEWIZ
  Progress p(sort_chr_names.size(), verbose);
#endif
  for(unsigned int i = 0; i < sort_chr_names.size(); i++) {
    // refID is reference ID as appears in BAM file; i is the nth chromosome as ordered in alpha order
    refID = chrs[i].refID;
    auto itChr = &chrName_vec_final[2].at(refID);
    int coverage = 0;
    bool covered = false;
    
    if (itChr->begin()->first == 0 && itChr->begin()->second > threshold) {
      covered = true;
    } else {
      // Write first coordinate
      *os << chrs[i].chr_name << "\t0\t";
    }
    for(auto it_pos = itChr->begin(); it_pos != itChr->end(); it_pos++) {
      coverage = it_pos->second;
      if(coverage > threshold) {
        if(!covered) {
          *os << it_pos->first << '\n';
          covered = true;
        }
      } else {
        if(covered) {
          *os << chrs[i].chr_name << "\t"
              << it_pos->first << "\t";
          covered = false;
        }
      }
    }
    // Write last entry
    if(!covered) {
      *os << chrs[i].chr_len << "\n";    
    }
#ifdef SPLICEWIZ
    p.increment(1);
#endif
  }
  return 0;
}
