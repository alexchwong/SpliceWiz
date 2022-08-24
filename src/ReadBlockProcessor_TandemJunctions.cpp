
/* ReadBlockProcessor_TandemJunctions.cpp Records tandem junctions

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

#include "ReadBlockProcessor_TandemJunctions.h"

// Constructors from strings:

TandemJunctions::TandemJunctions(std::string &refString) {
  std::istringstream inTandemJn;
  inTandemJn.str(refString);
  loadRef(inTandemJn);
}

void TandemJunctions::loadRef(std::istringstream &IN) {
  // ChrName, Start1, End1, Start2, End2
  std::string myLine;
  std::string myField;
  myLine.reserve(1000);
  myField.reserve(100);

  // Text inputs
  string s_chr;
  s_chr.reserve(100);
  
  uint32_t start1;
  uint32_t end1;
  uint32_t start2;
  uint32_t end2;
  
  string direction;
  
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
    
    getline(lineStream, s_chr, '\t');
    getline(lineStream, myField, '\t');
    start1 = stol(myField);
    getline(lineStream, myField, '\t');
    end1 = stol(myField);
    getline(lineStream, myField, '\t');
    start2 = stol(myField);
    getline(lineStream, myField, '\t');
    end2 = stol(myField);
    getline(lineStream, direction, '\t');

    if (direction == "-")  {
      (chrName_tandemJn)[s_chr][tandemJn(start1,end1,start2,end2)][2] += 1;
    }  else if (direction == "+") {
      (chrName_tandemJn)[s_chr][tandemJn(start1,end1,start2,end2)][2] += 2;
    }
  }
}

//chrName_junc_count holds the data structure -- ChrName(string) -> Junc Start/End -> count.
//chrID_junc_count holds the ChrID -> ...
//  where the ChrID is the ChrID relating to the appropriate ChrName, as understood by the currently processed BAM file.
void TandemJunctions::ChrMapUpdate(const std::vector<chr_entry> &chrmap) {
  chrID_tandemJn.resize(0);

  for (unsigned int i = 0; i < chrmap.size(); i++) {
    chrID_tandemJn.push_back( &(chrName_tandemJn)[chrmap.at(i).chr_name] );
  }
}

void TandemJunctions::ProcessBlocks(const FragmentBlocks &blocks) {
  for (int index = 0; index < blocks.readCount; index ++) {
    //Walk each *pair* of blocks. ie: ignore a read that is just a single block.
    for (unsigned int j = 2; j < blocks.rLens[index].size(); j++) {
      if (
        (blocks.rLens[index][j-2] >= 5) &&
        (blocks.rLens[index][j-1] >= 5) && 
        (blocks.rLens[index][j] >= 5)
      ) {
        // Check if this tandem junction exists
        tandemJn tmpJn(
          blocks.readStart[index] + blocks.rStarts[index][j-2] + blocks.rLens[index][j-2],
          blocks.readStart[index] + blocks.rStarts[index][j-1],
          blocks.readStart[index] + blocks.rStarts[index][j-1] + blocks.rLens[index][j-1],
          blocks.readStart[index] + blocks.rStarts[index][j]
        );
        (*chrID_tandemJn[blocks.chr_id])[tmpJn][blocks.direction]++;
      }
    }
  }
}

void TandemJunctions::Combine(const TandemJunctions &child) {
  for(unsigned int j = 0; j < 2; j++) {
    for(unsigned int i = 0; i < chrName_tandemJn.size(); i++) {
      auto itChr1=child.chrName_tandemJn.begin();
      for(unsigned int k = 0; k < i; k++) itChr1++;
      for (auto itPos = itChr1->second.begin(); itPos != itChr1->second.end(); itPos++) {
        // itPos is iterator to std::map<tandemJn, unsigned int[3]> object
        // Insert missing entries from child:
        chrName_tandemJn.at(itChr1->first)[itPos->first][j] += itPos->second[j];
      }
    }
  }
}

int TandemJunctions::WriteOutput(std::string& output, std::string& QC) const {
  std::ostringstream oss; std::ostringstream oss_qc; 

  for (auto itChr=chrName_tandemJn.begin(); itChr!=chrName_tandemJn.end(); itChr++) {
    string chr = itChr->first;
    for (auto itJuncs=itChr->second.begin(); itJuncs!=itChr->second.end(); ++itJuncs) {
      oss << chr << "\t" << 
        itJuncs->first.start1 << "\t" << 
        itJuncs->first.end1 << "\t" << 
        itJuncs->first.start2 << "\t" << 
        itJuncs->first.end2 << "\t" << 
        ( (itJuncs->second)[2] & 1 ? "-" : (itJuncs->second)[2] & 2 ? "+" : "." ) << "\t" << 
        ((itJuncs->second)[1] + (itJuncs->second)[0]) << "\t" << 
        (itJuncs->second)[1] << "\t" << (itJuncs->second)[0] << "\n";
    }
  }
  oss_qc << "";
  
  output = oss.str();
  QC.append(oss_qc.str());
  return 0;
}

