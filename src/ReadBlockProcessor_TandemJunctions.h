/* ReadBlockProcessor_TandemJunctions.h Lists tandem junctions in output

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

#ifndef CODE_READBLOCKPROCESSOR_TANDEMJUNCTION
#define CODE_READBLOCKPROCESSOR_TANDEMJUNCTION

#include "ReadBlockProcessor.h"

// tandemJn: a class that stores start1, end1, start2, end2
class tandemJn {
  public:
    uint32_t start1;
    uint32_t end1;
    uint32_t start2;
    uint32_t end2;
    
    tandemJn(
      uint32_t b, uint32_t c,
      uint32_t d, uint32_t e
    ) {
      start1 = b;
      end1 = c;
      start2 = d;
      end2 = e;
    };
};

// Sort a vector of tandemJn by the following
inline bool operator< (const tandemJn& lhs, const tandemJn& rhs){
  if(
    lhs.start1 == rhs.start1 &&
    lhs.end1 == rhs.end1 &&
    lhs.start2 == rhs.start2
  ) {
    return lhs.end2 < rhs.end2;
  } else if(
    lhs.start1 == rhs.start1 &&
    lhs.end1 == rhs.end1
  ) {
    return lhs.start2 < rhs.start2;
  } else if(
    lhs.start1 == rhs.start1
  ) {
    return lhs.end1 < rhs.end1;
  } else {
    return lhs.start1 < rhs.start1;
  }
}

class TandemJunctions : public ReadBlockProcessor {
	private:
		std::map<string, std::map<tandemJn, unsigned int[3]>> chrName_tandemJn;
		std::vector<std::map<tandemJn, unsigned int[3]>*> chrID_tandemJn;
		//unsigned int[3] - 0, neg strand count; 1, pos strand count; 2 = expected direction from ref: 0=unknown, 1=neg, 2=pos.
	public:
    TandemJunctions(std::string &refString);
    void Reset();
		void Combine(const TandemJunctions &child);
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		int WriteOutput(std::string& output, std::string& QC) const;
		void loadRef(std::istringstream &IN);
};

#endif
