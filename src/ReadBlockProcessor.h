/* ReadBlockProcessor.h Reads blocks

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

#ifndef CODE_READBLOCKPROCESSOR
#define CODE_READBLOCKPROCESSOR

#include "includedefine.h"
#include "FragmentBlocks.h"

/*
The code can be finished faster if we force a requirement that all input files are coordinate sorted by the start of each block.
ie: sort -k2,2n (for BED files).
Chromosome sorted or not won't matter, as these get split into different vectors in all cases.
*/



class ReadBlockProcessor {
	public:
    virtual ~ReadBlockProcessor() {}; // do nothing
		virtual void ProcessBlocks(const FragmentBlocks &fragblock) = 0;
		virtual void ChrMapUpdate(const std::vector<chr_entry> &chrmap) = 0; //Maybe some of these funcs shouldn't be pure virtual - overloadable if needed, but default often ok.
};

class JunctionCount : public ReadBlockProcessor {
	private:
		std::map<string, std::map<std::pair<unsigned int,unsigned int>,unsigned int[3]>> chrName_junc_count;
		std::vector<std::map<std::pair<unsigned int,unsigned int>,unsigned int[3]>*> chrID_junc_count;
		//unsigned int[3] - 0, neg strand count; 1, pos strand count; 2 = expected direction from ref: 0=unknown, 1=neg, 2=pos.

		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncLeft_count;
		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncLeft_count;

		std::map<string, std::map<unsigned int,unsigned int[2]>> chrName_juncRight_count;
		std::vector<std::map<unsigned int,unsigned int[2]>*> chrID_juncRight_count;
		  //chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:
    JunctionCount(std::string &refString);
		void Combine(const JunctionCount &child);
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		int WriteOutput(std::string& output, std::string& QC) const;
		void loadRef(std::istringstream &IN); //loadRef is optional, it allows directional detection to determine not just non-dir vs dir, but also which direction.

		int Directional(std::string& output) const;
		
		unsigned int lookup(std::string ChrName, unsigned int left, unsigned int right, bool direction) const;
		unsigned int lookup(std::string ChrName, unsigned int left, unsigned int right) const;
		unsigned int lookupLeft(std::string ChrName, unsigned int left, bool direction) const;
		unsigned int lookupLeft(std::string ChrName, unsigned int left) const;
		unsigned int lookupRight(std::string ChrName, unsigned int right, bool direction) const;
		unsigned int lookupRight(std::string ChrName, unsigned int right) const;

// Ideally we would read the XS junction strand attribute from the BAM if we want to count junctions from non-directional sequencing.
//   that will require BAM2blocks to be informed it should read the optional attributes looking for that attrib in that case.
// -- or we can just ignore direction -- the splice start/end information effectively determines the XS info (by ref to the reference)
};


class SpansPoint : public ReadBlockProcessor {
	private:
		std::map<string, std::vector<unsigned int>> chrName_pos;
		std::map<string, std::vector<unsigned int>> chrName_count[2];
		std::vector<std::vector<unsigned int>*> chrID_pos;
		std::vector<std::vector<unsigned int>*> chrID_count[2];
		char overhangLeft;
		char overhangRight;
		char overhangTotal;
		//chrID_... stores a fast access pointer to the appropriate structure in chrName_... 
	public:
    SpansPoint(std::string &refString);
		void Combine(const SpansPoint &child);
		void setSpanLength(unsigned int overhang_left, unsigned int overhang_right);
		void loadRef(std::istringstream &IN);
		void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		//void SetOutputStream(std::ostream *os);
		int WriteOutput(std::string& output, std::string& QC) const;
		unsigned int lookup(std::string ChrName, unsigned int pos, bool direction) const;
		unsigned int lookup(std::string ChrName, unsigned int pos) const;
};

class FragmentsInChr : public ReadBlockProcessor {
	// Counts the number of fragments in each Chromosome. (for both + & - strands).
	private:
		std::map<string, std::vector<unsigned int>> chrName_count; //only expecting 2 items in our vector.
		std::vector<std::vector<unsigned int>*> chrID_count;
	public:
		void Combine(const FragmentsInChr &child);
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		int WriteOutput(std::string& output, std::string& QC) const;		
};


class FragmentsInROI : public ReadBlockProcessor {
	// Counts the number of fragments fully contained within a ROI.
	//   the ROIs may not overlap. Direction ignored for overlap detect.
	private:
		std::map<string, unsigned long> RegionID_counter[2];
 
		std::map<string, std::vector<std::pair<unsigned int,unsigned int>>> chrName_ROI;
		std::map<string, std::vector<unsigned long*>> chrName_count[2];

		std::vector<std::vector<std::pair<unsigned int,unsigned int>>*> chrID_ROI;
		std::vector<std::vector<unsigned long*>*> chrID_count[2];

		// Perhaps we want to store some text relating to each record too? Easy to do if the input is pre-sorted (at least within each Chr).
		//   if pre-sorted, it may be easier to check for no overlapping blocks on read .. or can do this immediately after read with a single nested-walk.
		std::map<string, std::vector<string>> chrName_ROI_text;
	public:
    FragmentsInROI(std::string &refString);
		void Combine(const FragmentsInROI &child);
		void ProcessBlocks(const FragmentBlocks &blocks);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		void loadRef(std::istringstream &IN);
		int WriteOutput(std::string& output, std::string& QC) const;		
};


/*
class CoverageBlocks : public ReadBlockProcessor { ... }
// In it's own file -- bigger code.
*/

#endif
