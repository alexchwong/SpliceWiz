/* ReadBlockProcessor_CoverageBlocks.h Reads coverage blocks

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

#ifndef CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS
#define CODE_READBLOCKPROCESSOR_COVERAGEBLOCKS

#include "ReadBlockProcessor_FragmentsMap.h"


struct BEDrecord {
	std::string chrName;
	std::string name;
	unsigned int start;
	unsigned int end;
	bool direction;
	
	std::vector<std::pair<unsigned int,unsigned int>> blocks;
};

class CoverageBlocks : public ReadBlockProcessor {
	//Store the Blocked BED record for each ROI/intron. This won't be referred to again until the end.
	//XX Create the temporary vectors (per Chr) which simply list the blocks sequentially as read.
	//XX Sort the temporary vectors
	//XX Build the final vectors of "blocks of interest"
	//xx Delete the temporary vectors
	//xx Create the parallel vectors with counter objects. (do these as a batch at the end, once vector size is known - for best memory layout)
	//xx Process fragments against the counter structure. (have I already written a class/object for this?)
	
	//Produce summary statistical output for each Blocked BED record, using the counter structure.

	private:

  
	protected:
		std::vector<BEDrecord> BEDrecords;

	public:
		void Reset();
    void ProcessBlocks(const FragmentBlocks &fragblock);
		void ChrMapUpdate(const std::vector<chr_entry> &chrmap);
		void loadRef(std::istringstream &IN);
		int WriteOutput(std::string& output, const FragmentsMap &FM) const;
		
	  void fillHist(std::map<unsigned int,unsigned int> &hist, const unsigned int &refID, const std::vector<std::pair<unsigned int,unsigned int>> &blocks, const FragmentsMap &FM, bool debug = false) const;
		void fillHist(std::map<unsigned int,unsigned int> &hist, const unsigned int &refID, const std::vector<std::pair<unsigned int,unsigned int>> &blocks, bool direction, const FragmentsMap &FM, bool debug = false) const;

		double meanFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double coverageFromHist(const std::map<unsigned int,unsigned int> &hist) const;
		double percentileFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int percentile) const;
		double trimmedMeanFromHist(const std::map<unsigned int,unsigned int> &hist, unsigned int centerPercent, bool debug = false) const;

    vector<chr_entry> chrs;
};

class CoverageBlocksIRFinder : public CoverageBlocks {
	public:
    CoverageBlocksIRFinder();
    CoverageBlocksIRFinder(std::string &refString);
    void initialize(std::string &refString);
		void Combine(CoverageBlocksIRFinder &child);
		int WriteOutput(std::string& output, std::string& QC, const JunctionCount &JC, const SpansPoint &SP, const FragmentsMap &FM, int n_threads = 1, int directionality = 0) const;
};


#endif
