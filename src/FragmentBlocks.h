/* FragmentBlocks.h FragmentBlocks

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

#ifndef CODE_FRAGMENTBLOCKS
#define CODE_FRAGMENTBLOCKS

#include "includedefine.h"

/* A class to store up to 2 reads belonging to a single fragment.
 * It is a storage class, almost a struct, it does not perform processing itself.
 * Read1 is always valid.
 * Read2 is only valid if readCount == 2.
 *
 * There may only be a single read if:
 *  - the sequencing is single end rather than paired end..
 *  - the sequencing is paired end, but the two reads overlapped and have been combined
 *    into a single synthetic read / block of coverage.
 */
class FragmentBlocks {
	private:
		static const int initial_alloc = 100;
		static const int max_read_name = 300;
		std::vector<std::string> chr_names; //TODO - this is currently unused??
	public:
		FragmentBlocks();
		const std::string chrName() const;
		void ChrMapUpdate(const std::vector<std::string>& chrmap);

		std::string readName;
		std::vector<int> rStarts[2];
		std::vector<int> rLens[2];
		unsigned int readStart[2];
		unsigned int readEnd[2];
		int readCount;
		unsigned int chr_id; // Assumption that both r1 & r2 are on the same chromosome?
					//   if they aren't we shouldn't process them as a single fragment.
					//   perhaps a sanity check in pairing, only treat them as a pair
					//   if the name of the reads matches and the Chr matches.
		bool direction;
};

#endif