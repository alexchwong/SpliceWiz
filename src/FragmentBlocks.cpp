/* FragmentBlocks.cpp FragmentBlocks

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

#include "FragmentBlocks.h"
#include "includedefine.h"
// using namespace std;


// This class is an information storage container only -- pretty much a struct.
// It allows all the relevant information relating to an interpreted fragment to be passed
// to the variety of callback watchers that require fragment blocks to update their stats.

FragmentBlocks::FragmentBlocks() {
	rStarts[0].reserve(initial_alloc);
	rLens[0].reserve(initial_alloc);
	rStarts[1].reserve(initial_alloc);
	rLens[1].reserve(initial_alloc);
	readName.reserve(max_read_name);
	readCount = 0;
}

// Return a string representation of the Chromosome name.
const std::string FragmentBlocks::chrName() const {
	return chr_names.at(chr_id);
}

// Update the internal data structure with a new mapping between Chromosome ID# and Chromosome name (string).
void FragmentBlocks::ChrMapUpdate(const std::vector<string> &chrmap) {
	chr_names = chrmap;
}
