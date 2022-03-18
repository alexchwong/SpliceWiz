/* FastaReader.h Reads FASTA files

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

#ifndef CODE_FASTAREADER
#define CODE_FASTAREADER

#include "includedefine.h"

class FastaReader {
  private:
    istream * IN;
    bool FirstSeq;
  public:
    FastaReader();
  
    std::string seqname;
    std::string sequence;
    
    void SetInputHandle(std::istream *in_stream);
    void Profile();
    bool ReadSeq();

    // tab or space terminated chromosome names.
    std::vector<std::string> chr_names;
    // length of each chromosome
    // - not used when reading, used if optionally outputting an altered BAM file)
    std::vector<int32_t> chr_lens;
    size_t total_size;
};

#endif