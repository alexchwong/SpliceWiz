/* FastaReader.cpp Reads FASTA files

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

#include "FastaReader.h"

FastaReader::FastaReader() {
  total_size = 0;
  IN = NULL;
  FirstSeq = true;
}

void FastaReader::SetInputHandle(std::istream *in_stream) {
  IN = in_stream;
  FirstSeq = true;
}

void FastaReader::Profile() {
  // Profiles the FASTQ file to know the chr_names and chr_lens
  chr_names.clear();
  chr_lens.clear();
  IN->seekg (0, std::ios_base::beg);
  total_size = 0;
  
  while(!IN->eof() && !IN->fail()) {
    ReadSeq();
    chr_names.push_back(seqname);
    chr_lens.push_back((int32_t)sequence.length());
    total_size += sequence.length();
  }
  IN->clear();
  IN->seekg (0, std::ios_base::beg);
  FirstSeq = true;
}

bool FastaReader::ReadSeq() {
  std::string myLine;
  std::string sequence_raw;
  std::string line;
  std::string subline;
  std::string subline2;
  
  sequence.clear();
  if(FirstSeq) {
    std::getline(*IN, myLine, '>');
    FirstSeq = false;
  }
  std::getline(*IN, myLine, '\n');
  // ensure tabs and spaces are removed
  std::istringstream iss;
  iss.str(myLine);
  std::getline(iss, myLine, '\r');  // Remove Win-based \r from \r\n
  std::istringstream iss2;
  iss2.str(myLine);
  std::getline(iss2, myLine, '\t');  // chromosome name is tab-terminated
  std::istringstream iss3;
  iss3.str(myLine);
  std::getline(iss3, seqname, ' ');  // chromosome name is space-terminated
  
  std::getline(*IN, sequence_raw, '>');

  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), ' '), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\r'), sequence_raw.end() );
  sequence_raw.erase( std::remove(sequence_raw.begin(), sequence_raw.end(), '\n'), sequence_raw.end() );
  sequence.append(sequence_raw);
  return(true);
}