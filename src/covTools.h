/* covTools.h Reads and writes COV files (BAM coverage)

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


#ifndef _CODE_COVTOOLS
#define _CODE_COVTOOLS

#include "includedefine.h"
#include "SpliceWiz.h"

#include <zlib.h>
#include <zconf.h>

#include "pbam_defs.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

union stream_uint64 {
  char c[8];
  uint64_t u;
};
union stream_uint32 {
  char c[4];
  uint32_t u;
};
union stream_int32 {
  char c[4];
  int32_t i;
};
union stream_uint16 {
  char c[2];
  uint16_t u;
};

static const unsigned int BGZF_max = 65536 - 18 - 8;

class covReader {
	private:
    // Buffers
    char * compressed_buffer;
    char * buffer;

    unsigned long bufferPos;  // Position of decompressed buffer
    unsigned long bufferMax;  // Size of decompressed buffer
    
    uint32_t index_begin;     // File position of first byte of COV index
    uint32_t body_begin;      // File position of first byte of COV body
    
    istream * IN;
    
    int IS_EOF;               // Set to 1 if istream hits eof()
    int IS_FAIL;              // Set to 1 if istream hits fail()
    
    size_t IS_LENGTH;         // Size of COV file
    size_t EOF_POS;           // Position of first byte of BGZF EOF block
       
    std::vector<std::string> chr_names;   // seqnames
    std::vector<uint32_t> chr_lens;       // chromosome lengths

  public:
    covReader();
    ~covReader();
    void SetInputHandle(std::istream *in_stream);
    
    int ReadBuffer();
    int read(char * dest, unsigned int len);
    int ignore(unsigned int len);
    bool eof();
    bool fail();

    // Input functions
    
    int ReadHeader();
    int GetChrs(std::vector<chr_entry> &chrs);
    int FetchPos(const std::string seqname, const uint32_t start, const int strand,
      uint64_t * file_offset, uint32_t * block_start);
    int FetchRLE(const std::string seqname, 
      const uint32_t start, const uint32_t end, const int strand,
      std::vector<int> * values, std::vector<unsigned int> * lengths
    );
};


// Body buffer for CovWriter
class buffer_out_chunk {
  private:
    char * buffer = NULL;
    char * compressed_buffer = NULL;
    
    // current write position of buffer
    unsigned int buffer_pos = 0;
    // number of bytes that need to be compressed
    unsigned int buffer_size = 0;   
    // number of bytes compressed that need to be written out
    unsigned int compressed_size = 0;   
    
  public:
    buffer_out_chunk();
    ~buffer_out_chunk();


    unsigned int write(char * src, unsigned int len);
    int WriteToFile(ostream * OUT);
    int Compress();
    
    unsigned int getBGZFSize() { return(compressed_size); };

    unsigned int SetPos(unsigned int pos) {
      if(pos >= BGZF_max) return(buffer_pos);
      buffer_pos = pos;
      if(pos > buffer_size) {
        buffer_size = pos;
      }
      return(pos);
    };

    unsigned int GetPos() { return(buffer_pos); };
    
    unsigned int write_to_pos(char * src, unsigned int len, unsigned int pos) {
      if(len + pos > BGZF_max) return(0);
      SetPos(pos);
      return(write(src, len));
    };
    
    bool IsAtCap(unsigned int len) {
      if(len + buffer_pos >= BGZF_max) return(true);
      return(false);
    }
};

class covWriter {
  private:
    ostream * OUT;
    
    // When chrs is set, initialize these:
    std::vector<chr_entry> chrs;
    
    // The buffers
    std::vector< std::vector<buffer_out_chunk> > body;          
    // The start coords of each bgzf
    std::vector< std::vector<uint32_t> > block_coord_starts;    

    int WriteEmptyEntry(unsigned int refID);
    int WriteHeaderToFile();
    int WriteIndexToFile();
  public:
    covWriter();
    ~covWriter();
    
    void SetOutputHandle(std::ostream *out_stream);
    
    int InitializeCOV(std::vector<chr_entry> chrs_to_copy);
  
    int WriteFragmentsMap(
      std::vector< std::pair<unsigned int, int> > * vec, 
      unsigned int chrID, unsigned int strand,
      unsigned int n_threads_to_use = 1
    );
    
    int WriteToFile();

};

#endif