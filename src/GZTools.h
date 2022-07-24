/* GZTools.h Reads Gzipped files

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

#ifndef CODE_GZREADER
#define CODE_GZREADER

#include "includedefine.h"
#include "SpliceWiz.h"

#include <zlib.h>
#include <zconf.h>

#define CHUNK_gz 262144

class GZReader {
  private:
    gzFile gz_in = NULL;
    int GetBuffer();

    char * buffer = NULL;            // stores data
    unsigned long bufferLen = 0;  // amount read from file
    unsigned long bufferPos = 0;  // amount read using getline, read or ignore
    
    bool loaded = false;
    bool lazy = false;
    bool streamed = false;
  public:
    GZReader();
    ~GZReader();
    int LoadGZ(std::string s_filename, bool asStream = false, bool lazymode = false);
    int getline(std::string & s_myLine, const char delim);
    int closeGZ();
    
    void read(char * dest, const unsigned long len);
    void ignore(const unsigned long len);
    bool eof();
    
    // public access to istringstream version of buffer
    std::istringstream iss; 
};

class GZWriter {
private:
  ostream * OUT;
  
  char compressed_buffer[CHUNK_gz];
  
  char buffer[CHUNK_gz];
  unsigned int bufferPos = 0;
  
public:
  void SetOutputHandle(std::ostream *out_stream);

  int writebuffer(const char * src, unsigned int len);
  int writeline(const std::string& s_src);
  int writestring(const std::string& s_src);
  int flush(bool final = false);
};

#endif