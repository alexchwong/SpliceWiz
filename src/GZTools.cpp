/* GZTools.cpp Reads and writes gzipped files

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

#include "GZTools.h"
  
GZReader::GZReader() {
  bufferLen = 0;
  bufferPos = 0;
  buffer = NULL;
  
  loaded = false; lazy = false; streamed = false;
}

GZReader::~GZReader() {
  if(buffer != NULL) {
    free(buffer);
  }
}

// Only allowed for lazy = true
// Returns -1 if error; 0 if success; 1 if success and EOF
int GZReader::getline(std::string & s_myLine, const char delim) {
  if(!lazy || !loaded || streamed) return(-1);
  int ret = 0;
  unsigned long i = bufferPos;
  while(ret != 1) {
    if(i == bufferLen) {
      ret = GetBuffer();
    }
    while(i < bufferLen) {
      if(buffer[i] == delim) {
        break;
      }
      i++;
    }
    if(i < bufferLen || ret == 1) {
      s_myLine.clear();
      if(i > bufferPos) {
        char * str_ptr = new char[i - bufferPos + 1];
        memcpy(str_ptr, &buffer[bufferPos], i - bufferPos);
        str_ptr[i - bufferPos] = '\0';
        s_myLine.assign(str_ptr, i - bufferPos + 1);
        delete[] str_ptr;
      }
      bufferPos = i + 1;
			return(ret);
    }
  }
	return(ret);
}

// gets a single chunk of data and appends to buffer. Only for lazy = TRUE
int GZReader::GetBuffer() {
  unsigned char *data = NULL;
  int data_alloc = 0;
  int curpos = 0;
  
  int err;
  int bytes_read;
  unsigned char *data_tmp;
  
  data = (unsigned char *)realloc((data_tmp = data), data_alloc += CHUNK_gz - 1);
  bytes_read = gzread (gz_in, data + curpos, CHUNK_gz - 1);
  curpos += bytes_read;
  
  if (bytes_read < CHUNK_gz - 1) {
    if (gzeof (gz_in)) {
      data = (unsigned char *)realloc((data_tmp = data), data_alloc -= (CHUNK_gz - 1) - bytes_read );
    } else {
      const char * error_string;
      error_string = gzerror (gz_in, & err);
      if (err) {
        cout << "Exception during zlib decompression: (" << err << ") " << error_string;
        free(data);
				return(err);
      }
    }
  }

  char *buffer_tmp;
  buffer = (char*)realloc(buffer_tmp = buffer, bufferLen + curpos);
  memcpy(&buffer[bufferLen], data, curpos);

  bufferLen += curpos;
  free(data);
  if (gzeof (gz_in)) {
    return(1);
  } else {
    return(0);
  }
}

// Loads a file
//   Options:
//   - lazy = FALSE: opens file as well as reads entire file into memory
//   - asStream = TRUE: copies read data into istringstream object
int GZReader::LoadGZ(std::string s_filename, bool asStream, bool lazymode) {
  gz_in = gzopen(s_filename.c_str(), "r");
  
  if(lazymode == false) {
    unsigned char *data = NULL;
    int data_alloc = 0;
    int curpos = 0;
    
    while(true) {
      int err = 0;
      int bytes_read = 0;
      unsigned char *data_tmp = NULL;
      
      data = (unsigned char *)realloc((data_tmp = data), data_alloc += CHUNK_gz - 1);
      bytes_read = gzread (gz_in, data + curpos, CHUNK_gz - 1);
      curpos += bytes_read;
      
      if (bytes_read < CHUNK_gz - 1) {
        if (gzeof (gz_in)) {
          data = (unsigned char *)realloc((data_tmp = data), data_alloc -= (CHUNK_gz - 2) - bytes_read );
          *(data + data_alloc - 1) = 0;
          break;
        }
        else {
          const char * error_string;
          error_string = gzerror (gz_in, & err);
          if (err) {
            cout << "Exception during zlib decompression: (" << err << ") " << error_string;
            free(data);
						return(err);
          }
        }
      }
    }
    if(asStream) {
      iss.str((char*)data);
      loaded = true; streamed = true; lazy = false;
    } else {
      char *buffer_tmp = NULL;
      buffer = (char*)realloc(buffer_tmp = buffer, curpos);
      memcpy(buffer, data, curpos);
      bufferLen = curpos;
      loaded = true; streamed = false; lazy = false;
    }
    gzclose(gz_in);
    free(data);
  } else {
    loaded = true; streamed = false; lazy = true;
  }
	return(0);
}

// Operations to read gzip as if it was ifstream object. Not used in SpliceWiz

void GZReader::read(char * dest, const unsigned long len) {
  memcpy(dest, &buffer[bufferPos], len);
  bufferPos += len;
}
void GZReader::ignore(const unsigned long len) {
  bufferPos += len;
}

bool GZReader::eof() {
  return(gzeof(gz_in) && bufferPos == bufferLen);
}

// Only required for lazy mode
int GZReader::closeGZ() {
	gzclose(gz_in);
	return(0);
}


void GZWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}

// Writes a line to a gzipped file, including '\n'
int GZWriter::writeline(const std::string& s_src) {
  unsigned int s_size = s_src.size() + 1;
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size - 1);
  line[s_size - 1] = '\n';
  int ret = writebuffer(line, s_size);
  delete[] line;
  return(ret);
}

// Writes a string to a gzipped file, excluding '\n'
int GZWriter::writestring(const std::string& s_src) {
  unsigned int s_size = s_src.size();
  char * line = new char[s_size];
  memcpy(line, s_src.data(), s_size);
  int ret = writebuffer(line, s_size);
  delete[] line;
  return(ret);
}

// Writes from given char* src buffer of given length len. Used internally in SpliceWiz
int GZWriter::writebuffer(const char * src, unsigned int len) {
  unsigned int bytesremaining = len;
  unsigned int srcpos = 0;
  int ret;
  if(bufferPos >= CHUNK_gz) {
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }  
  while (bytesremaining + bufferPos > CHUNK_gz) {
    memcpy(&buffer[bufferPos], &src[srcpos], CHUNK_gz - bufferPos);
    srcpos += CHUNK_gz - bufferPos;
    bytesremaining -= CHUNK_gz - bufferPos;
    bufferPos = CHUNK_gz;
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }
  memcpy(&buffer[bufferPos], &src[srcpos], bytesremaining);
  bufferPos += bytesremaining;
  bytesremaining = 0;
  if(bufferPos >= CHUNK_gz) {
    ret = flush(0);
    if(ret != Z_OK) return(ret);
  }  
  return(0);
}

// Writes from memory to file via ostream
// Returns Z_OK if success, and error message otherwise
int GZWriter::flush(bool final) {
  if(bufferPos > 0) {
    int ret;
    unsigned int have;
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    
    ret = deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY);
    if (ret != Z_OK) {
      cout << "Exception during zlib initialization: (" << ret << ") "  << strm.msg;
			return(ret);
    }

    strm.avail_in = bufferPos;
    strm.next_in = (Bytef*)buffer;
    strm.avail_out = CHUNK_gz;
    strm.next_out = (Bytef*)compressed_buffer;
    
    ret = deflate(&strm, Z_FINISH);  

    if (ret != Z_OK && ret != Z_STREAM_END) {
        cout << "Exception during zlib deflate: (" << ret << ") " << strm.msg;
				return(ret);
    }

    have = strm.total_out;
  
    OUT->write(compressed_buffer, have);
    if(final) {   
      OUT->flush();
    }
    
    deflateEnd(&strm);
    bufferPos=0;
  }
  return(Z_OK);
}