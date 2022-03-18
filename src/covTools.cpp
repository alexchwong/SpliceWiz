/* covTools.cpp Reads and writes COV files (BAM coverage)

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

#include "covTools.h"

// Constructor
covReader::covReader() {
  bufferPos = 0;
  bufferMax = 0;
  index_begin = 0;
  body_begin = 0;

  compressed_buffer = (char*)malloc(65536);
  buffer = (char*)malloc(65536);
}

// Destructor
covReader::~covReader() {
  free(buffer);
  free(compressed_buffer);
}

void covReader::SetInputHandle(std::istream *in_stream) {
  IS_EOF = 0;
  IS_FAIL = 0;
  IS_LENGTH = 0;

  IN = in_stream;

  // Identify EOF
  IN->seekg (0, std::ios_base::end);
  IS_LENGTH = IN->tellg();
  
  // Check EOF bit
  IN->seekg (-bamEOFlength, std::ios_base::end);
  
  char check_eof_buffer[bamEOFlength+1];
  IN->read(check_eof_buffer, bamEOFlength);
       
  if(strncmp(check_eof_buffer, bamEOF, bamEOFlength) == 0) {
    EOF_POS = IS_LENGTH - bamEOFlength;
  } else {
    // cout << "EOF bit not detected\n";
    EOF_POS = 0;
    IS_EOF = 1;
    IS_FAIL = 1;				
  }
  IN->seekg (0, std::ios_base::beg);
}

int covReader::ReadBuffer() {
  // read compressed buffer
  if((size_t)IN->tellg() >= EOF_POS) {
    IS_EOF = 1;
    return(Z_STREAM_END);
  } else if(fail()) {
    return(Z_STREAM_ERROR);
  }

  stream_uint16 u16;
  int ret = 0;
  
  char GzipCheck[bamGzipHeadLength];
  IN->read(GzipCheck, bamGzipHeadLength);

   if(strncmp(bamGzipHead, GzipCheck, bamGzipHeadLength) != 0) {
    cout << "Exception during BAM decompression - BGZF header corrupt: (at " 
      << IN->tellg() << " bytes) ";
    return(Z_BUF_ERROR);
  }

  IN->read(u16.c, 2);
  IN->read(compressed_buffer, u16.u + 1 - 2  - bamGzipHeadLength);

  bufferMax = 65536;
  z_stream zs;
  zs.zalloc = NULL;
  zs.zfree = NULL;
  zs.msg = NULL;
  zs.next_in = (Bytef*)compressed_buffer;
  zs.avail_in = u16.u + 1 - 2  - bamGzipHeadLength;
  zs.next_out = (Bytef*)buffer;
  zs.avail_out = bufferMax;

  stream_uint32 u32;
  memcpy(u32.c, &compressed_buffer[u16.u + 1 - 2 - bamGzipHeadLength - 8],4);

  ret = inflateInit2(&zs, -15);
  if(ret != Z_OK) {
    cout << "Exception during BAM decompression - inflateInit2() fail: (" 
      << ret << ") ";
    return(ret);
  }
  ret = inflate(&zs, Z_FINISH);
  if(ret != Z_OK && ret != Z_STREAM_END) {
    cout << "Exception during BAM decompression - inflate() fail: (" 
      << ret << ") ";
    return(ret);
  }
  ret = inflateEnd(&zs);
  
  bufferMax -= zs.avail_out;
  
  // check CRC
  uint32_t crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, bufferMax);
  // CRC check:
  if(u32.u != crc) {
    std::ostringstream oss;
    oss << "CRC fail during BAM decompression: (at " << IN->tellg() << " bytes) ";
    return(ret);
  }
  bufferPos = 0;
  
  return(ret);
}

int covReader::read(char * dest, unsigned int len) {
  unsigned int remaining_bytes = 0;
  unsigned int dest_pos = 0;
  int ret = 0;
  // Read next block if buffer empty or if pos is at end of buffer
  if(bufferMax == 0 || bufferPos == bufferMax) {
    ret = ReadBuffer();
    if(ret != Z_OK) return(ret);
  }
  
  if (len <= bufferMax - bufferPos) {
    memcpy(&dest[0], &buffer[bufferPos], len);
    bufferPos += len;
    return(Z_OK);
  } else {
    memcpy(&dest[dest_pos], &buffer[bufferPos], bufferMax - bufferPos);
    remaining_bytes = len - (bufferMax - bufferPos);
    dest_pos += bufferMax - bufferPos;
    bufferMax = 0;
    bufferPos = 0;
    ret = ReadBuffer();
    if(ret != Z_OK) return(ret);

    while(remaining_bytes > bufferMax) {
      memcpy(&dest[dest_pos], &buffer[0], bufferMax);
      remaining_bytes -= bufferMax;
      dest_pos += bufferMax;
      bufferMax = 0;
      bufferPos = 0;
      ret = ReadBuffer();
      if(ret != Z_OK) return(ret);
    }
    
    memcpy(&dest[dest_pos], &buffer[bufferPos], remaining_bytes);
    bufferPos += remaining_bytes;
    dest_pos += remaining_bytes;
  }
  return(Z_OK);
}

int covReader::ignore(unsigned int len) {
  // Essentially copy read() but without memcpy etc.
  unsigned int remaining_bytes = len;

  if(bufferMax == 0 || bufferPos == bufferMax) {
    ReadBuffer();        
  }
  
  if (len <= bufferMax - bufferPos) {
    bufferPos += len;
    return(Z_OK);
  } else {
    remaining_bytes = len - (bufferMax - bufferPos);

    bufferMax = 0;
    bufferPos = 0;
    ReadBuffer();

    while(remaining_bytes > bufferMax) {
      remaining_bytes -= bufferMax;
      bufferMax = 0;
      bufferPos = 0;
      ReadBuffer();
    }
    
    bufferPos += remaining_bytes;
    // Note bufferPos can equal bufferMax. IN->tellg() will return the position of the next bgzf block
  }
  return(0);
}

bool covReader::eof() {
  if(IS_EOF == 1) {
    return (true);
  } else {
    if(IN->eof()) {
      IS_EOF = 1;
      return (true);
    } else {
      return (false);
    }
  }
}

bool covReader::fail() {
  if(IS_FAIL == 1) {
    return (true);
  } else {
    if(IN->fail()) {
      IS_FAIL = 1;
      return (true);
    } else {
      return (false);
    }
  }
}

int covReader::ReadHeader() {
  IN->seekg (0, std::ios_base::beg);    
  chr_names.clear();
  chr_lens.clear();
  bufferPos = 0;
  bufferMax = 0;    
  
  char cov_header[4];
  int ret = read(cov_header,4);
  if(ret != Z_OK) {
    cout << "File is not BGZF compressed; unlikely to be COV file\n";
    return(-1);
  }
  std::string s_cov_header = "COV\x01";
  if(strncmp(cov_header, s_cov_header.c_str(), 4) != 0) {
    cout << "COV file has incorrect header!\n";
    return(-1);
  }
  
  stream_uint32 n_ref;
  std::string chrName;

  read(n_ref.c, 4);
  for(unsigned int i = 0; i < n_ref.u; i++) {
    stream_uint32 l_name;
    read(l_name.c, 4);
    
    char * c_name = new char[l_name.u];
    read(c_name, l_name.u);
    chrName = string(c_name, l_name.u - 1);
    chr_names.push_back(chrName);

    stream_uint32 l_ref;
    read(l_ref.c, 4);
    chr_lens.push_back(l_ref.u);
    
    delete[] c_name;
  }
  
  // should be the start point of bgzf block containing index
  index_begin = IN->tellg();      
  
  bufferPos = 0;
  bufferMax = 0;    

  // keep profiling to identify where body_begin is
  stream_uint32 u32;
  for(unsigned int j = 0; j < 3; j++) {
    for(unsigned int i = 0; i < chr_names.size(); i++) {
      read(u32.c, 4);
      ignore(u32.u);            
    }
  }
  body_begin = IN->tellg();

  return(n_ref.u);
}

int covReader::GetChrs(std::vector<chr_entry> &chrs) {
  if(chr_names.size() > 0) {
    for(unsigned int i = 0; i < chr_names.size(); i++) {
      chrs.push_back(chr_entry(i, chr_names.at(i), chr_lens.at(i)));
    }
  }
  return(0);
}

int covReader::FetchPos(
  const std::string seqname, const uint32_t start, const int strand,
  uint64_t * file_offset, uint32_t * block_start
) {

  // Inputs seqname, and desired start position of query
  // Alters file_offset to point to the compressed offset position of bgzf block to start reading
  // Alters block_start to display the actual coordinate start of the first entry of block
  // Returns:
      // 0 = success
      // -1 = fail
  if(strand < 0 || strand > 2) {
    return -1;
  }        
     
  if(index_begin == 0) {
    ReadHeader();
    if(index_begin == 0) {
      return -1;
    }
  }
  
  int ref_index;
  auto it_chr = std::find(chr_names.begin(), chr_names.end(), seqname);
  if(it_chr == chr_names.end()) {
    return -1;
  } else {
    ref_index = distance(chr_names.begin(), it_chr);
    ref_index += strand * chr_names.size();
  }

  int i = 0;
  stream_uint32 u32;
  IN->seekg(index_begin, std::ios_base::beg);    
  bufferPos = 0;
  bufferMax = 0;    

  while(i < ref_index) {
    read(u32.c, 4);
    ignore(u32.u);
    i += 1;
  }
  
  // Now we are at the start of the ref_index
  stream_uint32 chr_block_size;
  stream_uint32 cur_block_start;
  stream_uint64 cur_offset;
  uint32_t prev_block_start;
  uint64_t prev_offset = 0;

  uint32_t block_counter = 0;
  prev_block_start = 0;
  
  read(chr_block_size.c, 4);     // use this to know when to stop
  while(block_counter < chr_block_size.u) {
    read(cur_block_start.c, 4);
    read(cur_offset.c, 8);
    block_counter += 12;
    if(cur_block_start.u > start) { 
      // first entry should always equal zero so should not be immediately called
      break;
    } else {
      prev_block_start = cur_block_start.u;  
      prev_offset = cur_offset.u;
    }
  }

  *file_offset = prev_offset + body_begin;
  *block_start = prev_block_start;
  return 0;
}

int covReader::FetchRLE(
    const std::string seqname, 
    const uint32_t start, const uint32_t end, 
    const int strand,
    std::vector<int> * values, 
    std::vector<unsigned int> * lengths
) {
 
  stream_int32 i32;
  stream_uint32 u32;
  
  uint64_t file_offset = 0;
  uint32_t block_start = 0;
  
  auto it_chr = std::find(chr_names.begin(), chr_names.end(), seqname);
  if(it_chr == chr_names.end()) {
    return -1;
  } else {
    int ref_index = distance(chr_names.begin(), it_chr);
    if(end > chr_lens[ref_index]) {
      return -1;
    }
  }
  
  if(FetchPos(seqname, start, strand, &file_offset, &block_start) != 0) {
    return -1;
  }
  
  IN->seekg(file_offset, std::ios_base::beg);
  bufferPos = 0;
  bufferMax = 0;    
  
  uint32_t prev_start = block_start;
  
  do {
    read(i32.c, 4);
    read(u32.c, 4);
    prev_start += u32.u;
  } while(prev_start < start);
  
  if(prev_start > start) {
    if(prev_start >= end) {
      (*values).push_back(i32.i);
      (*lengths).push_back((unsigned int)(end - start));
    } else {
      (*values).push_back(i32.i);
      (*lengths).push_back((unsigned int)(prev_start - start));         
    }
  }
  // main loop
  while(prev_start < end) {
    read(i32.c, 4);
    read(u32.c, 4);
    if(prev_start + u32.u >= end) {
      (*values).push_back(i32.i);
      (*lengths).push_back((unsigned int)(end - prev_start));
      break;
    } else {
      (*values).push_back(i32.i);
      (*lengths).push_back(u32.u);       
    }
    prev_start += u32.u;
  };
  return(0);
}


// ######################### COV WRITER ########################################


buffer_out_chunk::buffer_out_chunk() {
  buffer = (char*)malloc(65536);
  // leave compressed buffer uninitialized until needed
}

// Destructor
buffer_out_chunk::~buffer_out_chunk() {
  if(buffer) free(buffer);
  if(compressed_buffer) free(compressed_buffer);
}

// Writes from src to buffer at current buffer_pos
unsigned int buffer_out_chunk::write(char * src, unsigned int len) {
  if(buffer_pos + len > BGZF_max) return(0);
  
  memcpy(buffer + buffer_pos, src, len);
  buffer_pos += len;
  if(buffer_pos > buffer_size) buffer_size = buffer_pos;
  return(len);
}

// Writes from compressed_buffer to ostream
int buffer_out_chunk::WriteToFile(ostream * OUT) {
  if(compressed_size == 0) return(Z_DATA_ERROR);
  OUT->write(compressed_buffer, compressed_size);
  free(compressed_buffer);
  compressed_size = 0;
  compressed_buffer = NULL;
  return(0);
}

// Compresses buffer, 
int buffer_out_chunk::Compress() {
  if(buffer_size < 1) return(Z_DATA_ERROR);
  if(buffer_size > BGZF_max) return(Z_DATA_ERROR);
  
  stream_uint16 u16;
  stream_uint32 u32;
  uint32_t crc;
  z_stream zs;
  
  char * temp_comp_buffer;
  temp_comp_buffer = (char*)malloc(65536);

  zs.zalloc = NULL; zs.zfree = NULL;
  zs.msg = NULL;
  zs.next_in  = (Bytef*)buffer;
  zs.avail_in = buffer_size;
  zs.next_out = (Bytef*)temp_comp_buffer;
  zs.avail_out = 65536 - 18 - 8;

  // -15 to disable zlib header/footer
  int ret = deflateInit2(&zs, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY); 
  if(ret != Z_OK) {
    cout << "Exception during BAM decompression - deflateInit2() fail: (" 
      << ret << ") ";
    return(ret);
  }
  
  ret = deflate(&zs, Z_FINISH);
  if(ret != Z_OK && ret != Z_STREAM_END) {
    cout << "Exception during BAM decompression - deflate() fail: (" 
      << ret << ") ";
    return(ret);
  }
  
  ret = deflateEnd(&zs);
  if(ret != Z_OK) {
    cout << "Exception during BAM decompression - deflateEnd() fail: (" 
      << ret << ") ";
    return(ret);
  }
  
  int block_len = zs.total_out + 18 + 8;
  
  // Initialize compressed buffer
  compressed_buffer = (char*)malloc(block_len + 1);

  memcpy(compressed_buffer, bamGzipHead, 16);

  u16.u = block_len - 1;
  memcpy(compressed_buffer + 16, u16.c, 2);
  
  memcpy(compressed_buffer + 18, temp_comp_buffer, zs.total_out);

  crc = crc32(crc32(0L, NULL, 0L), (Bytef*)buffer, buffer_size);
  u32.u = crc;
  memcpy(compressed_buffer + 18 + zs.total_out, u32.c, 4);

  u32.u = buffer_size;
  memcpy(compressed_buffer + 18 + zs.total_out + 4, u32.c, 4);

  // Now that compressed buffer is done, remove buffer to save memory:
  free(buffer);
  buffer = NULL;
  
  compressed_size = block_len;

  free(temp_comp_buffer);
  return(ret);
}

covWriter::covWriter() {
  // do nothing for now
}

covWriter::~covWriter() {
  // also do nothing
}

// Internal function
// If there are no alignments to the given refID, 
// - write a single entry (0, chr_len)
// - means zero depth, for length of chromosome
int covWriter::WriteEmptyEntry(unsigned int refID) {
  if(chrs.size() == 0) {
    cout << "ERROR: COV header missing\n";
    return(-1);
  }
  if(refID >= 3 * chrs.size()) {
    cout << "ERROR: Invalid chrID parsed to covWriter\n";
    return(-1);
  }

  unsigned int chrID = refID;
  while(chrID > chrs.size()) chrID -= chrs.size();

  // Resize body and block_coord_starts to 1
  body.at(refID).resize(1);
  block_coord_starts.at(refID).resize(1);
  
  // Initialize first pos = 0
  block_coord_starts.at(refID).at(0) = 0;

  stream_int32 i32;
  stream_uint32 u32;

  i32.i = 0;
  body.at(refID).at(0).write(i32.c, 4);

  u32.u = chrs.at(chrID).chr_len;
  body.at(refID).at(0).write(u32.c, 4);
  
  body.at(refID).at(0).Compress();
  
  return(0);
}

// Writes the header to file. Internal
// - writes 'COV\1', then the chrom count
// - then l_chr_name, chr_name, chr_len
int covWriter::WriteHeaderToFile() {
  // Write the header to file:
  char zero = '\0';
  char wh_buffer[1000];
  std::string header_str = "COV\x01";
  stream_uint32 u32;
  
  buffer_out_chunk * header = new buffer_out_chunk;
  strncpy(wh_buffer, header_str.c_str(), 4);
  header->write(wh_buffer, 4);
  
  u32.u = chrs.size();    // number of chroms
  header->write(u32.c ,4);
  
  for(unsigned int i = 0; i < chrs.size(); i++) {
    unsigned int chr_buf_len = 8 + 1 + chrs.at(i).chr_name.length();
    
    // This is very unlikely to run
    if(header->IsAtCap(chr_buf_len)) {
      header->Compress();
      header->WriteToFile(OUT);
      delete header;
      
      header = new buffer_out_chunk;
    }
    
    u32.u = chrs.at(i).chr_name.length() + 1;
    header->write(u32.c, 4);
    
    strncpy(wh_buffer, chrs.at(i).chr_name.c_str(), chrs.at(i).chr_name.length());
    header->write(wh_buffer, chrs.at(i).chr_name.length());
    header->write(&zero, 1);

    u32.u = chrs.at(i).chr_len;
    header->write(u32.c ,4);
  }

  header->Compress();
  header->WriteToFile(OUT);
  delete header;
  
  return(0);
}

// Writes the index to file. Internal
int covWriter::WriteIndexToFile() {
  stream_uint32 u32;
  stream_uint64 u64;
  
  std::vector< buffer_out_chunk > index_buffer;
  
  uint32_t index_size = 0;
  uint64_t body_pos = 0;
  unsigned int cur_buffer = 0;
  
  for(unsigned int i = 0; i < 3 * chrs.size(); i++) {
    if(block_coord_starts.at(i).size() == 0 || body.at(i).size() == 0) 
      WriteEmptyEntry(i);

    index_size = 0;   // Resets to zero for every refID
    index_buffer.resize(1);
    index_buffer.at(cur_buffer).SetPos(4); // Write the index size at the very end

    for(unsigned int j = 0; j < body.at(i).size(); j++) {
      // Increase vector by 1 if new BGZF block needed
      if(index_buffer.at(cur_buffer).IsAtCap(12)) {
        index_buffer.resize(index_buffer.size() + 1);
        cur_buffer++;
      }
      // Write block starting POS
      u32.u = block_coord_starts.at(i).at(j);
      index_buffer.at(cur_buffer).write(u32.c, 4);

      // Write current body byte-wise displacement from start of COV body
      u64.u = body_pos;
      index_buffer.at(cur_buffer).write(u64.c, 8);
      
      // Displace body_pos by the size of compressed BGZF block
      body_pos += body.at(i).at(j).getBGZFSize();
      
      // Index entry consumes 12 bytes
      index_size += 12;
    }
    
    // Write uncompressed size of index entry to the beginning of index block
    u32.u = index_size;
    index_buffer.at(0).write_to_pos(u32.c, 4, 0);
    
    // Write all index buffers of current refID to file
    for(unsigned int j = 0; j < index_buffer.size(); j++) {
      index_buffer.at(j).Compress();
      index_buffer.at(j).WriteToFile(OUT);
    }
    index_buffer.clear();
  }

  return(0);
}

void covWriter::SetOutputHandle(std::ostream *out_stream) {
  OUT = out_stream;
}

// Initializes nested vectors according to the number of chromosomes
int covWriter::InitializeCOV(std::vector<chr_entry> chrs_to_copy) {
  for(auto chr : chrs_to_copy) {
    chrs.push_back(chr);
  }

  // Now initialise nested vectors:
  block_coord_starts.resize(chrs.size() * 3);
  body.resize(chrs.size() * 3);

  // Make sure stuff within nested vectors is empty:
  for(unsigned int i = 0; i < chrs.size() * 3; i++) {
    block_coord_starts.at(i).resize(0);
    body.at(i).resize(0);
  }
  
  return 0;
}

// Given a vector of pairs, chrID, and strand, write these to COV body
// - Called from FragmentsMap::WriteBinary(covWriter *os, bool verbose)
// - vec is an RLE of loci and depth
// - output is depth (NOT depth increment) and length
int covWriter::WriteFragmentsMap(
    std::vector< std::pair<unsigned int, int> > * vec, 
    unsigned int chrID, unsigned int strand,
    unsigned int n_threads_to_use
) {
  if(chrs.size() == 0) {
    cout << "ERROR: COV header missing\n";
    return(-1);
  }
  if(chrID >= chrs.size()) {
    cout << "ERROR: Invalid chrID parsed to covWriter\n";
    return(-1);
  }
  // Initialize the vector depending on vector size
  unsigned int vec_cap = (BGZF_max / 8);
  
  unsigned int vec_size = vec->size();
  unsigned int job_size = vec_size / vec_cap;
  if(job_size * vec_cap < vec_size) job_size++;
  // job_size is the number of BGZF blocks this refID will occupy
  
  unsigned int refID = chrID + chrs.size() * strand;
  body.at(refID).resize(job_size);
  block_coord_starts.at(refID).resize(job_size);
  
#ifdef _OPENMP
  #pragma omp parallel for num_threads(n_threads_to_use) schedule(static,1)
#endif
  for(unsigned int i = 0; i < job_size; i++) {
    stream_int32 i32;
    stream_uint32 u32;
    // Start coordinate for this bgzf block
    block_coord_starts.at(refID).at(i) = (uint32_t)vec->at(i * vec_cap).first;
    
    unsigned int cur_coord = vec->at(i * vec_cap).first;
    
    for(unsigned int j = i * vec_cap; j < (i+1) * vec_cap && j < vec_size; j++) {
      
      // if last entry, assume it extends till end of chromosome
      if(j == vec_size - 1) {
        if((unsigned int)chrs.at(chrID).chr_len > cur_coord) {
          i32.i = vec->at(j).second;
          body.at(refID).at(i).write(i32.c, 4);
          
          u32.u = (unsigned int)chrs.at(chrID).chr_len - cur_coord;
          body.at(refID).at(i).write(u32.c, 4);
        }
        cur_coord = chrs.at(chrID).chr_len;   // This step is probably pointless
      } else {
        // Avoid writing zero-length RLE entries
        if(vec->at(j + 1).first > cur_coord) {
          i32.i = vec->at(j).second;
          body.at(refID).at(i).write(i32.c, 4);
          
          u32.u = vec->at(j + 1).first - cur_coord;
          body.at(refID).at(i).write(u32.c, 4);
          cur_coord = vec->at(j + 1).first;
        }
      }
    }
    body.at(refID).at(i).Compress();
  }
  return(0);
}

// Writes everything to file
int covWriter::WriteToFile() {
  if(!OUT) {
    cout << "No COV file set to write to";
    return(-1);
  }
  if(chrs.size() == 0) {
    cout << "ERROR: COV header missing\n";
    return(-1);
  }
  
  WriteHeaderToFile();
  WriteIndexToFile();
  for(unsigned int i = 0; i < 3 * chrs.size(); i++) {
    for(unsigned int j = 0; j < body.at(i).size(); j++) {
      body.at(i).at(j).WriteToFile(OUT);
    }
  }

  OUT->write(bamEOF, bamEOFlength);
  OUT->flush();
  
  return(0);
}