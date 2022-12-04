/* synthReadGenerator.cpp Generates synthetic mappability reads

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

#include "synthReadGenerator.h"

synthReadGenerator::synthReadGenerator() {
  read_len = 70;
  error_pos = 35;
}

synthReadGenerator::synthReadGenerator(
  unsigned int read_length,
  int error_position
){
  read_len = read_length;
  error_pos = error_position;
}

// Replicate old PERL script; return true if N's constitute less than half of length
bool synthReadGenerator::checkDNA(char * input_read) {
  unsigned int numN = 0;
  for(unsigned int i = 0; i < read_len; i++) {
    if(input_read[i]!='A' && input_read[i]!='T' && input_read[i]!='G' && input_read[i]!='C' &&
      input_read[i]!='a' && input_read[i]!='t' && input_read[i]!='g' && input_read[i]!='c') {
      numN++;
    }
  }
  return(numN < read_len / 2);
}

std::string synthReadGenerator::GenerateReadError(
    char * input_read, 
    const unsigned int direction, 
    const size_t error_seed
) {
  
  // Copy https://github.com/williamritchie/IRFinder/blob/master/bin/util/generateReadsError.pl

  char * new_read_inter = new char[read_len + 1];
  new_read_inter[read_len] = '\0';
  memcpy(&new_read_inter[0], input_read, read_len);  

  char error_nuc = '\0';  // set this as something to avoid warning at compile
  if(error_seed % 3 == 0) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'G'; break;
      case 'C':
        error_nuc = 'A'; break;
      case 'G':
        error_nuc = 'T'; break;
      case 'T':
        error_nuc = 'C'; break;
      case 'a':
        error_nuc = 'g'; break;
      case 'c':
        error_nuc = 'a'; break;
      case 'g':
        error_nuc = 't'; break;
      case 't':
        error_nuc = 'c'; break;
      default:
        error_nuc = 'N';
    }
  } else if(error_seed % 3 == 1) {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'T'; break;
      case 'C':
        error_nuc = 'G'; break;
      case 'G':
        error_nuc = 'C'; break;
      case 'T':
        error_nuc = 'A'; break;
      case 'a':
        error_nuc = 't'; break;
      case 'c':
        error_nuc = 'g'; break;
      case 'g':
        error_nuc = 'c'; break;
      case 't':
        error_nuc = 'a'; break;
      default:
        error_nuc = 'N';
    }
  } else {
    switch(new_read_inter[error_pos - 1]) {
      case 'A':
        error_nuc = 'C'; break;
      case 'C':
        error_nuc = 'T'; break;
      case 'G':
        error_nuc = 'A'; break;
      case 'T':
        error_nuc = 'G'; break;
      case 'a':
        error_nuc = 'c'; break;
      case 'c':
        error_nuc = 't'; break;
      case 'g':
        error_nuc = 'a'; break;
      case 't':
        error_nuc = 'g'; break;
      default:
        error_nuc = 'N';
    }
  }
  
  memcpy(&new_read_inter[error_pos - 1], &error_nuc, 1);
  
  char * new_read = new char[read_len + 1];
  new_read[read_len] = '\0';
  if(direction == 0) {
    memcpy(&new_read[0], new_read_inter, read_len);  
  } else {
    for(unsigned int i = 0; i < read_len; i++) {
      switch(new_read_inter[i]) {   
        case 'A':
          new_read[read_len - i - 1] = 'T'; break;
        case 'T':
          new_read[read_len - i - 1] = 'A'; break;
        case 'G':
          new_read[read_len - i - 1] = 'C'; break;
        case 'C':
          new_read[read_len - i - 1] = 'G'; break;
        case 'a':
          new_read[read_len - i - 1] = 't'; break;
        case 't':
          new_read[read_len - i - 1] = 'a'; break;
        case 'g':
          new_read[read_len - i - 1] = 'c'; break;
        case 'c':
          new_read[read_len - i - 1] = 'g'; break;
        default :
          new_read[read_len - i - 1] = 'N';
      }         
    }
  }

  string return_str = string(new_read);

  delete[] new_read_inter;
  delete[] new_read;
  return(return_str);
}