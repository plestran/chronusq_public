/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2015 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#ifndef INCLUDED_FILEIO
#define INCLUDED_FILEIO
#define MAXBLOCK 1000
#include <global.h>
#include <cerr.h>
#include <tools.h>

/****************************/
/* Error Messages 1000-1999 */
/****************************/
namespace ChronusQ {
class FileIO {

  std::string  name_in;                 // name of the input file
  std::string  name_out;                // name of the output file
  std::string  name_scr;                // name of the scratch file
  std::string  name_bin;                // name of the binary file

public:

  fstream in;                    // file handler of the input file
  fstream out;                   // file handler of the output file
  fstream scr;                   // file handler of the scratch file
  fstream bin;                   // file handler of the binary file

  // constructor and destructor
  FileIO(const std::string);
  ~FileIO() {
    if(in.is_open()) in.close();
    if(scr.is_open()) scr.close();
    if(bin.is_open()) bin.close();
    if(scr.is_open()) scr.close();
    if(bin.is_open()) bin.close();
    if(out.is_open()) out.close();
  }

  // Python API
  void write(std::string);

};
} // namespace ChronusQ
#endif
