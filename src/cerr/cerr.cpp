/* 
 *  The Chronus Quantum (ChronusQ) software package is high-performace 
 *  computational chemistry software with a strong emphasis on explicitly 
 *  time-dependent and post-SCF quantum mechanical methods.
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
#include <cerr.h>

void ChronusQ::CErr(std::ostream & out){
  time_t currentTime;
  time(&currentTime); 
  if(getRank() == 0) {
    out << "\"Die Die Die\"" << endl;
    out << "Job terminated: "<<ctime(&currentTime)<<endl;
  }
  libint2::finalize();
  exit(EXIT_FAILURE);
}
void ChronusQ::CErr(std::string msg,std::ostream & out){
  time_t currentTime;
  time(&currentTime); 
  if(getRank() == 0) {
    out << msg << endl;
    out << "Job terminated: "<<ctime(&currentTime)<<endl;
  }
  libint2::finalize();
  exit(EXIT_FAILURE);
}
void ChronusQ::CErr(std::exception_ptr eptr,std::ostream & out) {
  time_t currentTime;
  time(&currentTime); 
  try {
    if(eptr) std::rethrow_exception(eptr);
  } catch(const std::exception & e) {
    if(getRank() == 0) {
      out << "Caught \"" << e.what() << "\"" << endl;
      out << "Job terminated: "<<ctime(&currentTime)<<endl;
    }
    libint2::finalize();
    exit(EXIT_FAILURE);
  }
}
void ChronusQ::CErr(std::exception_ptr eptr,std::string msg,std::ostream & out) {
  time_t currentTime;
  time(&currentTime); 
  try {
    if(eptr) std::rethrow_exception(eptr);
  } catch(const std::exception & e) {
    if(getRank() == 0) {
      out << msg + " caught \"" << e.what() << "\"" << endl;
      out << "Job terminated: "<<ctime(&currentTime)<<endl;
    }
    libint2::finalize();
    exit(EXIT_FAILURE);
  }
}
