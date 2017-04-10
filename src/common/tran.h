// File Location: 

//
// Copyright (c) 2016-2017 Geosim Ltd.
// 
// Written by Amit Henig
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

/******************************************************************************
*
*: Package Name: ifrmath_tran
*
*: Title:
*
******************************************************************************/

#ifndef __ifrmath_tran_H
#define __ifrmath_tran_H


/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/


  /******************************************************************************
  *                             EXPORTED CONSTANTS                              *
  ******************************************************************************/

  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

namespace std
{
  template<class _Ty> class complex;
}

namespace tpcl
{

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/

  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/



  /** Descrete Fourier transform using butterfly 
    * @param Xio_Data           both input data and output.
    * @param Xi_N               length of both input data and result. */
  bool DFT(unsigned int Xi_size, std::complex<float> *Xio_data, bool Xi_forward = true);

  bool DFT2D(unsigned int Xi_width, unsigned int Xi_height, std::complex<float> *Xio_data, bool Xi_forward=true);

  //bool DFT(unsigned int Xi_size, std::complex<double> *Xio_data, bool Xi_forward = true);

} // namespace tpcl

#endif
