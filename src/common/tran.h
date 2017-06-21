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

#ifndef __ifrmath_tran_H
#define __ifrmath_tran_H


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
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/
  
  //assumes dimensions are powers of 2.

  /** forward/backwards 1D Descrete Fourier transform using butterfly.
  * @param in_size            length of both input data and result.
  * @param io_data           both input data and output.
  * @param in_forward         if true - forward FFT. if false - backwards FFT. */
  bool DFT(unsigned int in_size, std::complex<float> *io_data, bool in_forward = true);

  /** forward/backwards 2D Descrete Fourier transform using butterfly.
  * @param in_width           width of both input data and result.
  * @param in_height          height of both input data and result.
  * @param io_data           both input data and output.
  * @param in_forward         if true - forward FFT. if false - backwards FFT. */
  bool DFT2D(unsigned int in_width, unsigned int in_height, std::complex<float> *io_data, bool in_forward = true);


  /** shift (per dimension) by half of it's size, in order to center origin of axes (0,0).
  * @param in_size            for 1D: length of both input data and result.
  * @param in_width           for 2D: width of both input data and result.
  * @param in_height          for 2D: height of both input data and result.
  * @param io_data           both input data and output. */
  void DFTshift0ToOrigin(std::complex<float> *io_data, unsigned int in_size);
  void DFTshift0ToOrigin(std::complex<float> *io_data, unsigned int in_width, unsigned int in_height);


  /** calculates phase correlation in the frequency domain (for both 1D and 2D signals).
  * @param in_DFT0            DFT of the first signal.
  * @param in_DFT1            DFT of the second signal.
  * @param out_PhCor           DFT of input signals' phase correlation.
  * @param in_in_width        for 1D/2D: size/width of the signals.
  * @param in_in_height       fill if 2D signals: height of the signals. */
  void PhaseCorrelation(std::complex<float>* in_DFT0, std::complex<float>* in_DFT1, std::complex<float>* out_PhCor, unsigned int in_width, unsigned int in_height = 1);


  /** calculates phase correlation in the frequency domain (for both 1D and 2D signals).
  *   in result size of each of the elements is either 1 or 0.
  * @param in_DFT0            DFT of the first signal.
  * @param in_DFT1            DFT of the second signal.
  * @param out_PhCor           DFT of input signals' phase correlation.
  * @param in_in_width        for 1D/2D: size/width of the signals.
  * @param in_in_height       fill if 2D signals: height of the signals. */
  void UnitPhaseCorrelation(std::complex<float>* in_DFT0, std::complex<float>* in_DFT1, std::complex<float>* out_PhCor, unsigned int in_width, unsigned int in_height = 1);


} // namespace tpcl

#endif
