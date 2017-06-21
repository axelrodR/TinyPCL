//
// Copyright (c) 2016-2017 Geosim Ltd.
// FFT transform is a copyright of LIBROW (see details there)
// FFT2D transform is a copyright of Paul Burke
//
// Written by Amit Henig
//
// Functions within written by other people are credited with
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

#include <complex>
#include "tran.h"

namespace tpcl
{

  /** Inplace version of rearrange function for the butterfly rearrangement of 1D FFT.
  * @param io_data           both input data and output.
  * @param in_size            length of both input data and result. */
  void Rearrange(std::complex<float>* io_data, unsigned int in_size)
  {
    //   Swap position
    unsigned int Target = 0;
    //   Process all positions of input signal
    for (unsigned int Position = 0; Position < in_size; ++Position)
    {
      //   Only for not yet swapped entries
      if (Target > Position)
      {
        //   Swap entries
        std::complex<float> Temp(io_data[Target]);
        io_data[Target] = io_data[Position];
        io_data[Position] = Temp;
      }
      //   Bit mask
      unsigned int Mask = in_size;
      //   While bit is set
      while (Target & (Mask >>= 1))
        //   Drop bit
        Target &= ~Mask;
      //   The current bit is 0 - set it
      Target |= Mask;
    }
  }


  //   FFT implementation
  /** 1D FFT implementation for the already butterfly rearranged input signal.
  * @param io_data           both input data and output.
  * @param in_size            length of both input data and result. */
  void Perform(std::complex<float>* io_data, unsigned int in_size,
    bool in_forward = true )
  {
    double pi = in_forward ? -3.14159265358979323846 : 3.14159265358979323846;
    //   Iteration through dyads, quadruples, octads and so on...
    for (unsigned int Step = 1; Step < in_size; Step <<= 1)
    {
      //   Jump to the next entry of the same transform factor
      unsigned int Jump = Step << 1;
      //   Angle increment
      double delta = pi / double(Step);
      //   Auxiliary sin(delta / 2)
      double Sine = sin(delta * .5);
      //   Multiplier for trigonometric recurrence
      std::complex<float> Multiplier(float(-2. * Sine * Sine), float(sin(delta)));
      //   Start value for transform factor, fi = 0
      std::complex<float> Factor(1.);
      //   Iteration through groups of different transform factor
      for (unsigned int Group = 0; Group < Step; ++Group)
      {
        //   Iteration within group 
        for (unsigned int Pair = Group; Pair < in_size; Pair += Jump)
        {
          //   Match position
          unsigned int Match = Pair + Step;
          //   Second term of two-point transform
          std::complex<float> Product(Factor * io_data[Match]);
          //   Transform for fi + pi
          io_data[Match] = io_data[Pair] - Product;
          //   Transform for fi
          io_data[Pair] += Product;
        }
        //   Successive transform factor via trigonometric recurrence
        Factor = Multiplier * Factor + Factor;
      }
    }

    if (!in_forward)
    {
      //scale:
      float invSize = 1 / float(in_size);
      for (unsigned int index = 0; index < in_size; index++)
      {
        io_data[index] *= invSize;
      }
    }
  }


  //   The FFT code is based on FFT code from LIBROW:
  //   see: http://www.librow.com/articles/article-10
  //   You can use it on your own
  //   When utilizing credit LIBROW site
  bool DFT(unsigned int in_size, std::complex<float> *io_data, bool in_forward)
  {
    //   Check input parameters
    if (!io_data || in_size < 1 || in_size & (in_size - 1))
      return false;
    //   Rearrange
    Rearrange(io_data, in_size);
    //   Call FFT implementation
    Perform(io_data, in_size, in_forward);
    //   Succeeded
    return true;
  }


  //        The FFT2D is based on http://paulbourke.net/miscellaneous/dft/
  bool DFT2D(unsigned int in_width, unsigned int in_height, std::complex<float> *io_data, bool in_forward)
  {
    //   Check input parameters
    if (!io_data || in_width < 1 || in_width & (in_width - 1) || in_height < 1 || in_height & (in_height - 1))
      return false;

    for (unsigned int row = 0; row < in_height; row++)
    {      
      if ((!DFT(in_width, io_data + row * in_width, in_forward)))
        return false;
    }


  std::complex<float>* DFTcols = new std::complex<float>[in_height];
  for (unsigned int col = 0; col < in_width; col++)
  {
    for (unsigned int row = 0; row < in_height; row++)
    {
      DFTcols[row] = io_data[(row * in_width) + col];
    }

    if ((!DFT(in_height, DFTcols, in_forward)))
    {
      delete[] DFTcols;
      return false;
    }

    for (unsigned int row = 0; row < in_height; row++)
    {
       io_data[(row * in_width) + col] = DFTcols[row];
    }
  }


  delete[] DFTcols;

  return true;
  }


  void DFTshift0ToOrigin(std::complex<float> *io_data, unsigned int in_size)
  {
    unsigned  int index0 = 0;
    unsigned int index1 = in_size >> 1;

    for (index0; index0 < (in_size >> 1); index0++, index1++)
    {
      std::swap(io_data[index0], io_data[index1]);
    }
  }

  void DFTshift0ToOrigin(std::complex<float> *io_data, unsigned int in_width, unsigned int in_height)
  {
    unsigned int stepH = in_height >> 1;
    unsigned int stepW = in_width >> 1;
    for (unsigned int row = 0; row < stepH; row++)
    {
        unsigned int index1Q = (row * in_width);
        unsigned int index2Q = index1Q + stepW;
        unsigned int index3Q = ((row + stepH) * in_width);
        unsigned int index4Q = index3Q + stepW;
        unsigned int stopAt = index2Q;

      for (index1Q; index1Q < stopAt; index1Q++, index2Q++, index3Q++, index4Q++)
      {
        std::swap(io_data[index1Q], io_data[index4Q]);
        std::swap(io_data[index2Q], io_data[index3Q]);
      }
    }
  }


  void PhaseCorrelation(std::complex<float>* in_DFT0, std::complex<float>* in_DFT1, std::complex<float>* out_PhCor, unsigned int in_width, unsigned int in_height)
  {
    std::complex<float> zeroComplex = 0;
    for (unsigned int index2 = 0; index2 < in_height; index2++)
    {
      for (unsigned int index1 = 0; index1 < in_width; index1++)
      {
        unsigned int index = (index2 * in_width) + index1;

        out_PhCor[index] = in_DFT0[index] * conj(in_DFT1[index]);
      }
    }
  }

  void UnitPhaseCorrelation(std::complex<float>* in_DFT0, std::complex<float>* in_DFT1, std::complex<float>* out_PhCor, unsigned int in_width, unsigned int in_height)
  {
    std::complex<float> zeroComplex = 0;
    for (unsigned int index2 = 0; index2 < in_height; index2++)
    {
      for (unsigned int index1 = 0; index1 < in_width; index1++)
      {
        unsigned int index = (index2 * in_width) + index1;

        out_PhCor[index] = in_DFT0[index] * conj(in_DFT1[index]);
        if (out_PhCor[index] != zeroComplex)
          out_PhCor[index] /= std::abs(out_PhCor[index]);
      }
    }
  }


} // namespace IfrMath

