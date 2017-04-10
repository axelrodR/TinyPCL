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
******************************************************************************/

#include <complex>
#include "tran.h"

namespace tpcl
{

  //   Inplace version of rearrange function
  void Rearrange(std::complex<float>* Xio_data, unsigned int Xi_size)
  {
    //   Swap position
    unsigned int Target = 0;
    //   Process all positions of input signal
    for (unsigned int Position = 0; Position < Xi_size; ++Position)
    {
      //   Only for not yet swapped entries
      if (Target > Position)
      {
        //   Swap entries
        std::complex<float> Temp(Xio_data[Target]);
        Xio_data[Target] = Xio_data[Position];
        Xio_data[Position] = Temp;
      }
      //   Bit mask
      unsigned int Mask = Xi_size;
      //   While bit is set
      while (Target & (Mask >>= 1))
        //   Drop bit
        Target &= ~Mask;
      //   The current bit is 0 - set it
      Target |= Mask;
    }
  }


  //   FFT implementation
  void Perform(std::complex<float>* Xio_data, unsigned int Xi_size,
    bool Inverse  = false )
  {
    double pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;
    //   Iteration through dyads, quadruples, octads and so on...
    for (unsigned int Step = 1; Step < Xi_size; Step <<= 1)
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
        for (unsigned int Pair = Group; Pair < Xi_size; Pair += Jump)
        {
          //   Match position
          unsigned int Match = Pair + Step;
          //   Second term of two-point transform
          std::complex<float> Product(Factor * Xio_data[Match]);
          //   Transform for fi + pi
          Xio_data[Match] = Xio_data[Pair] - Product;
          //   Transform for fi
          Xio_data[Pair] += Product;
        }
        //   Successive transform factor via trigonometric recurrence
        Factor = Multiplier * Factor + Factor;
      }
    }
  }

  bool DFT(unsigned int Xi_size, std::complex<float> *Xio_data, bool Xi_forward)
  {
    //   Check input parameters
    if (!Xio_data || Xi_size < 1 || Xi_size & (Xi_size - 1))
      return false;
    //   Rearrange
    Rearrange(Xio_data, Xi_size);
    //   Call FFT implementation
    Perform(Xio_data, Xi_size, Xi_forward);
    //   Succeeded
    return true;
  }



  bool DFT2D(unsigned int Xi_width, unsigned int Xi_height, std::complex<float> *Xio_data, bool Xi_forward)
  {
    //   Check input parameters
    if (!Xio_data || Xi_width < 1 || Xi_width & (Xi_width - 1) || Xi_height < 1 || Xi_height & (Xi_height - 1))
      return false;

    for (unsigned int row = 0; row < Xi_height; row++)
    {      
      if ((!DFT(Xi_width, Xio_data + row * Xi_width, Xi_forward)))
        return false;
    }


  std::complex<float>*DFTcols = new std::complex<float>[Xi_height];
  for (unsigned int col = 0; col < Xi_width; col++)
  {
    for (unsigned int row = 0; row < Xi_height; row++)
    {
      DFTcols[row] = Xio_data[(row * Xi_width) + col];
    }

    if ((!DFT(Xi_width, DFTcols, Xi_forward)))
    {
      delete[] DFTcols;
      return false;
    }

    for (unsigned int row = 0; row < Xi_height; row++)
    {
       Xio_data[(row * Xi_width) + col] = DFTcols[row];
    }
  }


  delete[] DFTcols;

  return true;
  }

} // namespace tpcl

