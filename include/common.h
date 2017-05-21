// File Location: 

//
// Copyright (c) 2016-2017 Geosim Ltd.
// 
// Written by Ramon Axelrod       ramon.axelrod@gmail.com
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
*: Package Name: common
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_common_H
#define __tpcl_common_H



namespace tpcl
{
  const double M_PI = 3.14159265358979323846264338327950288;

  template<class T> T MinT(T a, T b) { return ((a) < (b)) ? (a) : (b); }
  template<class T> T MaxT(T a, T b) { return ((a) > (b)) ? (a) : (b); }
  template<class T> void SwapT(T &a, T &b) { T c = a; a = b; b = c; }
  template<class T> inline T ClampT(T v, T mn, T mx) { return v < mn ? mn : (v > mx ? mx : v); }
  template<class T> inline T AbsT(T a) { return a < 0 ? -a : a; }


  //void Log(...) {}


} // namespace tpcl

#endif
