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
*: Package Name: tpcl
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_vec_H
#define __tpcl_vec_H

#include <math.h>

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#define DLL_Entry
#ifndef _LIB
#ifndef _LIB_LINK
#undef DLL_Entry

#ifdef TPCL_EXPORTS
#define DLL_Entry __declspec(dllexport)
#else
#define DLL_Entry __declspec(dllimport)
#endif
#endif
#endif


/******************************************************************************
*                             EXPORTED CONSTANTS                              *
******************************************************************************/

/******************************************************************************
*                        INCOMPLETE CLASS DECLARATIONS                        *
******************************************************************************/

namespace tpcl
{

  /******************************************************************************
  *                            EXPORTED CLASSES                                 *
  ******************************************************************************/

  /** 3D vector (floats or doubles) */
  template <typename T>
  struct TVec3
  {
    T x, y, z;

    // constructors
    TVec3() {}
    TVec3(T x_, T y_, T z_) { x=x_; y=y_; z=z_; }

    // addition, subtraction
    inline TVec3& operator+=(const TVec3& u);
    inline TVec3& operator-=(const TVec3& u);
    inline TVec3 operator+(const TVec3& rhs) const;
    inline TVec3 operator-(const TVec3& rhs) const;

    // multiplication, division
    inline TVec3& operator*=(T f);
    inline TVec3& operator/=(T f);
    inline TVec3 operator*(T rhs) const;
    inline TVec3 operator/(T rhs) const;

    // unary operators
    inline TVec3 operator-() const;

    // comparison
    inline bool operator==(const TVec3& ths) const;
    inline bool operator!=(const TVec3& rhs) const;

    // component access
    inline T& operator[](int i);
    inline T operator[](int i) const;

  };

  
  template <typename T> inline TVec3<T> operator*(T lhs, const TVec3<T>& rhs);


  /** 4D matrix (4x4) */
  template <typename T>
  struct TMat4
  {
    float m[4][4];

    // constructors
    TMat4() {}
    inline TMat4(T* Vals);

    // addition, subtraction
    inline TMat4& operator+=(const TMat4& u);
    inline TMat4& operator-=(const TMat4& u);
    inline TMat4 operator+(const TMat4& rhs);
    inline TMat4 operator-(const TMat4& rhs);

    // multiplication, division
    inline TMat4& operator*=(T f);
    inline TMat4& operator/=(T f);
    inline TMat4 operator*(T rhs);
    inline TMat4 operator/(T rhs);
    inline TMat4& operator*=(const TMat4& rhs);
    inline TMat4 operator*(const TMat4& rhs);

    // unary operator
    TMat4 operator-() const;

    // comparison
    inline bool operator==(const TMat4& ths) const;
    inline bool operator!=(const TMat4& rhs) const;
  };

  // unary operator
  template <typename TMat4> inline void MatrixIdentity(TMat4* rhs);
  template <typename TMat4> inline double MatrixDeterminant(const TMat4* rhs);
  template <typename TMat4> inline void Transpose(const TMat4& rhs, TMat4& res);
  template <typename TMat4, typename T> inline void Lengths(const TMat4& rhs, T& length_mat, T& length_vec);

  // multiplication in place
  template <typename T> inline void Multiply(const TMat4<T>& lhs, const TVec3<T>& rhs, TMat4<T>& res);
  template <typename TMat4, typename TVec3> inline void TransposeLeftMultiply(const TMat4& lhs, const TVec3& rhs, TMat4& res);
  template <typename TMat4, typename TVec3> inline void MultiplyVectorRightSide(const TMat4& lhs, const TVec3& rhs, TVec3& res);


  typedef TVec3<float> CVec3;
  typedef TMat4<float> CMat4;

  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/

  // length of Tvec3 by reference
  template <typename Tvec3>     inline float Length(const Tvec3& u);
  template <typename Tvec3>     inline float LengthSqr(const Tvec3& u);
  template <typename Tvec3>     inline float Length2D(const Tvec3& u);
  template <typename Tvec3>     inline float LengthSqr2D(const Tvec3& u);
  
  // distances between Tvec3 vectors
  template <typename Tvec3>     inline float Dist(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float DistSqr(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float Dist2D(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float DistSqr2D(const Tvec3& lhs, const Tvec3& rhs);

  // product, angles and normalzation
  template <typename Tvec3>     inline float DotProd(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float DotProd2D(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float CosAngle(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline float CosAngle2D(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline void Normalize(Tvec3& v);
  template <typename Tvec3>     inline Tvec3 CrossProd(const Tvec3& lhs, const Tvec3& rhs);

  // per component operations of Tvec3
  template <typename Tvec3>     inline Tvec3 Mul_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline Tvec3 Max_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline Tvec3 Min_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline bool Less_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline bool More_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline bool LessEq_ps(const Tvec3& lhs, const Tvec3& rhs);
  template <typename Tvec3>     inline bool MoreEq_ps(const Tvec3& lhs, const Tvec3& rhs);

  /** linear interploation */
  template <typename Tvec3>     inline Tvec3 Lerp(const Tvec3& lhs, const Tvec3& rhs, float t);

  bool IsValid(const TVec3<float>& Xi_v);
  bool IsValid(const TVec3<double>& Xi_v);




  

  /******************************************************************************
  *                 INLINE FUNCTIONS IMPLEMENTATION (3D vector)                 *
  ******************************************************************************/

  template <typename T>
  inline TVec3<T>& TVec3<T>::operator+=(const TVec3<T>& rhs)
  {
    x+=rhs.x; y+=rhs.y; z+=rhs.z; 
    return *this;
  }
  
  template <typename T>
  inline TVec3<T>& TVec3<T>::operator-=(const TVec3<T>& rhs)
  {
    x-=rhs.x; y-=rhs.y; z-=rhs.z; 
    return *this;
  }

  template <typename T>
  inline TVec3<T> TVec3<T>::operator+(const TVec3<T>& rhs) const
    { return TVec3<T>(x+rhs.x, y+rhs.y, z+rhs.z); }

  template <typename T>
  inline TVec3<T> TVec3<T>::operator-(const TVec3<T>& rhs) const
    { return TVec3<T>(x-rhs.x, y-rhs.y, z-rhs.z); }

  // multiplication, division
  template <typename T>
  inline TVec3<T>& TVec3<T>::operator*=(T f)
  {
    x*=f; y*=f; z*= f;
    return *this;
  }

  template <typename T>
  inline TVec3<T>& TVec3<T>::operator/=(T f)
  {
    (*this) *= 1.0f / f; 
    return *this;
  }

  template <typename T>
  inline TVec3<T> TVec3<T>::operator*(T rhs) const
    { return TVec3<T>(x*rhs, y*rhs, z*rhs); }

  template <typename T>
  inline TVec3<T> TVec3<T>::operator/(T rhs) const
  {
    rhs = ((T)1) / rhs;
    return TVec3<T>(x*rhs, y*rhs, z*rhs);
  }

  // unary operators
  template <typename T>
  inline TVec3<T> TVec3<T>::operator-() const
    { return TVec3<T>(-x,-y, -z); }

  template <typename T>
  inline bool TVec3<T>::operator==(const TVec3<T>& rhs) const
    { return (x==rhs.x && y==rhs.y && z==rhs.z); }

  template <typename T>
  inline bool TVec3<T>::operator!=(const TVec3<T>& rhs) const
    { return (x==rhs.x && y==rhs.y && z==rhs.z); }

  template <typename T>
  inline TVec3<T> operator*(T lhs, const TVec3<T>& rhs)
    { return TVec3<T>(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z); }


  template <typename T>
  inline T& TVec3<T>::operator[](int i)
    { return (&x)[i]; }

  template <typename T>
  inline T TVec3<T>::operator[](int i) const
    { return (&x)[i]; }


  // length of TVec3<T> by reference
  template <typename Tvec3>
  inline float Length(const Tvec3& u)
    { return sqrtf(u.x*u.x + u.y*u.y + u.z*u.z); }

  template <typename Tvec3>
  inline float LengthSqr(const Tvec3& u)
    { return u.x*u.x + u.y*u.y + u.z*u.z; }
  
  template <typename Tvec3>
  inline float Length2D(const Tvec3& u)
    { return sqrtf(u.x*u.x + u.y*u.y); }

  template <typename Tvec3>
  inline float LengthSqr2D(const Tvec3& u)
    { return u.x*u.x + u.y*u.y; }
  
  // distances between Tvec3 vectors
  template <typename Tvec3>
  inline float Dist(const Tvec3& lhs, const Tvec3& rhs)
    { return Length(lhs-rhs); }

  template <typename Tvec3>
  inline float DistSqr(const Tvec3& lhs, const Tvec3& rhs)
    { return LengthSqr(lhs-rhs); }
  
  template <typename Tvec3>
  inline float Dist2D(const Tvec3& lhs, const Tvec3& rhs)
    { return Length2D(lhs-rhs); }
  
  template <typename Tvec3>
  inline float DistSqr2D(const Tvec3& lhs, const Tvec3& rhs)
    { return LengthSqr2D(lhs-rhs); }

  // product, angles and normalzation
  template <typename Tvec3>
  inline float DotProd(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z; }

  template <typename Tvec3>
  inline float DotProd2D(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x*rhs.x + lhs.y*rhs.y; }

  template <typename Tvec3>
  inline float CosAngle(const Tvec3& lhs, const Tvec3& rhs)
  {
    float f = LengthSqr(lhs) * LengthSqr(rhs);
    f = sqrtf(f);
    return DotProd(lhs,rhs) * (1.0f/f);
  }
  
  template <typename Tvec3>
  inline float CosAngle2D(const Tvec3& lhs, const Tvec3& rhs)
  {
    float f = LengthSqr2D(lhs) * LengthSqr2D(rhs);
    f = sqrtf(f);
    return DotProd(lhs, rhs) * (1.0f/f);
  }

  template <typename Tvec3>
  inline void Normalize(Tvec3& v)
  {
    float f = Length(v); 
    v *= (1.0f/f);
  }

  template <typename Tvec3>
  inline Tvec3 CrossProd(const Tvec3& lhs, const Tvec3& rhs)
  {
    Tvec3 res;
    res.x = lhs.y*rhs.z - lhs.z*rhs.y;
    res.y = lhs.z*rhs.x - lhs.x*rhs.z;
    res.z = lhs.x*rhs.y - lhs.y*rhs.x;
    return res;
  }

  // per component operations of Tvec3
  template <typename Tvec3>
  inline Tvec3 Mul_ps(const Tvec3& lhs, const Tvec3& rhs)
  {
    Tvec3 u(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
    return u;
  }

  template <typename Tvec3>
  inline Tvec3 Max_ps(const Tvec3& lhs, const Tvec3& rhs)
  {
    Tvec3 u;
    u.x = lhs.x > rhs.x ? lhs.x : rhs.x;
    u.y = lhs.y > rhs.y ? lhs.y : rhs.y;
    u.z = lhs.z > rhs.z ? lhs.z : rhs.z;
    return u;
  }

  template <typename Tvec3>
  inline Tvec3 Min_ps(const Tvec3& lhs, const Tvec3& rhs)
  {
    Tvec3 u;
    u.x = lhs.x < rhs.x ? lhs.x : rhs.x;
    u.y = lhs.y < rhs.y ? lhs.y : rhs.y;
    u.z = lhs.z < rhs.z ? lhs.z : rhs.z;
    return u;
  }

  template <typename Tvec3>
  inline bool Less_ps(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x<rhs.x && lhs.y<rhs.y && lhs.z<rhs.z; }
  
  template <typename Tvec3>
  inline bool More_ps(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x>rhs.x && lhs.y>rhs.y && lhs.z>rhs.z; }

  template <typename Tvec3>
  inline bool LessEq_ps(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x<=rhs.x && lhs.y<=rhs.y && lhs.z<=rhs.z; }

  template <typename Tvec3>
  inline bool MoreEq_ps(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x>=rhs.x && lhs.y>=rhs.y && lhs.z>=rhs.z; }

  // linear interploation
  template <typename Tvec3>
  inline Tvec3 Lerp(const Tvec3& lhs, const Tvec3& rhs, float t)
    {return lhs*(1.0f-t) + rhs*t;}



  /******************************************************************************
  *                   INLINE FUNCTIONS (4x4 matrix)                              *
  ******************************************************************************/

  template <typename T>
  inline TMat4<T>& TMat4<T>::operator+=(const TMat4<T>& u)
  {
    m[0][0] += u.m[0][0]; m[0][1] += u.m[0][1]; m[0][2] += u.m[0][2]; m[0][3] += u.m[0][3];
    m[1][0] += u.m[1][0]; m[1][1] += u.m[1][1]; m[0][2] += u.m[1][2]; m[1][3] += u.m[1][3];
    m[2][0] += u.m[2][0]; m[2][1] += u.m[2][1]; m[0][2] += u.m[2][2]; m[2][3] += u.m[2][3];
    m[3][0] += u.m[3][0]; m[3][1] += u.m[3][1]; m[0][2] += u.m[3][2]; m[3][3] += u.m[3][3];
    return *this;
  }

  template <typename T>
  inline TMat4<T>& TMat4<T>::operator-=(const TMat4<T>& u)
  {
    (*this) += -u;
    return *this;
  }

  template <typename T>
  inline TMat4<T> TMat4<T>::operator+(const TMat4<T>& rhs)
  {
    TMat4 res(*this);
    res += rhs;
    return res;
  }

  template <typename T>
  inline TMat4<T> TMat4<T>::operator-(const TMat4<T>& rhs)
  {
    TMat4 res = -rhs;
    res += (*this);
    return res;
  }

  // multiplication, division
  template <typename T>
  inline TMat4<T>& TMat4<T>::operator*=(T f)
  {
    m[0][0] *= f; m[0][1] *= f; m[0][2] *= f; m[0][3] *= f;
    m[1][0] *= f; m[1][1] *= f; m[0][2] *= f; m[1][3] *= f;
    m[2][0] *= f; m[2][1] *= f; m[0][2] *= f; m[2][3] *= f;
    m[3][0] *= f; m[3][1] *= f; m[0][2] *= f; m[3][3] *= f;
  }

  template <typename T>
  inline TMat4<T>& TMat4<T>::operator/=(T f)
  {
    (*this) *= 1.0f / f;
    return *this;
  }

  template <typename T>
  inline TMat4<T> TMat4<T>::operator*(T f)
  {
    TMat4<T> res(*this);
    res *= f;
    return res;
  }

  template <typename T>
  inline TMat4<T> TMat4<T>::operator/(T f)
  {
    TMat4<T> res(*this);
    res *= 1.0f / f;
    return res;
  }

  // unary operators

  template <typename T>
  inline TMat4<T> TMat4<T>::operator-() const
  {
    TMat4 r;
    r.m[0][0] = -m[0][0]; r.m[0][1] = -m[0][1]; r.m[0][2] = -m[0][2]; r.m[0][3] = -m[0][3];
    r.m[1][0] = -m[1][0]; r.m[1][1] = -m[1][1]; r.m[0][2] = -m[1][2]; r.m[1][3] = -m[1][3];
    r.m[2][0] = -m[2][0]; r.m[2][1] = -m[2][1]; r.m[0][2] = -m[2][2]; r.m[2][3] = -m[2][3];
    r.m[3][0] = -m[3][0]; r.m[3][1] = -m[3][1]; r.m[0][2] = -m[3][2]; r.m[3][3] = -m[3][3];
    return r;
  }

  // matrix assignment
  template <typename T>
  inline TMat4<T>::TMat4(T* vals)
  {
    for (int r = 0; r<4; ++r)
    {
      for (int c = 0; c < 4; ++c)
        m[r][c] = vals[r * 4 + c];
    }
  }

  // Lengths of the matrix and the vector
  template <typename TMat4, typename T>
  inline void Lengths(const TMat4& rhs, T& length_mat, T& length_vec)
  {
    length_mat = length_vec = 0;
    for (int c = 0; c < 3; ++c)
    {
      length_vec += (T)(rhs.m[3][c]) * (T)(rhs.m[3][c]);
      for (int r = 0; r < 3; ++r)
        length_mat += (T)(rhs.m[r][c]) * (T)(rhs.m[r][c]);
    }

    length_mat = sqrt(length_mat);
    length_vec = sqrt(length_vec);
  }


  //set Matrix to Identity
  template <typename TMat4> 
  inline void MatrixIdentity(TMat4* rhs)
  {
    rhs->m[0][1] = rhs->m[0][2]  = rhs->m[0][3]  = 
    rhs->m[1][0] = rhs->m[1][2]  = rhs->m[1][3]  = 
    rhs->m[2][0] = rhs->m[2][1]  = rhs->m[1][3]  = 
    rhs->m[3][0] = rhs->m[3][1]  = rhs->m[3][2]  = 0;

    rhs->m[0][0] = rhs->m[1][1] = rhs->m[2][2] = rhs->m[3][3] = 1;
  }

  //Matrix Determinant
  template <typename TMat4> 
  inline double MatrixDeterminant(const TMat4* rhs)
  {
    return double(rhs->m[0][0] * rhs->m[1][1] * rhs->m[2][2] + 
                  rhs->m[0][1] * rhs->m[1][2] * rhs->m[2][0] + 
                  rhs->m[0][2] * rhs->m[1][0] * rhs->m[2][1] 
                - rhs->m[0][2] * rhs->m[1][1] * rhs->m[2][0] 
                - rhs->m[0][1] * rhs->m[1][0] * rhs->m[2][2]
                - rhs->m[0][0] * rhs->m[1][2] * rhs->m[2][1] );
  }


  // matrix Transpose
  template <typename TMat4>
  inline void Transpose(const TMat4& rhs, TMat4& res)
  {
    TMat4 tmp;
    for (int r = 0; r < 3; ++r)
    {
      for (int c = 0; c < 3; ++c)
        tmp.m[r][c] = rhs.m[c][r];
    }
    for (int c = 0; c<4; ++c)
      tmp.m[3][c] = rhs.m[3][c];
    for (int r = 0; r<3; ++r)
      tmp.m[r][3] = rhs.m[r][3];
    res = tmp;
  }

  // matrix multiplication
  template <typename T>
  inline void Multiply(const TMat4<T>& lhs, const TMat4<T>& rhs, TMat4<T>& res)
  {
    TMat4<T> tmp;
    for (int r=0; r<4; ++r)
    {
      for (int c=0; c<4; ++c)
        tmp.m[r][c] = lhs.m[r][0]*rhs.m[0][c] + lhs.m[r][1]*rhs.m[1][c] + lhs.m[r][2]*rhs.m[2][c] + lhs.m[r][3]*rhs.m[3][c];
    }
    res = tmp;
  }

  template <typename TMat4, typename TVec3>
  inline void TransposeLeftMultiply(const TMat4& lhs, const TVec3& rhs, TMat4& res)
  {
    TMat4 tmp;
    for (int r = 0; r<3; ++r)
    {
      for (int c = 0; c<4; ++c)
        tmp.m[r][c] = lhs.m[0][r] * rhs.m[0][c] + lhs.m[1][r] * rhs.m[1][c] + lhs.m[2][r] * rhs.m[2][c];
    }
    int r = 3;
    for (int c = 0; c<4; ++c)
      tmp.m[r][c] = lhs.m[r][0] * rhs.m[0][c] + lhs.m[r][1] * rhs.m[1][c] + lhs.m[r][2] * rhs.m[2][c] + lhs.m[r][3] * rhs.m[3][c];
    res = tmp;
  }

  template <typename TMat4, typename TVec3>
  inline void MultiplyVectorRightSide(const TMat4& lhs, const TVec3& rhs, TVec3& res)
  {
    TVec3 tmp;
    float fTmp[3] = { rhs.x , rhs.y, rhs.z };
    float fTmpR[3] = { 0 };
    for (int r = 0; r<3; ++r)
    {
      fTmpR[r] = lhs.m[r][0] * fTmp[0] + lhs.m[r][1] * fTmp[1] + lhs.m[r][2] * fTmp[2];
    }
    tmp.x = fTmpR[0];
    tmp.y = fTmpR[1];
    tmp.z = fTmpR[2];
    res = tmp;
  }

  // specialization for floats using SIMD
  template<> inline void Multiply(const TMat4<float>& lhs, const TVec3<float>& rhs, TMat4<float>& res);


  template <typename T>
  inline TMat4<T>& TMat4<T>::operator*=(const TMat4<T>& rhs)
  {
    Multiply(*this, rhs, *this);
    return *this;
  }
  
  template <typename T>
  inline TMat4<T> TMat4<T>::operator*(const TMat4<T>& rhs)
  {
    TMat4<T> res;
    Multiply(*this, rhs, res);
    return res;
  }


  template <typename T>
  inline bool TMat4<T>::operator==(const TMat4<T>& rhs) const
  {
    if (m[0][0] != rhs.m[0][0] || m[0][1] != rhs.m[0][1] || m[0][2] != rhs.m[0][2] || m[0][3] != rhs.m[0][3])
      return false;
    if (m[1][0] != rhs.m[1][0] || m[1][1] != rhs.m[1][1] || m[1][2] != rhs.m[1][2] || m[1][3] != rhs.m[1][3])
      return false;
    if (m[2][0] != rhs.m[2][0] || m[2][1] != rhs.m[2][1] || m[2][2] != rhs.m[2][2] || m[2][3] != rhs.m[2][3])
      return false;
    if (m[3][0] != rhs.m[3][0] || m[3][1] != rhs.m[3][1] || m[3][2] != rhs.m[3][2] || m[3][3] != rhs.m[3][3])
      return false;
    return true;
  }


  template <typename T>
  inline bool TMat4<T>::operator!=(const TMat4<T>& rhs) const
  {
    return !((*this) == rhs);
  }

  template <typename T>
  inline TMat4<T> operator*(T lhs, const TMat4<T>& rhs)
  {
    TMat4<T> res = rhs;
    res *= lhs;
    return res;
  }


} // namespace tpcl


#undef DLL_Entry
#endif // __tpcl_vec_H
