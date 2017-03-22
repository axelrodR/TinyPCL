// File Location: S:\gen\gengmtrx\gengmtrx_vec.h
/******************************************************************************
*
*: Package Name: gengmtrx
*
*: Title:
*
******************************************************************************/

#ifndef __gengmtrx_vec_H
#define __gengmtrx_vec_H

#include <math.h>

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#define DLL_Entry
#ifndef _LIB
#ifndef _LIB_LINK
#undef DLL_Entry

#ifdef GENGMTRX_EXPORTS
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
struct D3DXMATRIX;

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
    TVec3(float x_, float y_, float z_) { x=x_; y=y_; z=z_; }

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
  *                   INLINE FUNCTIONS (3D vector)                              *
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
    (*this) /= rhs;
    return (*this); 
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
  inline TVec3<T> operator*(float lhs, const TVec3<T>& rhs)
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
    { return lhs.x*rhs.x + lhs.y+rhs.y + lhs.z*rhs.z; }

  template <typename Tvec3>
  inline float DotProd2D(const Tvec3& lhs, const Tvec3& rhs)
    { return lhs.x*rhs.x + lhs.y+rhs.y; }

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
  inline TMat4<T> operator*(float lhs, const TMat4<T>& rhs)
  {
    TMat4<T> res = rhs;
    res *= lhs;
    return res;
  }


} // namespace tpcl


#undef DLL_Entry
#endif // __gengmtrx_vec_H
