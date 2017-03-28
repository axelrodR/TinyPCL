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
*: Package Name: plane
*
*: Title:
*
******************************************************************************/

#ifndef __gengmtrx_pln_H
#define __gengmtrx_pln_H


#include "../../include/vec.h"

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/


/******************************************************************************
*                        INCOMPLETE CLASS DECLARATIONS                        *
******************************************************************************/


/******************************************************************************
*                             EXPORTED CONSTANTS                              *
******************************************************************************/

namespace tpcl
{

/******************************************************************************
*                              EXPORTED CLASSES                               *
******************************************************************************/
  /******************************************************************************
  *
  *: Class name: CPlane
  *
  *: Abstract: 
  *
  ******************************************************************************/
  class CPlane
  {
  public:
    float m_A, m_B, m_C, m_D;       ///< plane parameters

    /// acccesors
    float& Get (int i)                             {return *(&m_A + i);}
    const float& Get (int i) const                 {return *(&m_A + i);}
    const CVec3& GetOrtho() const                  {return *(CVec3*)&m_A;} ///< get orthogonal vector
    //const D3DXVECTOR4& GetAsVec4() const           {return *(D3DXVECTOR4*)&m_A;} ///< get plane as 4-vector
    void Set(float* Xi_f)                          {m_A=Xi_f[0]; m_B=Xi_f[1]; m_C=Xi_f[2]; m_D=Xi_f[3];}  ///< set
    void Set(float Xi_A,float Xi_B,float Xi_C,float Xi_D)   {m_A=Xi_A; m_B=Xi_B; m_C=Xi_C; m_D=Xi_D;}     ///< set

    /** get the normal to the plane */
    CVec3 GetNormal() const                         { CVec3 v=GetOrtho(); tpcl::Normalize(v); return v; }

    /** convert plane to normal form (i.e. A,B,C are the normal vector) */
    void Normalize_m()
    {
      CVec3* v = (CVec3*)&m_A;
      float l = 1.0f / Length(*v);
      *v *= l;
      m_D *= l;
    }

    /** Get height at (x,y) point (m_C must be non zero) */
    float GetHeightAt(float x, float y) const      {return -(m_A*x + m_B*y + m_D) / m_C;}

    /** Get a point on the plane */
    CVec3 GetPointInside() const;

    /** Get the two directions on the plane */
    void GetDirsInside(CVec3& l_dir1, CVec3& l_dir2) const;


    /** @brief best fit a plane using mean-least-squares
     * 
     * @param Xi_normalize    normalize the plane (i.e. A,B,C are the normal vector)
     * @return fitness parameter (average distance of points from the plane),
     *         -1 if no plane can be fitted (i.e. all point lie on the same line)
     */
    float BestFit(int Xi_numPts, const CVec3* Xi_pts, bool Xi_normalize=true);


    /** @brief best fit plane using RanSac like algorithm 
     * @param Xi_thresholdDist      threshold distance for inliers
     * @param Xi_exitGrade          eaarly out will happen when the exit grade is reached
     * @return  the percentage of points within threshold distance
     */
    float RanSaC(int Xi_numPts, const CVec3* Xi_pts,
                 float Xi_thresholdDist, float Xi_exitGrade = 0.9f);


    /** @brief best fit plane using RanSac like algorithm
     * plane models are constrained to be within a certain angle of some constrainet normal
     * @param Xi_constraintAngle    cosine of constraint angle
     * @param Xi_contraint          normalized constraint normal
     * @return  the percentage of points within threshold distance
     */
    float RanSacConstrained(int Xi_numPts, const CVec3* Xi_pts, 
                            const CVec3& Xi_contraint, const float Xi_constraintAngle,
                            float Xi_thresholdDist, float Xi_exitGrade = 0.9f);

    /** distance of a point from (non normalized) plane 
     * @param Xi_sign   false=absolute disntace,  true=keep sign of distance compared to plane normal)
     */
    float Dist(const CVec3& Xi_pt, bool Xi_sign=false) const;

    /** Get closest point on the plane 
     * @return distance of Xi_pt from the plane (with sign compared to plane normal) */
    float GetClosestPt(const CVec3& Xi_pt,CVec3& Xo_closest) const;


    /** Filters points with distance below threshold from the plane
     * @param Xi_pt     points to filter
     * @param Xo_filteredPt     buffer for filtered points (pointer can be the same as Xi_pt)
     * @param Xi_AboveBelow     filter bitfield controls which points to keep: 0x1=close above, 0x2=close below, 0x4=far above, 0x8=far below
     * @return  number of points within threshold */
    int FilterDistTooHigh(int Xi_numPts, const CVec3* Xi_pt, CVec3* Xo_filteredPt,
                          float Xi_thresholdDist=0.0f, int Xi_AboveBelow=0, bool Xi_justCount=false) const;

    
    /** Intersction of two planes, as a line through two points
     * @return 0=parralel disjoint, -1=coincide, 1=line intersection */
    int Intersect(const CPlane& Xi_pl, CVec3& Xo_pt1, CVec3& Xo_pt2) const;


    /** Intersction of a plane and a segment/line
     * @param Xi_segment    true = treat the two points as segment, false = line through those point
     * @return 0=disjoint, -1=line lies in plane, 1=intersect at a single point */
    int Intersect(const CVec3& Xi_pt1, const CVec3& Xi_pt2, float& Xo_intersectPt, bool Xi_segment=true) const;


    // constrctors
    CPlane (){}
    CPlane (float a, float b, float c, float d)    {m_A=a; m_B=b; m_C=c; m_D=d;}
    CPlane (float* Xi_vec)                 {m_A=Xi_vec[0]; m_B=Xi_vec[1]; m_C=Xi_vec[2]; m_D=Xi_vec[3];}

    // destructor
    ~CPlane (){}
  };



  /******************************************************************************
  *
  *: Class name: CPlane
  *
  *: Abstract: simple structure to hold normal in compact form
  *
  ******************************************************************************/
  struct CCompactPlane
  {
    short x,y,z,w;                          ///< plane parameters
    static const int FIX_PT_NORMAL = 1024;  ///< normals and curvature are written as 10 bit fixed point
    static const int FIX_PT_D = 64;         ///< D parameter is is written as 6 bit fixed point

    /** set from plane */
    void Set(const CPlane& pl)
    {
      x = short(pl.m_A*FIX_PT_NORMAL), y = short(pl.m_B*FIX_PT_NORMAL); z = short(pl.m_C*FIX_PT_NORMAL);
      int h = int(pl.m_D*FIX_PT_D);
      if ( ((unsigned int)h) > 0xfff0)
        h = 0xffff;
      w = short(h);
    }

    /** get a plane from compact form
     * @return true = m_D was encoded, false=m_D was not encoded correctly
     */
    bool Get(CPlane& pl)
    {
      const float l_invFP1 = 1.0f / FIX_PT_NORMAL,l_invFP2 = 1.0f / FIX_PT_D;
      pl.m_A = float(x)*l_invFP1,    pl.m_B = float(y)*l_invFP1;
      pl.m_C = float(z)*l_invFP1,    pl.m_D = float(w)*l_invFP2;
      return (w != 0xffff);
    }

    /** constrcutor */
    CCompactPlane(){}
    CCompactPlane(float Xi_x, float Xi_y, float Xi_z, float Xi_d)
    {
      x = short(Xi_x*FIX_PT_NORMAL), y = short(Xi_y*FIX_PT_NORMAL);
      z = short(Xi_z*FIX_PT_NORMAL), w = short(Xi_d*FIX_PT_D);
    }
    CCompactPlane(const CPlane& Xi_pl) {Set(Xi_pl);}

    /** accessor */
    short& operator[](int i)                  {return *(&x + i);}
    const short& operator[](int i) const      {return *(&x + i);}
  };


/******************************************************************************
*                            EXPORTED FUNCTIONS                               *
******************************************************************************/

} //namespace tpcl

#endif
