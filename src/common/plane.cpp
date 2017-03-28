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
*: Package Name: gengmtrx_pln
*
******************************************************************************/

#include "plane.h"
#include "common.h"
#include "float.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

namespace tpcl
{

  /******************************************************************************
  *                             INTERNAL CONSTANTS                              *
  ******************************************************************************/

  const float MIN_SCALE = 0.0001f;   //   the minimal scale of interest is 0.1mm
  const float MIN_SCALE_SQR = MIN_SCALE * MIN_SCALE;

  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

  /******************************************************************************
  *                       FORWARD FUNCTION DECLARATIONS                         *
  ******************************************************************************/

  /******************************************************************************
  *                             STATIC VARIABLES                                *
  ******************************************************************************/

  /******************************************************************************
  *                      CLASS STATIC MEMBERS INITIALIZATION                    *
  ******************************************************************************/

  /******************************************************************************
  *                              INTERNAL CLASSES                               *
  ******************************************************************************/


  /** Fast random number generator 
   *
   * Marsaglia's xorshf. Characteristics that are good for RanSac
   * see: https://en.wikipedia.org/wiki/Xorshift
   */
  class Xorshift128
  {
    /* These state variables must be initialized so that they are not all zero. */
    unsigned int x, y, z, w;
  public:
    unsigned int rand()
    {
      unsigned int t = x;
      t ^= t << 11;
      t ^= t >> 8;
      x = y; y = z; z = w;
      w ^= w >> 19;
      w ^= t;
      return w;
    }


    /** random number in an interval [0,maxRange] */
    unsigned int rand(int Xi_maxRange)
    {
      unsigned int rI = rand() & ((unsigned int)0x0fffffff);  // 28 bit random integer
      double a = double(rI) * s_div28Bit * Xi_maxRange;
      return unsigned int(a);
    }
    
    Xorshift128()
    {
      x = 123456789;
      y = 362436069;
      z = 521288629;
      w = 175262959;
      s_div28Bit = 1.0 / (double(RAND_MAX_28bit) + 1.0);
    }

    Xorshift128(const unsigned int* Xi_seeds)
    {
      x = Xi_seeds[0];
      y = Xi_seeds[1];
      z = Xi_seeds[2];
      w = Xi_seeds[3];
      s_div28Bit = 1.0 / (double(RAND_MAX_28bit) + 1.0);
    }

    Xorshift128(const float* Xi_vecAsSeed)
    {
      float f[4];
      // for the seed take only 
      f[0] = floorf(Xi_vecAsSeed[0] * 10.0f) * 0.1f;
      f[1] = floorf(Xi_vecAsSeed[1] * 10.0f) * 0.1f;
      f[2] = floorf(Xi_vecAsSeed[2] * 10.0f) * 0.1f;
      f[3] = floorf(Xi_vecAsSeed[3] * 10.0f) * 0.1f;
      x = (*(unsigned int*)f) ^ 123456789;
      y = (*(unsigned int*)(f + 1)) ^ 362436069;
      z = (*(unsigned int*)(f + 2)) ^ 521288629;
      w = (*(unsigned int*)(f + 3)) ^ 175262959;
      s_div28Bit = 1.0 / (double(RAND_MAX_28bit) + 1.0);
    }

  protected:
    // for this application specificically we want only 28bits
    static const unsigned int RAND_MAX_28bit = 1<<28;
    static const unsigned int MASK_28bits = RAND_MAX_28bit - 1;
    static double s_div28Bit;
  };
  double Xorshift128::s_div28Bit;


  /** estimate the number of runs required for RANSAC to complete */
  int EstimateRanSacIters(float Xi_inliersRatio, float Xi_probSuccess, int Xi_sampleSize)
  {
    double l_numIterD = log(1.0 - Xi_probSuccess) / log(1.0 - pow(Xi_inliersRatio,Xi_sampleSize));
    int l_numIt = int(ceil(l_numIterD));
    return l_numIt;
   }


  template<class T> T MinT(const T& a, const T& b, const T& c) { return MinT(a, MinT(b,c)); }


  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CPlane
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: BestFit
  *
  ******************************************************************************/
  float CPlane::BestFit(int Xi_numPts, const CVec3* Xi_pts, bool Xi_normalize)
  {
    if (Xi_numPts < 3)
    {
      m_A=m_B=m_D=0; m_C=1.0f;
      return FLT_MAX;
    }

    // compute center and covariance matrix
    double xx = 0.0, xy = 0.0, xz = 0.0;  // covariance matrix (off-diagonal, excluding symmetries)
    double yy = 0.0, yz = 0.0, zz = 0.0;  // covariance matrix (diagonal terms)
    CVec3 l_cent(0,0,0);
    for (int j=0; j<Xi_numPts; ++j)
      l_cent += Xi_pts[j];
    l_cent *= (float)(1.0 / Xi_numPts);
    for (int j=0; j<Xi_numPts; ++j)
    {
      CVec3 v = Xi_pts[j] - l_cent;
      xx += v.x * v.x;  yy += v.y * v.y;  zz += v.z * v.z;
      xy += v.x * v.y;  xz += v.x * v.z;  yz += v.y * v.z;
    }

    // try fitting a plane in the Z direction (we prefer plane fitting in the z direction)
    double l_detZ = xx * yy - xy * xy;
    double l_detX = yy*zz - yz*yz;
    double l_detY = xx*zz - xz*xz;
    double l_maxDet = MaxT(MaxT(l_detX,l_detY),l_detZ);
    if (l_maxDet <= MIN_SCALE_SQR)
    {
      m_A = m_B = m_D = 0; m_C = 1.0f;
      return FLT_MAX;   // points are collinear
    }
    if (l_maxDet == l_detX)
    {
      double D = 1.0 / l_detX;
      m_A = 1.0f;
      m_B = float( (xz*yz - xy*zz) * D );
      m_C = float( (xy*yz - xz*yy) * D );
    } 
    else if (l_maxDet == l_detY)
    {
      double D = 1.0 / l_detY;
      m_A = float( (yz*xz - xy*zz) * D );
      m_B = 1.0f;
      m_C = float( (xy*xz - yz*xx) * D );
    } 
    else /** Z is still the best choice */
    {
      double D = 1.0 / l_detZ;
      m_A = float( (yz * xy - xz * yy) * D );
      m_B = float( (xy * xz - xx * yz) * D );
      m_C = 1.0f;
    }

    // plane offset
    CVec3 l_norm;
    float l_D;
    if (Xi_normalize)
    {
      Normalize_m();
      l_norm = GetOrtho();
      l_D = m_D = -DotProd(l_norm, l_cent);
    }
    else // compute a normalized plane locally
    {
      l_norm = GetOrtho();
      m_D = -DotProd(l_norm, l_cent);
      float f = (float)(1.0 / Length(l_norm));
      l_norm *= f;
      l_D = m_D *= f;
    }
    float fitness = 0.0f;
    for (int j=0; j<Xi_numPts; ++j)
      fitness += fabsf(DotProd(l_norm, Xi_pts[j]) + l_D);
    fitness /= Xi_numPts;
    return fitness;
  }


  float CPlane::RanSaC(int Xi_numPts, const CVec3* Xi_pts, float Xi_thresholdDist, float Xi_exitGrade)
  {
    if (Xi_numPts < 3)
      return 0;

    return RanSacConstrained(Xi_numPts, Xi_pts, CVec3(0,0,0), -1, Xi_thresholdDist, Xi_exitGrade); // run with no constraint
  }
    


  float CPlane::RanSacConstrained(int Xi_numPts, const CVec3* Xi_pts, 
                          const CVec3& Xi_contraint, const float Xi_constraintAngle,
                          float Xi_thresholdDist, float Xi_exitGrade)
  {
    if (Xi_numPts < 3)
      return 0;

    // constants
    const int MAX_ITERS = 100;
    const int RANDOM_SAMPLE_SIZE = 4;
    const float PROB_SUCCESS = 0.995f; // required success rate
    const float l_exitgrade = Xi_exitGrade * Xi_numPts;
    // factor for reducing score score with ditance from plane (within threshold distance). A deviation from the normal Ransac.
    const float iFactor = float(1.0 / (3.0 * Xi_thresholdDist));

    CVec3 l_sample[RANDOM_SAMPLE_SIZE];     // random samples holder
    //Xorshift128 l_rand();                       // seed is the same for all polygons.
    Xorshift128 l_rand((float*)Xi_pts);           // seed from the points array. I.e. changes every time
    float l_bestScore = 0;
    int l_numIter = MAX_ITERS;

    // random concensus iteration
    for (int i=0; i<l_numIter; i++)
    {
      // choose a random sample
      for (int i=0; i<RANDOM_SAMPLE_SIZE; i++)
      {
        unsigned int l_rId = l_rand.rand(Xi_numPts);
        l_sample[i] = Xi_pts[l_rId];
      }
      // fit a plane
      CPlane l_plane;
      float f = l_plane.BestFit(RANDOM_SAMPLE_SIZE, l_sample, true);
      if (f > 100)
        continue;
      const CVec3* l_norm = &(l_plane.GetOrtho());
      float l_cs = DotProd(Xi_contraint, *l_norm);
      if (fabsf(l_cs) < Xi_constraintAngle)
        continue;   // plane not in the right direction
      if (l_cs < 0)
        l_plane.Set(-l_plane.m_A,-l_plane.m_B,-l_plane.m_C,-l_plane.m_D);   // facing opposite direction so reverse normal

      // calculate the rank (based on how many points are within threshold of that plane)
      float l_score = 0;
      int l_numInliers = 0;
      for (int j=0; j<Xi_numPts; j++)
      {
        float l_dist = fabsf(DotProd(*l_norm, Xi_pts[j]) + l_plane.m_D);
        if (l_dist > Xi_thresholdDist)
          continue;
        l_numInliers++;
        l_score += 1.0f - (l_dist * iFactor);
      }

      // better plane than best plane
      if (l_bestScore < l_score)
      {
        l_bestScore = l_score;
        *this = l_plane;
        if (l_bestScore > l_exitgrade)
          break;    // good enough to stop

        // estimate the number of iterations needed based on the inliers ratio
        float l_inlierRatio = float(double(l_numInliers) / double(Xi_numPts));
        int l_numIt = EstimateRanSacIters(l_inlierRatio, PROB_SUCCESS, RANDOM_SAMPLE_SIZE);
        l_numIter = MinT(l_numIt, l_numIter, MAX_ITERS);
      }
    }
    return l_bestScore / Xi_numPts;
  }


  float CPlane::Dist(const CVec3& Xi_pt,bool Xi_sign) const
  {
    float norm = Length(GetOrtho());
    float dot = DotProd(Xi_pt, GetOrtho()) + m_D;
    float dist = dot * (1.0f / norm);
    if (!Xi_sign)
      return fabsf(dist);
    else
      return dist;
  }



  float CPlane::GetClosestPt(const CVec3& Xi_pt,CVec3& Xo_closest) const
  {
    float norm = Length(GetOrtho());
    float dot = DotProd(Xi_pt,GetOrtho()) + m_D;
    float dist = dot * (1.0f / norm);
    Xo_closest = Xi_pt - GetOrtho() * dist;
    return dist;
  }



  int CPlane::FilterDistTooHigh(int Xi_numPts, const CVec3* Xi_pt, CVec3* Xo_filteredPt,
                          float Xi_thresholdDist, int Xi_AboveBelow, bool Xi_justCount) const
  {
    const CVec3& l_ortho = GetOrtho();
    float l_norm = Length(l_ortho);
    float l_invNorm = 1.0f / l_norm;

    int l_n = 0;
    for (int i=0; i<Xi_numPts; ++i)
    {
      float l_dot = DotProd(Xi_pt[i], l_ortho) + m_D;
      float l_dist = l_dot * l_invNorm;
      if (fabsf(l_dist) > Xi_thresholdDist)
      {
        if (l_dist > 0) // above?
        {
          if ( !(Xi_AboveBelow & 0x4) )       // far above flag bit set?
            continue;
        }
        else if ( !(Xi_AboveBelow & 0x8) )    // far below flag bit set?
          continue;
      }
      else // near
      {
        if (l_dist >= 0) // above?
        {
          if ( !(Xi_AboveBelow & 0x1) )       // far above flag bit set?
            continue;
        }
        else if ( !(Xi_AboveBelow & 0x2) )    // far below flag bit set?
          continue;
      }
      if (!Xi_justCount)
        Xo_filteredPt[l_n] = Xi_pt[i];
      l_n++;
    }
    return l_n;
  }


  CVec3 CPlane::GetPointInside() const
  {
    if (fabsf(m_C) > 0.001f)
      return CVec3(0,0, -m_D/m_C);
    else if (fabsf(m_B) > 0.001f)
      return CVec3(0,-m_D/m_B,0);
    return CVec3(-m_D/m_A,0,0);
  }


  void CPlane::GetDirsInside(CVec3& l_dir1, CVec3& l_dir2) const
  {
    CVec3 l_pt;
    if (fabsf(m_C) > 0.001f)
    {
      l_pt = CVec3(0,0, -m_D/m_C);
      l_dir1 = CVec3(1,0, -m_D/m_C) - l_pt;
      l_dir2 = CVec3(0,1, -m_D/m_C) - l_pt;
    }
    else if (fabsf(m_B) > 0.001f)
    {
      l_pt = CVec3(0,-m_D/m_B,0);
      l_dir1 = CVec3(0,-m_D/m_B,1) - l_pt;
      l_dir2 = CVec3(1,-m_D/m_B,0) - l_pt;
    }
    else
    {
      l_pt = CVec3(-m_D/m_A,0,0);
      l_dir1 = CVec3(-m_D/m_A,1,0) - l_pt;
      l_dir2 = CVec3(-m_D/m_A,0,1) - l_pt;
    }
    Normalize(l_dir1);
    Normalize(l_dir2);
  }


  /** based on: http://geomalgorithms.com/a05-_intersect-1.html */
  int CPlane::Intersect(const CPlane& Xi_pl, CVec3& Xo_pt1, CVec3& Xo_pt2) const
  {
    CVec3 l_pl1N = GetOrtho();
    CVec3 l_pl2N = Xi_pl.GetOrtho();
    CVec3 l_cs = CrossProd(l_pl1N, l_pl2N);    // TODO: maybe make cross product with double for more precision.
    float ax = (l_cs.x >= 0 ? l_cs.x : -l_cs.x);
    float ay = (l_cs.y >= 0 ? l_cs.y : -l_cs.y);
    float az = (l_cs.z >= 0 ? l_cs.z : -l_cs.z);

    // test if the two planes are parallel
    if ((ax+ay+az) < 0.1f)
    {
      // get a point on the plane and test if it is also on the other plane
      CVec3 l_pt = GetPointInside() - Xi_pl.GetPointInside();
      if (DotProd(l_pl1N, l_pt) == 0)
        return -1; // planes coincide
      else 
        return 0; // planes are disjoint
    }

    // we have aline. We need to compute a point on the line
    float l_d1 = m_D, l_d2 = Xi_pl.m_D;
    float l_max = MaxT(ax,MaxT(ay,az));
    if (l_max == ax)
    {
      Xo_pt1.x = 0;
      Xo_pt1.y = (l_d2*l_pl1N.z - l_d1*l_pl2N.z) /  l_cs.x;
      Xo_pt1.z = (l_d1*l_pl2N.y - l_d2*l_pl1N.y) /  l_cs.x;
    }
    else if (l_max == ay)
    {
      Xo_pt1.x = (l_d1*l_pl2N.z - l_d2*l_pl1N.z) /  l_cs.y;
      Xo_pt1.y = 0;
      Xo_pt1.z = (l_d2*l_pl1N.x - l_d1*l_pl2N.x) /  l_cs.y;
    }
    else // (l_max == az)
    {
      Xo_pt1.x = (l_d2*l_pl1N.y - l_d1*l_pl2N.y) /  l_cs.z;
      Xo_pt1.y = (l_d1*l_pl2N.x - l_d2*l_pl1N.x) /  l_cs.z;
      Xo_pt1.z = 0;
    }

    Xo_pt2 = Xo_pt1 + l_cs;
    return 1;
  }


  /** based on: http://geomalgorithms.com/a05-_intersect-1.html */
  int CPlane::Intersect(const CVec3& Xo_pt1, const CVec3& Xo_pt2, float& Xo_intersectPt, bool Xi_segment) const
  {
    CVec3 u = Xo_pt2 - Xo_pt1;
    CVec3 w = Xo_pt1 - GetPointInside();

    CVec3 l_n = GetNormal();
    float D = DotProd(l_n, u);
    float N = -DotProd(l_n, w);

    if (fabsf(D) < 0.00001) {  // segment is parallel to plane
      if (N == 0)    // segment lies in plane
        return -1;
      else
        return 0;    // no intersection
    }
    // they are not parallel - compute intersect param
    Xo_intersectPt = N / D;
    if (!Xi_segment)
      return 1;
    // check if the intersection point is in the segment
    if (Xo_intersectPt < 0.0f || Xo_intersectPt > 1.0f)
      return 0;    // no intersection
    return 1;
  }



  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/

  /******************************************************************************
  *                              Private methods                                *
  ******************************************************************************/

  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/

  /******************************************************************************
  *                            INTERNAL FUNCTIONS                               *
  ******************************************************************************/

} // namespace tpcl
