#include "RegICP.h"
#include "SpatialHash.h"
#include "common.h"
#include <vector>
#include "../../include/vec.h"
#include "../include/ptCloud.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


namespace tpcl
{
  /** Try to find match in point cloud for another point.
  * Currently we match by finding closest point within a given threshold radius ("inliers")
  * @param in_pcl1    	    Spatial hashing of main point cloud.
  * @param in_p2            a point from 2nd point cloud.
  * @param in_normal        !!curently not used!! - normal of the point to which we want to find a match.
  * @param out_match         match found.
  * @param indist           inline distance. 
  * return                  true if a match was found, flase otherwise*/
  bool MatchPoint(const CSpatialHash2D& in_pcl1, const CVec3& in_p2, const CVec3& in_normal, const double in_distThreshold, CVec3& out_match, double& out_dist)
  {
    // go over retrieved list (if size = 0, return false), check if closer than in_distThreshold, update out_match if better normal match. return true.
    if (in_pcl1.FindNearest(in_p2, &out_match, float(in_distThreshold))) // search nearest neighbor
    {
      // check if it is an inlier
      out_dist = Dist(in_p2, out_match);
      if (out_dist < in_distThreshold)
        return true;
    }
    // if no match/neighbor found closer than indist, return false
    return false;
  }


  /** calculates pythagoras output = sqrt(a^2 + b^2) */
  double pythag(double a, double b)
  {
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
      return double(absa*sqrt(1.0 + pow(absb / absa, 2)));
    else
      return double((absb == 0.0 ? 0.0 : absb*sqrt(1.0 + pow(absa / absb, 2))));
  }


  // build bidiagonal form of using Householder reduction
  void BiDiag(double io_U[3][3], double out_W[3], double out_V[3][3], double rv1[4])
  {
    // Householder reduction to bidiagonal form.
    int i, j, k;
    double f, g = 0, h, s;

    for (i = 0; i<3; i++)
    {
      rv1[i] = g;
      g = s = 0.0;

      // act on columns
      for (k = i; k<3; k++)
        s += io_U[k][i] * io_U[k][i];
      if (s)
      {
        f = io_U[i][i];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        io_U[i][i] = f - g;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i; k<3; k++)
            s += io_U[k][i] * io_U[k][j];
          f = s / h;
          for (k = i; k<3; k++)
            io_U[k][j] += f*io_U[k][i];
        }
      }
      out_W[i] = g;

      // act on rows
      g = s = 0.0;
      for (k = i + 1; k < 3; k++)
        s += io_U[i][k] * io_U[i][k];
      if (s)
      {
        f = io_U[i][i + 1];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        io_U[i][i + 1] = f - g;
        for (k = i + 1; k<3; k++)
          rv1[k] = io_U[i][k] / h;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += io_U[j][k] * io_U[i][k];
          for (k = i + 1; k<3; k++)
            io_U[j][k] += s*rv1[k];
        }
      }
    }

    // Accumulation of right-hand transformations.
    out_V[2][2] = 1.0;
    g = rv1[2];
    for (i = 2; i >= 0; i--)
    {
      if (g)
      {
        for (j = i + 1; j<3; j++)
          out_V[j][i] = (io_U[i][j] / io_U[i][i + 1]) / g;   // Double division to avoid possible underflow.
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += io_U[i][k] * out_V[k][j];
          for (k = i + 1; k<3; k++)
            out_V[k][j] += s*out_V[k][i];
        }
      }
      for (j = i + 1; j<3; j++)
        out_V[i][j] = out_V[j][i] = 0.0;
      out_V[i][i] = 1.0;
      g = rv1[i];
    }

    // Accumulation of left-hand transformations.
    for (i = 2; i >= 0; i--)
    {
      g = out_W[i];
      for (j = i + 1; j<3; j++)
        io_U[i][j] = 0.0;
      if (g)
      {
        g = 1.0 / g;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += io_U[k][i] * io_U[k][j];
          f = (s / io_U[i][i])*g;
          for (k = i; k<3; k++)
            io_U[k][j] += f*io_U[k][i];
        }
        for (j = i; j<3; j++)
          io_U[j][i] *= g;
      }
      else
        for (j = i; j<3; j++)
          io_U[j][i] = 0.0;
      ++io_U[i][i];
    }
  }



  // compute the SVD of a 3x3 matrix
  bool svd3x3(const double in_M[3][3], double out_U[3][3], double out_W[3], double out_V[3][3])
  {
    double w[3];
    double rv1[3];

    int flag, i, its, j, jj, k, nm;
    double c, f, g, h, s, x, y, z;

    for (i = 0; i < 9; ++i)
      out_U[0][i] = in_M[0][i];

    // convert to bidiagonal form
    BiDiag(out_U, w, out_V, rv1);

    // claculate scale of stuff
    double anorm = fabs(w[0]) + fabs(rv1[0]);
    anorm = MaxT(anorm, fabs(w[1]) + fabs(rv1[1]));
    anorm = MaxT(anorm, fabs(w[2]) + fabs(rv1[2]));

    // Diagonalization of the bidiagonal form: Loop over singular values
    for (k = 2; k >= 0; k--)
    {
      // and over allowed iterations.
      for (its = 0; its<30; its++)
      {
        flag = 1;
        int l;
        // Test for splitting.
        for (l = k; l >= 0; l--)
        {
          nm = l - 1;
          if ((double)(fabs(rv1[l]) + anorm) == anorm)
          {
            flag = 0;
            break;
          }
          if ((double)(fabs(w[nm]) + anorm) == anorm)
          {
            break;
          }
        }
        if (flag)
        {
          c = 0.0; // Cancellation of rv1[l], if l > 1.
          s = 1.0;
          for (i = l; i <= k; i++)
          {
            f = s*rv1[i];
            rv1[i] = c*rv1[i];
            if ((double)(fabs(f) + anorm) == anorm)
              break;
            g = w[i];
            h = pythag(f, g);
            w[i] = h;
            h = 1.0 / h;
            c = g*h;
            s = -f*h;
            for (j = 0; j<3; j++)
            {
              y = out_U[j][nm];
              z = out_U[j][i];
              out_U[j][nm] = y*c + z*s;
              out_U[j][i] = z*c - y*s;
            }
          }
        }
        z = w[k];
        if (l == k)
        { // Convergence.
          if (z<0.0)
          { // Singular value is made nonnegative.
            w[k] = -z;
            for (j = 0; j<3; j++)
              out_V[j][k] = -out_V[j][k];
          }
          break;
        }
        if (its == 29)
          return false; // No convergence
        x = w[l]; // Shift from bottom 2-by-2 minor.
        nm = k - 1;
        y = w[nm];
        g = rv1[nm];
        h = rv1[k];
        f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.f*h*y);
        g = pythag(f, 1.0);
        f = ((x - z)*(x + z) + h*((y / (f + SIGN(g, f))) - h)) / x;
        c = s = 1.0; // Next QR transformation:
        for (j = l; j <= nm; j++)
        {
          i = j + 1;
          g = rv1[i];
          y = w[i];
          h = s*g;
          g = c*g;
          z = pythag(f, h);
          rv1[j] = z;
          c = f / z;
          s = h / z;
          f = x*c + g*s;
          g = g*c - x*s;
          h = y*s;
          y *= c;
          for (jj = 0; jj<3; jj++)
          {
            x = out_V[jj][j];
            z = out_V[jj][i];
            out_V[jj][j] = x*c + z*s;
            out_V[jj][i] = z*c - x*s;
          }
          z = pythag(f, h);
          w[j] = z; // Rotation can be arbitrary if z = 0.
          if (z)
          {
            z = 1.0 / z;
            c = f*z;
            s = h*z;
          }
          f = c*g + s*y;
          x = c*y - s*g;
          for (jj = 0; jj<3; jj++)
          {
            y = out_U[jj][j];
            z = out_U[jj][i];
            out_U[jj][j] = y*c + z*s;
            out_U[jj][i] = z*c - y*s;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }

    // sort singular values and corresponding columns of u and v
    // by decreasing magnitude. Also, signs of corresponding columns are
    // flipped so as to maximize the number of positive elements.
    int s2;
    double   sw;
    double su[3];
    double sv[3];
    for (i = 1; i<3; i++)
    {
      sw = w[i];
      for (k = 0; k<3; k++)
        su[k] = out_U[k][i];
      for (k = 0; k<3; k++)
        sv[k] = out_V[k][i];
      j = i;
      while (w[j - 1] < sw)
      {
        w[j] = w[j - 1];
        for (k = 0; k<3; k++)
          out_U[k][j] = out_U[k][j - 1];
        for (k = 0; k<3; k++)
          out_V[k][j] = out_V[k][j - 1];
        j -= 1;
        if (j < 1) break;
      }
      w[j] = sw;
      for (k = 0; k<3; k++)
        out_U[k][j] = su[k];
      for (k = 0; k<3; k++)
        out_V[k][j] = sv[k];
    }

    // flip signs
    for (k = 0; k<3; k++)
    {
      s2 = 0;
      for (i = 0; i<3; i++)
        if (out_U[i][k] < 0.0)
          s2++;
      for (j = 0; j<3; j++)
        if (out_V[j][k] < 0.0)
          s2++;
      if (s2 > 3)
      {
        for (i = 0; i<3; i++)
          out_U[i][k] = -out_U[i][k];
        for (j = 0; j<3; j++)
          out_V[j][k] = -out_V[j][k];
      }
    }

    // create vector and copy singular values
    for (int r = 0; r < 3; r++)
      out_W[r] = w[r];
    return true;
  }


  bool svd3x3(double* in_M, CMat4& out_U, CMat4& out_W, CMat4& out_V)
  {
    const int M1size = 3;
    double M[M1size][M1size] = { 0 };
    double U[M1size][M1size] = { 0 };
    double W[M1size] = { 0 };
    double V[M1size][M1size] = { 0 };

    for (int row = 0; row < M1size; row++)
    {
      for (int col = 0; col < M1size; col++)
      {
        M[row][col] = in_M[row*M1size + col];
      }
    }
    bool convergence = svd3x3(M, U, W, V);

    if (convergence)
    {
      for (int row = 0; row < M1size; row++)
      {
        for (int col = 0; col < M1size; col++)
        {
          out_U.m[row][col] = (float)U[row][col];
          out_W.m[row][col] = 0.f;
          out_V.m[row][col] = (float)V[row][col];
        }
        out_W.m[row][row] = (float)W[row];
      }

      for (int i = 0; i < M1size; ++i)
      {
        out_U.m[M1size][i] = out_U.m[i][M1size] = 0.f;
        out_W.m[M1size][i] = out_W.m[i][M1size] = 0.f;
        out_V.m[M1size][i] = out_V.m[i][M1size] = 0.f;
      }
      out_U.m[M1size][M1size] = out_W.m[3][3] = out_V.m[M1size][M1size] = 1.f;
    }

    return convergence;
  }


  typedef TVec3<double> CVec3D;

  void PerformIter(CSpatialHash2D& in_pcl1, const CPtCloud& in_pcl2, CMat4& io_Rt, const float in_regRes, double& out_transformationChange, double& out_PreviousFitnessScore)
  {
    //double l_distThreshold = 2 * in_regRes;
    double l_scoreDistThreshold = 2 * in_regRes;
    double l_regDistThreshold = 2 * in_regRes;   //l_regDistThreshold <= l_scoreDistThreshold
    double l_accError = 0;
    // extract matrix and translation vector
    CVec3D l_massCenter1(0, 0 ,0), l_massCenter2(0, 0, 0);  // center of masses for both clouds (as doubles)
    CVec3 l_massCenter1f, l_massCenter2f;   // center of masses in floats
    int matchSize = 0;
    int accErrorSize = 0;

    double H[9] = { 0 };

    #pragma omp parallel
    {
      CVec3D partialMC1(0, 0, 0), partialMC2(0, 0, 0);
      double partialH[9] = { 0 };
      CVec3* pts1Matched = new CVec3[in_pcl2.m_numPts];
      CVec3* pts2Matched = new CVec3[in_pcl2.m_numPts]; //and transformed.
      int numPts = 0;

      // establish correspondences
      #pragma omp for reduction(+:matchSize, accErrorSize, l_accError)
      for (int i = 0; i<in_pcl2.m_numPts; i++)
      {
        // transform point according to R|t
        MultiplyVectorRightSidePlusOffset(io_Rt, in_pcl2.m_pos[i], pts2Matched[numPts]);


        // search nearest neighbor
        CVec3 normal(0, 0, 1);
        double l_dist;
        if (!MatchPoint(in_pcl1, pts2Matched[numPts], normal, l_scoreDistThreshold, pts1Matched[numPts], l_dist))
          continue;   // no nearest point within radius
        
        l_accError += Dist(pts2Matched[numPts], pts1Matched[numPts]);
        accErrorSize++;

        if (!(l_dist < l_regDistThreshold))
          continue;

        partialMC1 += CVec3D(pts1Matched[numPts].x, pts1Matched[numPts].y, pts1Matched[numPts].z);
        partialMC2 += CVec3D(pts2Matched[numPts].x, pts2Matched[numPts].y, pts2Matched[numPts].z);

        matchSize++;
        numPts++;
      }

      #pragma omp critical
      {
        l_massCenter1 += partialMC1;
        l_massCenter2 += partialMC2;
      }

      #pragma omp barrier 
      // compute center of mass (average) from sums
      #pragma omp master
      {
        //TODO: see if matchSize == 0 -> Zero points were matched with current registration

        out_PreviousFitnessScore = l_accError / accErrorSize;
        l_massCenter1 /= (double)matchSize;
        l_massCenter2 /= (double)matchSize;
        l_massCenter1f = CVec3(float(l_massCenter1.x), float(l_massCenter1.y), float(l_massCenter1.z));
        l_massCenter2f = CVec3(float(l_massCenter2.x), float(l_massCenter2.y), float(l_massCenter2.z));
      }
      #pragma omp barrier

      // subtract center of mass
      //#pragma omp for 
      for (int i = 0; i < numPts; i++)
      {
        pts1Matched[i] -= l_massCenter1f;
        pts2Matched[i] -= l_massCenter2f;
      }

      // compute relative rotation matrix R and translation vector t
      for (int i = 0; i < numPts; i++)
      {
        partialH[0] += double(pts2Matched[i].x * pts1Matched[i].x); partialH[1] += double(pts2Matched[i].x * pts1Matched[i].y); partialH[2] += double(pts2Matched[i].x * pts1Matched[i].z);
        partialH[3] += double(pts2Matched[i].y * pts1Matched[i].x); partialH[4] += double(pts2Matched[i].y * pts1Matched[i].y); partialH[5] += double(pts2Matched[i].y * pts1Matched[i].z);
        partialH[6] += double(pts2Matched[i].z * pts1Matched[i].x); partialH[7] += double(pts2Matched[i].z * pts1Matched[i].y); partialH[8] += double(pts2Matched[i].z * pts1Matched[i].z);
      }

      #pragma omp critical
      {
        H[0] += partialH[0]; H[1] += partialH[1]; H[2] += partialH[2];
        H[3] += partialH[3]; H[4] += partialH[4]; H[5] += partialH[5];
        H[6] += partialH[6]; H[7] += partialH[7]; H[8] += partialH[8];
      }

      delete[] pts1Matched;
      delete[] pts2Matched;
    } // omp

    CMat4 U, W, V;
    svd3x3(H, U, W, V);

    CMat4 TranU; Transpose(U, TranU);
    CMat4 RChange = V * TranU;

    double det = MatrixDeterminant(&RChange);
    // fix improper matrix problem
    if (det<0) {
      CMat4 B; MatrixIdentity(&B);
      B.m[2][2] = float(det);
      RChange = V*B*TranU;
    }

    CVec3 R_mut; MultiplyVectorRightSide(RChange, l_massCenter2f, R_mut);
    CVec3 tChange = l_massCenter1f - R_mut;

    // compose: R|t = R_|t_ * R|t
    io_Rt = RChange * io_Rt;
    CVec3 t(io_Rt.m[3][0], io_Rt.m[3][1], io_Rt.m[3][2]);
    MultiplyVectorRightSide(RChange, t, R_mut);
    t = R_mut + tChange;
    io_Rt.m[3][0] = t.x;    io_Rt.m[3][1] = t.y;    io_Rt.m[3][2] = t.z; io_Rt.m[3][3] = 1;

    // return max delta in parameters
    RChange.m[3][0] = tChange.x;    RChange.m[3][1] = tChange.y;    RChange.m[3][2] = tChange.z;
    double length_mat, length_vec;
    Lengths(RChange, length_mat, length_vec);

    out_transformationChange = length_mat + length_vec;
  }


  double FinalError(CSpatialHash2D& in_pcl1, const CPtCloud& in_pcl2, const CMat4& in_Rt, const double in_scoreDistThreshold)
  {
    double l_accError = 0;
    int accErrorSize = 0;

    // extract matrix and translation vector
    for (int i = 0; i<in_pcl2.m_numPts; i++)
    {
      CVec3 transformedPt;
      CVec3 closestPt;

      // transform point according to R|t
      MultiplyVectorRightSidePlusOffset(in_Rt, in_pcl2.m_pos[i], transformedPt);
      // search for match
      CVec3 normal(0, 0, 1);
      double dist;
      if (!MatchPoint(in_pcl1, transformedPt, normal, in_scoreDistThreshold, closestPt, dist))
        continue;
      l_accError += Dist(transformedPt, closestPt);
      accErrorSize++;
    }
    if (accErrorSize != 0)
      l_accError /= accErrorSize;
    return l_accError;
  }



  /******************************************************************************
  *
  *: Class name: ICP
  *
  ******************************************************************************/
  ICP::ICP()
  {
    initMembers();
  }

  ICP::ICP(float in_regRes)
  {
    initMembers();
    m_regRes = in_regRes;
  }

  ICP::~ICP()
  {
    if (!m_outsourceMainPC)
      delete m_mainHashed;
  }

  void ICP::SetMainPtCloud(const CPtCloud& in_pcl, bool in_append)
  {
    if (m_outsourceMainPC)
    {
      m_mainHashed = new CSpatialHash2D(m_regRes);
      m_outsourceMainPC = false;
    }
    
    if (!in_append)
      m_mainHashed->Clear();

    for (int ptrIndex = 0; ptrIndex < in_pcl.m_numPts; ptrIndex++)
    {
      m_mainHashed->Add(in_pcl.m_pos[ptrIndex], (void*)(1));
    }
  }

  void ICP::SetMainPtCloud(CSpatialHash2D* in_mainHashed)
  {
    if (!m_outsourceMainPC)
    {
      delete m_mainHashed;
    }

    m_mainHashed = in_mainHashed;
    m_outsourceMainPC = true;
  }


  void* ICP::getMainHashedPtr()
  {
    return m_mainHashed;
  }


  void ICP::setRegistrationResolution(float in_regRes)
  {
    m_regRes = in_regRes;
  };


  float ICP::RegisterCloud(const CPtCloud& in_pcl, CMat4& out_registration, CMat4* in_estimatedOrient)
  {
    //TODO: see if in_pcl.m_numPts <5 -> ICP registration called with less than 5 points

    // initial guess of orientation
    if (in_estimatedOrient)
      out_registration = *in_estimatedOrient;
    else
      MatrixIdentity(&out_registration);


    double l_transformationEpsilon = 0.75 * m_regRes;
    double l_fitnessEpsilon = 0.2 * m_regRes;

    double l_transformationChange;
    double l_PreviousFitnessScore;

    PerformIter(*m_mainHashed, in_pcl, out_registration, m_regRes, l_transformationChange, l_PreviousFitnessScore);
    bool converged = (l_PreviousFitnessScore < l_fitnessEpsilon) || (l_transformationChange <= l_transformationEpsilon);

    int l_iterLeft = 150;
    while (!converged)
    {
      PerformIter(*m_mainHashed, in_pcl, out_registration, m_regRes, l_transformationChange, l_PreviousFitnessScore);
      l_iterLeft--;
      converged = (l_PreviousFitnessScore < l_fitnessEpsilon) || (l_transformationChange < l_transformationEpsilon) || (l_iterLeft == 0);
    }

    return float(FinalError(*m_mainHashed, in_pcl, out_registration, 5 * m_regRes));
  }


  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/

  void ICP::initMembers()
  {
    m_regRes = 0.5;
    m_mainHashed = new CSpatialHash2D(m_regRes);
    m_mainHashed->Clear();
    m_outsourceMainPC = false;
  }

} //namespace tpcl