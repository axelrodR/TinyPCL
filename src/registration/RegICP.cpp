/******************************************************************************
*
*: Package Name: sldrcr_ftr
*
******************************************************************************/
#include "RegICP.h"
#include "SpatialHash.h"
#include "common.h"
#include <vector>
#include "../../include/vec.h"


//#include <D3dx9core.h> // uncommenet if using DirectX
////#include <ifr/ifrgen/ifrgen_stnd.h>
//#include "sldrcr_icp.h"
//#include <gen/gengmtrx/gengmtrx_vec.h>
//#include <IFR\ifrlog\ifrlog_prfl.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif



namespace tpcl
{

  /******************************************************************************
  *                             INTERNAL CONSTANTS  / Functions                 *
  ******************************************************************************/

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
  /** Try to find match in point cloud for another point.
  * Currently we match by finding closest point within a given threshold radius ("inliers")
  * @param Xi_pcl1    	    Spatial hashing of main point cloud.
  * @param Xi_p2            a point from 2nd point cloud.
  * @param Xi_normal        !!curently not used!! - normal of the point to which we want to find a match.
  * @param Xo_match         match found.
  * @param indist           inline distance.
  * return                  true if a match was found, flase otherwise*/
  bool MatchPoint(const CSpatialHash2D& Xi_pcl1, const CVec3& Xi_p2, const CVec3& /*Xi_normal*/, const double indist, CVec3& Xo_match)
  {
    // int GetNear(const CVec3& Xi_pos, int xi_bufSize, void** Xo_buf, CVec3* Xo_pos=0, float Xi_max2DRadius=0.0f) const;
    // go over retrieved list (if size = 0, return false), check if closer than indist, update Xo_match if better normal match. return true.


    // search nearest neighbor
    if (Xi_pcl1.FindNearest(Xi_p2, &Xo_match, float(indist)))
    {
      // check if it is an inlier
      double dist = Dist(Xi_p2, Xo_match);
      if (dist < indist)
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
  void BiDiag(double Xio_U[3][3], double Xo_W[3], double Xo_V[3][3], double rv1[4])
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
        s += Xio_U[k][i] * Xio_U[k][i];
      if (s)
      {
        f = Xio_U[i][i];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        Xio_U[i][i] = f - g;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i; k<3; k++)
            s += Xio_U[k][i] * Xio_U[k][j];
          f = s / h;
          for (k = i; k<3; k++)
            Xio_U[k][j] += f*Xio_U[k][i];
        }
      }
      Xo_W[i] = g;

      // act on rows
      g = s = 0.0;
      for (k = i + 1; k < 3; k++)
        s += Xio_U[i][k] * Xio_U[i][k];
      if (s)
      {
        f = Xio_U[i][i + 1];
        g = -SIGN(sqrt(s), f);
        h = f*g - s;
        Xio_U[i][i + 1] = f - g;
        for (k = i + 1; k<3; k++)
          rv1[k] = Xio_U[i][k] / h;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += Xio_U[j][k] * Xio_U[i][k];
          for (k = i + 1; k<3; k++)
            Xio_U[j][k] += s*rv1[k];
        }
      }
    }

    // Accumulation of right-hand transformations.
    Xo_V[2][2] = 1.0;
    g = rv1[2];
    for (i = 2; i >= 0; i--)
    {
      if (g)
      {
        for (j = i + 1; j<3; j++)
          Xo_V[j][i] = (Xio_U[i][j] / Xio_U[i][i + 1]) / g;   // Double division to avoid possible underflow.
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += Xio_U[i][k] * Xo_V[k][j];
          for (k = i + 1; k<3; k++)
            Xo_V[k][j] += s*Xo_V[k][i];
        }
      }
      for (j = i + 1; j<3; j++)
        Xo_V[i][j] = Xo_V[j][i] = 0.0;
      Xo_V[i][i] = 1.0;
      g = rv1[i];
    }

    // Accumulation of left-hand transformations.
    for (i = 2; i >= 0; i--)
    {
      g = Xo_W[i];
      for (j = i + 1; j<3; j++)
        Xio_U[i][j] = 0.0;
      if (g)
      {
        g = 1.0 / g;
        for (j = i + 1; j<3; j++)
        {
          for (s = 0.0, k = i + 1; k<3; k++)
            s += Xio_U[k][i] * Xio_U[k][j];
          f = (s / Xio_U[i][i])*g;
          for (k = i; k<3; k++)
            Xio_U[k][j] += f*Xio_U[k][i];
        }
        for (j = i; j<3; j++)
          Xio_U[j][i] *= g;
      }
      else
        for (j = i; j<3; j++)
          Xio_U[j][i] = 0.0;
      ++Xio_U[i][i];
    }
  }



  // compute the SVD of a 3x3 matrix
  bool svd3x3(const double Xi_M[3][3], double Xo_U[3][3], double Xo_W[3], double Xo_V[3][3])
  {
    double w[3];
    double rv1[3];

    int flag, i, its, j, jj, k, nm;
    double c, f, g, h, s, x, y, z;

    for (i = 0; i < 9; ++i)
      Xo_U[0][i] = Xi_M[0][i];

    // convert to bidiagonal form
    BiDiag(Xo_U, w, Xo_V, rv1);

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
              y = Xo_U[j][nm];
              z = Xo_U[j][i];
              Xo_U[j][nm] = y*c + z*s;
              Xo_U[j][i] = z*c - y*s;
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
              Xo_V[j][k] = -Xo_V[j][k];
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
            x = Xo_V[jj][j];
            z = Xo_V[jj][i];
            Xo_V[jj][j] = x*c + z*s;
            Xo_V[jj][i] = z*c - x*s;
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
            y = Xo_U[jj][j];
            z = Xo_U[jj][i];
            Xo_U[jj][j] = y*c + z*s;
            Xo_U[jj][i] = z*c - y*s;
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
        su[k] = Xo_U[k][i];
      for (k = 0; k<3; k++)
        sv[k] = Xo_V[k][i];
      j = i;
      while (w[j - 1] < sw)
      {
        w[j] = w[j - 1];
        for (k = 0; k<3; k++)
          Xo_U[k][j] = Xo_U[k][j - 1];
        for (k = 0; k<3; k++)
          Xo_V[k][j] = Xo_V[k][j - 1];
        j -= 1;
        if (j < 1) break;
      }
      w[j] = sw;
      for (k = 0; k<3; k++)
        Xo_U[k][j] = su[k];
      for (k = 0; k<3; k++)
        Xo_V[k][j] = sv[k];
    }

    // flip signs
    for (k = 0; k<3; k++)
    {
      s2 = 0;
      for (i = 0; i<3; i++)
        if (Xo_U[i][k] < 0.0)
          s2++;
      for (j = 0; j<3; j++)
        if (Xo_V[j][k] < 0.0)
          s2++;
      if (s2 > 3)
      {
        for (i = 0; i<3; i++)
          Xo_U[i][k] = -Xo_U[i][k];
        for (j = 0; j<3; j++)
          Xo_V[j][k] = -Xo_V[j][k];
      }
    }

    // create vector and copy singular values
    for (int r = 0; r < 3; r++)
      Xo_W[r] = w[r];
    return true;
  }


  bool svd3x3(double* Xi_M, CMat4& Xo_U, CMat4& Xo_W, CMat4& Xo_V)
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
        M[row][col] = Xi_M[row*M1size + col];
      }
    }
    bool convergence = svd3x3(M, U, W, V);

    if (convergence)
    {
      for (int row = 0; row < M1size; row++)
      {
        for (int col = 0; col < M1size; col++)
        {
          Xo_U.m[row][col] = (float)U[row][col];
          Xo_W.m[row][col] = 0.f;
          Xo_V.m[row][col] = (float)V[row][col];
        }
        Xo_W.m[row][row] = (float)W[row];
      }

      for (int i = 0; i < M1size; ++i)
      {
        Xo_U.m[M1size][i] = Xo_U.m[i][M1size] = 0.f;
        Xo_W.m[M1size][i] = Xo_W.m[i][M1size] = 0.f;
        Xo_V.m[M1size][i] = Xo_V.m[i][M1size] = 0.f;
      }
      Xo_U.m[M1size][M1size] = Xo_W.m[3][3] = Xo_V.m[M1size][M1size] = 1.f;
    }

    return convergence;
  }


  typedef TVec3<double> CVec3D;

  double fitStep(CSpatialHash2D& Xi_pcl1, int Xi_pcl2size, CVec3* Xi_pcl2, CMat4& Xio_Rt, const double indist)
  {
    //PROFILE("fitStep");
    // extract matrix and translation vector
    double r00 = Xio_Rt.m[0][0]; double r01 = Xio_Rt.m[0][1]; double r02 = Xio_Rt.m[0][2];
    double r10 = Xio_Rt.m[1][0]; double r11 = Xio_Rt.m[1][1]; double r12 = Xio_Rt.m[1][2];
    double r20 = Xio_Rt.m[2][0]; double r21 = Xio_Rt.m[2][1]; double r22 = Xio_Rt.m[2][2];
    double t0 = Xio_Rt.m[3][0]; double t1 = Xio_Rt.m[3][1]; double t2 = Xio_Rt.m[3][2];
    CVec3D cm1d(0, 0, 0), cm2d(0, 0, 0);  // center of masses for both clouds (as doubles)
    CVec3 cm1, cm2;   // center of masses in floats
    int activeSize = 0;

    double H[9] = { 0 };

#pragma omp parallel
    {
      CVec3D partialCm1d(0, 0, 0), partialCm2d(0, 0, 0);
      double partialH[9] = { 0 };
      CVec3* pts1Matched = new CVec3[Xi_pcl2size];
      CVec3* pts2Matched = new CVec3[Xi_pcl2size]; //and transformed.
      int numPts = 0;

      // establish correspondences
#pragma omp for reduction(+:activeSize)
      for (int i = 0; i<Xi_pcl2size; i++)
      {
        // transform point according to R|t
        CVec3 pos = Xi_pcl2[i];
        pts2Matched[numPts].x = float(r00*pos.x + r01*pos.y + r02*pos.z + t0);
        pts2Matched[numPts].y = float(r10*pos.x + r11*pos.y + r12*pos.z + t1);
        pts2Matched[numPts].z = float(r20*pos.x + r21*pos.y + r22*pos.z + t2);

        // search nearest neighbor
        CVec3 normal(0, 0, 1);
        if (!MatchPoint(Xi_pcl1, pts2Matched[numPts], normal, indist, pts1Matched[numPts]))
          continue;   // no nearest point within radius

        partialCm1d += CVec3D(pts1Matched[numPts].x, pts1Matched[numPts].y, pts1Matched[numPts].z);
        partialCm2d += CVec3D(pts2Matched[numPts].x, pts2Matched[numPts].y, pts2Matched[numPts].z);

        activeSize++;
        numPts++;
      }

#pragma omp critical
      {
        cm1d += partialCm1d;
        cm2d += partialCm2d;
      }

#pragma omp barrier 
      // compute center of mass (average) from sums
#pragma omp master
      {
        cm1d /= (double)activeSize;
        cm2d /= (double)activeSize;
        cm1 = CVec3(float(cm1d.x), float(cm1d.y), float(cm1d.z));
        cm2 = CVec3(float(cm2d.x), float(cm2d.y), float(cm2d.z));
      }
#pragma omp barrier

      // subtract center of mass
      //#pragma omp for 
      for (int i = 0; i < numPts; i++)
      {
        pts1Matched[i] -= cm1;
        pts2Matched[i] -= cm2;
      }

      // compute relative rotation matrix R and translation vector t
      //H = q_t^T*q_m = sum(~q_t[i] 3x1 * q_m[i] 1x3) = 3x3.
      //#pragma omp /*for*/ reduction(+: H[0],H[1],H[2],H[3],H[4],H[5],H[6],H[7],H[8],H[9],H[10],H[11])
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
    CMat4 R_ = V * TranU;

    double det = MatrixDeterminant(&R_);
    // fix improper matrix problem
    if (det<0) {
      CMat4 B; MatrixIdentity(&B);
      B.m[2][2] = float(det);
      R_ = V*B*TranU;
    }

    CVec3 R_mut; MultiplyVectorRightSide(R_, cm2, R_mut);
    CVec3 t_ = cm1 - R_mut;

    // compose: R|t = R_|t_ * R|t
    Xio_Rt = R_ * Xio_Rt;
    CVec3 t(Xio_Rt.m[3][0], Xio_Rt.m[3][1], Xio_Rt.m[3][2]);
    MultiplyVectorRightSide(R_, t, R_mut);
    t = R_mut + t_;
    Xio_Rt.m[3][0] = t.x;    Xio_Rt.m[3][1] = t.y;    Xio_Rt.m[3][2] = t.z; Xio_Rt.m[3][3] = 1;

    // return max delta in parameters
    for (int i = 0; i < 3; i++)
      R_.m[i][i] -= 1;
    R_.m[3][0] = t_.x;    R_.m[3][1] = t_.y;    R_.m[3][2] = t_.z;
    double length_mat, length_vec;
    Lengths(R_, length_mat, length_vec);

    return MaxT(length_mat, length_vec);
  }




  double getResidual(CSpatialHash2D& Xi_pcl1, int Xi_pcl2size, CVec3* Xi_pcl2, const CMat4& Xi_Rt, const double Xi_indist)
  {
    double residual = 0;
    int inliers = 0;

    // extract matrix and translation vector
    double r00 = Xi_Rt.m[0][0]; double r01 = Xi_Rt.m[0][1]; double r02 = Xi_Rt.m[0][2];
    double r10 = Xi_Rt.m[1][0]; double r11 = Xi_Rt.m[1][1]; double r12 = Xi_Rt.m[1][2];
    double r20 = Xi_Rt.m[2][0]; double r21 = Xi_Rt.m[2][1]; double r22 = Xi_Rt.m[2][2];
    double t0 = Xi_Rt.m[3][0]; double t1 = Xi_Rt.m[3][1]; double t2 = Xi_Rt.m[3][2];
    for (int i = 0; i<Xi_pcl2size; i++)
    {
      CVec3 query;
      CVec3 result;

      // transform point according to R|t
      CVec3 Pos = Xi_pcl2[i];
      query.x = (float)(r00*Pos.x + r01*Pos.y + r02*Pos.z + t0);
      query.y = (float)(r10*Pos.x + r11*Pos.y + r12*Pos.z + t1);
      query.z = (float)(r20*Pos.x + r21*Pos.y + r22*Pos.z + t2);
      // search for match
      CVec3 normal(0, 0, 1);
      if (!MatchPoint(Xi_pcl1, query, normal, Xi_indist, result))
        continue;
      residual += Dist(query, result);
      inliers++;
    }
    if (inliers != 0)
      residual /= inliers;
    return residual;
  }


  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/

  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: SLDR_ICP_ICP
  *
  ******************************************************************************/
  ICP::ICP()
  {
    initMembers();
  }

  ICP::ICP(double Xi_regThresh)
  {
    initMembers();
    m_minDelta = Xi_regThresh;
  }

  /******************************************************************************
  *
  *: Method name: ~SLDR_ICP_ICP
  *
  ******************************************************************************/
  ICP::~ICP()
  {
    if (!m_outsourceMainPC)
      delete m_mainHashed;
  }

  /******************************************************************************
  *
  *: Method name: MainPointCloudUpdate
  *
  ******************************************************************************/
  void ICP::MainPointCloudUpdate(int Xi_numPts, const CVec3* Xi_pts, bool Xi_clean)
  {
    if (m_outsourceMainPC)
    {
      m_mainHashed = new CSpatialHash2D(m_voxelSize);
      m_outsourceMainPC = false;
    }

    if (Xi_clean)
      m_mainHashed->Clear();

    for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
    {
      m_mainHashed->Add(Xi_pts[ptrIndex], (void*)(1));
    }
  }

  void ICP::MainPointCloudUpdate(void* Xi_mainHashed)
  {
    if (!m_outsourceMainPC)
    {
      delete m_mainHashed;
    }

    m_mainHashed = (CSpatialHash2D*)Xi_mainHashed;
    m_outsourceMainPC = true;
  }


  /******************************************************************************
  *
  *: Method name: getMainHashedPtr
  *
  ******************************************************************************/
  void* ICP::getMainHashedPtr()
  {
    return m_mainHashed;
  }


  /******************************************************************************
  *
  *: Method name: setRegistrationThresh
  *
  ******************************************************************************/
  void ICP::setRegistrationThresh(double Xi_regThresh)
  {
    m_minDelta = Xi_regThresh;
  };



  /******************************************************************************
  *
  *: Method name: FillPointCloud
  *
  ******************************************************************************/
  //SoudRegistration(CMat4& Xo_registration, CVec3* Xi_pts, int Xi_numPts, int Xi_lineWidth = -1, CVec3* Xi_estimatedOrient = NULL);
  float ICP::SecondaryPointCloudRegistration(CMat4& Xo_registration, CVec3* Xi_pts, int Xi_numPts, int /*Xi_lineWidth*/, CMat4* Xi_estimatedOrient)
  {
    if (Xi_numPts<5) {
      //ERROR: Icp works only with at least 5 template points
      return 0;
    }

    // initial guess of orientation
    if (Xi_estimatedOrient)
      Xo_registration = *Xi_estimatedOrient;
    else
      MatrixIdentity(&Xo_registration);

    std::vector<int>	l_active;
    l_active.reserve(Xi_numPts);
    //l_active.resize(Xi_numPts);
    //for (int i = 0; i<Xi_numPts; i++)
    //  l_active[i] = i;

    //double m_inlier_ratio = 1;
    double l_inlierDist = m_voxelSize;

    double delta;
    const int max_iter = 200;

    for (int iter = 0; iter < max_iter; ++iter)
    {
      l_inlierDist = MaxT(l_inlierDist*0.9f, 0.05);

      delta = fitStep(*m_mainHashed, Xi_numPts, Xi_pts, Xo_registration, l_inlierDist);
      if (delta <= m_minDelta)
        break;
    }

    return (float)getResidual(*m_mainHashed, Xi_numPts, Xi_pts, Xo_registration, m_voxelSize);
  }


  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/

  void ICP::initMembers()
  {
    m_voxelSize = 0.5;
    m_mainHashed = new CSpatialHash2D(m_voxelSize);
    m_mainHashed->Clear();
    m_outsourceMainPC = false;
    m_minDelta = 1e-4f;
  }


  /******************************************************************************
  *                              Private methods                                *
  ******************************************************************************/


  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/

  /******************************************************************************
  *                            INTERNAL FUNCTIONS                               *
  ******************************************************************************/
  //////////////////////////////////////////////////////////////////

} //namespace tpcl