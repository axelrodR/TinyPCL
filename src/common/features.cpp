// File Location: S:\gen\gengmtrx\gengmtrx_grid.h

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
*: Package Name: features
*
******************************************************************************/
//#include <D3dx9core.h> // uncommenet if using DirectX
//#include <ifr/ifrgen/ifrgen_stnd.h>
#include "features.h"
////#include <gen/gengmtrx/gengmtrx_spat.h>
#include "plane.h"
//#include <vector>



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

  /** 2D median filter. assums image width/height are power of 2.
  *   edges are treated with cyclic indices.
  * @param Xi_lineWidth                image width.
  * @param Xi_numlines                 image height.
  * @param Xi_medFiltSize              median filter size.
  * @param Xi_pts                      input image.
  * @param Xo_ptsFiltered              2D median filtered image. assumes  Xo_ptsFiltered != Xi_pts*/
  void Median2DPowOf2(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize, float* Xi_pts, float* Xo_ptsFiltered)
  {
    int windowSize = Xi_medFiltSize*Xi_medFiltSize;
    int half1DWindow = Xi_medFiltSize >> 1;

    #pragma omp parallel
    {
      float* window = new float[windowSize];
      #pragma omp for
      for (int row = 0; row < Xi_numlines; row++)
      {
        for (int col = 0; col < Xi_lineWidth; col++)
        {
          for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
          {
            for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
            {
              int cyclicRow = ((unsigned int)(row + winRow)) & (Xi_numlines - 1);
              int cyclicCol = ((unsigned int)(col + winCol)) & (Xi_lineWidth - 1);
              int index = (cyclicRow * Xi_lineWidth) + cyclicCol;
              int winIndex = ((winRow + half1DWindow) * Xi_medFiltSize) + (winCol + half1DWindow);
              window[winIndex] = Xi_pts[index];
            }
          }
          std::nth_element(window, window + (windowSize >> 1), window + windowSize);
          Xo_ptsFiltered[row*Xi_lineWidth + col] = window[windowSize >> 1];
        }
      }
      delete[] window;
    }
  }


  /** 2D median filter.
  *   edges are treated with mirroring.
  * @param Xi_lineWidth                image width.
  * @param Xi_numlines                 image height.
  * @param Xi_medFiltSize              median filter size.
  * @param Xi_pts                      input image.
  * @param Xo_ptsFiltered              2D median filtered image. assumes  Xo_ptsFiltered != Xi_pts*/
  void Median2D(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize, float* Xi_pts, float* Xo_ptsFiltered)
  {
    if (Xi_medFiltSize < 2)
    {
      for (int index = 0; index < Xi_lineWidth * Xi_numlines; index++)
        Xo_ptsFiltered[index] = Xi_pts[index];
    }
    else
    {
      int windowSize = Xi_medFiltSize*Xi_medFiltSize;
      int half1DWindow = Xi_medFiltSize >> 1;

      #pragma omp parallel
      {
        //filter on center of image (without edges):
        float* window = new float[windowSize];
        #pragma omp for
        for (int row = half1DWindow; row < Xi_numlines - half1DWindow; row++)
        {
          for (int col = half1DWindow; col < Xi_lineWidth - half1DWindow; col++)
          {
            for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
            {
              for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
              {
                int index = ((row + winRow) * Xi_lineWidth) + (col + winCol);
                int winIndex = ((winRow + half1DWindow) * Xi_medFiltSize) + (winCol + half1DWindow);
                window[winIndex] = Xi_pts[index];
              }
            }
            std::nth_element(window, window + (windowSize >> 1), window + windowSize);
            Xo_ptsFiltered[row*Xi_lineWidth + col] = window[windowSize >> 1];
          }
        }
        delete[] window;
      }


      //filter on edges of image:
      float* window = new float[windowSize];
      for (int row = 0; row < Xi_numlines; row++)
      {
        bool center = (row >= half1DWindow) && (row < Xi_numlines - half1DWindow);

        for (int col = 0; col < Xi_lineWidth; col++)
        {
          if (center && (col >= half1DWindow)) //not working on the center again
          {
            col = Xi_lineWidth - half1DWindow;
            center = false;
          }

          for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
          {
            for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
            {
              int indexRow = row + winRow;
              indexRow = (indexRow < 0) ? -indexRow :
                ((indexRow >= Xi_numlines) ? (2 * (Xi_numlines - 1) - indexRow) :
                  indexRow);

              int indexCol = col + winCol;
              indexCol = (indexCol < 0) ? -indexCol :
                ((indexCol >= Xi_lineWidth) ? (2 * (Xi_lineWidth - 1) - indexCol) :
                  indexCol);

              int index = (indexRow * Xi_lineWidth) + indexCol;
              int winIndex = ((winRow + half1DWindow) * Xi_medFiltSize) + (winCol + half1DWindow);
              window[winIndex] = Xi_pts[index];
            }
          }
          std::nth_element(window, window + (windowSize >> 1), window + windowSize);
          Xo_ptsFiltered[row*Xi_lineWidth + col] = window[windowSize >> 1];
        }
      }
      delete[] window;
    }
  }
 


  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/

  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: Features
  *
  ******************************************************************************/
  Features::Features()
  {
  }
 
  /******************************************************************************
  *
  *: Method name: ~Features
  *
  ******************************************************************************/
  Features::~Features()
  {
  }


  /******************************************************************************
  *
  *: Method name: FillPointCloud
  *
  ******************************************************************************/
  void Features::FillPointCloud(int Xi_numPts, const CVec3* Xi_pts, int Xo_numPts, const CVec3* Xo_pts)
  {
  }


  /******************************************************************************
  *
  *: Method name: FindNormal
  *
  ******************************************************************************/
  void Features::FindNormal(int Xi_numPts, CVec3* Xio_pts, CVec3* Xo_Normals, float Xi_radius, GenGmtrx::CSpatialHash2D* Xi_globalHashed, bool Xi_fixZ)
  {
    const int bufSize = 100;
    bool gotHashed = Xi_globalHashed != NULL;

    //if didn't get global hashed point cloud - create one from the point cloud for whoch the normals are found:
    if (!gotHashed)
    {
      Xi_globalHashed = new GenGmtrx::CSpatialHash2D;
      for (int Index = 0; Index < Xi_numPts; Index++)
      {
        Xi_globalHashed->Add(Xio_pts[Index], NULL);
      }
    }

    #pragma omp parallel
    {
    void* unsuedBuf[bufSize];
      CVec3 closePts[bufSize];
    int numOfClose = 0;

      #pragma omp for
    for (int ptIndex = 0; ptIndex < Xi_numPts; ptIndex++)
    {
        float radius = Xi_radius * 2;
      //get close points:

        numOfClose = Xi_globalHashed->GetNear(Xio_pts[ptIndex], bufSize, unsuedBuf, closePts, radius);

      //find plane:
        GenGmtrx::CPlane approxPlane;
        if (approxPlane.RanSaC(numOfClose, closePts, 1.0f) == 0)
          Xo_Normals[ptIndex] = CVec3(0, 0, 1); //TODO: desice what to do if not enoguh points.
        else
        {

      //update point's z:
          if (Xi_fixZ)
      Xio_pts[ptIndex].z = approxPlane.GetHeightAt(Xio_pts[ptIndex].x, Xio_pts[ptIndex].y);

      //get normal and normalize it, up:
      Xo_Normals[ptIndex] = approxPlane.GetNormal();
      if (Xo_Normals[ptIndex].z < 0)
        Xo_Normals[ptIndex].z = -Xo_Normals[ptIndex].z;
      Normalize(Xo_Normals[ptIndex]);
    }
  }
    }

    if (!gotHashed)
    {
      delete Xi_globalHashed;
    }
  }



  /******************************************************************************
  *
  *: Method name: DenoiseRangeOfOrderedPointCloud
  *
  ******************************************************************************/
  int Features::DenoiseRangeOfOrderedPointCloud(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, CVec3* Xi_pts, CVec3* Xo_ptsDenoised)
  {
    PROFILE("DenoiseRangeOfOrderedPointCloud");
    int totalSize = Xi_lineWidth * Xi_numlines;

    //convert to range image (each pixle in ptsLocal has it's xyz).
    float* rangeImage = new float[totalSize];
#pragma omp parallel for
    for (int ptrIndex = 0; ptrIndex < totalSize; ptrIndex++)
    {
      rangeImage[ptrIndex] = Length(Xi_pts[ptrIndex]);
    }

    //find diferent points:
    float* distFiltered = new float[totalSize];
    float* threshFiltered = new float[totalSize];

    Median2D(Xi_lineWidth, Xi_numlines, Xi_medFiltSize0, rangeImage, distFiltered);

#pragma omp parallel for
    for (int row = 0; row < Xi_numlines; row++)
    {
      for (int col = 0; col < Xi_lineWidth; col++)
      {
        int index = (row * Xi_lineWidth) + col;
        distFiltered[index] = (abs(distFiltered[index] - rangeImage[index]) < Xi_distFromMedianThresh) ? 1.0f : 0.0f;
      }
    }

    Median2D(Xi_lineWidth, Xi_numlines, Xi_medFiltSize1, distFiltered, threshFiltered);

    int outputSize = 0;
    for (int row = 0; row < Xi_numlines; row++)
    {
      for (int col = 0; col < Xi_lineWidth; col++)
      {
        int index = (row * Xi_lineWidth) + col;

        if ((distFiltered[index] == 1) && (threshFiltered[index] == 1))
        {
          Xo_ptsDenoised[outputSize] = Xi_pts[index];
          outputSize++;
        }
      }
    }

    delete[] rangeImage;
    delete[] distFiltered;
    delete[] threshFiltered;

    return outputSize;
  }




  /******************************************************************************
  *
  *: Method name: DenoiseRangeOfPointCloud
  *
  ******************************************************************************/
  int Features::DenoiseRangeOfPointCloud(float Xi_res, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, int Xi_numPts, CVec3* Xi_pts, CVec3* Xo_ptsDenoised)
  {
    PROFILE("DenoiseRangeOfPointCloud");
    //temp: input = outpout
    for (int index = 0; index < Xi_numPts; index++)
      Xo_ptsDenoised[index] = Xi_pts[index];
    return 0;
  }



  /******************************************************************************
  *
  *: Method name: DownSamplePointCloud
  *
  ******************************************************************************/
  void Features::DownSamplePointCloud(float Xi_voxelSize, int Xi_numPts, const CVec3* Xi_pts, int& Xo_numPts, CVec3* Xo_pts)
  {
    PROFILE("DownSamplePointCloud");
    float InvVoxelSize = 1.0f / Xi_voxelSize;

    //find min/max x/y/z:
    CVec3 l_bbox[2] = { Xi_pts[0], Xi_pts[0] };
    for (int ptrIndex = 1; ptrIndex < Xi_numPts; ptrIndex++)
    {
      l_bbox[0] = Min_ps(l_bbox[0], Xi_pts[ptrIndex]);
      l_bbox[1] = Max_ps(l_bbox[1], Xi_pts[ptrIndex]);
    }

    //create downsampling vector:
    int Mx = int(ceil((l_bbox[1].x - l_bbox[0].x) * InvVoxelSize));
    int My = int(ceil((l_bbox[1].y - l_bbox[0].y) * InvVoxelSize));
    int Mz = int(ceil((l_bbox[1].z - l_bbox[0].z) * InvVoxelSize));
    int Mxy = Mx*My;
    unsigned int Mxyz = Mxy * Mz;

    bool* VisitedVoxel = new bool[Mxyz];
    memset(VisitedVoxel, false, Mxyz * sizeof(bool));
    std::vector<CVec3> DownSampledPts;

    for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
    {
      float x = Xi_pts[ptrIndex].x;
      float y = Xi_pts[ptrIndex].y;
      float z = Xi_pts[ptrIndex].z;

      //find point's voxel index:
      int xInd = int(floor((x - l_bbox[0].x) * InvVoxelSize));
      int yInd = int(floor((y - l_bbox[0].y) * InvVoxelSize));
      int zInd = int(floor((z - l_bbox[0].z) * InvVoxelSize));

      if (xInd == Mx) xInd--;
      if (yInd == My) yInd--;
      if (zInd == Mz) zInd--;

      unsigned int index = zInd*Mxy + yInd*Mx + xInd;

      //if we haven't filled this voxel with a point yet, add current point:
      if (!VisitedVoxel[index])
      {
        DownSampledPts.push_back(Xi_pts[ptrIndex]);
        VisitedVoxel[index] = true;
      }
    }
    delete[] VisitedVoxel;


    //save selected points:
    for (int VecIndex = 0; VecIndex < DownSampledPts.size(); VecIndex++)
    {
      Xo_pts[VecIndex] = DownSampledPts[VecIndex];
    }
    Xo_numPts = int(DownSampledPts.size());
  }




  /******************************************************************************
  *
  *: Method name: RMSEofRegistration
  *
  ******************************************************************************/
  float Features::RMSEofRegistration(float hashRes, GenGmtrx::CSpatialHash2D* Xi_pcl1, int Xi_pcl2size, CVec3* Xi_pcl2, CMat4& Xi_Rt)
  {
    PROFILE("RMSEofRegistration");
    double RMSE = 0;
    float max2DRadius = hashRes * 2;
    float outPenalty = max2DRadius*max2DRadius*1.5;

    // extract matrix and translation vector
    double r00 = Xi_Rt.m[0][0]; double r01 = Xi_Rt.m[0][1]; double r02 = Xi_Rt.m[0][2];
    double r10 = Xi_Rt.m[1][0]; double r11 = Xi_Rt.m[1][1]; double r12 = Xi_Rt.m[1][2];
    double r20 = Xi_Rt.m[2][0]; double r21 = Xi_Rt.m[2][1]; double r22 = Xi_Rt.m[2][2];
    double t0  = Xi_Rt.m[3][0]; double t1  = Xi_Rt.m[3][1]; double t2  = Xi_Rt.m[3][2];
    #pragma omp parallel for reduction(+:RMSE)
    for (int i = 0; i<Xi_pcl2size; i++)
    {
      CVec3 query;
      CVec3 result;

      // transform point according to R|t
      CVec3 Pos = Xi_pcl2[i];
      query.x = r00*Pos.x + r01*Pos.y + r02*Pos.z + t0;
      query.y = r10*Pos.x + r11*Pos.y + r12*Pos.z + t1;
      query.z = r20*Pos.x + r21*Pos.y + r22*Pos.z + t2;
      // search for match
      float dist = outPenalty;
      if (Xi_pcl1->FindNearest(query, &result, max2DRadius))
      {
        dist = DistSqr(query, result);
      }

      RMSE += dist;
    }
    
    RMSE /= Xi_pcl2size;
    return sqrt(RMSE);
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


} //namespace SLDR