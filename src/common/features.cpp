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
#include "features.h"
#include <algorithm>
#include "plane.h"
#include "../../include/vec.h"



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
  *: Method name: FindNormal
  *
  ******************************************************************************/
  void Features::FindNormal(CPtCloud& Xio_pcl, float Xi_radius, CSpatialHash2D* Xi_globalHashed, bool Xi_fixZ)
  {
    const int bufSize = 100;
    bool gotHashed = Xi_globalHashed != NULL;

    //if didn't get global hashed point cloud - create one from the point cloud for whoch the normals are found:
    if (!gotHashed)
    {
      Xi_globalHashed = new CSpatialHash2D;
      for (int Index = 0; Index < Xio_pcl.m_numPts; Index++)
      {
        Xi_globalHashed->Add(Xio_pcl.m_pos[Index], NULL);
      }
    }

    #pragma omp parallel
    {
    void* unsuedBuf[bufSize];
      CVec3 closePts[bufSize];
    int numOfClose = 0;

      #pragma omp for
      for (int ptIndex = 0; ptIndex < Xio_pcl.m_numPts; ptIndex++)
      {
        float radius = Xi_radius * 2;
        //get close points:

        numOfClose = Xi_globalHashed->GetNear(Xio_pcl.m_pos[ptIndex], bufSize, unsuedBuf, closePts, radius);

        //find plane:
        CPlane approxPlane;
        if (approxPlane.RanSaC(numOfClose, closePts, 1.0f) == 0)
          Xio_pcl.m_normal[ptIndex] = CVec3(0, 0, 1); //TODO: desice what to do if not enoguh points.
        else
        {

          //update point's z:
          if (Xi_fixZ)
            Xio_pcl.m_pos[ptIndex].z = approxPlane.GetHeightAt(Xio_pcl.m_pos[ptIndex].x, Xio_pcl.m_pos[ptIndex].y);

          //get normal and normalize it, up:
          Xio_pcl.m_normal[ptIndex] = approxPlane.GetNormal();
          if (Xio_pcl.m_normal[ptIndex].z < 0)
            Xio_pcl.m_normal[ptIndex].z = -Xio_pcl.m_normal[ptIndex].z;
          Normalize(Xio_pcl.m_normal[ptIndex]);
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
  void Features::DenoiseRangeOfOrderedPointCloud(const CPtCloud& Xi_pcl, CPtCloud& Xo_pcl, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh)
  {
    int numlines = Xi_pcl.m_numPts / Xi_pcl.m_lineWidth;

    //convert to range image (each pixle in ptsLocal has it's xyz).
    float* rangeImage = new float[Xi_pcl.m_numPts];
    #pragma omp parallel for
    for (int ptrIndex = 0; ptrIndex < Xi_pcl.m_numPts; ptrIndex++)
    {
      rangeImage[ptrIndex] = Length(Xi_pcl.m_pos[ptrIndex]);
    }

    //find diferent points:
    float* distFiltered = new float[Xi_pcl.m_numPts];
    float* threshFiltered = new float[Xi_pcl.m_numPts];

    Median2D(Xi_pcl.m_lineWidth, numlines, Xi_medFiltSize0, rangeImage, distFiltered);

#pragma omp parallel for
    for (int row = 0; row < numlines; row++)
    {
      for (int col = 0; col < Xi_pcl.m_lineWidth; col++)
      {
        int index = (row * Xi_pcl.m_lineWidth) + col;
        distFiltered[index] = (abs(distFiltered[index] - rangeImage[index]) < Xi_distFromMedianThresh) ? 1.0f : 0.0f;
      }
    }

    Median2D(Xi_pcl.m_lineWidth, numlines, Xi_medFiltSize1, distFiltered, threshFiltered);

    int outputSize = 0;
    for (int row = 0; row < numlines; row++)
    {
      for (int col = 0; col < Xi_pcl.m_lineWidth; col++)
      {
        int index = (row * Xi_pcl.m_lineWidth) + col;

        if ((distFiltered[index] == 1) && (threshFiltered[index] == 1))
        {
          Xo_pcl.m_pos[outputSize] = Xi_pcl.m_pos[index];
          if ((Xi_pcl.m_color != 0) && (Xo_pcl.m_color != 0))
            Xo_pcl.m_color[outputSize] = Xi_pcl.m_color[index];
          if ((Xi_pcl.m_normal != 0) && (Xo_pcl.m_normal != 0))
            Xo_pcl.m_normal[outputSize] = Xi_pcl.m_normal[index];
          outputSize++;
        }
      }
    }

    Xo_pcl.m_numPts = outputSize;

    delete[] rangeImage;
    delete[] distFiltered;
    delete[] threshFiltered;
  }




  /******************************************************************************
  *
  *: Method name: DenoiseRangeOfPointCloud
  *
  ******************************************************************************/
  void Features::DenoiseRangeOfPointCloud()
  {

  }



  /******************************************************************************
  *
  *: Method name: DownSamplePointCloud
  *
  ******************************************************************************/
  void Features::DownSamplePointCloud(const CPtCloud& Xi_pcl, CPtCloud& Xo_pcl, float Xi_voxelSize)
  {
    float InvVoxelSize = 1.0f / Xi_voxelSize;

    //find min/max x/y/z:
    CVec3 l_bbox[2] = { Xi_pcl.m_pos[0], Xi_pcl.m_pos[0] };
    for (int ptrIndex = 1; ptrIndex < Xi_pcl.m_numPts; ptrIndex++)
    {
      l_bbox[0] = Min_ps(l_bbox[0], Xi_pcl.m_pos[ptrIndex]);
      l_bbox[1] = Max_ps(l_bbox[1], Xi_pcl.m_pos[ptrIndex]);
    }

    //create downsampling vector:
    int Mx = int(ceil((l_bbox[1].x - l_bbox[0].x) * InvVoxelSize));
    int My = int(ceil((l_bbox[1].y - l_bbox[0].y) * InvVoxelSize));
    int Mz = int(ceil((l_bbox[1].z - l_bbox[0].z) * InvVoxelSize));
    int Mxy = Mx*My;
    unsigned int Mxyz = Mxy * Mz;

    bool* VisitedVoxel = new bool[Mxyz];
    memset(VisitedVoxel, false, Mxyz * sizeof(bool));

    int outputSize = 0;
    for (int ptrIndex = 0; ptrIndex < Xi_pcl.m_numPts; ptrIndex++)
    {
      float x = Xi_pcl.m_pos[ptrIndex].x;
      float y = Xi_pcl.m_pos[ptrIndex].y;
      float z = Xi_pcl.m_pos[ptrIndex].z;

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
        Xo_pcl.m_pos[outputSize] = Xi_pcl.m_pos[ptrIndex];
        if ((Xi_pcl.m_color != 0) && (Xo_pcl.m_color != 0))
          Xo_pcl.m_color[outputSize] = Xi_pcl.m_color[ptrIndex];
        if ((Xi_pcl.m_normal != 0) && (Xo_pcl.m_normal != 0))
          Xo_pcl.m_normal[outputSize] = Xi_pcl.m_normal[ptrIndex];
        outputSize++;

        VisitedVoxel[index] = true;
      }
    }
    delete[] VisitedVoxel;

    Xo_pcl.m_numPts = outputSize;
  }




  /******************************************************************************
  *
  *: Method name: RMSEofRegistration
  *
  ******************************************************************************/
  float Features::RMSEofRegistration(CSpatialHash2D* Xi_pcl1, const CPtCloud& Xi_pcl2, float Xi_max2DRadius, const CMat4& Xi_Rt)
  {
    double RMSE = 0;
    float outPenalty = Xi_max2DRadius*Xi_max2DRadius*1.5f;

    #pragma omp parallel for reduction(+:RMSE)
    for (int i = 0; i<Xi_pcl2.m_numPts; i++)
    {
      CVec3 transformedPt;
      CVec3 closestPt;

      // transform point according to R|t
      MultiplyVectorRightSidePlusOffset(Xi_Rt, Xi_pcl2.m_pos[i], transformedPt);

      // search for match
      float dist = outPenalty;
      if (Xi_pcl1->FindNearest(transformedPt, &closestPt, Xi_max2DRadius))
      {
        dist = DistSqr(transformedPt, closestPt);
      }

      RMSE += dist;
    }
    
    RMSE /= Xi_pcl2.m_numPts;
    return float(sqrt(RMSE));
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