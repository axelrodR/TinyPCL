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

#include "features.h"
#include <algorithm>
#include "SpatialHash.h"
#include "plane.h"
#include "../include/ptCloud.h"


namespace tpcl
{
  /** 2D median filter. assums image width/height are power of 2.
  *   edges are treated with cyclic indices.
  * @param in_lineWidth                image width.
  * @param in_numlines                 image height.
  * @param in_medFiltSize              median filter size.
  * @param in_pts                      input image.
  * @param out_ptsFiltered              2D median filtered image. assumes  out_ptsFiltered != in_pts*/
  void Median2DPowOf2(int in_lineWidth, int in_numlines, int in_medFiltSize, float* in_pts, float* out_ptsFiltered)
  {
    int windowSize = in_medFiltSize*in_medFiltSize;
    int half1DWindow = in_medFiltSize >> 1;

    #pragma omp parallel
    {
      float* window = new float[windowSize];
      #pragma omp for
      for (int row = 0; row < in_numlines; row++)
      {
        for (int col = 0; col < in_lineWidth; col++)
        {
          for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
          {
            for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
            {
              int cyclicRow = ((unsigned int)(row + winRow)) & (in_numlines - 1);
              int cyclicCol = ((unsigned int)(col + winCol)) & (in_lineWidth - 1);
              int index = (cyclicRow * in_lineWidth) + cyclicCol;
              int winIndex = ((winRow + half1DWindow) * in_medFiltSize) + (winCol + half1DWindow);
              window[winIndex] = in_pts[index];
            }
          }
          std::nth_element(window, window + (windowSize >> 1), window + windowSize);
          out_ptsFiltered[row*in_lineWidth + col] = window[windowSize >> 1];
        }
      }
      delete[] window;
    }
  }


  /** 2D median filter.
  *   edges are treated with mirroring.
  * @param in_lineWidth                image width.
  * @param in_numlines                 image height.
  * @param in_medFiltSize              median filter size.
  * @param in_pts                      input image.
  * @param out_ptsFiltered              2D median filtered image. assumes  out_ptsFiltered != in_pts*/
  void Median2D(int in_lineWidth, int in_numlines, int in_medFiltSize, float* in_pts, float* out_ptsFiltered)
  {
    if (in_medFiltSize < 2)
    {
      for (int index = 0; index < in_lineWidth * in_numlines; index++)
        out_ptsFiltered[index] = in_pts[index];
    }
    else
    {
      int windowSize = in_medFiltSize*in_medFiltSize;
      int half1DWindow = in_medFiltSize >> 1;

      #pragma omp parallel
      {
        //filter on center of image (without edges):
        float* window = new float[windowSize];
        #pragma omp for
        for (int row = half1DWindow; row < in_numlines - half1DWindow; row++)
        {
          for (int col = half1DWindow; col < in_lineWidth - half1DWindow; col++)
          {
            for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
            {
              for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
              {
                int index = ((row + winRow) * in_lineWidth) + (col + winCol);
                int winIndex = ((winRow + half1DWindow) * in_medFiltSize) + (winCol + half1DWindow);
                window[winIndex] = in_pts[index];
              }
            }
            std::nth_element(window, window + (windowSize >> 1), window + windowSize);
            out_ptsFiltered[row*in_lineWidth + col] = window[windowSize >> 1];
          }
        }
        delete[] window;
      }


      //filter on edges of image:
      float* window = new float[windowSize];
      for (int row = 0; row < in_numlines; row++)
      {
        bool center = (row >= half1DWindow) && (row < in_numlines - half1DWindow);

        for (int col = 0; col < in_lineWidth; col++)
        {
          if (center && (col >= half1DWindow)) //not working on the center again
          {
            col = in_lineWidth - half1DWindow;
            center = false;
          }

          for (int winRow = -half1DWindow; winRow <= half1DWindow; winRow++)
          {
            for (int winCol = -half1DWindow; winCol <= half1DWindow; winCol++)
            {
              int indexRow = row + winRow;
              indexRow = (indexRow < 0) ? -indexRow :
                ((indexRow >= in_numlines) ? (2 * (in_numlines - 1) - indexRow) :
                  indexRow);

              int indexCol = col + winCol;
              indexCol = (indexCol < 0) ? -indexCol :
                ((indexCol >= in_lineWidth) ? (2 * (in_lineWidth - 1) - indexCol) :
                  indexCol);

              int index = (indexRow * in_lineWidth) + indexCol;
              int winIndex = ((winRow + half1DWindow) * in_medFiltSize) + (winCol + half1DWindow);
              window[winIndex] = in_pts[index];
            }
          }
          std::nth_element(window, window + (windowSize >> 1), window + windowSize);
          out_ptsFiltered[row*in_lineWidth + col] = window[windowSize >> 1];
        }
      }
      delete[] window;
    }
  }



  /******************************************************************************
  *
  *: Class name: Features
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
  void Features::FillNormals(CPtCloud& io_pcl, float in_radius, CSpatialHash2D* in_pclHash, bool in_fixZ)
  {
    const int bufSize = 100;
    bool gotHashed = in_pclHash != NULL;

    //if didn't get global hashed point cloud - create one from the point cloud for whoch the normals are found:
    if (!gotHashed)
    {
      in_pclHash = new CSpatialHash2D;
      for (int Index = 0; Index < io_pcl.m_numPts; Index++)
      {
        in_pclHash->Add(io_pcl.m_pos[Index], NULL);
      }
    }

    #pragma omp parallel
    {
    void* unsuedBuf[bufSize];
      CVec3 closePts[bufSize];
    int numOfClose = 0;

      #pragma omp for
      for (int ptIndex = 0; ptIndex < io_pcl.m_numPts; ptIndex++)
      {
        float radius = in_radius * 2;
        //get close points:

        numOfClose = in_pclHash->GetNear(io_pcl.m_pos[ptIndex], bufSize, unsuedBuf, closePts, radius);

        //find plane:
        CPlane approxPlane;
        if (approxPlane.RanSaC(numOfClose, closePts, 0.1f) == 0)
          io_pcl.m_normal[ptIndex] = CVec3(0, 0, 1); //TODO: desice what to do if not enoguh points.
        else
        {

          //update point's z:
          if (in_fixZ)
            io_pcl.m_pos[ptIndex].z = approxPlane.GetHeightAt(io_pcl.m_pos[ptIndex].x, io_pcl.m_pos[ptIndex].y);

          //get normal and normalize it, up:
          io_pcl.m_normal[ptIndex] = approxPlane.GetNormal();
          if (io_pcl.m_normal[ptIndex].z < 0)
            io_pcl.m_normal[ptIndex].z = -io_pcl.m_normal[ptIndex].z;
          Normalize(io_pcl.m_normal[ptIndex]);
        }
      }
    }

    if (!gotHashed)
    {
      delete in_pclHash;
    }
  }


  void Features::DenoiseRange(const CPtCloud& in_pcl, CPtCloud& out_pcl, int in_windowSize, float in_noiseTh)
  {
    int medFiltSize0 = in_windowSize;
    int medFiltSize1 = std::max(medFiltSize0 - 2, 1);
    int numlines = in_pcl.m_numPts / in_pcl.m_lineWidth;

    //convert to range image (each pixle in ptsLocal has it's xyz).
    float* rangeImage = new float[in_pcl.m_numPts];
    #pragma omp parallel for
    for (int ptrIndex = 0; ptrIndex < in_pcl.m_numPts; ptrIndex++)
    {
      rangeImage[ptrIndex] = Length(in_pcl.m_pos[ptrIndex]);
    }

    //find diferent points:
    float* distFiltered = new float[in_pcl.m_numPts];
    float* threshFiltered = new float[in_pcl.m_numPts];

    Median2D(in_pcl.m_lineWidth, numlines, medFiltSize0, rangeImage, distFiltered);

#pragma omp parallel for
    for (int row = 0; row < numlines; row++)
    {
      for (int col = 0; col < in_pcl.m_lineWidth; col++)
      {
        int index = (row * in_pcl.m_lineWidth) + col;
        distFiltered[index] = (abs(distFiltered[index] - rangeImage[index]) < in_noiseTh) ? 1.0f : 0.0f;
      }
    }

    Median2D(in_pcl.m_lineWidth, numlines, medFiltSize1, distFiltered, threshFiltered);

    int outputSize = 0;
    for (int row = 0; row < numlines; row++)
    {
      for (int col = 0; col < in_pcl.m_lineWidth; col++)
      {
        int index = (row * in_pcl.m_lineWidth) + col;

        if ((distFiltered[index] == 1) && (threshFiltered[index] == 1))
        {
          out_pcl.m_pos[outputSize] = in_pcl.m_pos[index];
          if ((in_pcl.m_color != 0) && (out_pcl.m_color != 0))
            out_pcl.m_color[outputSize] = in_pcl.m_color[index];
          if ((in_pcl.m_normal != 0) && (out_pcl.m_normal != 0))
            out_pcl.m_normal[outputSize] = in_pcl.m_normal[index];
          outputSize++;
        }
      }
    }

    out_pcl.m_numPts = outputSize;

    delete[] rangeImage;
    delete[] distFiltered;
    delete[] threshFiltered;
  }


  void Features::DenoiseRangeOfPointCloud()
  {

  }


  void Features::DownSample(const CPtCloud& in_pcl, CPtCloud& out_pcl, float in_voxelSize)
  {
    float InvVoxelSize = 1.0f / in_voxelSize;

    //find min/max x/y/z:
    CVec3 l_bbox[2] = { in_pcl.m_pos[0], in_pcl.m_pos[0] };
    for (int ptrIndex = 1; ptrIndex < in_pcl.m_numPts; ptrIndex++)
    {
      l_bbox[0] = Min_ps(l_bbox[0], in_pcl.m_pos[ptrIndex]);
      l_bbox[1] = Max_ps(l_bbox[1], in_pcl.m_pos[ptrIndex]);
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
    for (int ptrIndex = 0; ptrIndex < in_pcl.m_numPts; ptrIndex++)
    {
      float x = in_pcl.m_pos[ptrIndex].x;
      float y = in_pcl.m_pos[ptrIndex].y;
      float z = in_pcl.m_pos[ptrIndex].z;

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
        out_pcl.m_pos[outputSize] = in_pcl.m_pos[ptrIndex];
        if ((in_pcl.m_color != 0) && (out_pcl.m_color != 0))
          out_pcl.m_color[outputSize] = in_pcl.m_color[ptrIndex];
        if ((in_pcl.m_normal != 0) && (out_pcl.m_normal != 0))
          out_pcl.m_normal[outputSize] = in_pcl.m_normal[ptrIndex];
        outputSize++;

        VisitedVoxel[index] = true;
      }
    }
    delete[] VisitedVoxel;

    out_pcl.m_numPts = outputSize;
  }


  float Features::RMSEofRegistration(CSpatialHash2D* in_pcl1, const CPtCloud& in_pcl2, float in_max2DRadius, const CMat4& in_Rt)
  {
    double RMSE = 0;
    float outPenalty = in_max2DRadius*in_max2DRadius*1.5f;

    #pragma omp parallel for reduction(+:RMSE)
    for (int i = 0; i<in_pcl2.m_numPts; i++)
    {
      CVec3 transformedPt;
      CVec3 closestPt;

      // transform point according to R|t
      MultiplyVectorRightSidePlusOffset(in_Rt, in_pcl2.m_pos[i], transformedPt);

      // search for match
      float dist = outPenalty;
      if (in_pcl1->FindNearest(transformedPt, &closestPt, in_max2DRadius))
      {
        dist = DistSqr(transformedPt, closestPt);
      }

      RMSE += dist;
    }
    
    RMSE /= in_pcl2.m_numPts;
    return float(sqrt(RMSE));
  }


  void Features::CalcRotateMatZaxisToNormal(const CVec3& Xi_Normal, CMat4& Xo_RotateMat, const CVec3& Xi_Pos)
  {
    //find refFrame matrix:
    //zVec = Xi_Normal
    CVec3 xVec(1, 0, 0);
    xVec = xVec - DotProd(xVec, Xi_Normal)*Xi_Normal;
    Normalize(xVec);
    CVec3 yVec = CrossProd(Xi_Normal, xVec);
    Normalize(yVec);

    //update transformation matrix:
    //float OrientVals[] = { xVec.x   , yVec.x  , Xi_Normal.x, 0.0f,
    //                       xVec.y   , yVec.y  , Xi_Normal.y, 0.0f,
    //                       xVec.z   , yVec.z  , Xi_Normal.z, 0.0f,
    //                       Xi_Pos.x , Xi_Pos.y, Xi_Pos.z   , 1.0f };
    //TODO: temp logic fix to rotate matrix direction ambiguity, waiting for David's reply.
    float OrientVals[] = { xVec.x     , xVec.y     , xVec.z     , 0.0f,
      yVec.x     , yVec.y     , yVec.z     , 0.0f,
      Xi_Normal.x, Xi_Normal.y, Xi_Normal.z, 0.0f,
      Xi_Pos.x   , Xi_Pos.y   , Xi_Pos.z   , 1.0f };
    Xo_RotateMat = CMat4(OrientVals);
  }


} //namespace SLDR