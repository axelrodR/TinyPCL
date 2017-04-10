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
*: Package Name: sldrcr_sp
*
******************************************************************************/
#include "registration.h"
#include "features.h"
#include "SpatialHash.h"
#include "tran.h"
#include "common.h"
#include <vector>
#include <complex>
#include <algorithm>


//#ifdef _DEBUG
//#define new_file DEBUG_NEW
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

  struct CRegOptions
  {
    //parameters:
    float m_voxelSizeGlobal;      // if >0, use voxel grid downsampling.
    float m_voxelSizeLocal;       // if >0, use voxel grid downsampling.
    float m_d_grid;               // grid resolution.
    float m_d_sensor;             // dist from ground for grid point.
    int m_lineWidth;              // width of the polar depth map
    int m_numlines;               // height of the polar depth map 
    float m_searchRange;          // range of grid locations to check. if not using GPS this needs to be set to inf.
    int m_medFiltSize0;           // deniseing median filter size for the range imgae.
    int m_medFiltSize1;           // deniseing median filter size for the thresh image.
    float m_distFromMedianThresh; // max distance between point and median filter's result.
    float m_r_max;                // maximum distance from grid point for descriptor creation.
    float m_r_min;                // minimum distance from grid point for descriptor creation.

    void SetDefaults();
    CRegOptions() { SetDefaults(); }
  };



  class CRegDictionary
  {
  public:
    CMat4* m_Orient;      // 4x4 matrix, camera location & orientations.
    int m_size;
  };



  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CRegDictionary
  //
  ///////////////////////////////////////////////////////////////////////////////
  class CRegDictionaryE : public CRegDictionary
  {
  public:
    float** m_descriptors;          // descriptors array     (per orientation).
    std::complex<float>** m_descriptorsDFT;       // descriptors DFT array (per orientation).

    CSpatialHash2D m_globalHashed;  //hashed global point cloud.

    /** destructor */
    CRegDictionaryE::CRegDictionaryE()
    {
      m_Orient = NULL;
      m_descriptors = NULL;
      m_descriptorsDFT = NULL;
      m_size = -1;
    }

    /** destructor */
    CRegDictionaryE::~CRegDictionaryE()
    {
      DeleteDictionary();
    }

    /** free allocated memory of the dictionary. */
    void DeleteDictionary()
    {
      for (int dicIndex = 0; dicIndex < m_size; dicIndex++)
      {
        delete[] m_descriptors[dicIndex];
        delete[] m_descriptorsDFT[dicIndex];

      }
      delete[] m_descriptors;
      delete[] m_descriptorsDFT;
      m_descriptors = NULL;
      m_descriptorsDFT = NULL;

      delete[] m_Orient;
      m_Orient = NULL;

      m_size = -1;

      m_globalHashed.Clear();
    }


    /** builds the hashed point cloud of the input point cloud, in the member m_globalHashed.
    * @param Xi_pts           global point cloud */
    void HashPointCloud(int Xi_numPts, const CVec3* Xi_pts)
    {
      for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
      {
        m_globalHashed.Add(Xi_pts[ptrIndex], NULL);
      }
    }

    /** deletes previous dictionary and fills the dictionary with the grid - sets of location & rotation per entry.
    * @param Xi_d_grid           Resolution of the grid.
    * @param Xi_d_sensor         dist from ground for grid point.
    * @param Xi_minNumNBRS       Min number of neighbors in global for a grid point.
    * @param Xi_pts           global point cloud */
    void ViewpointGridCreation(float Xi_d_grid, float Xi_d_sensor, int Xi_numPts, const CVec3* Xi_pts)
    {
      DeleteDictionary();
      float invGridRes = 1.0f / Xi_d_grid;

      //hash global cloud:
      HashPointCloud(Xi_numPts, Xi_pts);

      //find grid's boundaries:
      CVec3 bbox[2] = { Xi_pts[0] , Xi_pts[0] };
      for (int i = 1; i < Xi_numPts; ++i)
      {
        bbox[0] = Min_ps(bbox[0], Xi_pts[i]);
        bbox[1] = Max_ps(bbox[1], Xi_pts[i]);
      }

      int GridWidth = int(ceil((bbox[1].x - bbox[0].x) * invGridRes));
      int GridHeight = int(ceil((bbox[1].y - bbox[0].y) * invGridRes));
      int totalGridPoint = GridHeight*GridWidth;

      //create grid's locations (taking z value from the closest point):
      CVec3* gridPositions = new CVec3[totalGridPoint];
      for (int yGrid = 0; yGrid < GridHeight; yGrid++)
      {
        for (int xGrid = 0; xGrid < GridWidth; xGrid++)
        {
          CVec3 pos(bbox[0].x + xGrid*Xi_d_grid, bbox[0].y + yGrid*Xi_d_grid, 0);
          CVec3 closest;
          m_globalHashed.FindNearest(pos, &closest);
          pos.z = closest.z;
          gridPositions[yGrid*GridWidth + xGrid] = pos;
        }
      }

      //find normals:
      const float maxDistForPlane = 10;
      CVec3* Normals = new CVec3[totalGridPoint];
      //TODO:      Features::FillPointCloud();
      Features::FindNormal(maxDistForPlane, m_globalHashed, totalGridPoint, gridPositions, Normals);

      //set viewpoints height above ground (in the normal vector direction above the plane that was found):
      for (int yGrid = 0; yGrid < GridHeight; yGrid++)
      {
        for (int xGrid = 0; xGrid < GridWidth; xGrid++)
        {
          int index = yGrid*GridWidth + xGrid;
          gridPositions[index] = gridPositions[index] + (Xi_d_sensor * Normals[index]);
        }
      }

      m_size = totalGridPoint;
      m_Orient = new CMat4[totalGridPoint];

      //creating final grid's transformation matrix:
      for (int yGrid = 0; yGrid < GridHeight; yGrid++)
      {
        for (int xGrid = 0; xGrid < GridWidth; xGrid++)
        {
          //find refFrame matrix:
          int index = yGrid*GridWidth + xGrid;
          //zVec = Normals[index]
          CVec3 xVec(1, 0, 0);
          xVec = xVec - DotProd(xVec, Normals[index])*Normals[index];
          Normalize(xVec);
          CVec3 yVec = CrossProd(Normals[index], xVec);
          Normalize(yVec);

          //update transformation matrix:
          float OrientVals[] = { xVec.x, yVec.x, Normals[index].x, 0.0f,
                                 xVec.y, yVec.y, Normals[index].y, 0.0f,
                                 xVec.z, yVec.z, Normals[index].z, 0.0f,
                                 gridPositions[index].x, gridPositions[index].y, gridPositions[index].z, 1.0f };
          m_Orient[index] = CMat4(OrientVals);
        }
      }

      //release memory:
      delete[] gridPositions;
      delete[] Normals;

      //create vectors for the descriptors:
      m_descriptors    = new float*[m_size];
      m_descriptorsDFT = new std::complex<float>*[m_size];

    }



    /** Convert a point cloud into the descriptor - range image, a polar map of distance around a point
    * @param Xi_lineWidth      width of the polar depth map.
    * @param Xi_numlines       height of the polar depth map.
    * @param Xi_pts            list of points in the cloud.
    * @param Xi_origin         the center of the polar map.
    * @param Xo_RangeImage   polar depth map.  */
    void ConvertPCL2descriptor(int Xi_lineWidth, int Xi_numlines, int Xi_numPts, const CVec3* Xi_pts, float* Xo_RangeImage)
    {
      float azimuthRes = 2 * float(M_PI) / Xi_lineWidth;
      float elevationRes = float(M_PI) / Xi_numlines;
      //fill range image with default value - 0.
      for (int elevation = 0; elevation < Xi_numlines; elevation++)
      {
        for (int azimuth = 0; azimuth < Xi_lineWidth; azimuth++)
        {
          Xo_RangeImage[elevation*Xi_lineWidth + azimuth] = 0;
        }
      }

      //put closest range per azimuth & elevation in the range image.
      for (int i = 0; i < Xi_numPts; i++)
      {
        float x = Xi_pts[i].x;
        float y = Xi_pts[i].y;
        float z = Xi_pts[i].z;
        float azimuth = atan2(y, x);
        float elevation = atan2(z, sqrt(x*x + y*y));
        float r = sqrt(x*x + y*y + z*z);


        int azimuthInds = int(floor((azimuth + M_PI) / azimuthRes));
        int elevationInds = int(floor((elevation + M_PI / 2) / elevationRes));

        int index = elevationInds*Xi_lineWidth + azimuthInds;
        if (r < Xo_RangeImage[index])
          Xo_RangeImage[index] = r;
      }

    }


    /** find FFT of a 2D descriptor.
    * @param Xi_desc1Dsize      size of the descriptor's first dimension (width of the range image - lineWidth)
    * @param Xi_desc2Dsize      size of the descriptor's second dimension (height of the range image - numlines).
    * @param Xi_Descriptor      2D descriptor - range image. 
    * @param Xo_DescriptorFFT   descriptor's 2D FFT.  */
    void ConvertDescriptor2FFT(int Xi_desc1Dsize, int Xi_desc2Dsize, float* Xi_Descriptor, std::complex<float>* Xo_DescriptorFFT)
    {
      for (int index2 = 0; index2 < Xi_desc2Dsize; index2++)
      {
        for (int index1 = 0; index1 < Xi_desc1Dsize; index1++)
        {
          int index = (index2 * Xi_desc1Dsize) + index1;
          Xo_DescriptorFFT[index] = Xi_Descriptor[index];
        }
      }

      DFT2D(unsigned int(Xi_desc1Dsize), unsigned int(Xi_desc2Dsize), Xo_DescriptorFFT);
    }


    /** fills the dictionary entries with their descriptors.
    * @param Xi_r_max        maximum distance from grid point for descriptor creation.
    * @param Xi_r_min        minimum distance from grid point for descriptor creation.
    * @param Xi_pts          global point cloud
    * @param Xi_desc1Dsize   size of the descriptor's first dimension (width of the range image - lineWidth)
    * @param Xi_desc2Dsize   size of the descriptor's second dimension (height of the range image - numlines).  */
    void DescriptorsCreation(float Xi_r_max, float Xi_r_min, int Xi_numPts, const CVec3* Xi_pts, int Xi_desc1Dsize, int Xi_desc2Dsize = 0)
    {
      const int maxInRangePoints = 256;
      int nearGotten = 0;
      CVec3 ptsInRange[maxInRangePoints];
      CVec3 ptsTransformed[maxInRangePoints];
      void* UnusedBuffer[maxInRangePoints];

      
      //go over grid points and create descriptors for each point:
      for (int gridIndex = 0; gridIndex < m_size; gridIndex++)
      {
        CMat4 orients = m_Orient[gridIndex];

        CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);
        nearGotten = m_globalHashed.GetNear(Pos, maxInRangePoints, UnusedBuffer, ptsInRange, Xi_r_max);

        int TransformedSize = 0;

        for (int nearIndex = 0; nearIndex < nearGotten; nearIndex++)
        {
          //check if point too close:
          float dist = Dist(Pos, ptsInRange[nearIndex]);
          if (dist < Xi_r_min)
            continue;

          //transform point:
          CVec3 PosShifted = ptsInRange[nearIndex] - Pos;
          //DO: make sure this is the mult he meant.
          CVec3 MatRow0 = CVec3(orients.m[0][0], orients.m[0][1], orients.m[0][2]);
          CVec3 MatRow1 = CVec3(orients.m[1][0], orients.m[1][1], orients.m[1][2]);
          CVec3 MatRow2 = CVec3(orients.m[2][0], orients.m[2][1], orients.m[2][2]);
          ptsTransformed[TransformedSize].x = DotProd(PosShifted, MatRow0);
          ptsTransformed[TransformedSize].y = DotProd(PosShifted, MatRow1);
          ptsTransformed[TransformedSize].z = DotProd(PosShifted, MatRow2);

          TransformedSize++;
        }
        
        //create descriptor - range image, FFT
        delete[] m_descriptors[gridIndex];
        m_descriptors[gridIndex] = new float[Xi_desc1Dsize*Xi_desc2Dsize];
        ConvertPCL2descriptor(Xi_desc1Dsize, Xi_desc2Dsize, TransformedSize, ptsTransformed, m_descriptors[gridIndex]);

        delete[] m_descriptorsDFT[gridIndex];
        m_descriptorsDFT[gridIndex] = new std::complex<float>[Xi_desc1Dsize*Xi_desc2Dsize];
        ConvertDescriptor2FFT(Xi_desc1Dsize, Xi_desc2Dsize, m_descriptors[gridIndex], m_descriptorsDFT[gridIndex]);
      }

    }

    /** calculate best phase correlation.
    * @param Xi_descriptorFFT0             input fisrt descriptor.
    * @param Xi_descriptorFFT1             input second descriptor.
    * @param Xo_bestRow                    row of best score.
    * @param Xo_bestCol                    column of best score.
    * @param Xo_bestScore                  list of the candidates grades.
    * @param Xi_desc1Dsize                 size of the descriptor's first dimension (width of the range image - lineWidth).
    * @param Xi_desc2Dsize                 size of the descriptor's second dimension (height of the range image - numlines). */
    void BestPhaseCorr(std::complex<float>* Xi_descriptorFFT0, std::complex<float>* Xi_descriptorFFT1, int& Xo_bestRow, int& Xo_bestCol, float& Xo_bestScore, int Xi_desc1Dsize, int Xi_desc2Dsize = 0)
    {
      Xo_bestScore = FLT_MIN;
      Xo_bestRow = 0;
      Xo_bestCol = 0;

      std::complex<float>* PhCor = new std::complex<float>[Xi_desc1Dsize * Xi_desc2Dsize];
      float lengthSum = 0;
      for (int index2 = 0; index2 < Xi_desc2Dsize; index2++)
      {
        for (int index1 = 0; index1 < Xi_desc1Dsize; index1++)
        {
          int index = (index2 * Xi_desc1Dsize) + index1;

          PhCor[index] = Xi_descriptorFFT0[index] * conj(Xi_descriptorFFT1[index]);
          lengthSum += std::abs(PhCor[index]);
        }
      }

      if (lengthSum == 0)
        return;

      //normalize:
      for (int index = 0; index < Xi_desc1Dsize*Xi_desc2Dsize; index++)
      {
        PhCor[index] /= lengthSum;
      }

      DFT2D(Xi_desc1Dsize, Xi_desc2Dsize, PhCor, false);

      //TODO: check if shift is needed. sent same image twice and see if in the result the delta is in the center or not.

      //find max:
      for (int index2 = 0; index2 < Xi_desc2Dsize; index2++)
      {
        for (int index1 = 0; index1 < Xi_desc1Dsize; index1++)
        {
          int index = (index2 * Xi_desc1Dsize) + index1;

          if (PhCor[index].real() > Xo_bestScore)
          {
            Xo_bestScore = PhCor[index].real();
            Xo_bestRow = index2;
            Xo_bestCol = index1;
          }
        }
      }
    }




    /** return candidates suitable for the input descriptor.
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_descriptor        input local descriptor.
    * @param Xo_candidates        list of the candidates.
    * @param Xo_grades            list of the candidates grades. 
    * @param Xo_rotations         list of the candidates rotation.
    * @param Xi_desc1Dsize        size of the descriptor's first dimension (width of the range image - lineWidth)
    * @param Xi_desc2Dsize        size of the descriptor's second dimension (height of the range image - numlines).  
    * @return                     number of candidates.*/
    int SearchDictionary(int Xi_maxCandidates, float Xi_searchRadius, std::complex<float>* Xi_descriptorFFT, int* Xo_candidates, float* Xo_grades, CMat4* Xo_rotations, int Xi_desc1Dsize, int Xi_desc2Dsize = 0, CVec3* Xi_estimatePos = NULL)
    {
      float finalMax = FLT_MIN;
      int finalIndex = -1;
      int finalrow = -1;
      int finalcol = -1;


      int NumOfCandidates = 0;
      int minIndex = 0;
      Xo_grades[minIndex] = FLT_MAX;

      int* bestCols = new int[Xi_maxCandidates];

      int gridIndex = 0;
      //take first Xi_maxCandidates entries from dictionary which are closer than Xi_searchRadius from GPS estimation:
      for (gridIndex; gridIndex < m_size; gridIndex++)
      {
        CMat4 orients = m_Orient[gridIndex];
        CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);

        //ignore entries which are too far from the GPS guess:
        if (Xi_estimatePos)
        {
          if (Dist(*Xi_estimatePos, Pos) > Xi_searchRadius)
            continue;
        }

        std::complex<float>* gridDescDFT = m_descriptorsDFT[gridIndex];
        float bestMax = FLT_MIN;
        int bestRow = -1;
        int bestCol = -1;
        BestPhaseCorr(gridDescDFT, Xi_descriptorFFT, bestRow, bestCol, bestMax, Xi_desc1Dsize, Xi_desc2Dsize);

        Xo_candidates[NumOfCandidates] = gridIndex;
        Xo_grades[NumOfCandidates] = bestMax;
        bestCols[NumOfCandidates] = bestCol;

        if (bestMax < Xo_grades[minIndex])
        {
          minIndex = NumOfCandidates;
        }

        NumOfCandidates++;

        if (NumOfCandidates == Xi_maxCandidates)
          break;
      }


      //continue combing dictionary (remaining with best Xi_maxCandidates Candidates:
      for (gridIndex; gridIndex < m_size; gridIndex++)
      {
        CMat4 orients = m_Orient[gridIndex];
        CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);
        //ignore entries which are too far from the GPS guess:
        if (Xi_estimatePos)
        {
          if (Dist(*Xi_estimatePos, Pos) > Xi_searchRadius)
            continue;
        }

        std::complex<float>* gridDescDFT = m_descriptorsDFT[gridIndex];
        float bestMax = FLT_MIN;
        int bestRow = -1;
        int bestCol = -1;
        BestPhaseCorr(gridDescDFT, Xi_descriptorFFT, bestRow, bestCol, bestMax, Xi_desc1Dsize, Xi_desc2Dsize);
        if (bestMax > Xo_grades[minIndex])
        {
          Xo_grades[minIndex] = bestMax;
          Xo_candidates[minIndex] = gridIndex;
          bestCols[minIndex] = bestCol;


          for (int candIndex = 0; candIndex < NumOfCandidates; candIndex++)
          {
            if (Xo_grades[candIndex] < Xo_grades[minIndex])
              minIndex = candIndex;
          }
        }
      }


      //compute rotation matrix for candidates:
      float azimuthRes = 2 * float(M_PI) / Xi_desc1Dsize;

      for (int candIndex = 0; candIndex < NumOfCandidates; candIndex++)
      {
        CMat4 orients = m_Orient[Xo_candidates[candIndex]];
        CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);
        
        //calc beast azimuth in radians.
        float centerShift = ((Xi_desc1Dsize & 1) == 0) ? 0.5f : 0.0f;
        float peak_azimuth = float(bestCols[candIndex] + centerShift) * azimuthRes - float(M_PI);
        
        float OrientVals[] = { cos(peak_azimuth), -sin(peak_azimuth), 0.0f, 0.0f,
                               sin(peak_azimuth),  cos(peak_azimuth), 0.0f, 0.0f,
                               0.0f             ,  0.0f             , 1.0f, 0.0f,
                               0.0f             ,  0.0f             , 0.0f, 1.0f };


        TransposeLeftMultiply(orients, CMat4(OrientVals), Xo_rotations[candIndex]); //TODO: confirm.

        Xo_rotations[candIndex].m[3][0] = Pos.x;
        Xo_rotations[candIndex].m[3][1] = Pos.y;
        Xo_rotations[candIndex].m[3][2] = Pos.z;
        Xo_rotations[candIndex].m[0][3] = 0;
        Xo_rotations[candIndex].m[1][3] = 0;
        Xo_rotations[candIndex].m[2][3] = 0;
        Xo_rotations[candIndex].m[3][3] = 1;

      }

      delete[] bestCols;
      return NumOfCandidates;
    }
  };





  /** 2D median filter.
  * @param Xi_lineWidth                image width.
  * @param Xi_numlines                 image height.
  * @param Xi_medFiltSize              median filter size.
  * @param Xi_pts                      input image.
  * @param Xo_ptsFiltered              2D median filtered image. assumes  Xo_ptsFiltered != Xi_pts*/
  void Median2DPowOf2(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize, float* Xi_pts, float* Xo_ptsFiltered)
  {
    int windowSize = Xi_medFiltSize*Xi_medFiltSize;
    int half1DWindow = Xi_medFiltSize >> 1;
    float* window = new float[windowSize];


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


  /** denose a xyz image and return a denoised point cloud.
  * @param Xi_lineWidth                width of the polar depth map.
  * @param Xi_numlines                 height of the polar depth map
  * @param Xi_medFiltSize0             median filter size for the range imgae.
  * @param Xi_medFiltSize1             median filter size for the thresh image.
  * @param Xi_distFromMedianThresh     max distance between point and median filter's result.
  * @param Xi_pts                      input xyz image. assuming row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
  * @param Xo_ptsDenoised              denoised point cloud. assumes size at least as Xi_pts size.
  * @return                            size of point cloud.  */
  int DenoiseOrderedPointCloud(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, CVec3* Xi_pts, CVec3* Xo_ptsDenoised)
  {
    int totalSize = Xi_lineWidth * Xi_numlines;

    //convert to range image (each pixle in ptsLocal has it's xyz).
    float* rangeImage = new float[totalSize];
    for (int ptrIndex = 0; ptrIndex < totalSize; ptrIndex++)
    {
      rangeImage[ptrIndex] = Length(Xi_pts[ptrIndex]);
    }

    //find diferent points:
    float* distFiltered = new float[totalSize];
    float* threshFiltered = new float[totalSize];

    Median2DPowOf2(Xi_lineWidth, Xi_numlines, Xi_medFiltSize0, rangeImage, distFiltered);

    for (int row = 0; row < Xi_numlines; row++)
    {
      for (int col = 0; col < Xi_lineWidth; col++)
      {
        int index = (row * Xi_lineWidth) + col;
        distFiltered[index] = (abs(distFiltered[index] - rangeImage[index]) < Xi_distFromMedianThresh) ? 1.0f : 0.0f;
      }
    }

    Median2DPowOf2(Xi_lineWidth, Xi_numlines, Xi_medFiltSize1, distFiltered, threshFiltered);

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



  /** denose a point cloud.
  * @param Xi_lineWidth                width of the polar depth map.
  * @param Xi_numlines                 height of the polar depth map
  * @param Xi_medFiltSize0             median filter size for the range imgae.
  * @param Xi_medFiltSize1             median filter size for the thresh image.
  * @param Xi_distFromMedianThresh     max distance between point and median filter's result.
  * @param Xi_pts                      input xyz image. assuming row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
  * @param Xo_ptsDenoised              denoised point cloud. assumes size at least as Xi_pts size.
  * @return                            size of point cloud.  */
  int DenoisePointCloud(float Xi_res, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, int Xi_numPts, CVec3* Xi_pts, CVec3* Xo_ptsDenoised)
  {
    //float res = float(2 * M_PI) / (128 * 5); //TODO

    //CSpatialHash2D hashedCloud(Xi_res);
    //float* range = new float[Xi_numPts];
    //bool* underThresh = new bool[Xi_numPts];
    //for (int ptIndex = 0; ptIndex < Xi_numPts; ptIndex++)
    //{
    //  underThresh[ptIndex] = false;
    //  range[ptIndex] = Length(Xi_pts[ptIndex]);
    //  hashedCloud.Add(Xi_pts[ptIndex], range + ptIndex);
    //}

    //int bufSize = Xi_medFiltSize0 * Xi_medFiltSize0;
    //void** buf = new void* [bufSize];
    //CVec3* nearPts = new CVec3[bufSize];
    //float max2DRadius = (Xi_medFiltSize0 >> 1) * res;
    //for (int ptIndex = 0; ptIndex < Xi_numPts; ptIndex++)
    //{
    //  hashedCloud.GetNear(Xi_pts[ptIndex], bufSize, buf, nearPts, max2DRadius);

    //  std::nth_element((float*)(*buf), (float*)(*buf) + (bufSize>>1), (float*)(*buf) + bufSize);
    //}
    //
    return 0;
  }



  /** downsample a point cloud. Divides to grid from minXYZ (of pts) to max XYZ, of size m_voxelSize.
   *  If more than one point in same grid index, takes only first one.
   * @param Xi_voxelSize            size of a voxel in grid. if equal zero - NO downsampling.
   * @param Xio_pts                 in: point cloud,  out: downsampled point cloud. */
  void DownSamplePointCloud(float Xi_voxelSize, int& Xio_numPts, CVec3* Xio_pts)
  {
    float InvVoxelSize = 1.0f / Xi_voxelSize;

    //find min/max x/y/z:
    float xMin = Xio_pts[0].x; float xMax = xMin;
    float yMin = Xio_pts[0].y; float yMax = yMin;
    float zMin = Xio_pts[0].z; float zMax = zMin;

    for (int ptrIndex = 1; ptrIndex < Xio_numPts; ptrIndex++)
    {
      float x = Xio_pts[ptrIndex].x;
      float y = Xio_pts[ptrIndex].y;
      float z = Xio_pts[ptrIndex].z;

      xMin = MinT(xMin, x);    yMin = MinT(yMin, y);    zMin = MinT(zMin, z);
      xMax = MaxT(xMax, x);    yMax = MaxT(yMax, y);    zMax = MaxT(zMax, z);
    }

    //create downsampling vector:
    int Mx = int(ceil((xMax - xMin) * InvVoxelSize));
    int My = int(ceil((yMax - yMin) * InvVoxelSize));
    int Mz = int(ceil((zMax - zMin) * InvVoxelSize));
    int Mxy = Mx*My;

    bool* VisitedVoxel = new bool[Mxy*Mz];
    memset(VisitedVoxel, false, Mxy*Mz * sizeof(bool));
    std::vector<CVec3> DownSampledPts;

    for (int ptrIndex = 0; ptrIndex < Xio_numPts; ptrIndex++)
    {
      float x = Xio_pts[ptrIndex].x;
      float y = Xio_pts[ptrIndex].y;
      float z = Xio_pts[ptrIndex].z;

      int xInd = int(floor((x - xMin) * InvVoxelSize));
      int yInd = int(floor((y - yMin) * InvVoxelSize));
      int zInd = int(floor((z - zMin) * InvVoxelSize));

      if (xInd == Mx) xInd--;
      if (yInd == My) yInd--;
      if (zInd == Mz) zInd--;

      int index = zInd*Mxy + yInd*Mx + xInd;
      if (!VisitedVoxel[index])
      {
        DownSampledPts.push_back(Xio_pts[ptrIndex]);
        VisitedVoxel[index] = true;
      }
    }
    delete[] VisitedVoxel;


    //save selected points:
    for (int VecIndex = 0; VecIndex < DownSampledPts.size(); VecIndex++)
    {
      Xio_pts[VecIndex] = DownSampledPts[VecIndex];
    }
    Xio_numPts = int(DownSampledPts.size());
  }



  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CCoarseRegister
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: tpclCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::CCoarseRegister()
  {
    m_opts = new CRegOptions;
    //CRegOptions* optsP = (CRegOptions*)m_opts;
    m_dictionary = new CRegDictionaryE;
  }

  /******************************************************************************
  *
  *: Method name: ~tpclCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::~CCoarseRegister()
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    delete optsP;

    CRegDictionaryE* dictionaryP = (CRegDictionaryE*)m_dictionary;
    delete dictionaryP;

    delete[] m_ptsGlobal;
  }


  /******************************************************************************
  *
  *: Method name: SetGlobalCloud
  *
  ******************************************************************************/
  void CCoarseRegister::BuildDictionary(int Xi_numPts, const CVec3* Xi_pts)
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionaryE* dictionaryP = (CRegDictionaryE*)m_dictionary;

    //copy global cloud to local variable:
    delete[] m_ptsGlobal;
    m_ptsGlobal = new CVec3[Xi_numPts];
    for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
    {
      m_ptsGlobal[ptrIndex] = Xi_pts[ptrIndex];
    }
    m_numPtsGlobal = Xi_numPts;

    //preprocess global cloud:
    DownSamplePointCloud(optsP->m_voxelSizeGlobal, m_numPtsGlobal, m_ptsGlobal);

    //create grid's position and orientation: 
    dictionaryP->ViewpointGridCreation(optsP->m_d_grid, optsP->m_d_sensor, Xi_numPts, Xi_pts);

    //fill the dictionary entries with their descriptors:
    dictionaryP->DescriptorsCreation(optsP->m_r_max, optsP->m_r_min, Xi_numPts, Xi_pts, optsP->m_lineWidth, optsP->m_numlines);
  }


  /******************************************************************************
  *
  *: Method name: GetLocalRegistrationCandidates
  *
  ******************************************************************************/
  int CCoarseRegister::GetLocalRegistrationCandidates(int Xi_maxCandidates, CVec3* Xi_pts, CVec3 Xi_originApprox, int* Xo_candidates, float* Xo_grades, CMat4* Xo_rotations, int Xi_numPts, int Xi_numlines, CVec3* Xi_GPS)
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionaryE* dictionaryP = (CRegDictionaryE*)m_dictionary;
    CVec3* ptsPrePro = new CVec3[Xi_numPts];
    int totalPixels;

    //preprocess local cloud:
    if (Xi_numlines > 0)
    {
      int lineWidth = Xi_numPts / Xi_numlines;
      totalPixels = DenoiseOrderedPointCloud(lineWidth, Xi_numlines, optsP->m_medFiltSize0, optsP->m_medFiltSize1, optsP->m_distFromMedianThresh, Xi_pts, ptsPrePro);
    }
    else
    {
      float res = float(2 * M_PI) / (128 * 5);
      totalPixels = DenoisePointCloud(res, optsP->m_medFiltSize0, optsP->m_medFiltSize1, optsP->m_distFromMedianThresh, Xi_numPts, Xi_pts, ptsPrePro);
    }

    DownSamplePointCloud(optsP->m_voxelSizeLocal, totalPixels, ptsPrePro);

    ////remove points that too close to sensor:
    int postProSize = 0;
    for (int ptrIndex = 0; ptrIndex < totalPixels; ptrIndex++)
    {
      float len = Length(ptsPrePro[ptrIndex]);
      if (len < optsP->m_r_min)
        continue;
      ptsPrePro[postProSize] = ptsPrePro[ptrIndex];
      postProSize++;
    }
    totalPixels = postProSize;

    //create range image:
    float* descriptor = new float[optsP->m_lineWidth * optsP->m_numlines];
    dictionaryP->ConvertPCL2descriptor(optsP->m_lineWidth, optsP->m_numlines, totalPixels, ptsPrePro, descriptor);

    //create range image's 2D FFT:
    std::complex<float>* descriptorFFT = new std::complex<float>[totalPixels];
    dictionaryP->ConvertDescriptor2FFT(optsP->m_lineWidth, optsP->m_numlines, descriptor, descriptorFFT);
    
    
    int NumOfCandidates = dictionaryP->SearchDictionary(Xi_maxCandidates, optsP->m_searchRange, descriptorFFT, Xo_candidates, Xo_grades, Xo_rotations, optsP->m_lineWidth, optsP->m_numlines, Xi_GPS);

    delete[] ptsPrePro;
    delete[] descriptor;
    delete[] descriptorFFT;

    return NumOfCandidates;
  }



  /******************************************************************************
  *
  *: Method name: GetLocalRegistration
  *
  ******************************************************************************/
  void CCoarseRegister::GetLocalRegistration(int Xi_NumOfCandidates, int* Xi_candidates, float* Xi_grades, CMat4* Xi_rotations, CMat4& Xo_best)
  {
    int maxIndex = 0;
    float maxGrade = Xi_grades[0];
    //search for candidate with highest grade:
    for (int canIndex = 1; canIndex < Xi_NumOfCandidates; canIndex++)
    {
      if (Xi_grades[canIndex] > maxGrade)
      {
        maxIndex = canIndex;
        maxGrade = Xi_grades[canIndex];
      }
    }

    Xo_best = Xi_rotations[maxIndex];
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


  void CRegOptions::SetDefaults()
  {
    m_voxelSizeGlobal = 0.5;       // if >0, use voxel grid downsampling.
    m_voxelSizeLocal = 0.25;       // if >0, use voxel grid downsampling.
    m_d_grid = 3;                  // grid resolution.
    m_d_sensor = 2;                // dist from ground for grid point.
    m_lineWidth = 128;             // width of the polar depth map
    m_numlines = 64;               // height of the polar depth map 
    m_searchRange = 30;            // range of grid locations to check. if not using GPS this needs to be set to inf.
    m_medFiltSize0 = 7;            // deniseing median filter size for the range imgae.
    m_medFiltSize1 = 5;            // deniseing median filter size for the thresh image.
    m_distFromMedianThresh = 0.03f; // max distance between point and median filter's result.
    m_r_max = 60;                  // maximum distance from grid point for descriptor creation.
    m_r_min = 2;                   // minimum distance from grid point for descriptor creation.
  }


} //namespace tpcl