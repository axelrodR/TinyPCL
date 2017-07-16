
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


#include "OrientDict.h"
#include "features.h"
#include "SpatialHash.h"
#include "common.h"
#include "tran.h"
#include <complex>
#include <vector>


//#define DEBUG_LOCAL_RANGE_IMAGE
#ifdef DEBUG_LOCAL_RANGE_IMAGE
#include "dbg.h"
#endif


namespace tpcl
{
  /** resize an array, retaining previous values and set new ellements to zero.
  * @param io_ptr          input: array of size in_size. output: array of size in_newSize with all ellements from input array.
  * @param in_size          size of input array.
  * @param in_newSize       size of output array. */
  template <typename T>
  static void resizeArray(T* &io_ptr, const int& in_size, const int& in_newSize)
  {
    T* ptrExtended = new T[in_newSize];
    memset(ptrExtended, 0, in_newSize * sizeof(T));
    memcpy(ptrExtended, io_ptr, (MinT(in_size, in_newSize)) * sizeof(T));
    delete[] io_ptr;
    io_ptr = ptrExtended;
  };




  /******************************************************************************
  *
  *: Class name: COrientedGrid
  *
  ******************************************************************************/


  COrientedGrid::COrientedGrid()
  {
    m_pclMain.m_color = NULL;    m_pclMain.m_normal = NULL;    m_pclMain.m_type = PCL_TYPE_FUSED;
    initMembers();
  }


  COrientedGrid::COrientedGrid(float in_voxelSize)
  {
    m_pclMain.m_color = NULL;    m_pclMain.m_normal = NULL;    m_pclMain.m_type = PCL_TYPE_FUSED;
    initMembers();
    m_voxelSize = in_voxelSize;
    m_mainHashed = new CSpatialHash2D(m_voxelSize);
    ((CSpatialHash2D*)(m_mainHashed))->Clear();
  }


  COrientedGrid::~COrientedGrid()
  {
    DeleteGrid();
  }


  float COrientedGrid::getVoxelSize()
  {
    return m_voxelSize;
  }


  void* COrientedGrid::getMainHashedPtr()
  {
    return m_mainHashed;
  }


  void COrientedGrid::getPclMainPtr(CPtCloud* &out_pclMain)
  {
    out_pclMain = &m_pclMain;
  }


  int COrientedGrid::getGrid(CMat4* &out_Orient)
  {
    out_Orient = m_Orient;
    return m_size;
  }


  void COrientedGrid::getBBox(CVec3 &out_minBBox, CVec3 &out_maxBBox)
  {
    out_minBBox = m_minBBox;
    out_maxBBox = m_maxBBox;
  }


  void COrientedGrid::DeleteGrid()
  {
    delete m_mainHashed;
    delete[] m_pclMain.m_pos;
    delete[] m_Orient;
  }


  void COrientedGrid::ResetGrid()
  {
    DeleteGrid();

    m_mainHashed = new CSpatialHash2D(m_voxelSize);
    ((CSpatialHash2D*)(m_mainHashed))->Clear();

    m_pclMain.m_pos = NULL;
    m_Orient = NULL;

    m_pclMain.m_numPts = 0;
    m_size = 0;
    m_minBBox = CVec3(0, 0, 0);
    m_maxBBox = CVec3(0, 0, 0);
  }


  void COrientedGrid::DeleteAndSetVoxelSize(float in_voxelSize)
  {
    m_voxelSize = in_voxelSize;
    ResetGrid();
  }


  void COrientedGrid::PointCloudUpdate(const CPtCloud& in_pcl, CVec3& out_minBox, CVec3& out_maxBox)
  {
    CSpatialHash2D& mainHashed = *((CSpatialHash2D*)(m_mainHashed));

    int totalPts = m_pclMain.m_numPts + in_pcl.m_numPts;
    resizeArray(m_pclMain.m_pos, m_pclMain.m_numPts, totalPts);

    out_minBox = out_maxBox = in_pcl.m_pos[0];
    for (int ptrIndex = 0; ptrIndex < in_pcl.m_numPts; ptrIndex++)
    {
      m_pclMain.m_pos[m_pclMain.m_numPts] = in_pcl.m_pos[ptrIndex];
      m_pclMain.m_numPts++;

      mainHashed.Add(in_pcl.m_pos[ptrIndex], (void*)(1));
      out_minBox = Min_ps(out_minBox, in_pcl.m_pos[ptrIndex]);
      out_maxBox = Max_ps(out_maxBox, in_pcl.m_pos[ptrIndex]);
    }
    //at end of for loop: m_pclMain.m_numPts = totalPts;

    m_minBBox = Min_ps(out_minBox, m_minBBox);;
    m_maxBBox = Max_ps(out_maxBox, m_maxBBox);;
  }


  int COrientedGrid::ViewpointGridUpdate(float in_d_grid, float in_d_sensor, CVec3& in_minBox, CVec3& in_maxBox)
  {
    Features feat;
    CSpatialHash2D& mainHashed = *((CSpatialHash2D*)(m_mainHashed));

    float invGridRes = 1.0f / in_d_grid;

    //create grid's locations (taking z value from the closest point):
    int GridWidth = int(ceil((in_maxBox.x - in_minBox.x) * invGridRes));
    int GridHeight = int(ceil((in_maxBox.y - in_minBox.y) * invGridRes));

    CPtCloud gridPositions;  gridPositions.m_type = PCL_TYPE_FUSED;   gridPositions.m_color = 0;
    gridPositions.m_numPts = GridHeight*GridWidth;
    gridPositions.m_pos = new CVec3[gridPositions.m_numPts];
    gridPositions.m_normal = new CVec3[gridPositions.m_numPts];

    #pragma omp parallel for
    for (int yGrid = 0; yGrid < GridHeight; yGrid++)
    {
      for (int xGrid = 0; xGrid < GridWidth; xGrid++)
      {
        CVec3 pos(in_minBox.x + xGrid*in_d_grid, in_minBox.y + yGrid*in_d_grid, 0);
        CVec3 closest;
        if (mainHashed.FindNearest(pos, &closest, in_d_grid))
          pos.z = closest.z;
        else
          pos.z = in_minBox.z;
        gridPositions.m_pos[yGrid*GridWidth + xGrid] = pos;
      }
    }

    //find normals:
    const float maxDistForPlane = 4;// m_voxelSize * 2;
    feat.FillNormals(gridPositions, maxDistForPlane, &mainHashed, true);

    int preSize = m_size;
    m_size += gridPositions.m_numPts;
    resizeArray(m_Orient, preSize, m_size);

    //creating final grid's transformation matrix:
    #pragma omp parallel for //private(xGrid) collapse(2)
    for (int index = 0; index < gridPositions.m_numPts; index++)
    {
      //set viewpoints height above ground (in the normal vector direction above the plane that was found):
      gridPositions.m_pos[index] = gridPositions.m_pos[index] + (in_d_sensor * gridPositions.m_normal[index]);

      //find refFrame matrix:
      feat.CalcRotateMatZaxisToNormal(gridPositions.m_normal[index], m_Orient[preSize + index], gridPositions.m_pos[index]);
    }

    //release memory:
    delete[] gridPositions.m_pos;
    delete[] gridPositions.m_normal;

    return preSize;
  }


  int COrientedGrid::PointCloudAndGridUpdate(const CPtCloud& in_pcl, float in_d_grid, float in_d_sensor)
  {
    CVec3 minBox, maxBox;

    PointCloudUpdate(in_pcl, minBox, maxBox);

    return ViewpointGridUpdate(in_d_grid, in_d_sensor, minBox, maxBox);
  }


  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/
  void COrientedGrid::initMembers()
  {
    m_pclMain.m_numPts = 0;
    m_pclMain.m_pos = NULL;
    m_voxelSize = 0.5;
    m_Orient = NULL;
    m_size = 0;
    m_mainHashed = NULL;
    m_minBBox = CVec3(0, 0, 0);
    m_maxBBox = CVec3(0, 0, 0);
  }



  /******************************************************************************
  *
  *: Class name: CRegDictionary
  *
  ******************************************************************************/

  /** constructor */
  CRegDictionary::CRegDictionary()
  {
    m_descriptors = NULL;
    m_descriptorsDFT = NULL;

    m_r_max = m_r_min = 0;
    m_descWidth = m_descHeight = 1;
  }


  CRegDictionary::CRegDictionary(float in_voxelSize, float in_r_max, float in_r_min, int in_descWidth, int in_descHeight) : COrientedGrid(in_voxelSize)
  {
    m_descriptors = NULL;
    m_descriptorsDFT = NULL;

    m_r_max = in_r_max;
    m_r_min = in_r_min;
    m_descWidth = in_descWidth;
    m_descHeight = in_descHeight;
  }


  CRegDictionary::~CRegDictionary()
  {
    DeleteDescriptors();
  }


  void CRegDictionary::setParameters(float in_r_max, float in_r_min, int in_descWidth, int in_descHeight)
  {
    m_r_max = in_r_max;
    m_r_min = in_r_min;
    m_descWidth = in_descWidth;
    m_descHeight = in_descHeight;
  }


  void CRegDictionary::getParameters(float& out_r_max, float& out_r_min, int& out_descWidth, int& out_descHeight)
  {
    out_r_max = m_r_max;
    out_r_min = m_r_min;
    out_descWidth = m_descWidth;
    out_descHeight = m_descHeight;
  }


  void CRegDictionary::ResetDictionary()
  {
    DeleteDescriptors();

    ResetGrid();
  }


  void CRegDictionary::DeleteDescriptors()
  {
    if (m_descriptors != NULL)
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
    }
  }


  void CRegDictionary::DictionaryUpdate(const CPtCloud& in_pcl, float in_d_grid, float in_d_sensor)
  {
    int preSize = PointCloudAndGridUpdate(in_pcl, in_d_grid, in_d_sensor);

    //create vectors for the descriptors:
    resizeArray(m_descriptors, preSize, m_size);
    resizeArray(m_descriptorsDFT, preSize, m_size);
  }


  void CRegDictionary::PCL2descriptor(const CPtCloud& in_pcl, float* out_RangeImage)
  {
    int totalDescSize = m_descHeight * m_descWidth;
    float azimuthRes = m_descWidth / (2.0f * float(M_PI));
    float elevationRes = m_descHeight / float(M_PI);
    memset(out_RangeImage, 0, totalDescSize * sizeof(float));

    //put closest range per azimuth & elevation in the range image.
    #pragma omp parallel
    {
      float* partialRangeImage = new float[totalDescSize];
      memset(partialRangeImage, 0, totalDescSize * sizeof(float));

      #pragma omp for
      for (int i = 0; i < in_pcl.m_numPts; i++)
      {
        float x = in_pcl.m_pos[i].x;
        float y = in_pcl.m_pos[i].y;
        float z = in_pcl.m_pos[i].z;
        float azimuth = atan2(y, x);
        float elevation = atan2(z, sqrt(x*x + y*y));
        float r = sqrt(x*x + y*y + z*z);

        int azimuthInds = int(floor((azimuth + M_PI) * azimuthRes));
        int elevationInds = m_descHeight - 1 - int(floor((elevation + M_PI / 2) * elevationRes));

        int index = elevationInds*m_descWidth + azimuthInds;

        if ((partialRangeImage[index] == 0) || (r < partialRangeImage[index]))
          partialRangeImage[index] = r;
      }

      #pragma omp critical
      {
        for (int row=0; row<m_descHeight; row++)
          for (int col = 0; col < m_descWidth; col++)
          {
            int index = row*m_descWidth + col;

            if (
              (out_RangeImage[index] == 0)
              ||
              ((partialRangeImage[index] != 0) && (partialRangeImage[index] < out_RangeImage[index]))
              )
              out_RangeImage[index] = partialRangeImage[index];
          }
      }
    }
  }


  void CRegDictionary::Descriptor2DFT(float* in_Descriptor, std::complex<float>* out_DescriptorDFT)
  {
    for (int index2 = 0; index2 < m_descHeight; index2++)
    {
      for (int index1 = 0; index1 < m_descWidth; index1++)
      {
        int index = (index2 * m_descWidth) + index1;
        out_DescriptorDFT[index] = in_Descriptor[index];
      }
    }

    tpcl::DFT2D(unsigned int(m_descWidth), unsigned int(m_descHeight), out_DescriptorDFT);
  }


  std::complex<float>* CRegDictionary::GetEntryDescriptorDFT(int in_entryIndex)
  {
    if (m_descriptors[in_entryIndex] == NULL)
    {
      CPtCloud ptsTran;
      ptsTran.m_pos = new CVec3[m_pclMain.m_numPts];
      ptsTran.m_numPts = 0;

      CMat4 orients = m_Orient[in_entryIndex];
      CVec3 Pos     = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);

      //transform main point cloud to grid point's orientation:
      //#pragma omp parallel for
      float rMinSqr = m_r_min*m_r_min;
      float rMaxSqr = m_r_max*m_r_max;

      for (int Index = 0; Index < m_pclMain.m_numPts; Index++)
      {
        //transform point:
        CVec3 posShifted = m_pclMain.m_pos[Index] - Pos;
        float rQsr = LengthSqr(posShifted);
        if (rQsr < rMinSqr)
          continue;
        if (rQsr > rMaxSqr)
          continue;

        MultiplyVectorRightSide(orients, posShifted, ptsTran.m_pos[ptsTran.m_numPts]);

        ptsTran.m_numPts++;
      }

      //create descriptor - range image, DFT
      int totalDescSize = m_descHeight * m_descWidth;

      m_descriptors[in_entryIndex] = new float[totalDescSize];
      PCL2descriptor(ptsTran, m_descriptors[in_entryIndex]);

      m_descriptorsDFT[in_entryIndex] = new std::complex<float>[totalDescSize];
      Descriptor2DFT(m_descriptors[in_entryIndex], m_descriptorsDFT[in_entryIndex]);

      delete[] ptsTran.m_pos;
    }

    return m_descriptorsDFT[in_entryIndex];
  }


  void CRegDictionary::BestPhaseCorr(std::complex<float>* in_descriptorDFT0, std::complex<float>* in_descriptorDFT1, int& out_bestRow, int& out_bestCol, float& out_bestScore)
  {
    out_bestScore = FLT_MIN;
    out_bestRow = 0;
    out_bestCol = 0;

    std::complex<float>* PhCor = new std::complex<float>[m_descWidth * m_descHeight];

    tpcl::UnitPhaseCorrelation(in_descriptorDFT0, in_descriptorDFT1, PhCor, unsigned int(m_descWidth), unsigned int(m_descHeight));
    tpcl::DFT2D(unsigned int(m_descWidth), unsigned int(m_descHeight), PhCor, false);
    tpcl::DFTshift0ToOrigin(PhCor, unsigned int(m_descWidth), unsigned int(m_descHeight));

    //find max:
    for (int index2 = 0; index2 < m_descHeight; index2++)
    {
      for (int index1 = 0; index1 < m_descWidth; index1++)
      {
        int index = (index2 * m_descWidth) + index1;

        if (PhCor[index].real() > out_bestScore)
        {
          out_bestScore = PhCor[index].real();
          out_bestRow = index2;
          out_bestCol = index1;
        }
      }
    }

    delete[] PhCor;
  }


  int CRegDictionary::SearchDictionary(int in_maxCandidates, float in_searchRadius, std::complex<float>* in_descriptorDFT, int* out_candidates, float* out_grades, CMat4* out_orientations, const CVec3& in_estimatePos)
  {
    int NumOfCandidates = 0;
    int minIndex = 0;
    out_grades[minIndex] = FLT_MAX;

    int* bestCols = new int[in_maxCandidates];

    std::vector<int> inSearchRadius;


    for (int gridIndex = 0; gridIndex < m_size; gridIndex++)
    {
      CMat4 orients = m_Orient[gridIndex];
      CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);

      //ignore entries which are too far from the estimation:
      if (Dist(in_estimatePos, Pos) <= in_searchRadius)
        inSearchRadius.push_back(gridIndex);
    }

    //TODO: see if inSearchRadius.size() == 0 -> Estimated position not in range of any grid point

    #pragma omp parallel for
    for (int inSrIndex = 0; inSrIndex < int(inSearchRadius.size()); inSrIndex++)
      GetEntryDescriptorDFT(inSearchRadius[inSrIndex]);


    int inSrIndex = 0;
    //take first in_maxCandidates entries from dictionary which were considered close enough to estimation:
    for (inSrIndex; inSrIndex < int(inSearchRadius.size()); inSrIndex++)
    {
      int gridIndex = inSearchRadius[inSrIndex];

      std::complex<float>* gridDescDFT = GetEntryDescriptorDFT(gridIndex);
      float bestMax = FLT_MIN;
      int bestRow = -1;
      int bestCol = -1;
      BestPhaseCorr(gridDescDFT, in_descriptorDFT, bestRow, bestCol, bestMax);

      out_candidates[NumOfCandidates] = gridIndex;
      out_grades[NumOfCandidates] = bestMax;
      bestCols[NumOfCandidates] = bestCol;

      if (bestMax < out_grades[minIndex])
      {
        minIndex = NumOfCandidates;
      }

      NumOfCandidates++;

      if (NumOfCandidates == in_maxCandidates)
        break;
    }


    //continue combing dictionary (remaining with best in_maxCandidates Candidates:
    for (inSrIndex; inSrIndex < int(inSearchRadius.size()); inSrIndex++)
    {
      int gridIndex = inSearchRadius[inSrIndex];

      std::complex<float>* gridDescDFT = GetEntryDescriptorDFT(gridIndex);
      float bestMax = FLT_MIN;
      int bestRow = -1;
      int bestCol = -1;
      BestPhaseCorr(gridDescDFT, in_descriptorDFT, bestRow, bestCol, bestMax);
      if (bestMax > out_grades[minIndex])
      {
        out_grades[minIndex] = bestMax;
        out_candidates[minIndex] = gridIndex;
        bestCols[minIndex] = bestCol;


        for (int candIndex = 0; candIndex < NumOfCandidates; candIndex++)
        {
          if (out_grades[candIndex] < out_grades[minIndex])
            minIndex = candIndex;
        }
      }
    }


    //compute rotation matrix for candidates:
    float azimuthRes = 2 * float(M_PI) / m_descWidth;
    float centerShift = ((m_descWidth & 1) == 0) ? 0.5f : 0.0f;

    #pragma omp parallel for
    for (int candIndex = 0; candIndex < NumOfCandidates; candIndex++)
    {
      CMat4 orients = m_Orient[out_candidates[candIndex]];
      CVec3 Pos = CVec3(orients.m[3][0], orients.m[3][1], orients.m[3][2]);

      //calc beast azimuth in radians.
      float peak_azimuth = float(bestCols[candIndex] + centerShift) * azimuthRes - float(M_PI);

      //create corresponding rotation matrix and calc final orientation:
      float OrientVals[] = { cos(peak_azimuth), -sin(peak_azimuth), 0.0f, 0.0f,
        sin(peak_azimuth),  cos(peak_azimuth), 0.0f, 0.0f,
        0.0f             ,  0.0f             , 1.0f, 0.0f,
        0.0f             ,  0.0f             , 0.0f, 1.0f };


      TransposeLeftMultiply(orients, CMat4(OrientVals), out_orientations[candIndex]);
      out_orientations[candIndex].m[3][0] = Pos.x;
      out_orientations[candIndex].m[3][1] = Pos.y;
      out_orientations[candIndex].m[3][2] = Pos.z;
    }

    delete[] bestCols;

    #ifdef DEBUG_LOCAL_RANGE_IMAGE //DEBUG - SearchDictionary print candidates

    int bestIndex = 0;
    float bestGrade = out_grades[0];
    for (int index = 1; index < NumOfCandidates; index++)
    {
      if (out_grades[index] > bestGrade)
      {
        bestIndex = index;
        bestGrade = out_grades[index];
      }
    }
    for (int index = 0; index < NumOfCandidates; index++)
    {
      CRegDebug debugDic("L:\\code\\SLDR\\sldrcr\\Debug");

      char fileName[100] = "Candidate_";
      char indexC[10];  itoa(index, indexC, 10);
      char candC[10];   itoa(out_candidates[index], candC, 10);
      char gradeC[10];  itoa(int(out_grades[index] * 100), gradeC, 10);

      if (index == bestIndex)
      {
        strcat(fileName, indexC); strcat(fileName, "_best_");
      }
      else
      {
        strcat(fileName, indexC); strcat(fileName, "_");
      }
      strcat(fileName, candC);  strcat(fileName, "_");
      strcat(fileName, gradeC); strcat(fileName, ".bmp");

      debugDic.SaveAsBmp(fileName, m_descriptors[out_candidates[index]], m_descWidth, m_descHeight, 2, 60);//in_r_min, in_r_max); //SearchDictionary debug
    }
    #endif


    return NumOfCandidates;
  }


} //namespace SLDR