
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

//#include "registration.h"
#include "RegistDict.h"
#include "OrientDict.h"
#include "features.h"
#include "RegICP.h"
#include "SpatialHash.h"
#include <complex>
#include "tran.h"
#include "common.h"


//#define DEBUG_LOCAL_RANGE_IMAGE
#ifdef DEBUG_LOCAL_RANGE_IMAGE
#include "sldrcr_dbg.h"
#endif


namespace tpcl
{
  struct CRegOptions
  {
    //parameters:
    float m_voxelSizeGlobal;      // if >0, use voxel grid downsampling.
    float m_voxelSizeLocal;       // if >0, use voxel grid downsampling.
    float m_d_grid;               // grid resolution.
    float m_d_sensor;             // dist from ground for grid point.
    int m_lineWidth;              // width of the polar depth map
    int m_numlines;               // height of the polar depth map 
    float m_searchRange;          // range of grid locations to check. if not using GPS this parameter is considered as inf.
    int m_medFiltSize0;           // deniseing median filter size for the range imgae.
    int m_medFiltSize1;           // deniseing median filter size for the thresh image.
    float m_distFromMedianThresh; // max distance between point and median filter's result.
    float m_r_max;                // maximum distance from grid point to be included in the descriptor creation.
    float m_r_min;                // minimum distance from grid point to be included in the descriptor creation.

    CRegOptions() { SetDefaults(); }
    
    void SetDefaults()
    {
      float m_scale = 1.0f;
      m_voxelSizeGlobal = 0.5;
      m_voxelSizeLocal = 0.25;
      m_d_grid = 3;
      m_d_sensor = 2;
      m_lineWidth = 128;
      m_numlines = 64;
      m_searchRange = 50 * m_scale;
      m_medFiltSize0 = 7;
      m_medFiltSize1 = 5;
      m_distFromMedianThresh = 0.03f;
      m_r_max = 60;
      m_r_min = 2;
    }
  };

  /******************************************************************************
  *
  *: Class name: CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::CCoarseRegister()
  {
    m_opts = new CRegOptions;
    CRegOptions* optsP = (CRegOptions*)m_opts;

    m_dictionary = new CRegDictionary(optsP->m_voxelSizeGlobal, optsP->m_r_max, optsP->m_r_min, optsP->m_lineWidth, optsP->m_numlines);
  }


  CCoarseRegister::~CCoarseRegister()
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    delete optsP;

    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;
    delete dictionaryP;
  }


  float CCoarseRegister::RangeNeeded()
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    return (optsP->m_searchRange + optsP->m_r_max);
  }


  void CCoarseRegister::SetMainPtCloud(const CPtCloud& in_pcl, bool in_append)
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;
    Features feat;

    //preprocess main cloud:
    CPtCloud ptsMain; ptsMain.m_type = PCL_TYPE_FUSED; ptsMain.m_color = NULL; ptsMain.m_normal = NULL;
    ptsMain.m_numPts = in_pcl.m_numPts; ptsMain.m_pos = new CVec3[in_pcl.m_numPts];

    feat.DownSample(in_pcl, ptsMain, optsP->m_voxelSizeGlobal);

    if (!in_append)
      dictionaryP->ResetDictionary();

    //create grid's position and orientation: 
    dictionaryP->DictionaryUpdate(ptsMain, optsP->m_d_grid, optsP->m_d_sensor);

    delete[] ptsMain.m_pos;
  }


  void* CCoarseRegister::getMainHashedPtr()
  {
    return ((CRegDictionary*)m_dictionary)->getMainHashedPtr();
  }


  float CCoarseRegister::RegisterCloud(const CPtCloud& in_pcl, CMat4& out_registration, CMat4* in_estimatedOrient)
  {
    const int maxCandidates = 10;
    Features feat;

    float grades[maxCandidates];
    CMat4 candRegistrations[maxCandidates];

    CRegOptions* optsP = (CRegOptions*)m_opts;
    CPtCloud ptsPrePro; ptsPrePro.m_type = PCL_TYPE_SINGLE_ORIGIN; ptsPrePro.m_color = NULL; ptsPrePro.m_normal = NULL;
    ptsPrePro.m_numPts = in_pcl.m_numPts; ptsPrePro.m_pos = new CVec3[in_pcl.m_numPts];

    //preprocess local cloud:
    //if (in_pcl.m_type == PCL_TYPE_SINGLE_ORIGIN_SCAN)
    //{
      //int lineHeight = in_pcl.m_numPts / in_pcl.m_lineWidth;
      feat.DenoiseRange(in_pcl, ptsPrePro, optsP->m_medFiltSize0, optsP->m_distFromMedianThresh);
    //}
    //else
    //{
    //  float res = float(2 * M_PI) / (128 * 5);
    //  Features::DenoiseRangeOfPointCloud(in_pcl, ptsPrePro, optsP->m_medFiltSize0, optsP->m_medFiltSize1, optsP->m_distFromMedianThresh, res);
    //}

    feat.DownSample(ptsPrePro, ptsPrePro, optsP->m_voxelSizeLocal);

    //get registration candidates from dictinary:
    int NumOfCandidates = SecondaryPointCloudRegistrationCandidates(ptsPrePro, maxCandidates, grades, candRegistrations, in_estimatedOrient);

    feat.DownSample(ptsPrePro, ptsPrePro, 2);

    //find final registration:
    float bestGrade = GetRegistrationFromListOfCandidates(NumOfCandidates, ptsPrePro, candRegistrations, out_registration);

    delete[] ptsPrePro.m_pos;

    return bestGrade;
  }
 
  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/
  int CCoarseRegister::SecondaryPointCloudRegistrationCandidates(const CPtCloud& in_pcl, int in_maxCandidates, float* out_grades, CMat4* out_rotations, CMat4* in_estimatedOrient)
  {

    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;
    int* candidates = new int[in_maxCandidates];

    //pick points in range:
    CPtCloud ptsInRange;
    ptsInRange.m_pos = new CVec3[in_pcl.m_numPts];
    int num = 0;
    float rMinSqr = optsP->m_r_min*optsP->m_r_min;
    float rMaxSqr = optsP->m_r_max*optsP->m_r_max;

    for (int Index = 0; Index < in_pcl.m_numPts; Index++)
    {
      float rQsr = LengthSqr(in_pcl.m_pos[Index]);
      if (rQsr < rMinSqr)
        continue;
      if (rQsr > rMaxSqr)
        continue;

      ptsInRange.m_pos[num] = in_pcl.m_pos[Index];
      num++;
    }
    ptsInRange.m_numPts = num;

    //create range image:
    float* descriptor = new float[optsP->m_lineWidth * optsP->m_numlines];
    dictionaryP->PCL2descriptor(ptsInRange, descriptor);

    
    #ifdef DEBUG_LOCAL_RANGE_IMAGE //DEBUG
    CRegDebug debugDic("L:\\code\\SLDR\\sldrcr\\Debug");
    debugDic.SaveAsBmp("local_rangeImage.bmp", descriptor, optsP->m_lineWidth, optsP->m_numlines, optsP->m_r_min, optsP->m_r_max);
    #endif

    //create range image's 2D DFT:
    std::complex<float>* descriptorDFT = new std::complex<float>[optsP->m_lineWidth * optsP->m_numlines];
    dictionaryP->Descriptor2DFT(descriptor, descriptorDFT);

    CVec3 estimatedOrient;
    float searchRange = optsP->m_searchRange;
    if (in_estimatedOrient != NULL)
    {
      estimatedOrient = CVec3(in_estimatedOrient->m[3][0], in_estimatedOrient->m[3][1], in_estimatedOrient->m[3][2]);
    }
    else
    {
      CVec3 dicMinBBox, dicMaxBBox;
      dictionaryP->getBBox(dicMinBBox, dicMaxBBox);
      estimatedOrient = (dicMinBBox + dicMaxBBox) / 2;
      searchRange = Dist2D(dicMinBBox, dicMaxBBox);
    }

    int numOfCandidates = dictionaryP->SearchDictionary(in_maxCandidates, optsP->m_searchRange, descriptorDFT, candidates, out_grades, out_rotations, estimatedOrient);

    delete[] descriptor;
    delete[] descriptorDFT;
    delete[] candidates;

    return numOfCandidates;
  }


  float CCoarseRegister::GetRegistrationFromListOfCandidates(int in_NumOfCandidates, const CPtCloud& in_pcl, CMat4* in_registrations, CMat4& out_registration)
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    Features feat;

    //stay with candidates of minimum RMSE:
    float* CandRMSEs = new float[in_NumOfCandidates];
    #pragma omp parallel for
    for (int cand = 0; cand < in_NumOfCandidates; cand++)
    {
      CandRMSEs[cand] = feat.RMSEofRegistration((CSpatialHash2D*)(getMainHashedPtr()), in_pcl, 4 * optsP->m_voxelSizeGlobal, in_registrations[cand]);
    }

    //TODO: finish RMSE candidate filter
    const int fNumOfCandWanted = 10;// 3; //Final Num Of Candidates wanted
    int finalCandidates[fNumOfCandWanted] = { 0 };
    int fNumOfCand = MinT(fNumOfCandWanted, in_NumOfCandidates);
    for (int fCand = 0; fCand < fNumOfCand; fCand++)
    {
      for (int cand = 0; cand < in_NumOfCandidates; cand++)
      {
        if (CandRMSEs[cand] < CandRMSEs[finalCandidates[fCand]])
          finalCandidates[fCand] = cand;
      }
      CandRMSEs[finalCandidates[fCand]] = FLT_MAX;
    }



    ////select best registration of candidates according to ICP registration:
    ICP icpRegistration(1.5f * optsP->m_voxelSizeGlobal);
    icpRegistration.SetMainPtCloud((CSpatialHash2D*)getMainHashedPtr());
    double bestGrade = DBL_MAX;
    for (int fCand = 0; fCand < fNumOfCand; fCand++)
    {
      int cand = finalCandidates[fCand];
      CMat4 l_icpReg;
      double grade = icpRegistration.RegisterCloud(in_pcl, l_icpReg, in_registrations + cand);
      if (grade < bestGrade)
      {
        bestGrade = grade;
        out_registration = l_icpReg;
      }
    }

    icpRegistration.setRegistrationResolution(0.5f*optsP->m_voxelSizeGlobal);
    icpRegistration.RegisterCloud(in_pcl, out_registration, &out_registration);

    delete[] CandRMSEs;

    return float(bestGrade);
  }

} //namespace SLDR