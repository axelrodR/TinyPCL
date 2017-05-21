
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

/******************************************************************************
*
*: Package Name: sldrcr_sp
*
******************************************************************************/
#include <D3dx9core.h> // uncommenet if using DirectX
#include <ifr/ifrgen/ifrgen_stnd.h>
#include "sldrcr_sp.h"
#include "sldrcr_rdi.h"
#include "sldrcr_ftr.h"
#include "sldrcr_icp.h"
#include <gen/gengmtrx/gengmtrx_spat.h>
#include <gen/gengmtrx/gengmtrx_vec.h>
#include <vector>
#include <complex>
#include <IFR/ifrmath/ifrmath_tran.h>
#include <algorithm>
//#include <IFR\ifrlog\ifrlog.h>
#include <IFR\ifrlog\ifrlog_prfl.h>
#include <omp.h>

//#define DEBUG_LOCAL_RANGE_IMAGE
#ifdef DEBUG_LOCAL_RANGE_IMAGE
#include "sldrcr_dbg.h"
#endif


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


using namespace GenGmtrx;

namespace SLDR
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
    float m_searchRange;          // range of grid locations to check. if not using GPS this parameter is considered as inf.
    int m_medFiltSize0;           // deniseing median filter size for the range imgae.
    int m_medFiltSize1;           // deniseing median filter size for the thresh image.
    float m_distFromMedianThresh; // max distance between point and median filter's result.
    float m_r_max;                // maximum distance from grid point to be included in the descriptor creation.
    float m_r_min;                // minimum distance from grid point to be included in the descriptor creation.

    CRegOptions() { SetDefaults(); }

    void SetDefaults();
  };



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
  *: Method name: SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::CCoarseRegister()
    {
    m_opts = new CRegOptions;
    CRegOptions* optsP = (CRegOptions*)m_opts;

    m_dictionary = new CRegDictionary(optsP->m_voxelSizeGlobal, optsP->m_r_max, optsP->m_r_min, optsP->m_lineWidth, optsP->m_numlines);
    }

  /******************************************************************************
  *
  *: Method name: ~SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::~CCoarseRegister()
      {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    delete optsP;

    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;
    delete dictionaryP;
    }



  /******************************************************************************
  *
  *: Method name: SearchRange
  *
  ******************************************************************************/
  float CCoarseRegister::RangeNeeded()
        {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    return (optsP->m_searchRange + optsP->m_r_max);
      }



  /******************************************************************************
  *
  *: Method name: MainPointCloudUpdate
  *
  ******************************************************************************/
  void CCoarseRegister::MainPointCloudUpdate(int Xi_numPts, const D3DXVECTOR3* Xi_pts, bool Xi_clean)
        {
    PROFILE("MainPointCloudUpdate");

    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;

    //preprocess main cloud:
    D3DXVECTOR3* ptsGlobal = new D3DXVECTOR3[Xi_numPts];
    int numPtsGlobal = Xi_numPts;

    Features::DownSamplePointCloud(optsP->m_voxelSizeGlobal, Xi_numPts, Xi_pts, numPtsGlobal, ptsGlobal);

    if (Xi_clean)
      dictionaryP->ResetDictionary();

    //create grid's position and orientation: 
    dictionaryP->DictionaryUpdate(numPtsGlobal, ptsGlobal, optsP->m_d_grid, optsP->m_d_sensor);

    delete[] ptsGlobal;
    }


  /******************************************************************************
  *
  *: Method name: getMainHashedPtr
  *
  ******************************************************************************/
  void* CCoarseRegister::getMainHashedPtr()
        {
    return ((CRegDictionary*)m_dictionary)->getMainHashedPtr();
    }


  /******************************************************************************
  *
  *: Method name: SecondaryPointCloudRegistration
  *
  ******************************************************************************/
  float CCoarseRegister::SecondaryPointCloudRegistration(D3DXMATRIX& Xo_registration, D3DXVECTOR3* Xi_pts, int Xi_numPts, int Xi_lineWidth, D3DXMATRIX* Xi_estimatedOrient)
    {
    PROFILE("SecondaryPointCloudRegistration");
    const int maxCandidates = 10;
      
    float grades[maxCandidates];
    D3DXMATRIX candRegistrations[maxCandidates];

    CRegOptions* optsP = (CRegOptions*)m_opts;
    D3DXVECTOR3* ptsPrePro = new D3DXVECTOR3[Xi_numPts];

    int totalPixels;

    //preprocess local cloud:
    if (Xi_lineWidth > 0)
        {
      int lineHeight = Xi_numPts / Xi_lineWidth;
      totalPixels = Features::DenoiseRangeOfOrderedPointCloud(Xi_lineWidth, lineHeight, optsP->m_medFiltSize0, optsP->m_medFiltSize1, optsP->m_distFromMedianThresh, Xi_pts, ptsPrePro);
    }
    else
        {
      float res = float(2 * M_PI) / (128 * 5);
      totalPixels = Features::DenoiseRangeOfPointCloud(res, optsP->m_medFiltSize0, optsP->m_medFiltSize1, optsP->m_distFromMedianThresh, Xi_numPts, Xi_pts, ptsPrePro);
      }

    Features::DownSamplePointCloud(optsP->m_voxelSizeLocal, totalPixels, ptsPrePro, totalPixels, ptsPrePro);


    //get registration candidates from dictinary:
    int NumOfCandidates = SecondaryPointCloudRegistrationCandidates(maxCandidates, ptsPrePro, grades, candRegistrations, totalPixels, Xi_estimatedOrient);

    //find final registration:
    float bestGrade = GetRegistrationFromListOfCandidates(NumOfCandidates, totalPixels, ptsPrePro, candRegistrations, Xo_registration);

    delete[] ptsPrePro;

    return bestGrade;
      }

  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/


  /******************************************************************************
  *
  *: Method name: SecondaryPointCloudRegistrationCandidates
  *
  ******************************************************************************/
  int CCoarseRegister::SecondaryPointCloudRegistrationCandidates(int Xi_maxCandidates, D3DXVECTOR3* Xi_pts, float* Xo_grades, D3DXMATRIX* Xo_rotations, int Xi_numPts, D3DXMATRIX* Xi_estimatedOrient)
        {
    PROFILE("SecondaryPointCloudRegistrationCandidates");

    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;
    int* candidates = new int[Xi_maxCandidates];

    //create range image:
    float* descriptor = new float[optsP->m_lineWidth * optsP->m_numlines];
    dictionaryP->PCL2descriptor(Xi_numPts, Xi_pts, descriptor);


    #ifdef DEBUG_LOCAL_RANGE_IMAGE //DEBUG
    CRegDebug debugDic("L:\\code\\SLDR\\sldrcr\\Debug");
    debugDic.SaveAsBmp("local_rangeImage.bmp", descriptor, optsP->m_lineWidth, optsP->m_numlines, optsP->m_r_min, optsP->m_r_max);
    #endif

    //create range image's 2D DFT:
    std::complex<float>* descriptorDFT = new std::complex<float>[optsP->m_lineWidth * optsP->m_numlines];
    dictionaryP->Descriptor2DFT(descriptor, descriptorDFT);

    D3DXVECTOR3 estimatedOrient;
    float searchRange = optsP->m_searchRange;
    if (Xi_estimatedOrient != NULL)
        {
      estimatedOrient = D3DXVECTOR3(Xi_estimatedOrient->m[3][0], Xi_estimatedOrient->m[3][1], Xi_estimatedOrient->m[3][2]);
        }
    else
          {
      D3DXVECTOR3 dicMinBBox, dicMaxBBox;
      dictionaryP->getBBox(dicMinBBox, dicMaxBBox);
      estimatedOrient = (dicMinBBox + dicMaxBBox) / 2;
      searchRange = Dist2D(dicMinBBox, dicMaxBBox);
      }

    int NumOfCandidates = dictionaryP->SearchDictionary(Xi_maxCandidates, optsP->m_searchRange, descriptorDFT, candidates, Xo_grades, Xo_rotations, estimatedOrient);

    delete[] descriptor;
    delete[] descriptorDFT;
    delete[] candidates;

      return NumOfCandidates;
    }



  /******************************************************************************
  *
  *: Method name: GetRegistrationFromListOfCandidates
  *
  ******************************************************************************/
  float CCoarseRegister::GetRegistrationFromListOfCandidates(int Xi_NumOfCandidates, int Xi_numPts, D3DXVECTOR3* Xi_pts, D3DXMATRIX* Xi_registrations, D3DXMATRIX& Xo_registration)
  {
    PROFILE("GetRegistrationFromListOfCandidates");
    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionary* dictionaryP = (CRegDictionary*)m_dictionary;


    //stay with candidates of minimum RMSE:
    float* CandRMSEs = new float[Xi_NumOfCandidates];
    #pragma omp parallel for
    for (int cand = 0; cand < Xi_NumOfCandidates; cand++)
    {
      CandRMSEs[cand] = Features::RMSEofRegistration(optsP->m_voxelSizeGlobal, (CSpatialHash2D*)(getMainHashedPtr()), Xi_numPts, Xi_pts, Xi_registrations[cand]);
  }

    const int fNumOfCandWanted = 3; //Final Num Of Candidates wanted
    int finalCandidates[fNumOfCandWanted] = { 0 };
    int fNumOfCand = min(fNumOfCandWanted, Xi_NumOfCandidates);
    for (int fCand = 0; fCand < fNumOfCand; fCand++)
  {
      for (int cand = 0; cand < Xi_NumOfCandidates; cand++)
    {
        if (CandRMSEs[cand] < CandRMSEs[finalCandidates[fCand]])
          finalCandidates[fCand] = cand;
    }
      CandRMSEs[finalCandidates[fCand]] = FLT_MAX;
    }


    
    ////select best registration of candidates according to ICP registration:
    ICP icpRegistration(2.5 * optsP->m_voxelSizeGlobal);
    icpRegistration.MainPointCloudUpdate(getMainHashedPtr());
    double bestGrade = DBL_MAX;
    for (int fCand = 0; fCand < fNumOfCand; fCand++)
    {
      int cand = finalCandidates[fCand];
      D3DXMATRIX l_icpReg;
      double grade = icpRegistration.SecondaryPointCloudRegistration(l_icpReg, Xi_pts, Xi_numPts, 1, Xi_registrations + cand);
      if (grade < bestGrade)
      {
        bestGrade = grade;
        Xo_registration = l_icpReg;
      }
    }

    icpRegistration.setRegistrationThresh(optsP->m_voxelSizeGlobal);
    icpRegistration.SecondaryPointCloudRegistration(Xo_registration, Xi_pts, Xi_numPts, 1, &Xo_registration);

    delete[] CandRMSEs;
 
    return bestGrade;
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


  void CRegOptions::SetDefaults()
  {
    m_voxelSizeGlobal = 2;// 0.5;
    m_voxelSizeLocal = 2;// 0.25;
    m_d_grid = 3;                   
    m_d_sensor = 2;                 
    m_lineWidth = 128;              
    m_numlines = 64;                
    m_searchRange = 30;             
    m_medFiltSize0 = 7;             
    m_medFiltSize1 = 5;             
    m_distFromMedianThresh = 0.03f; 
    m_r_max = 60;                   
    m_r_min = 2;                    
  }


} //namespace SLDR