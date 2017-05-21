// File Location: 

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
*: Title:
*
******************************************************************************/

#ifndef __sldrcr_sp_H
#define __sldrcr_sp_H

#include "vec.h"


/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/


  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

namespace tpcl
{
  class  CRegDictionary;         // registration Dictionary

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: TPCL_CCoarseRegister
  *
  *: Abstract: Coarse registeration vased on dictionary
  *            Derived from the thesis of David Avidar (Malaach's group)
  *            see: <FILL PAPER REF HERE> 
  *
  ******************************************************************************/

  class IRegister
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor */
    IRegister();

    /** destructor */
    virtual ~IRegister();


    /** Build dictionary from global point cloud. Xi_pts IS being altered.
    * @param Xi_pts           point cloud */
    virtual void BuildDictionary(int Xi_numPts, const CVec3* Xi_pts);


    /** Create and get list of registration matches for a local point cloud
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_pts               local point cloud - arranged by azimuth and latitude. i.e. row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
    * @param Xi_originApprox      measured/approximate origin.
    * @param Xo_candidates        list of the candidates (indices of the candidates in the dictionary).
    * @param Xo_grades            list of the candidates grades.
    * @param Xo_rotations         list of the candidates rotation.
    * @param Xi_numPts            total points in the point cloud.
    * @param Xi_numlines          if it's an order point cloud then this is the height of the local point cloud.    
    * @param Xi_GPS               GPS estimation of registration location, if NULL then no estimation.
    * @return                     number of candidates.*/
    virtual int GetLocalRegistrationCandidates(int Xi_maxCandidates, CVec3* Xi_pts, CVec3 Xi_originApprox, int* Xo_candidates, float* Xo_grades, CMat4* Xo_rotations, int Xi_numPts, int Xi_numlines = -1, CVec3* Xi_GPS = NULL);


    /** Get best registration match from created candidates list.
    * @param Xi_NumOfCandidates   number of candidates.
    * @param Xi_candidates        list of the candidates (indices of the candidates in the dictionary).
    * @param Xi_grades            list of the candidates grades.
    * @param Xi_rotations         list of the candidates rotation.
    * @param Xo_best              best registration from candidates */
    virtual void GetLocalRegistration(int Xi_NumOfCandidates, int* Xi_candidates, float* Xi_grades, CMat4* Xi_rotations, CMat4& Xo_best);



  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/

    int m_numPtsGlobal;
    CVec3* m_ptsGlobal;       ///< a copyu of the original point cloud
    CRegDictionary* m_dictionary;   ///< internal data used to store features of global cloud
    void* m_opts;                   ///< implementation specific options

    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/

  };




} // namespace tpcl

#endif