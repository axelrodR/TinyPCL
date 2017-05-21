
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

#ifndef __tpcl_register_dict_H
#define __tpcl_register_dict_H


#include "../include/registration.h"


  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
  struct CVec3;
  struct CMat4;

namespace tpcl
{
  class  CRegDictionary;         // registration Dictionary


  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: SLDR_SP_CCoarseRegister
  *
  *: Abstract: Coarse registeration based on dictionary
  *            Derived from the thesis of David Avidar (Malaach's group)
  *            see: <FILL PAPER REF HERE> 
  *
  ******************************************************************************/

  class CCoarseRegister : public IRegister
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor */
    CCoarseRegister();

    /** destructor */
    virtual ~CCoarseRegister(); 


    /** Get the range of main point cloud needed for dictionary registration, from PGS estimation.
    /   this is equal to the search range + maximum range taken for entry's descriptor creation.
    * @return     range needed.*/
    float RangeNeeded();

    /** update the main point cloud in which registration is searched for, for a secondary point cloud.
    *   adds the dictionary's entries created from the input point cloud.
    * @param Xi_pts           point cloud to add to existing main point cloud.
    * @param Xi_clean         if true deletes all previous information of main point cloud. otherwise adds the new point cloud to the previous. */
    void MainPointCloudUpdate(int Xi_numPts, const CVec3* Xi_pts, bool Xi_clean = false);

    /** Get hashed main point cloud.
    * @param return         pointer to hashed main point cloud. */
    void* getMainHashedPtr();

    /** Get best dictionary registration for a secondary point cloud, with ICP registration refinements.
    * @param Xo_registration      best registration found.
    * @param Xi_pts               secondary point cloud - arranged by azimuth and latitude. i.e. row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
    * @param Xi_numPts            total points in the secondary point cloud.
    * @param Xi_lineWidth         if it's an order point cloud then this is the width of the local point cloud.
    * @param Xi_estimatedOrient   estimation of registration location (only uses location vector). if NULL then compare to all dictionary's entries.
    * return                      registration's score (ICP based) - the lower the better.
    * !!!!currently not supporting unordered point cloud!!!! */
    float SecondaryPointCloudRegistration(CMat4& Xo_registration, CVec3* Xi_pts, int Xi_numPts, int Xi_lineWidth = -1, CMat4* Xi_estimatedOrient = 0);



  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/
    CRegDictionary* m_dictionary;   ///< internal data used to store features of main point cloud.
    void* m_opts;                   ///< implementation specific options.

    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/
    /** Get list of top registration matches for a secondary point cloud and thier corresponding grades of the matches.
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_pts               secondary point cloud.
    * @param Xo_grades            list of the candidates' grades - the higher the better.
    * @param Xo_rotations         list of the candidates' registration results.
    * @param Xi_numPts            total points in the secondary point cloud.
    * @param Xi_estimatedOrient   estimation of registration location (only uses location vector). if NULL then compare to all dictionary's entries.
    * @return                     number of candidates found.*/
    int SecondaryPointCloudRegistrationCandidates(int Xi_maxCandidates, CVec3* Xi_pts, float* Xo_grades, CMat4* Xo_rotations, int Xi_numPts, CMat4* Xi_estimatedOrient = NULL);


    /** returns final registration from list of candidates, using RMSE to reduce list of candidates and ICP for final selection.
    * @param Xi_NumOfCandidates   number of candidates.
    * @param Xi_numPts            total points in the secondary point cloud.
    * @param Xi_pts               secondary point cloud.
    * @param Xi_registrations     list of the candidates' registration.
    * @param Xo_registration      best registration from candidates.
    * @return                     grade of final registration (ICP grade - the lower the better). */
    float GetRegistrationFromListOfCandidates(int Xi_NumOfCandidates, int Xi_numPts, CVec3* Xi_pts, CMat4* Xi_registrations, CMat4& Xo_registration);

  };



} // namespace SLDR

#endif
