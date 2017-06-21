
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

#ifndef __tpcl_register_dict_H
#define __tpcl_register_dict_H


#include "../include/registration.h"
#include "../include/vec.h"

namespace tpcl
{
  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
  class  CRegDictionary;         // registration Dictionary


  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: CCoarseRegister
  *
  *: Abstract: Coarse registeration based on dictionary, with ICP registration refinements.
  *            Derived from the thesis of David Avidar (Malaach's group)
  *            see: <FILL PAPER REF HERE> 
  *
  ******************************************************************************/

  class CCoarseRegister : public IRegister
  {
  public:
    /** Constructor */
    CCoarseRegister();

    /** destructor */
    virtual ~CCoarseRegister(); 


    /** Get the range of main point cloud needed for dictionary registration, from PGS estimation.
    /   this is equal to the search range + maximum range taken for entry's descriptor creation.
    * @return     range needed.*/
    float RangeNeeded();

    /** Set main cloud point.
    * Registration of secondary cloud points are done against this cloud using RegisterCloud()
    * @param in_pcl           point cloud.
    * @param in_append        if true, append points to the existing cloud
    */
    void SetMainPtCloud(const CPtCloud& in_pcl, bool in_append = false);

    /** Get hashed main point cloud.
    * @param return         pointer to hashed main point cloud. */
    void* getMainHashedPtr();

    /** Get registration for a secondary point cloud against the main cloud
    * The second cloud is not stored
    * //TODO/// !!!!currently not supporting unordered point cloud!!!!
    * @param in_pcl               secondary point cloud. //TODO/// - arranged by azimuth and latitude. i.e. row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
    * @param out_registration      best registration found.
    * @param in_estimatedOrient   estimation of registration, if 0 then estimation is identity.
    * @return                     registration's grade/error - the lower the better.
    */
    float RegisterCloud(const CPtCloud& in_pcl, CMat4& out_registration, CMat4* in_estimatedOrient = 0);



  protected:
    CRegDictionary* m_dictionary;   ///< internal data used to store features of main point cloud.
    void* m_opts;                   ///< implementation specific options.


    /** Get list of top registration matches for a secondary point cloud and thier corresponding grades of the matches.
    * @param in_maxCandidates     maximum number of matches to return.
    * @param in_pts               secondary point cloud.
    * @param out_grades            list of the candidates' grades - the higher the better.
    * @param out_rotations         list of the candidates' registration results.
    * @param in_numPts            total points in the secondary point cloud.
    * @param in_estimatedOrient   estimation of registration location (only uses location vector). if NULL then compare to all dictionary's entries.
    * @return                     number of candidates found.*/
    int SecondaryPointCloudRegistrationCandidates(const CPtCloud& in_pcl, int in_maxCandidates, float* out_grades, CMat4* out_rotations, CMat4* in_estimatedOrient = NULL);


    /** returns final registration from list of candidates, using RMSE to reduce list of candidates and ICP for final selection.
    * @param in_NumOfCandidates   number of candidates.
    * @param in_numPts            total points in the secondary point cloud.
    * @param in_pts               secondary point cloud.
    * @param in_registrations     list of the candidates' registration.
    * @param out_registration      best registration from candidates.
    * @return                     grade of final registration (ICP grade - the lower the better). */
    float GetRegistrationFromListOfCandidates(int in_NumOfCandidates, const CPtCloud& in_pcl, CMat4* in_registrations, CMat4& out_registration);
  };



} // namespace SLDR

#endif
