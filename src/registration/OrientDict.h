
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

#ifndef __tpcl_orient_dict_H
#define __tpcl_orient_dict_H

#include "../../include/vec.h"
#include "../include/ptCloud.h"

  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
  struct CVec3;
  struct CMat4;

  namespace std
  {
    template<class _Ty> class complex;
  }

namespace tpcl
{
  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: COrientedGrid
  *
  *: Abstract: creating an oriented 2D grid of 3D normals from a 3D point cloud.
  *            the main cloud is the union of all point cloud added from last grid reset.
  *
  ******************************************************************************/

  class COrientedGrid
  {
  public:
    /** Constructor 
    * @param in_voxelSize   the voxel size parameter of the hashed main point cloud. */
    COrientedGrid();
    COrientedGrid(float in_voxelSize);

    /** destructor */
    ~COrientedGrid(); 

    /** Get the voxel size parameter of the hashed main point cloud.
    * @return     voxel size. */
    float getVoxelSize();

    /** Get pointer to the hashed main point cloud.
    * @return     pointer to the hashed main point cloud. */
    void* getMainHashedPtr();

    /** Get pointer to the main point cloud.
    * @param out_ptsMain   pointer to the main point cloud.
    * @return               size of the main point cloud. */
    void getPclMainPtr(CPtCloud* &out_pclMain);

    /** Get pointer to the grid (location and normal orientation per grid point).
    * @param out_ptsMain   pointer to the grid.
    * @return             size of the grid. */
    int getGrid(CMat4* &out_Orient);

    /** Get BBox of main point cloud.
    * @param out_minBBox   minimum bounding box of main point cloud.
    * @param out_maxBBox   maximum bounding box of main point cloud. */
    void getBBox(CVec3 &out_minBBox, CVec3 &out_maxBBox);

    /** delete the grid, releasing all allocated memory !!making the COrientedGrid object unusable!!.*/
    void DeleteGrid();

    /** clearing all data from the grid.*/
    void ResetGrid();

    /** clearing all data from the grid, and updating the voxel size parameter of the hashed main point cloud.
    * @in_voxelSize   the voxel size parameter of the hashed main point cloud. */
    void DeleteAndSetVoxelSize(float in_voxelSize);

    /** saves the input point cloud, adds it to the hashed main cloud and finds bounding box.
    * @param in_pcl           input point cloud.
    * @param out_minBox        minimum of the PC's bounding box.
    * @param out_maxBox        maximum of the PC's bounding box. */
    void PointCloudUpdate(const CPtCloud& in_pcl, CVec3& out_minBox, CVec3& out_maxBox);

    /** adds grid points (location and normal to ground for that location according to the main cloud) to the grid.
    *   the 2D location of grid points added is derived only from the bouding box and the distance set between the grid points.
    * @param in_d_grid        1D distance between (the 2D) grid's points.
    * @param in_d_sensor      final location of the grid points is set to be in_d_sensor above ground detected from the main cloud.
    * @param out_minBox        minimum of bounding box of the area to which we wish to add grid points.
    * @param out_maxBox        maximum of bounding box of the area to which we wish to add grid points. 
    * return                  the previous size of the grid. */
    int ViewpointGridUpdate(float in_d_grid, float in_d_sensor, CVec3& in_minBox, CVec3& in_maxBox);

    
    /** basically serial call to the functions PointCloudUpdate and ViewpointGridUpdate:
    *   saves the input point cloud, adds it to the hashed main cloud and adds grid points to the grid according to the
    *   bounding box of the input point cloud.
    * @param in_pts           input point cloud.
    * @param in_d_grid        1D distance between (the 2D) grid's points.
    * @param in_d_sensor      final location of the grid points is set to be in_d_sensor above ground detected from the main cloud. 
    * return                  the previous size of the grid. */
    int PointCloudAndGridUpdate(const CPtCloud& in_pcl, float in_d_grid, float in_d_sensor);

  protected:
    CPtCloud m_pclMain;       ///< main point cloud.
    void* m_mainHashed;         ///< a hashed copy of the original point cloud.

    int m_size;                 ///< number of entries (grid points) in the grid.
    float m_voxelSize;          ///< the voxel size parameter of the hashed main point cloud.
    CMat4* m_Orient;       ///< location and normal orientation per grid point (created from the main point cloud).
    CVec3 m_minBBox;        ///< minimum of boounding box of accumulated main point cloud.
    CVec3 m_maxBBox;        ///< maximum of boounding box of accumulated main point cloud.

    /** Set default values to members. */
    void initMembers();

  };





  /******************************************************************************
  *
  *: Class name: SLDR_RDI_CRegDictionary
  *
  *: Abstract: Dictionary who's entries represents a 2D grid, each has it's location, normal according to a main point cloud, 
  *            a range image descriptor as seen from the entry's location&orientation and its corresponding DFT.
  *            the main cloud is the union of all point cloud added to the dictinary from last dictionary reset.
  *            Derived from the thesis of David Avidar (Malaach's group)
  *            see: <FILL PAPER REF HERE>
  *
  ******************************************************************************/

  class CRegDictionary : public COrientedGrid
  {
  public:
    /** Constructor 
    * @param in_voxelSize   the voxel size parameter of the hashed main point cloud.
    * @param in_r_max       maximum distance from grid point for descriptor creation.
    * @param in_r_min       minimum distance from grid point for descriptor creation.
    * @param in_descWidth   descriptor's width.
    * @param in_descHeight  descriptor's Height. */
    CRegDictionary();
    CRegDictionary(float in_voxelSize, float in_r_max, float in_r_min, int in_descWidth, int in_descHeight);

    /** destructor */
    ~CRegDictionary();

    /** set dictionary's parameters.
    * @param in_r_max       maximum distance from grid point for descriptor creation.
    * @param in_r_min       minimum distance from grid point for descriptor creation.
    * @param in_descWidth   descriptor's width.
    * @param in_descHeight  descriptor's Height. */
    void setParameters(float in_r_max, float in_r_min, int in_descWidth, int in_descHeight);

    /** get dictionary's parameters.
    * @param out_r_max       maximum distance from grid point for descriptor creation.
    * @param out_r_min       minimum distance from grid point for descriptor creation.
    * @param out_descWidth   descriptor's width.
    * @param out_descHeight  descriptor's Height. */
    void getParameters(float& out_r_max, float& out_r_min, int& out_descWidth, int& out_descHeight);

    /* clear ALL data from dictionary. */
    void ResetDictionary();

    /* deletes all already created entries' descriptors. entries are not deleted, descriptors will need to be calculated again next time they are needed. */
    void DeleteDescriptors();

    /** adds dictionary entries according to the input point cloud:
    *   adds input point cloud to main cloud. 
    *   updates grid points (location and normal to ground for that location according to the main cloud) to the dictionary's grid (entries).
    *   the 2D location of grid points added is derived only from the bouding box of the input PC and the distance set between the grid points.
    *   rest of dictionary's entries' parameters are calculated first time they are needed.
    * @param in_pcl           input point cloud.
    * @param in_d_grid        1D distance between (the 2D) grid's points.
    * @param in_d_sensor      final location of the grid points is set to be in_d_sensor above ground detected from the main cloud. */
    void DictionaryUpdate(const CPtCloud& in_pcl, float in_d_grid, float in_d_sensor);


    /** Convert a point cloud into the entry's descriptor - range image, a polar map of distance around the entry's grid point.
    * @param in_pcl          point cloud.
    * @param out_RangeImage   polar depth map. */
    void PCL2descriptor(const CPtCloud& in_pcl, float* out_RangeImage);


    /** find DFT of a 2D descriptor.
    * @param in_Descriptor      2D descriptor - range image.
    * @param out_DescriptorDFT   descriptor's 2D DFT.  */
    void Descriptor2DFT(float* in_Descriptor, std::complex<float>* out_DescriptorDFT);


    /** gets an entry's DFT descriptor (if it doesn't exist yet, it is made).
    * @param in_entryIndex     entry's index to which the descriptor DFT will be returned.  */
    std::complex<float>* GetEntryDescriptorDFT(int in_entryIndex);



    /** calculate best phase correlation between 2 descriptors DFTs.
    * @param in_descriptorDFT0             input fisrt descriptor DFT.
    * @param in_descriptorDFT1             input second descriptor DFT.
    * @param out_bestRow                    row of best score.
    * @param out_bestCol                    column of best score.
    * @param out_bestScore                  grade of the best result. */
    void BestPhaseCorr(std::complex<float>* in_descriptorDFT0, std::complex<float>* in_descriptorDFT1, int& out_bestRow, int& out_bestCol, float& out_bestScore);


    /** returns candidates most suitable for the input descriptor.
    * @param in_maxCandidates     maximum number of matches to return.
    * @param in_searchRadius      range of grid locations to check from estimation. if no input estimation, in_searchRadius will be considered as inf.
    * @param in_descriptorDFT     input descriptor to search for.
    * @param out_candidates        list of the candidates (entries' indices) suitable for the input descriptor.
    * @param out_grades            list of the candidates' grades.
    * @param out_orientations      list of the candidates' orientations.
    * @param in_estimatePos       optional: orientation estimation of the point cloud relevant to the input descriptor. if NULL compare to all dictionary's entries.
    * @return                     number of candidates found.*/
    int SearchDictionary(int in_maxCandidates, float in_searchRadius, std::complex<float>* in_descriptorDFT, int* out_candidates, float* out_grades, CMat4* out_orientations, const CVec3& in_estimatePos);


  protected:
    float m_r_max, m_r_min;                   // maximum/minimum distance from grid point to be included in the descriptor creation.
    int m_descWidth, m_descHeight;            // width/height of the descriptors/DFTs.

    float** m_descriptors;                    // descriptors array     (per entry).
    std::complex<float>** m_descriptorsDFT;   // descriptors DFT array (per entry).

    using COrientedGrid::DeleteGrid;
    using COrientedGrid::ResetGrid;
    using COrientedGrid::DeleteAndSetVoxelSize;
    using COrientedGrid::PointCloudUpdate;
    using COrientedGrid::ViewpointGridUpdate;
  };





} // namespace tpcl

#endif