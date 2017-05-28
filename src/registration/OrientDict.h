
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
*: Package Name: sldrcr_rdi
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_orient_dict_H
#define __tpcl_orient_dict_H



/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#include "../../include/vec.h"
#include "pcl.h"

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
  *: Class name: SLDR_RDI_COrientedGrid
  *
  *: Abstract: creating an oriented 2D grid of 3D normals from a 3D point cloud.
  *            the main cloud is the union of all point cloud added from last grid reset.
  *
  ******************************************************************************/

  class COrientedGrid
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor 
    * @param Xi_voxelSize   the voxel size parameter of the hashed main point cloud. */
    COrientedGrid();
    COrientedGrid(float Xi_voxelSize);

    /** destructor */
    ~COrientedGrid(); 

    /** Get the voxel size parameter of the hashed main point cloud.
    * @return     voxel size. */
    float getVoxelSize();

    /** Get pointer to the hashed main point cloud.
    * @return     pointer to the hashed main point cloud. */
    void* getMainHashedPtr();

    /** Get pointer to the main point cloud.
    * @param Xo_ptsMain   pointer to the main point cloud.
    * @return               size of the main point cloud. */
    void getPclMainPtr(CPtCloud* &Xo_pclMain);

    /** Get pointer to the grid (location and normal orientation per grid point).
    * @param Xo_ptsMain   pointer to the grid.
    * @return             size of the grid. */
    int getGrid(CMat4* &Xo_Orient);

    /** Get BBox of main point cloud.
    * @param Xo_minBBox   minimum bounding box of main point cloud.
    * @param Xo_maxBBox   maximum bounding box of main point cloud. */
    void getBBox(CVec3 &Xo_minBBox, CVec3 &Xo_maxBBox);

    /** delete the grid, releasing all allocated memory !!making the COrientedGrid object unusable!!.*/
    void DeleteGrid();

    /** clearing all data from the grid.*/
    void ResetGrid();

    /** clearing all data from the grid, and updating the voxel size parameter of the hashed main point cloud.
    * @Xi_voxelSize   the voxel size parameter of the hashed main point cloud. */
    void DeleteAndSetVoxelSize(float Xi_voxelSize);

    /** saves the input point cloud, adds it to the hashed main cloud and finds bounding box.
    * @param Xi_pcl           input point cloud.
    * @param Xo_minBox        minimum of the PC's bounding box.
    * @param Xo_maxBox        maximum of the PC's bounding box. */
    void PointCloudUpdate(const CPtCloud& Xi_pcl, CVec3& Xo_minBox, CVec3& Xo_maxBox);

    /** adds grid points (location and normal to ground for that location according to the main cloud) to the grid.
    *   the 2D location of grid points added is derived only from the bouding box and the distance set between the grid points.
    * @param Xi_d_grid        1D distance between (the 2D) grid's points.
    * @param Xi_d_sensor      final location of the grid points is set to be Xi_d_sensor above ground detected from the main cloud.
    * @param Xo_minBox        minimum of bounding box of the area to which we wish to add grid points.
    * @param Xo_maxBox        maximum of bounding box of the area to which we wish to add grid points. 
    * return                  the previous size of the grid. */
    int ViewpointGridUpdate(float Xi_d_grid, float Xi_d_sensor, CVec3& Xi_minBox, CVec3& Xi_maxBox);

    
    /** basically serial call to the functions PointCloudUpdate and ViewpointGridUpdate:
    *   saves the input point cloud, adds it to the hashed main cloud and adds grid points to the grid according to the
    *   bounding box of the input point cloud.
    * @param Xi_pts           input point cloud.
    * @param Xi_d_grid        1D distance between (the 2D) grid's points.
    * @param Xi_d_sensor      final location of the grid points is set to be Xi_d_sensor above ground detected from the main cloud. 
    * return                  the previous size of the grid. */
    int PointCloudAndGridUpdate(const CPtCloud& Xi_pcl, float Xi_d_grid, float Xi_d_sensor);

  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/
    CPtCloud m_pclMain;       ///< main point cloud.
    void* m_mainHashed;         ///< a hashed copy of the original point cloud.

    int m_size;                 ///< number of entries (grid points) in the grid.
    float m_voxelSize;          ///< the voxel size parameter of the hashed main point cloud.
    CMat4* m_Orient;       ///< location and normal orientation per grid point (created from the main point cloud).
    CVec3 m_minBBox;        ///< minimum of boounding box of accumulated main point cloud.
    CVec3 m_maxBBox;        ///< maximum of boounding box of accumulated main point cloud.


    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/
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
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor 
    * @param Xi_voxelSize   the voxel size parameter of the hashed main point cloud.
    * @param Xi_r_max       maximum distance from grid point for descriptor creation.
    * @param Xi_r_min       minimum distance from grid point for descriptor creation.
    * @param Xi_descWidth   descriptor's width.
    * @param Xi_descHeight  descriptor's Height. */
    CRegDictionary();
    CRegDictionary(float Xi_voxelSize, float Xi_r_max, float Xi_r_min, int Xi_descWidth, int Xi_descHeight);

    /** destructor */
    ~CRegDictionary();

    /** set dictionary's parameters.
    * @param Xi_r_max       maximum distance from grid point for descriptor creation.
    * @param Xi_r_min       minimum distance from grid point for descriptor creation.
    * @param Xi_descWidth   descriptor's width.
    * @param Xi_descHeight  descriptor's Height. */
    void setParameters(float Xi_r_max, float Xi_r_min, int Xi_descWidth, int Xi_descHeight);

    /** get dictionary's parameters.
    * @param Xo_r_max       maximum distance from grid point for descriptor creation.
    * @param Xo_r_min       minimum distance from grid point for descriptor creation.
    * @param Xo_descWidth   descriptor's width.
    * @param Xo_descHeight  descriptor's Height. */
    void getParameters(float& Xo_r_max, float& Xo_r_min, int& Xo_descWidth, int& Xo_descHeight);

    /* clear ALL data from dictionary. */
    void ResetDictionary();

    /* deletes all already created entries' descriptors. entries are not deleted, descriptors will need to be calculated again next time they are needed. */
    void DeleteDescriptors();

    /** adds dictionary entries according to the input point cloud:
    *   adds input point cloud to main cloud. 
    *   updates grid points (location and normal to ground for that location according to the main cloud) to the dictionary's grid (entries).
    *   the 2D location of grid points added is derived only from the bouding box of the input PC and the distance set between the grid points.
    *   rest of dictionary's entries' parameters are calculated first time they are needed.
    * @param Xi_pcl           input point cloud.
    * @param Xi_d_grid        1D distance between (the 2D) grid's points.
    * @param Xi_d_sensor      final location of the grid points is set to be Xi_d_sensor above ground detected from the main cloud. */
    void DictionaryUpdate(const CPtCloud& Xi_pcl, float Xi_d_grid, float Xi_d_sensor);


    /** Convert a point cloud into the entry's descriptor - range image, a polar map of distance around the entry's grid point.
    * @param Xi_pcl          point cloud.
    * @param Xo_RangeImage   polar depth map. */
    void PCL2descriptor(const CPtCloud& Xi_pcl, float* Xo_RangeImage);


    /** find DFT of a 2D descriptor.
    * @param Xi_Descriptor      2D descriptor - range image.
    * @param Xo_DescriptorDFT   descriptor's 2D DFT.  */
    void Descriptor2DFT(float* Xi_Descriptor, std::complex<float>* Xo_DescriptorDFT);


    /** gets an entry's DFT descriptor (if it doesn't exist yet, it is made).
    * @param Xi_entryIndex     entry's index to which the descriptor DFT will be returned.  */
    std::complex<float>* GetEntryDescriptorDFT(int Xi_entryIndex);



    /** calculate best phase correlation between 2 descriptors DFTs.
    * @param Xi_descriptorDFT0             input fisrt descriptor DFT.
    * @param Xi_descriptorDFT1             input second descriptor DFT.
    * @param Xo_bestRow                    row of best score.
    * @param Xo_bestCol                    column of best score.
    * @param Xo_bestScore                  grade of the best result. */
    void BestPhaseCorr(std::complex<float>* Xi_descriptorDFT0, std::complex<float>* Xi_descriptorDFT1, int& Xo_bestRow, int& Xo_bestCol, float& Xo_bestScore);


    /** returns candidates most suitable for the input descriptor.
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_searchRadius      range of grid locations to check from estimation. if no input estimation, Xi_searchRadius will be considered as inf.
    * @param Xi_descriptorDFT     input descriptor to search for.
    * @param Xo_candidates        list of the candidates (entries' indices) suitable for the input descriptor.
    * @param Xo_grades            list of the candidates' grades.
    * @param Xo_orientations      list of the candidates' orientations.
    * @param Xi_estimatePos       optional: orientation estimation of the point cloud relevant to the input descriptor. if NULL compare to all dictionary's entries.
    * @return                     number of candidates found.*/
    int SearchDictionary(int Xi_maxCandidates, float Xi_searchRadius, std::complex<float>* Xi_descriptorDFT, int* Xo_candidates, float* Xo_grades, CMat4* Xo_orientations, const CVec3& Xi_estimatePos);


  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/
    float m_r_max, m_r_min;                   // maximum/minimum distance from grid point to be included in the descriptor creation.
    int m_descWidth, m_descHeight;            // width/height of the descriptors/DFTs.

    float** m_descriptors;                    // descriptors array     (per entry).
    std::complex<float>** m_descriptorsDFT;   // descriptors DFT array (per entry).



    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/
    using COrientedGrid::DeleteGrid;
    using COrientedGrid::ResetGrid;
    using COrientedGrid::DeleteAndSetVoxelSize;
    using COrientedGrid::PointCloudUpdate;
    using COrientedGrid::ViewpointGridUpdate;
  };





} // namespace tpcl

#endif