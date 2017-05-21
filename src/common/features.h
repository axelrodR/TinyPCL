// File Location: 

// File Location: S:\gen\gengmtrx\gengmtrx_grid.h

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
*: Package Name: features
*
*: Title:
*
******************************************************************************/

#ifndef __sldrcr_ftr_H
#define __sldrcr_ftr_H

#include "SpatialHash.h"
#include "../../include/vec.h"

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


  /** basic features calculations on point cloud. */
  class Features
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor */
    Features();

    /** destructor */
    virtual ~Features();



    /** TODO: !!empty function!!
    *fills a point cloud.
    * @param Xi_pts           input point cloud.
    * @param Xo_pts           filled point cloud. */
    static void FillPointCloud(int Xi_numPts, const CVec3* Xi_pts, int Xo_numPts, const CVec3* Xo_pts);


    /** finds the normalized normals for the input points from the global hashed point cloud.
    *   optional: fix z value of input points to z in the point's xy from the estimated plane around that point.
    *   also, updates the z value of the points to be on the found plane.
    * @param Xio_pts            input: points where to look for the normal. output: point on the plane with the same x,y.
    * @param Xo_Normals         normals per input point.
    * @param Xi_radius          radius around input point to get global points for the plane estimation (from which the also the normals are calculated).
    * @param Xi_globalHashed    input hashed global point cloud. if NULL input points are hashed and considered the global point cloud.
    * @param Xi_fixZ            if true z value of input points are fixed according to the plane estimated around them. if false output points = input points. */
    static void FindNormal(int Xi_numPts, CVec3* Xio_pts, CVec3* Xo_Normals, float Xi_radius=10.0f, CSpatialHash2D* Xi_globalHashed = NULL, bool Xi_fixZ = false);


    /** denose by range an xyz image and return a denoised point cloud.
    * @param Xi_lineWidth                width of the polar depth map.
    * @param Xi_numlines                 height of the polar depth map
    * @param Xi_medFiltSize0             median filter size for the range imgae.
    * @param Xi_medFiltSize1             median filter size for the thresh image.
    * @param Xi_distFromMedianThresh     max distance between point and median filter's result.
    * @param Xi_pts                      input xyz image. assuming row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
    * @param Xo_ptsDenoised              denoised point cloud. assumes buffer size at least as Xi_pts's size.
    * @return                            size of point cloud.  */
    static int DenoiseRangeOfOrderedPointCloud(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, CVec3* Xi_pts, CVec3* Xo_ptsDenoised);


    /** TODO: !!empty function!!
    denose by range a point cloud.
    * @param Xi_lineWidth                width of the polar depth map.
    * @param Xi_numlines                 height of the polar depth map
    * @param Xi_medFiltSize0             median filter size for the range imgae.
    * @param Xi_medFiltSize1             median filter size for the thresh image.
    * @param Xi_distFromMedianThresh     max distance between point and median filter's result.
    * @param Xi_pts                      input xyz image. assuming row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j.
    * @param Xo_ptsDenoised              denoised point cloud. assumes size at least as Xi_pts size.
    * @return                            size of point cloud.  */
    static int DenoiseRangeOfPointCloud(float Xi_res, int Xi_medFiltSize0, int Xi_medFiltSize1, float Xi_distFromMedianThresh, int Xi_numPts, CVec3* Xi_pts, D3DXVECTOR3* Xo_ptsDenoised);


    /** downsample a point cloud. Divides to grid from minXYZ (of pts) to max XYZ, of size m_voxelSize.
    *  If more than one point in same grid index, takes first one encountered.
    *  supports Xi_pts = Xo_pts. assumes size of Xo_pts >= Xi_numPts.
    * @param Xi_voxelSize           size of a voxel in grid. assums bigger than 0.
    * @param Xi_pts                 point cloud.
    * @param Xo_pts                 downsampled point cloud. */
    static void DownSamplePointCloud(float Xi_voxelSize, int Xi_numPts, const CVec3* Xi_pts, int& Xo_numPts, CVec3* Xo_pts);


    /** calculates the RMSE of a registration.
    *   !! points further away than 2*hashRes will be considered as 2*hashRes*sqrt(1.5) away.
    * @param Xi_pcl1        main hashed point cloud.
    * @param Xi_pcl2size    2nd point cloud size.
    * @param Xi_pcl2        2nd point cloud.
    * @param Xi_Rt         registration from 2nd point cloud to main. 
    * return                RMSE of registration. */
    static float RMSEofRegistration(float hashRes, CSpatialHash2D* Xi_pcl1, int Xi_pcl2size, CVec3* Xi_pcl2, CMat4& Xi_Rt);

  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/

    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/

  };
 

} // namespace tpcl

#endif
