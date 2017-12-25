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


#ifndef __tpcl_ftr_H
#define __tpcl_ftr_H


#include "../../include/iFeatures.h"
#include "../../include/vec.h"

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
  *: Class name: Features
  *
  *: Abstract: a bundle of point cloud processing functions.
  *
  ******************************************************************************/

  /** basic features calculations on point cloud. */
  class Features : public IFeatures
  {
  public:
    /** Constructor */
    Features();

    /** destructor */
    virtual ~Features();


    /** finds the normalized normals for the input points from the global hashed point cloud.
    *   optional: fix z value of input points to z in the point's xy from the estimated plane around that point.
    *   also, updates the z value of the points to be on the found plane.
    * @param io_pcl            input: points where to look for the normal. output: filled with the normals. pos updated to be on plane with the same x,y.
    * @param in_radius         radius around input point to get global points for the plane estimation (from which the also the normals are calculated).
    * @param in_pclHash        input hashed global point cloud. if NULL input points are hashed and considered the global point cloud.
    * @param in_fixZ           if true z value of input points are fixed according to the plane estimated around them. if false output points = input points. */
    virtual void FillNormals(CPtCloud& io_pcl, float in_radius =10.0f, CSpatialHash2D* in_pclHash = 0, bool in_fixZ = false);


    /** denose by range an xyz image and return a denoised point cloud.
    * @param in_pcl                      input point cloud of type PCL_TYPE_SINGLE_ORIGIN_SCAN. (assuming row_i > row_j -> latitude_i > latitude_j. col_i > col_j -> azimuth_i > azimuth_j).
    * @param out_pcl                      denoised point cloud. assumes buffer size at least as in_pts's size.
    * @param in_windowSize             median filter size for the range imgae.
    * @param in_noiseTh     max distance between point and median filter's result. */
    virtual void DenoiseRange(const CPtCloud& in_pcl, CPtCloud& out_pcl, int in_windowSize, float in_noiseTh);


    /** TODO: !!empty function!!
    denose by range a point cloud. */
    virtual void DenoiseRangeOfPointCloud();


    /** downsample a point cloud. Divides to grid from minXYZ (of pts) to max XYZ, of size m_voxelSize.
    *  If more than one point in same grid index, takes first one encountered.
    *  supports in_pts = out_pts. assumes size of out_pts >= in_numPts.
    * @param in_pcl                 input point cloud.
    * @param out_pcl                 downsampled point cloud.
    * @param in_voxelSize           size of a voxel in grid. assums bigger than 0. */
    virtual void DownSample(const CPtCloud& in_pcl, CPtCloud& out_pcl, float in_voxelSize);


    /** calculates the RMSE of a registration.
    *   !! points further away than max2DRadius will be considered as max2DRadius*sqrt(1.5) away.
    * @param in_max2DRadius  maximum 2D radius to looks for matches of projected pcl1 in pcl2.
    * @param in_pcl1         main hashed point cloud.
    * @param in_pcl2size     2nd point cloud size.
    * @param in_pcl2         2nd point cloud.
    * @param in_Rt           registration from 2nd point cloud to main. 
    * return                 RMSE of registration. */
    virtual float RMSEofRegistration(CSpatialHash2D* in_pcl1, const CPtCloud& in_pcl2, float in_max2DRadius, const CMat4& in_Rt);


    /** finds the rotation matrix so that the new z axis will be in the normal direction. to be used: x_rotated = R * x.
    * @param Xi_Normal          input normal.
    * @param Xo_RotateMat       output rotation matrix.
    * @param Xi_Pos             optional: input vector to replace the default vector in output CMat4, which is (0, 0, 0).*/
    virtual void CalcRotateMatZaxisToNormal(const CVec3& Xi_Normal, CMat4& Xo_RotateMat, const CVec3& Xi_Pos = CVec3(0, 0, 0));
  };
 

} // namespace tpcl

#endif
