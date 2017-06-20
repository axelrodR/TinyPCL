// File Location: 

// File Location: S:\gen\gengmtrx\gengmtrx_grid.h

//
// Copyright (c) 2016-2017 Geosim Ltd.
// 
// Written by Amit Henig 
//            Ramon Axelrod
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



#ifndef __tpcl_iftr_H
#define __tpcl_iftr_H


namespace tpcl
{

  class  CRegDictionary;         // registration Dictionary
  class  CSpatialHash2D;
  struct CPtCloud;


  /** Basic features calculations on point cloud. */
  class IFeatures
  {
  public:

    /** Estiamtes the normals for the input points
     *   optional: attfix z value of input points to z in the point's xy from the estimated plane around that point.
     *   also, updates the z value of the points to be on the found plane.
     * @param io_pcl             point cloud to which normals are calculated
     * @param in_radius          radius around each position to use for normal estimation
     * @param in_pclHash         optional spatial hashing of point cloud (used for NN search).
     *                           can include more points than io_pcl
     * @param in_fixZ            true: points are attached to the tangent plane
     */
    virtual void FillNormals(CPtCloud& io_pcl, float in_radius=10.0f,
                             CSpatialHash2D* in_pclHash = 0, bool in_fixZ = false) = 0;


    /** denoise by range a point cloud. 
     * @param in_pcl            input point cloud
     * @param out_pcl           denoised point cloud.
     * @param in_windowSize     "averaging" window size
     * @param in_noiseTh        scale of noise. Example (for median filter): if distance from averge
     *                          is below in_noiseTh - point is snapped to average)
     */
    virtual void DenoiseRange(const CPtCloud& in_pcl, CPtCloud& out_pcl,
                              float in_windowSize, float in_noiseTh) = 0;


    /** downsample a point cloud. Divides to grid from minXYZ (of pts) to max XYZ, of size m_voxelSize.
     *  If more than one point in same grid index, takes the first one encountered.
     * @param in_pcl             input point cloud.
     * @param out_pcl            downsampled point cloud. Can be the same as in_pcl
     * @param in_voxelSize       size of a voxel in grid. assums bigger than 0
     */
    virtual void DownSample(const CPtCloud& in_pcl, CPtCloud& out_pcl, float in_voxelSize) = 0;

  protected:
    /** destructor */
    virtual ~IFeatures() = 0;
    
    friend class CTinyPCL;
  };
 

} // namespace tpcl

#endif
