/******************************************************************************
*
*: Package Name: 
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_ptCloud_H
#define __tpcl_ptcloud_H

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/
#include "vec.h"



  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
namespace tpcl
{
  typedef unsigned int CColor;
  enum EPtCloudType;

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: PtCloud
  *
  *: Abstract: Point cloud
  *
  ******************************************************************************/

  struct  CPtCloud
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/

    /** positions of each cloud point (number of points in m_numPts */
    CVec3* m_pos;

    /** optional: color associated with each point (0 = no color information) */
    CColor* m_color;

    /** normals associated with each point (0 = no color information) */
    CVec3* m_normal;

    /** number of points in the cloud */
    int m_numPts;

    /** type of point cloud see EPtCloudType */
    EPtCloudType m_type;

    /** the number of point in a line for ordered 2d grid (usually polar coords)
     * relevant when m_tpye == PCL_TYPE_SINGLE_ORIGIN_SCAN */
    int m_lineWidth;


    // default constructor
    CPtCloud()
    { 
      m_numPts = 0;
      m_pos = m_normal = 0;
      m_color = 0; 
    }

    //// constructor with size
    //CPtCloud(int size, int components) 
    //{ 
    //  m_numPts = size; 
    //  m_pos = new CVec3[size];
    //  if (components & 1)
    //    m_normal = new CVec3[size];
    //  if (components & 2)
    //    m_color = new CColor[size];
    //}

    //// 
    //~CPtCloud()
    //{
    //  delete[] m_pos;
    //  delete[] m_color;
    //  delete[] m_normal;
    //}
  };


  /** Point cloud types */
  enum EPtCloudType
  {
    PCL_TYPE_SINGLE_ORIGIN = 1, ///< cloud created from a single point of view (succeptible to occlusions)
    PCL_TYPE_FUSED         = 2, ///< cloud created from using multiple sources (less occlusions, different resolutions, etc.)
    PCL_TYPE_SINGLE_ORIGIN_SCAN = 3, ///< single point of view and data is ordered
  };






} // namespace SLDR

#endif
