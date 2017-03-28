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


  /******************************************************************************
  *
  *: Class name: GenGmtrxCR_SP_CCoarseRegister
  *
  *: Abstract: basic features calculation on point cloud.
  *            see: <FILL PAPER REF HERE> 
  *
  ******************************************************************************/

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


    /** TODO: fills a point cloud.
    * @param Xi_pts           input point cloud.
    * @param Xo_pts           filled point cloud. */
    static void FillPointCloud(int Xi_numPts, const CVec3* Xi_pts, int Xo_numPts, const CVec3* Xo_pts);


    /** finds the normalized normals for a set of points in a global hashed point cloud.
    *   also, updates the z value of the points to be on the found plane.
    * @param Xi_globalHashed    input hashed global point cloud.
    * @param Xio_pts            input: points where to look for the normal. output: point on the plane with the same x,y.
    * @param Xo_Normals         filled point cloud. */
    static void FindNormal(int Xi_radius, CSpatialHash2D& Xi_globalHashed, int Xi_numPts, CVec3* Xio_pts, CVec3* Xo_Normals);


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
