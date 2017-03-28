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
******************************************************************************/
//#include <D3dx9core.h> // uncommenet if using DirectX
//#include <ifr/ifrgen/ifrgen_stnd.h>
#include "features.h"
////#include <gen/gengmtrx/gengmtrx_spat.h>
#include "plane.h"
//#include <vector>


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


namespace tpcl
{

  /******************************************************************************
  *                             INTERNAL CONSTANTS  / Functions                 *
  ******************************************************************************/
  
  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

  /******************************************************************************
  *                       FORWARD FUNCTION DECLARATIONS                         *
  ******************************************************************************/

  /******************************************************************************
  *                             STATIC VARIABLES                                *
  ******************************************************************************/

  /******************************************************************************
  *                      CLASS STATIC MEMBERS INITIALIZATION                    *
  ******************************************************************************/

  /******************************************************************************
  *                              INTERNAL CLASSES                               *
  ******************************************************************************/

 


  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CCoarseRegister
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: tpcl_SP_Features
  *
  ******************************************************************************/
  Features::Features()
  {
  }
 
  /******************************************************************************
  *
  *: Method name: ~tpcl_SP_Features
  *
  ******************************************************************************/
  Features::~Features()
  {
  }


  /******************************************************************************
  *
  *: Method name: FillPointCloud
  *
  ******************************************************************************/
  void Features::FillPointCloud(int Xi_numPts, const CVec3* Xi_pts, int Xo_numPts, const CVec3* Xo_pts)
  {
  }


  /******************************************************************************
  *
  *: Method name: FillPointCloud
  *
  ******************************************************************************/
  void Features::FindNormal(int Xi_radius, CSpatialHash2D& Xi_globalHashed, int Xi_numPts, CVec3* Xio_pts, CVec3* Xo_Normals)
  {
    const int bufSize = 100;
    void* unsuedBuf[bufSize];
    CVec3 closePts[bufSize];
    int numOfClose = 0;

    for (int ptIndex = 0; ptIndex < Xi_numPts; ptIndex++)
    {
      //get close points:
      numOfClose = Xi_globalHashed.GetNear(Xio_pts[ptIndex], bufSize, unsuedBuf, closePts, Xi_radius);

      //find plane:
      CPlane approxPlane;
      approxPlane.RanSaC(numOfClose, closePts, 1.0f);

      //update point's z:
      Xio_pts[ptIndex].z = approxPlane.GetHeightAt(Xio_pts[ptIndex].x, Xio_pts[ptIndex].y);

      //get normal and normalize it, up:
      Xo_Normals[ptIndex] = approxPlane.GetNormal();
      if (Xo_Normals[ptIndex].z < 0)
        Xo_Normals[ptIndex].z = -Xo_Normals[ptIndex].z;
      Normalize(Xo_Normals[ptIndex]);
    }
  }


  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/




  /******************************************************************************
  *                              Private methods                                *
  ******************************************************************************/


  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/

  /******************************************************************************
  *                            INTERNAL FUNCTIONS                               *
  ******************************************************************************/


} //namespace tpcl