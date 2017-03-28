//
// Copyright (c) 2016-2017 Geosim Ltd.
// 
// Written by Ramon Axelrod       ramon.axelrod@gmail.com
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
*: Package Name: SpatialHash
*
*: Description: defines a simple spatial hashing geared toward 2.5D data
*
******************************************************************************/


#ifndef __SPATIAL_HASH_H
#define __SPATIAL_HASH_H

#include "../../include/vec.h"

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#define DLL_Entry 
#ifndef _LIB 
#ifndef _LIB_LINK
#undef DLL_Entry 
#ifdef GENGMTRX_EXPORTS
#define DLL_Entry __declspec(dllexport)
#else
#define DLL_Entry __declspec(dllimport)
#endif
#endif
#endif

/******************************************************************************
*                              EXPORTED CLASSES                               *
******************************************************************************/

namespace tpcl
{

  /**************************************************************************//**
  *
  * Spatial hashing for nearest neighbor search
  *
  * Meant for 2.5D searches.
  * Projected 2D space is divided into bins (i.e. bins ignore z). Nearwest neighbor
  * searches each cell in a radius.
  *
  ******************************************************************************/
  class DLL_Entry CSpatialHash2D
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/

    /** add an object+position pair to the spatial hash */
    void Add(const CVec3& Xi_pos, void* Xi_obj);

    /** 
     * Find nearest object (using its associated point)
     * @param Xo_pMinPt (optional) closest point
     * @param max2DRadius maximum 2D radius to search in  */
    void* FindNearest(const CVec3& Xi_pos, CVec3* Xo_pMinPt=0, float Xi_max2DRadius=0.0f) const;


    /** Get all objects in 2D radius 
     * @param Xo_buf        buffer to fill with objects
     * @param max2DRadius   maximum 2D radius to search in  
     * @return    nuumber of objects*/
    int GetNear(const CVec3& Xi_pos, int xi_bufSize, void** Xo_buf, CVec3* Xo_pos=0, float Xi_max2DRadius=0.0f) const;

    /** Clear data */
    void Clear();

    /** constructor */
    CSpatialHash2D (float res=0.1f);

    /** destructor */
    virtual ~CSpatialHash2D ();

  protected:
  /******************************************************************************
  *                             Protected members                               *
  ******************************************************************************/

    /** Each cell is hashed to an array of object+position pairs */
    void* m_data; 
  
#pragma warning (disable : 4251)
    float m_res, m_resInv;            ///< // internal "grid" resolution
    CVec3 m_minBox, m_maxBox;   ///< bounding box size
    CVec3 m_pivot;              ///< used internally to shift everything to be around (0,0,0)
#pragma warning (default : 4251)
  };

  /******************************************************************************
  *             a template wrapper to the Spatial hash 2D class                 *
  ******************************************************************************/
  
  /** a template wrapper to the Spatial hash 2D class */
  template <typename T> class TSpatialHash2D : public CSpatialHash2D
  {
  public:
    void Add(const CVec3& Xi_pos, T& Xi_obj)
      {CSpatialHash2D::Add(Xi_pos, (void*)&Xi_obj);}
    
    T* FindNearest(const CVec3& Xi_pos, CVec3* Xo_pMinPt=0, float Xi_max2DRadius=0.0f) const
      {return (T*)CSpatialHash2D::FindNearest(Xi_pos, Xo_pMinPt, Xi_max2DRadius);}
    
    int GetNear(const CVec3& Xi_pos, int Xi_bufSize, T** Xo_buf, CVec3* Xo_pos=0, float Xi_max2DRadius=0.0f) const
    {return CSpatialHash2D::GetNear(Xi_pos, Xi_bufSize, (void**)Xo_buf, Xo_pos, Xi_max2DRadius);}

  TSpatialHash2D(float res=0.1f) : CSpatialHash2D(res){}
  };



}// namespace __SPATIAL_HASH_H

/******************************************************************************
*                            EXPORTED FUNCTIONS                               *
******************************************************************************/

#undef DLL_Entry
#endif

