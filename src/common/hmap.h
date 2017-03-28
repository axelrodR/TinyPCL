// File Location: 

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
*: Package Name: tpcl_hmap
*
*: Title: Height maps
*
*  We define several types of height maps: 
*    - CSimpleHgtMap: single surface height map (2B per cell)
*    - CDynamicHeightMap: multiple surfaces, dynamic height map with extra
*                         info/color per cell.
*    - CCompactHeightMap: multiple surfaces, static with small foot-print with
                          extra info/color per cell and connectivity to adjacent cells
                          to represent continous surfaces
******************************************************************************/

#ifndef __tpcl_hmap_H
#define __tpcl_hmap_H

#include "grid.h"        // for CGrid2D
#include "mesh.h"        // for IRasterizable

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#define DLL_Entry 
#ifndef _LIB 
#ifndef _LIB_LINK
#undef DLL_Entry 
#ifdef TPCL_EXPORTS
#define DLL_Entry __declspec(dllexport)
#else
#define DLL_Entry __declspec(dllimport)
#endif
#endif
#endif

namespace tpcl
{
  /******************************************************************************
  *                             EXPORTED CONSTANTS                              *
  ******************************************************************************/

  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
  struct CSpanPool;                 // in this file
  struct CVertexUV;                 // in mesh.h
  struct CCompactSpan;              // in this file

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/


  ///////////////////////////////////////////////////////////////////////////////
  /// Simple height map.
  /// 
  /// Single (top) surface, 2B per cell
  ///////////////////////////////////////////////////////////////////////////////
  class DLL_Entry CSimpleHgtMap : public CGrid2D<unsigned short>, public tpcl::IRasterizable
  {
  public:
    // height resolution
    float GetHeightRes() const                  {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const               {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** Get height */
    float GetZ(int Xi_x, int Xi_y)  const       {return Get(Xi_x,Xi_y) * m_resH + m_bbMin.z;}
    float GetZ(CVec3& Xi_pos) const       {return Get(Xi_pos) * m_resH + m_bbMin.z;}
    void SetZ(int Xi_x,int Xi_y,float Xi_h)     {Get(Xi_x,Xi_y) = unsigned short((Xi_h-m_bbMin.z)*m_invResH);}      ///< set height

    /** get the world position associated with a cell coordinates (corner of a cell) */
    inline CVec3 GetPos(int Xi_x, int Xi_y) const;

    /** Rasterize a triangle height to the map
     * @param Xi_info  not used */
    bool U_RasterizeTri(const CVec3& Xi_v0, const CVec3& Xi_v1, const CVec3& Xi_v2, int Xi_info);

    /** Rasterize a triangle height and info/color
     * @param Xi_info  not used */
    bool U_RasterizeTri(const CVertexUV& Xi_v0, const CVertexUV& Xi_v1, const CVertexUV& Xi_v2,
                        int Xi_imgWidth, int Xi_imgHeight, const unsigned int* Xi_img);

    // inhertied from IRasterizable
    /** Get vertices with texture */
    virtual const CMesh* U_GetMesh();
    virtual void U_GetBBox(CVec3& Xo_min, CVec3& Xo_max) const;
    virtual const unsigned int* U_GetTexture(int& Xo_texWidth, int& Xo_texHeight);

    /** constructors */
    CSimpleHgtMap(int Xi_width, int Xi_height, const CVec3& Xi_bbmin, float Xi_res=1.0f);
    CSimpleHgtMap(const CVec3& Xi_bbMin, const CVec3& Xi_bbMax, float Xi_res=1.0f);
    
    /** destructor */
    virtual ~CSimpleHgtMap();

  protected:
    float m_resH, m_invResH;    ///< height resolution: the size of each cell (m_invRes=1.0/m_res)
  };



  ///////////////////////////////////////////////////////////////////////////////
  /// Dynamic height map.
  ///
  /// multiple surfaces, dynamic height map with extra info/color per cell.
  /// a variation of Recast/Detour's rcHeightfield.
  /// see: https://github.com/recastnavigation/recastnavigation
  ///////////////////////////////////////////////////////////////////////////////
  class DLL_Entry CDynamicHeightMap : public CGrid2dBase
  {
  public:
    static const int SPAN_HEIGHT_BITS = 16; // number of bits used for height
    static const unsigned int NULL_INFO  = 0;

    /** height span: a colum of filled space */
    struct CSpan
    {
      unsigned short smin;  ///< The lower limit of the span. [Limit: < #smax]
      unsigned short smax;  ///< The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT]
      unsigned int info : 30;///< 30bit information
      unsigned int facing: 2;///< back facing polygons = 1, front facing = 0
      CSpan* m_next;        ///< The next span higher up in column.
    };

    /** Convert height into world coordinates */
    float GetSpanHeight(CSpan& Xi_span)             {return Xi_span.smax * m_resH + m_bbMin.z;}

    // wrapper to access spans
    CSpan* Get(int Xi_x, int Xi_y)                  {return *(CSpan**)CGrid2dBase::Get(Xi_x,Xi_y);}    ///< the first span using cell coordinates
    CSpan* Get(int Xi_ind)                          {return *(CSpan**)CGrid2dBase::Get(Xi_ind);}       ///< the first span using an index
    CSpan* Get(const CVec3& Xi_pos)           {return *(CSpan**)CGrid2dBase::Get(Xi_pos);}       ///< get first span using a world pos
    const CSpan* Get(int Xi_x, int Xi_y) const      {return *(CSpan**)CGrid2dBase::Get(Xi_x,Xi_y);}    ///< the first span using cell coordinates
    const CSpan* Get(int Xi_ind) const              {return *(CSpan**)CGrid2dBase::Get(Xi_ind);}       ///< the first span using an index
    const CSpan* Get(const CVec3& Xi_pos) const {return *(CSpan**)CGrid2dBase::Get(Xi_pos);}       ///< get first span using a world pos

    // height resolution
    float GetHeightRes() const                      {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const                   {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                    {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** get the total number of spans used */
    int GetSpanCount()  const                       {return m_spanCount;} 

    /** Rasterize a triangle into a dynamic height map (rasterization buffer) */
    bool U_RasterizeTri(const CVec3& Xi_v0, const CVec3& Xi_v1, const CVec3& Xi_v2, int Xi_area);

    /** @brief Rasterize a triangle with texture into a dynamic height map (rasterization buffer)
     * @param Xi_img default: RGB + area (byte) in the alpha channel */
    bool U_RasterizeTri(const CVertexUV& Xi_v0, const CVertexUV& Xi_v1, const CVertexUV& Xi_v2,
                        int Xi_imgWidth, int Xi_imgHeight, const unsigned int* Xi_img);

    /** add a span (a height range that is blocked */
    bool AddSpan(int Xi_x, int Xi_y, int Xi_hBegin, int Xi_hEnd, int Xi_area, int facing, int Xi_mergeThr);

    /** constructors */
    CDynamicHeightMap(int Xi_width, int Xi_height, const CVec3& Xi_bbmin, float Xi_res=1.0f);
    CDynamicHeightMap(const CVec3& Xi_bbMin, const CVec3& Xi_bbMax, float Xi_res=1.0f);
    
    /** destructor */
    virtual ~CDynamicHeightMap();

    bool HasColorInfo() const                       {return m_hasColorInfo;}

  protected:
    /** allocate a span from the pool */
    CSpan* AllocSpan();   
    /** free a span to the pool*/
    void FreeSpan(CSpan* Xi_pSpan);

    float m_resH, m_invResH;    ///< height resolution: the size of each cell (m_invRes=1.0/m_res)
    CSpan* m_freeList;          ///< linked list of free spans
    CSpanPool* m_pools;         ///< memory pools
    int m_spanCount;            ///< number of used spans
    bool m_hasColorInfo;        ///< true=information is color, false=information is area (0..255)
  };



  /// Provides information on the content of a cell column in a compact heightfield. 
  struct CCompactCell
  {
	  unsigned int index : 24;    ///< Index to the first span in the column.
	  unsigned int count : 8;     ///< Number of spans in the column.
  };


  ///////////////////////////////////////////////////////////////////////////////
  /// Compact height map.
  /// 
  /// a variation of Recast/Detour's rcCompactHeightfield.
  /// Holds information about "voxels" in column major format (i.e. multiple arrays
  /// and an index to access them all). 
  /// Each cell points to an arrray of compact spans that are above that cell.
  /// Each span holds information about "connected" neighbors (see CCompactSpan)
  /// All spans are stored in one array (the cell uses an index)  
  /// see: https://github.com/recastnavigation/recastnavigation
  ///////////////////////////////////////////////////////////////////////////////
  class DLL_Entry CCompactHeightMap : public CGrid2D<CCompactCell>
  {
  public:

    // height resolution
    float GetHeightRes() const                      {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const                   {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                    {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** get the total number of spans used */
    int GetSpanCount()  const                       {return m_spanCount;}   

    /** Fill spans connection information to form continous surfaces 
     * @param Xi_walkableHeight   open height above surface to be considered "seperate"
     * @param Xi_walkableClimb    height difference below which surface are connected */
    void ComputeSurfaces(float Xi_walkableHeight, float Xi_walkableClimb);

    /** Fill buffer with all spans height above an x-y point (ordered from high to low)
     * @param Xi_bufSize    size of floats buffer to fill (note: 1 gets top height)
     * @return number of spans
     */
    int GetHeightAt(const CVec3& Xi_pos, int Xi_bufSize, float* Xo_heights) const;

    /** Get map height by span index (syntactic sugar) */
    inline float GetZ(int Xi_spanId) const;

    /** stntactic sugar to get the index of the top span */
    int GetTopSpanId(const CVec3& Xi_pos) const;

    /** construct from a dynamic height map
    * @param Xi_walkableHeight   open height above surface to be considered "seperate"
    * @param Xi_walkableClimb    height difference below which surface are connected (-1 do not connect)
    * @param Xi_backCull         Remove back facing traingles */
    CCompactHeightMap(const CDynamicHeightMap& Xi_dynHmap, float Xi_walkableHeight, float Xi_walkableClimb, bool Xi_backCull=false);

    /** destructor */
    virtual ~CCompactHeightMap();

    int m_borderSize;             ///< The AABB border size used during the build of the field. (See: rcConfig::borderSize)
    
    // extra data per span
    unsigned short* m_colors;     ///< color information
    unsigned char* m_areas;       ///< Array containing area ids
  
    static const unsigned char  NULL_AREA  = 0;
    static const unsigned char  BACK_FACE = 0x80; // back facing area mask
    static const unsigned short BORDER_REG = 0x8000;

  /*protected:*/
    CCompactSpan* m_spans;        ///< compact spans array
    float m_resH, m_invResH;      ///< height resolution: the size of each cell (m_invRes=1.0/m_res)
    int m_spanCount;              ///< number of used spans
  };


  /// Represents a span of unobstructed space within a compact heightfield.
  struct CCompactSpan
  {
    unsigned short y;           ///< The lower extent of the span. (Measured from the heightfield's base.)
    //unsigned short reg;       ///< The id of the region the span belongs to. (Or zero if not in a region.)
    unsigned int con : 24;      ///< Packed neighbor connection data.
    unsigned int h : 8;         ///< The height of free space above the span.  (Measured from #y.)
    unsigned int hv : 8;        ///< the height of the filled space above the y

    /** Gets index of neighbor for the specified direction 
     * @param[in]	dir	The index of the neighbor span (0-3)
     */
    int GetCon(int dir) const
    {
	    const unsigned int shift = (unsigned int)dir*6;
	    return (con >> shift) & 0x3f;
    }

    /** Sets the neighbor connection data in CompactSpan.con for the specified direction.
     * @param[in]	dir	The index of the neighbor span (0-3)
     */
    inline void SetCon(int dir, int i)
    {
	    const unsigned int shift = (unsigned int)dir*6;
	    unsigned int tmp = con;
	    con = (tmp & ~(0x3f << shift)) | (((unsigned int)i & 0x3f) << shift);
    }

    /** The value returned by CompactSpan::GetCon if the specified direction is not connected
     * to another span. (Has no neighbor.) */
    static const int NOT_CONNECTED = 0x3f;
  };



/******************************************************************************
*                        INLINE FUNCTIONS IMPLEMENTATION                      *
******************************************************************************/


  inline float CCompactHeightMap::GetZ(int Xi_spanId) const
  {
    return m_spans[Xi_spanId].y * m_resH + m_bbMin.z;
  }

  inline CVec3 CSimpleHgtMap::GetPos(int Xi_x, int Xi_y) const
  {
    return m_bbMin + CVec3(float(Xi_x)*m_res, float(Xi_y)*m_res, Get(Xi_x, Xi_y) * m_resH);
  }


} // namespace tpcl


#undef DLL_Entry
#endif


/******************************************************************************
*                    OLD STYLE TYPEDEFS                                       *
******************************************************************************/

typedef tpcl::CSimpleHgtMap TPCL_HMAP_CSimpleHgtMap;
typedef tpcl::CDynamicHeightMap TPCL_HMAP_CDynamicHeightMap;
typedef tpcl::CCompactHeightMap TPCL_HMAP_CCompactHeightMap;
