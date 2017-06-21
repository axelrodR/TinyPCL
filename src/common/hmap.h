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
  class CSimpleHgtMap : public CGrid2D<unsigned short>, public tpcl::IRasterizable
  {
  public:
    // height resolution
    float GetHeightRes() const                  {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const               {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** Get height */
    float GetZ(int in_x, int in_y)  const       {return Get(in_x,in_y) * m_resH + m_bbMin.z;}
    float GetZ(CVec3& in_pos) const       {return Get(in_pos) * m_resH + m_bbMin.z;}
    void SetZ(int in_x,int in_y,float in_h)     {Get(in_x,in_y) = unsigned short((in_h-m_bbMin.z)*m_invResH);}      ///< set height

    /** get the world position associated with a cell coordinates (corner of a cell) */
    inline CVec3 GetPos(int in_x, int in_y) const;

    /** Rasterize a triangle height to the map
     * @param in_info  not used */
    bool U_RasterizeTri(const CVec3& in_v0, const CVec3& in_v1, const CVec3& in_v2, int in_info);

    /** Rasterize a triangle height and info/color
     * @param in_info  not used */
    bool U_RasterizeTri(const CVertexUV& in_v0, const CVertexUV& in_v1, const CVertexUV& in_v2,
                        int in_imgWidth, int in_imgHeight, const unsigned int* in_img);

    // inhertied from IRasterizable
    /** Get vertices with texture */
    virtual const CMesh* U_GetMesh();
    virtual void U_GetBBox(CVec3& out_min, CVec3& out_max) const;
    virtual const unsigned int* U_GetTexture(int& out_texWidth, int& out_texHeight);

    /** constructors */
    CSimpleHgtMap(int in_width, int in_height, const CVec3& in_bbmin, float in_res=1.0f);
    CSimpleHgtMap(const CVec3& in_bbMin, const CVec3& in_bbMax, float in_res=1.0f);
    
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
  class CDynamicHeightMap : public CGrid2dBase
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
    float GetSpanHeight(CSpan& in_span)             {return in_span.smax * m_resH + m_bbMin.z;}

    // wrapper to access spans
    CSpan* Get(int in_x, int in_y)                  {return *(CSpan**)CGrid2dBase::Get(in_x,in_y);}    ///< the first span using cell coordinates
    CSpan* Get(int in_ind)                          {return *(CSpan**)CGrid2dBase::Get(in_ind);}       ///< the first span using an index
    CSpan* Get(const CVec3& in_pos)           {return *(CSpan**)CGrid2dBase::Get(in_pos);}       ///< get first span using a world pos
    const CSpan* Get(int in_x, int in_y) const      {return *(CSpan**)CGrid2dBase::Get(in_x,in_y);}    ///< the first span using cell coordinates
    const CSpan* Get(int in_ind) const              {return *(CSpan**)CGrid2dBase::Get(in_ind);}       ///< the first span using an index
    const CSpan* Get(const CVec3& in_pos) const {return *(CSpan**)CGrid2dBase::Get(in_pos);}       ///< get first span using a world pos

    // height resolution
    float GetHeightRes() const                      {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const                   {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                    {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** get the total number of spans used */
    int GetSpanCount()  const                       {return m_spanCount;} 

    /** Rasterize a triangle into a dynamic height map (rasterization buffer) */
    bool U_RasterizeTri(const CVec3& in_v0, const CVec3& in_v1, const CVec3& in_v2, int in_area);

    /** @brief Rasterize a triangle with texture into a dynamic height map (rasterization buffer)
     * @param in_img default: RGB + area (byte) in the alpha channel */
    bool U_RasterizeTri(const CVertexUV& in_v0, const CVertexUV& in_v1, const CVertexUV& in_v2,
                        int in_imgWidth, int in_imgHeight, const unsigned int* in_img);

    /** add a span (a height range that is blocked */
    bool AddSpan(int in_x, int in_y, int in_hBegin, int in_hEnd, int in_area, int facing, int in_mergeThr);

    /** constructors */
    CDynamicHeightMap(int in_width, int in_height, const CVec3& in_bbmin, float in_res=1.0f);
    CDynamicHeightMap(const CVec3& in_bbMin, const CVec3& in_bbMax, float in_res=1.0f);
    
    /** destructor */
    virtual ~CDynamicHeightMap();

    bool HasColorInfo() const                       {return m_hasColorInfo;}

  protected:
    /** allocate a span from the pool */
    CSpan* AllocSpan();   
    /** free a span to the pool*/
    void FreeSpan(CSpan* in_pSpan);

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
  class CCompactHeightMap : public CGrid2D<CCompactCell>
  {
  public:

    // height resolution
    float GetHeightRes() const                      {return m_resH;}       ///< get resolution
    float GetInvHeightRes() const                   {return m_invResH;}    ///< get 1/resolution
    void SetHeightRes(float res)                    {m_resH = res; m_invRes=1.0f/res;} ///< set resolution

    /** get the total number of spans used */
    int GetSpanCount()  const                       {return m_spanCount;}   

    /** Fill spans connection information to form continous surfaces 
     * @param in_walkableHeight   open height above surface to be considered "seperate"
     * @param in_walkableClimb    height difference below which surface are connected */
    void ComputeSurfaces(float in_walkableHeight, float in_walkableClimb);

    /** Fill buffer with all spans height above an x-y point (ordered from high to low)
     * @param in_bufSize    size of floats buffer to fill (note: 1 gets top height)
     * @return number of spans
     */
    int GetHeightAt(const CVec3& in_pos, int in_bufSize, float* out_heights) const;

    /** Get map height by span index (syntactic sugar) */
    inline float GetZ(int in_spanId) const;

    /** stntactic sugar to get the index of the top span */
    int GetTopSpanId(const CVec3& in_pos) const;

    /** construct from a dynamic height map
    * @param in_walkableHeight   open height above surface to be considered "seperate"
    * @param in_walkableClimb    height difference below which surface are connected (-1 do not connect)
    * @param in_backCull         Remove back facing traingles */
    CCompactHeightMap(const CDynamicHeightMap& in_dynHmap, float in_walkableHeight, float in_walkableClimb, bool in_backCull=false);

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


  inline float CCompactHeightMap::GetZ(int in_spanId) const
  {
    return m_spans[in_spanId].y * m_resH + m_bbMin.z;
  }

  inline CVec3 CSimpleHgtMap::GetPos(int in_x, int in_y) const
  {
    return m_bbMin + CVec3(float(in_x)*m_res, float(in_y)*m_res, Get(in_x, in_y) * m_resH);
  }


} // namespace tpcl

#endif
