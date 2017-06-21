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
******************************************************************************/
//#include <d3dx9core.h>
//#include <ifr\ifrgen\ifrgen_stnd.h>
#include "hmap.h"
#include "mesh.h"
#include "grid.h"
//#include <ifr\ifrlog\ifrlog.h>
#include "vec.h"
#include "common.h"
#include <string.h>
#include <stdio.h>

#ifdef _DEBUG
#define new_file DEBUG_NEW
#endif

#define Log printf

namespace tpcl 
{

/******************************************************************************
*                             INTERNAL CONSTANTS                              *
******************************************************************************/

/// Defines the maximum value for rcSpan::smin and rcSpan::smax.
static const int SPAN_MAX_HEIGHT = (1 << CDynamicHeightMap::SPAN_HEIGHT_BITS) - 1;

/// The number of spans allocated per span spool.
static const int SPANS_PER_POOL = 2048;

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

struct CSpanPool
{
  CSpanPool* m_next;        ///< The next span pool.
  CDynamicHeightMap::CSpan m_items[SPANS_PER_POOL]; ///< Array of spans in the pool.
  CSpanPool(){}
};


bool overlapBounds(const CVec3& amin, const CVec3& amax, const CVec3& bmin, const CVec3& bmax)
{
  bool overlap = true;
  overlap = (amin.x > bmax.x|| amax.x < bmin.x) ? false : overlap;
  overlap = (amin.y > bmax.y || amax.y < bmin.y) ? false : overlap;
  overlap = (amin.z > bmax.z || amax.z < bmin.z) ? false : overlap;
  return overlap;
}


bool overlapInterval(unsigned short amin, unsigned short amax,
                     unsigned short bmin, unsigned short bmax)
{
  if (amax < bmin) return false;
  if (amin > bmax) return false;
  return true;
}


/// clap value to a range
template<class T> inline T Clamp(T v, T mn, T mx) { return v < mn ? mn : (v > mx ? mx : v); }

/** absolute of variable */
template<class T> inline T Abs(T a) { return a < 0 ? -a : a; }

/** absolute per component of 4 char in an int*/
inline unsigned int Max4C(unsigned int a , unsigned  int b) 
{ 
  unsigned char* ca = (unsigned char*)&a;
  unsigned char* cb = (unsigned char*)&b;
  int res;
  unsigned char* cr = (unsigned char*)&res;
  cr[0] = ca[0] > cb[0] ? ca[0] : cb[0]; 
  cr[1] = ca[1] > cb[1] ? ca[1] : cb[1]; 
  cr[2] = ca[2] > cb[2] ? ca[2] : cb[2]; 
  cr[3] = ca[3] > cb[3] ? ca[3] : cb[3]; 
  return res;
}


// divides a convex polygons into two convex polygons on both sides of a line
static void DividePoly(const CVec3* in, int nin,
                       CVec3* out1, int* nout1,
                       CVec3* out2, int* nout2,
                       float x, int axis)
{
  float d[12];
  for (int i = 0; i < nin; ++i)
    d[i] = x - in[i][axis];

  int m = 0, n = 0;
  for (int i = 0, j = nin-1; i < nin; j=i, ++i)
  {
    bool ina = d[j] >= 0;
    bool inb = d[i] >= 0;
    if (ina != inb)
    {
      float s = d[j] / (d[j] - d[i]);
      out1[m].x = in[j].x + (in[i].x - in[j].x)*s;
      out1[m].y = in[j].y + (in[i].y - in[j].y)*s;
      out1[m].z = in[j].z + (in[i].z - in[j].z)*s;
      out2[n] = out1[m];
      m++;
      n++;
      // add the i'th point to the right polygon. Do NOT add points that are on the dividing line
      // since these were already added above
      if (d[i] > 0)
      {
        out1[m] = in[i];
        m++;
      }
      else if (d[i] < 0)
      {
        out2[n] = in[i];
        n++;
      }
    }
    else // same side
    {
      // add the i'th point to the right polygon. Addition is done even for points on the dividing line
      if (d[i] >= 0)
      {
        out1[m] = in[i];
        m++;
        if (d[i] != 0)
          continue;
      }
      out2[n] = in[i];
      n++;
    }
  }

  *nout1 = m;
  *nout2 = n;
}



// divides a convex polygon with UVs into two convex polygons with UVs on both sides of a line
static void DividePoly(const CVertexUV* in, int nin,
                       CVertexUV* out1, int* nout1,
                       CVertexUV* out2, int* nout2,
                       float x, int axis)
{
  float d[12];
  for (int i = 0; i < nin; ++i)
    d[i] = x - in[i].m_vtx[axis];

  int m = 0, n = 0;
  for (int i = 0, j = nin-1; i < nin; j=i, ++i)
  {
    bool ina = d[j] >= 0;
    bool inb = d[i] >= 0;
    if (ina != inb)
    {
      float s = d[j] / (d[j] - d[i]);
      out1[m].m_vtx.x = in[j].m_vtx.x + (in[i].m_vtx.x - in[j].m_vtx.x)*s;
      out1[m].m_vtx.y = in[j].m_vtx.y + (in[i].m_vtx.y - in[j].m_vtx.y)*s;
      out1[m].m_vtx.z = in[j].m_vtx.z + (in[i].m_vtx.z - in[j].m_vtx.z)*s;
      out1[m].m_uv.x = in[j].m_uv.x + (in[i].m_uv.x - in[j].m_uv.x)*s;
      out1[m].m_uv.y = in[j].m_uv.y + (in[i].m_uv.y - in[j].m_uv.y)*s;
      out2[n] = out1[m];
      m++;
      n++;
      // add the i'th point to the right polygon. Do NOT add points that are on the dividing line
      // since these were already added above
      if (d[i] > 0)
      {
        out1[m] = in[i];
        m++;
      }
      else if (d[i] < 0)
      {
        out2[n] = in[i];
        n++;
      }
    }
    else // same side
    {
      // add the i'th point to the right polygon. Addition is done even for points on the dividing line
      if (d[i] >= 0)
      {
        out1[m] = in[i];
        m++;
        if (d[i] != 0)
          continue;
      }
      out2[n] = in[i];
      n++;
    }
  }

  *nout1 = m;
  *nout2 = n;
}


/** find the best height resolution to fit the box and spatial resolution */
float FindBestHeightResolution(float xi_minZ, float in_maxZ, float in_xyRes)
{
  float l_diff = in_maxZ - xi_minZ;   // the entire height range should use 16bits
  float l_stdDiff = in_xyRes * 20.0f;   // heights several times larger than the x-y resolution should settle nicely in 8bits
  for (float l_scale = 0.001f ; l_scale<1000.0f; l_scale*=1.5f)
  {
    float l_maxDiff = l_scale * (float)((1<<CDynamicHeightMap::SPAN_HEIGHT_BITS) -1);
    float l_minStdDiff = l_scale * (float)0xff;       
    if ( (l_diff*2 < l_maxDiff) && (l_stdDiff < l_minStdDiff))
      return l_scale;
  }
  return 1.0f;
}

/******************************************************************************
*                           EXPORTED CLASS METHODS                            *
******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//
//                           tpcl::CSimpleHgtMap
//
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
*                               Public methods                                *
******************************************************************************/


CSimpleHgtMap::CSimpleHgtMap(int in_width, int in_height, const CVec3& in_bbmin, float in_res)
  : CGrid2D<unsigned short>(in_width, in_height, in_bbmin, in_res)
{
  m_resH = 0.01f;
  m_invResH = 1.0f/m_resH;
  // creation of parent failed?
  if (m_data != 0)
    memset(m_data, 0, sizeof(unsigned short)*m_width*m_height);
}


CSimpleHgtMap::CSimpleHgtMap(const CVec3& in_bbMin, const CVec3& in_bbMax, float in_res)
  : CGrid2D<unsigned short>(in_bbMin, in_bbMax, in_res)
{
  // some automatic setup to fit the scale of the problem. The user should use SetHeightRes() to override
  m_resH = FindBestHeightResolution(in_bbMin.z, in_bbMax.z, in_res);
  m_invResH = 1.0f/m_resH;
  // creation of parent failed?
  if (m_data != 0)
    memset(m_data, 0, sizeof(unsigned short)*m_width*m_height);
}

CSimpleHgtMap::~CSimpleHgtMap()
{
}

    
bool CSimpleHgtMap::U_RasterizeTri(const CVec3& in_v0, const CVec3& in_v1, const CVec3& in_v2, int in_info)
{
  CVertexUV l_v0(in_v0, CVec3(0,0, 0));
  CVertexUV l_v1(in_v1, CVec3(1,0, 0));
  CVertexUV l_v2(in_v2, CVec3(1,1, 0));
  return U_RasterizeTri(l_v0, l_v1, l_v2, 1, 1, (unsigned int*)&in_info);
}


bool  CSimpleHgtMap::U_RasterizeTri(const CVertexUV& in_v0, const CVertexUV& in_v1, const CVertexUV& in_v2,
                                    int /*in_imgWidth*/, int /*in_imgHeight*/, const unsigned int* /*in_img*/)
{
  const int w = GetWidth();
  const int h = GetHeight();
  const float cs = GetRes();
  const float ics = 1.0f / cs;
  const float ch = GetHeightRes();
  const float ich = GetInvHeightRes();
  unsigned short* cells = (unsigned short*)m_data;
  CVec3 tmin, tmax;
  const float bz = m_bbMax.z - m_bbMin.z;

  // Calculate the bounding box of the triangle.
  tmin = in_v0.m_vtx;
  tmax = in_v0.m_vtx;
  tmin = Min_ps(tmin, in_v1.m_vtx);
  tmin = Min_ps(tmin, in_v2.m_vtx);
  tmax = Max_ps(tmax, in_v1.m_vtx);
  tmax = Max_ps(tmax, in_v2.m_vtx);

  // If the triangle does not touch the bbox of the heightfield, skip the triagle.
  if (!overlapBounds(m_bbMin, m_bbMax, tmin, tmax))
    return true;

  // Calculate the footprint of the triangle on the grid's y-axis
  int y0 = (int)((tmin.y - m_bbMin.y)*ics);
  int y1 = (int)((tmax.y - m_bbMin.y)*ics);
  y0 = Clamp(y0, 0, h-1);
  y1 = Clamp(y1, 0, h-1);

  // Clip the triangle into all grid cells it touches.
  CVec3 buf[7*4];
  CVec3 *in = buf, *inrow = buf+7, *p1 = inrow+7, *p2 = p1+7;
  
  in[0] = in_v0.m_vtx;
  in[1] = in_v1.m_vtx;
  in[2] = in_v2.m_vtx;
  int nvrow, nvIn = 3;

  for (int y = y0; y <= y1; ++y)
  {
    // Clip polygon to row. Store the remaining polygon as well
    const float cz = m_bbMin[1] + y*cs;
    DividePoly(in, nvIn, inrow, &nvrow, p1, &nvIn, cz+cs, 1);
    SwapT(in, p1);
    if (nvrow < 3) continue;

    // find the horizontal bounds in the row
    float minX = inrow[0].x, maxX = inrow[0].x;
    for (int i=1; i<nvrow; ++i)
    {
      if (minX > inrow[i].x)	minX = inrow[i].x;
      if (maxX < inrow[i].x)	maxX = inrow[i].x;
    }
    int x0 = (int)((minX - m_bbMin[0])*ics);
    int x1 = (int)((maxX - m_bbMin[0])*ics);
    x0 = Clamp(x0, 0, w-1);
    x1 = Clamp(x1, 0, w-1);

    int nv, nv2 = nvrow;

    for (int x = x0; x <= x1; ++x)
    {
      // Clip polygon to column. store the remaining polygon as well
      const float cx = m_bbMin[0] + x*cs;
      DividePoly(inrow, nv2, p1, &nv, p2, &nv2, cx+cs, 0);
      SwapT(inrow, p2);
      if (nv < 3) continue;

      // Calculate min and max of the span.
      float smin = p1[0].z, smax = p1[0].z;
      for (int i = 1; i < nv; ++i)
      {
        smin = MinT(smin, p1[i].z);
        smax = MaxT(smax, p1[i].z);
      }
      smin -= m_bbMin[2];
      smax -= m_bbMin[2];
      // Skip the span if it is outside the heightfield bbox
      if (smax < 0.0f) continue;
      if (smin > bz) continue;
      // Clamp the span to the heightfield bbox.
      if (smin < 0.0f) smin = 0;
      if (smax > bz) smax = bz;

      // Snap the span to the heightfield height grid.
      unsigned int ismax = (unsigned short)Clamp((int)ceilf(smax * ich), 0, 0xffff);
      if (ismax > cells[x+y*w])
        cells[x+y*w] = unsigned short(ismax);
    }
  }

  return true;
}



const CMesh* CSimpleHgtMap::U_GetMesh()
{
  int h = m_height, w = m_width;
  int m_numVtx = (h+1)*(w+1);
  int m_numFaces = h*w*2;
  CMesh* l_mesh = new CMesh(m_numVtx,m_numFaces, false);
  // generate vertices
  int iVtx = 0;
  CVec3 l_pos = m_bbMin;
  for (int iy=0; iy<=h; ++iy)
  {
    l_pos.x = m_bbMin.x;
    for (int ix=0; ix<=w; ++ix)
    {
      l_pos.z = GetZ(ix, iy);
      l_mesh->m_verts[iVtx++] = l_pos;
      l_pos.x += m_res;
    }
    l_pos.y += m_res;
  }
  // generate faces
  int iFace = 0;
  for (int iy=0; iy<h; ++iy)
  {
    l_pos = CVec3(m_bbMin.x, l_pos.y+m_res, 0);
    for (int ix=0; ix<w; ++ix)
    {
      l_mesh->m_faces[iFace*3]   = ix + iy*(w+1);
      l_mesh->m_faces[iFace*3+1] = (ix+1) + iy*(w+1);
      l_mesh->m_faces[iFace*3+2] = (ix+1) + (iy+1)*(w+1);
      iFace++;
      l_mesh->m_faces[iFace*3]   = ix + iy*(w+1);
      l_mesh->m_faces[iFace*3+1] = (ix+1) + (iy+1)*(w+1);
      l_mesh->m_faces[iFace*3+1] = ix + (iy+1)*(w+1);
      iFace++;
    }
  }
  return l_mesh;
}



void CSimpleHgtMap::U_GetBBox(CVec3& out_min, CVec3& out_max) const
{
  out_min = m_bbMin;
  out_max = m_bbMax;
}


static unsigned int s_singlePixel = 0xffffffff;
const unsigned int* CSimpleHgtMap::U_GetTexture(int& out_texWidth, int& out_texHeight)
{
  out_texWidth = out_texHeight = 1;
  return &s_singlePixel;
}




///////////////////////////////////////////////////////////////////////////////
//
//                            tpcl::CDynamicHeightMap
//
///////////////////////////////////////////////////////////////////////////////

/******************************************************************************
*                               Public methods                                *
******************************************************************************/
/******************************************************************************
*
*: Method name: CDynamicHeightMap
*
******************************************************************************/
CDynamicHeightMap::CDynamicHeightMap(int in_width, int in_height, const CVec3& in_bbmin, float in_res)
  : CGrid2dBase(in_width, in_height, in_bbmin, in_res, sizeof(CSpan*))
{
  m_freeList = 0;
  m_pools = 0;
  m_resH = 0.01f;
  m_invResH = 1.0f/m_resH;
  m_hasColorInfo = false;
  // creation of parent failed?
  if (m_data != 0)
    memset(m_data, 0, sizeof(CSpan*)*m_width*m_height);
}


CDynamicHeightMap::CDynamicHeightMap(const CVec3& in_bbMin, const CVec3& in_bbMax, float in_res)
  : CGrid2dBase(in_bbMin, in_bbMax, in_res, sizeof(CSpan*))
{
  m_freeList = 0;
  m_pools = 0;
  m_spanCount = 0;
  
  // some automatic setup to fit the scale of the problem. THe user should use SetHeightRes() to override
  m_resH = FindBestHeightResolution(in_bbMin.z, in_bbMax.z, in_res);
  m_invResH = 1.0f/m_resH;
  m_hasColorInfo = false;

  // creation of parent failed?
  if (m_data != 0)
    memset(m_data, 0, sizeof(CSpan*)*m_width*m_height);
}

/******************************************************************************
*
*: Method name: ~CDynamicHeightMap
*
******************************************************************************/
CDynamicHeightMap::~CDynamicHeightMap()
{
  // Delete span pools.
  while (m_pools)
  {
    CSpanPool* next = m_pools->m_next;
    delete[] m_pools;
    m_pools = next;
  }
}



CDynamicHeightMap::CSpan* CDynamicHeightMap::AllocSpan()
{
	// If running out of memory, allocate new page and update the freelist.
	if (!m_freeList || !m_freeList->m_next)
	{
    // Create new page.
    // Allocate memory for the new pool.
    CSpanPool* pool = new CSpanPool();
    if (!pool)
      return 0;

    // Add the pool into the list of pools.
    pool->m_next = m_pools;
    m_pools = pool;
    // Add new items to the free list.
    //CSpan* freelist = m_freeList;
    CSpan* head = &pool->m_items[0];
    CSpan* it = &pool->m_items[SPANS_PER_POOL];
    do
    {
      --it;
      it->m_next = m_freeList;
      m_freeList = it;
    }
    while (it != head);
    m_freeList = it;
  }

  // Pop item from in front of the free list.
  CSpan* it = m_freeList;
  m_freeList = m_freeList->m_next;
  return it;
}

    
void CDynamicHeightMap::FreeSpan(CSpan* in_pSpan)
{
  if (!in_pSpan) return;
  // Add the node in front of the free list.
  in_pSpan->m_next = m_freeList;
  m_freeList = in_pSpan;
}


// average feilds in a bitfelid
#define AVG(x, y, mask)    ( ((x&mask)+(y&mask)) >> 1 ) & mask

bool CDynamicHeightMap::AddSpan(int in_x, int in_y, int in_hBegin, int in_hEnd, int in_area, int facing, int in_mergeThr)
{
  int idx = GetIndex(in_x,in_y);
  
  CSpan* l_s = AllocSpan();
  if (!l_s)
    return false;
  l_s->smin = (unsigned short)in_hBegin;
  l_s->smax = (unsigned short)in_hEnd;
  l_s->info = in_area;
  l_s->facing = facing;
  l_s->m_next = 0;

  // Empty cell, add the first span.
  CSpan** l_ppCell = (CSpan**)CGrid2dBase::Get(idx);
  if (!*l_ppCell)
  {
    *l_ppCell = l_s;
    m_spanCount++;
    return true;
  }
  CSpan* l_prev = 0;
  CSpan* l_cur = *l_ppCell;

  // Insert and merge spans.
  while (l_cur)
  {
    if (l_cur->smin > l_s->smax)
    {
      // Current span is further than the new span, break.
      break;
    }
    else if (l_cur->smax < l_s->smin)
    {
      // Current span is before the new span advance.
      l_prev = l_cur;
      l_cur = l_cur->m_next;
    }
    else
    {
      // Merge spans.
      if (l_cur->smin < l_s->smin)
        l_s->smin = l_cur->smin;
      if (l_cur->smax > l_s->smax)
      {
        if (l_cur->smax - l_s->smax <= in_mergeThr)
        {
          if (!m_hasColorInfo)
            l_s->info = MaxT(l_s->info, l_cur->info);
          else  // average color
          {
            int b = AVG(l_s->info, l_cur->info, 0xff);
            int g = AVG(l_s->info, l_cur->info, 0xff<<8);
            int r = AVG(l_s->info, l_cur->info, 0xff<<16);
            int a = MaxT(l_s->info & (0xff<24), l_cur->info & (0xff<24));
            l_s->info = b | g | r | a;
          }
        }
        else
          l_s->info = l_cur->info;
        l_s->facing = l_cur->facing;
        l_s->smax = l_cur->smax;
      }
      else // l_cur->smax <= l_s->smax
      {
        if (l_s->smax - l_cur->smax <= in_mergeThr)
        {
          if (!m_hasColorInfo)
            l_s->info = MaxT(l_s->info, l_cur->info);
          else // average color
          {
            int b = AVG(l_s->info, l_cur->info, 0xff);
            int g = AVG(l_s->info, l_cur->info, 0xff<<8);
            int r = AVG(l_s->info, l_cur->info, 0xff<<16);
            int a = MaxT(l_s->info & (0xff<24), l_cur->info & (0xff<24));
            l_s->info = b | g | r | a;
          }
        }
      }

      // Remove current span.
      CSpan* l_next = l_cur->m_next;
      FreeSpan(l_cur);
      m_spanCount--;
      if (l_prev)
        l_prev->m_next = l_next;
      else
        *l_ppCell = l_next;
      l_cur = l_next;
    }
  }

  // Insert new span.
  if (l_prev)
  {
    l_s->m_next = l_prev->m_next;
    l_prev->m_next = l_s;
  }
  else
  {
    l_s->m_next = *l_ppCell;
    *l_ppCell = l_s;
  }

  m_spanCount++;
  return true;
}





/******************************************************************************
*
*: Method name: RasterizeTri
*
* a variation of Recast/Detour's rcHeightfield.
* see: https://github.com/recastnavigation/recastnavigation
******************************************************************************/
bool CDynamicHeightMap::U_RasterizeTri(const CVec3& in_v0, const CVec3& in_v1, const CVec3& in_v2, int in_area)
{
  CVertexUV l_v0(in_v0,CVec3(0,0, 0));
  CVertexUV l_v1(in_v1,CVec3(1,0, 0));
  CVertexUV l_v2(in_v2,CVec3(1,1, 0));
  return U_RasterizeTri(l_v0, l_v1, l_v2, 1, 1, (unsigned int*)&in_area);
}


/******************************************************************************
*
*: Method name: RasterizeTri
*
* a variation of Recast/Detour's rcHeightfield.
* see: https://github.com/recastnavigation/recastnavigation
******************************************************************************/
bool CDynamicHeightMap::U_RasterizeTri(const CVertexUV& in_v0, const CVertexUV& in_v1, const CVertexUV& in_v2,
                                     int in_imgWidth, int in_imgHeight, const unsigned int* in_img)
{
  const int w = GetWidth();
  const int h = GetHeight();
  const float cs = GetRes();
  const float ics = 1.0f / cs;
  const float ch = GetHeightRes();
  const float ich = GetInvHeightRes();
  CVec3 tmin, tmax;
  const float bz = m_bbMax.z - m_bbMin.z;
  const float l_imgW = (float)in_imgWidth;
  const float l_imgH = (float)in_imgHeight;
  m_hasColorInfo = true;

  // Calculate the bounding box of the triangle.
  tmin = in_v0.m_vtx;
  tmax = in_v0.m_vtx;
  tmin = Min_ps(tmin, in_v1.m_vtx);
  tmin = Min_ps(tmin, in_v2.m_vtx);
  tmax = Max_ps(tmax, in_v1.m_vtx);
  tmax = Max_ps(tmax, in_v2.m_vtx);

  // If the triangle does not touch the bbox of the heightfield, skip the triagle.
  if (!overlapBounds(m_bbMin, m_bbMax, tmin, tmax))
    return true;

  // claculate the front facing / back facing
  CVec3 e1 = in_v1.m_vtx - in_v0.m_vtx;
  CVec3 e2 = in_v2.m_vtx - in_v0.m_vtx;
  float crossZ = e1.x*e2.y - e2.x*e1.y;
  int backFace = crossZ < 0.0f ? 1 : 0;

  // Calculate the footprint of the triangle on the grid's y-axis
  int y0 = (int)((tmin.y - m_bbMin.y)*ics);
  int y1 = (int)((tmax.y - m_bbMin.y)*ics);
  y0 = Clamp(y0, 0, h-1);
  y1 = Clamp(y1, 0, h-1);

  // Clip the triangle into all grid cells it touches.
  CVertexUV buf[7*4];
  CVertexUV *in = buf, *inrow = buf+7, *p1 = inrow+7, *p2 = p1+7;
  
  in[0] = in_v0;
  in[1] = in_v1;
  in[2] = in_v2;
  int nvrow, nvIn = 3;

  for (int y = y0; y <= y1; ++y)
  {
    // Clip polygon to row. Store the remaining polygon as well
    const float cz = m_bbMin[1] + y*cs;
    DividePoly(in, nvIn, inrow, &nvrow, p1, &nvIn, cz+cs, 1);
    SwapT(in, p1);
    if (nvrow < 3) continue;

    // find the horizontal bounds in the row
    float minX = inrow[0].m_vtx.x, maxX = inrow[0].m_vtx.x;
    for (int i=1; i<nvrow; ++i)
    {
      if (minX > inrow[i].m_vtx.x)	minX = inrow[i].m_vtx.x;
      if (maxX < inrow[i].m_vtx.x)	maxX = inrow[i].m_vtx.x;
    }
    int x0 = (int)((minX - m_bbMin[0])*ics);
    int x1 = (int)((maxX - m_bbMin[0])*ics);
    x0 = Clamp(x0, 0, w-1);
    x1 = Clamp(x1, 0, w-1);

    int nv, nv2 = nvrow;

    for (int x = x0; x <= x1; ++x)
    {
      // Clip polygon to column. store the remaining polygon as well
      const float cx = m_bbMin[0] + x*cs;
      DividePoly(inrow, nv2, p1, &nv, p2, &nv2, cx+cs, 0);
      SwapT(inrow, p2);
      if (nv < 3) continue;

      // Calculate min and max of the span.
      float smin = p1[0].m_vtx.z, smax = p1[0].m_vtx.z;
      for (int i = 1; i < nv; ++i)
      {
        smin = MinT(smin, p1[i].m_vtx.z);
        smax = MaxT(smax, p1[i].m_vtx.z);
      }
      smin -= m_bbMin[2];
      smax -= m_bbMin[2];
      // Skip the span if it is outside the heightfield bbox
      if (smax < 0.0f) continue;
      if (smin > bz) continue;
      // Clamp the span to the heightfield bbox.
      if (smin < 0.0f) smin = 0;
      if (smax > bz) smax = bz;

      // Snap the span to the heightfield height grid.
      unsigned int ismin = (unsigned short)Clamp((int)floorf(smin * ich), 0, SPAN_MAX_HEIGHT);
      unsigned int ismax = (unsigned short)Clamp((int)ceilf(smax * ich), (int)ismin+1, SPAN_MAX_HEIGHT);

      // color from the center of the cell
      const float l_dcx = (float)(1.0/ 3.0);
      CVec3 l_cen = p1[0].m_uv; // (in_v0.m_uv + in_v1.m_uv + in_v2.m_uv) * l_dcx;
      int l_cenx = int(l_cen.x * l_imgW);
      int l_ceny = int(l_cen.y * l_imgH);
      int l_area = in_img[l_cenx + l_ceny * in_imgWidth];
      if (l_area == NULL_INFO)
        l_area = 1;

      if (!AddSpan(x, y, ismin, ismax, l_area, backFace, 0))
        return false;
    }
  }

  return true;
}



///////////////////////////////////////////////////////////////////////////////
//
//                           tpcl::CCompactHeightMap
//
///////////////////////////////////////////////////////////////////////////////



/** construct from a dynamic height map */
CCompactHeightMap::CCompactHeightMap(const CDynamicHeightMap& in_dynHmap, float in_walkableHeight, float in_walkableClimb, bool in_backCull)
  : CGrid2D<CCompactCell>(in_dynHmap.GetWidth(), in_dynHmap.GetHeight(), in_dynHmap.GetBBoxMin(), in_dynHmap.GetRes())
{
  // creation of parent failed?
  if (m_data == 0)
  {
    m_spans = 0;
    m_areas = 0;
    return;
  }

  const int w = m_width;
  const int h = m_height;
  const int NOT_CONNECTED = CCompactSpan::NOT_CONNECTED;
  m_resH = in_dynHmap.GetHeightRes();
  m_invResH = 1.0f / m_resH;

  m_bbMax = in_dynHmap.GetBBoxMax();
  m_spanCount = in_dynHmap.GetSpanCount();
  m_spans = new CCompactSpan[m_spanCount];
  m_areas = new unsigned char[m_spanCount];
  m_colors = in_dynHmap.HasColorInfo() ? new unsigned short[m_spanCount] : 0;
  bool l_colorCreationFailed = in_dynHmap.HasColorInfo() ? (m_colors==0) : false;
  if ( !m_spans || !m_areas || l_colorCreationFailed )
  {
    Log("Could not create compact height map size %d x %d. Out of memory?\n", m_width, m_height);
    return;
  }
  CCompactCell* cells = (CCompactCell*)m_data;
  memset(cells, 0, sizeof(CCompactCell)*m_width*m_height);
  memset(m_spans, 0, sizeof(CCompactSpan)*m_spanCount);
  memset(m_areas, NULL_AREA, sizeof(unsigned char)*m_spanCount);
  if (m_colors)
    memset(m_colors, 0, sizeof(unsigned short)*m_spanCount);

  const int MAX_HEIGHT = 0xffff;

  // Fill in cells and spans.
  int idx = 0;
  for (int y = 0; y < h; ++y)
  {
    for (int x = 0; x < w; ++x)
    {
      const CDynamicHeightMap::CSpan* s = in_dynHmap.Get(x,y);
      // If there are no spans at this cell, just leave the data to index=0, count=0.
      if (!s) continue;
      CCompactCell& c = cells[x+y*w];
      c.index = idx;
      c.count = 0;
      int lastTop = 0;
      for ( ; s!=0 ; s = s->m_next)
      {
        if (s->info == CDynamicHeightMap::NULL_INFO)
        {
          lastTop = (int)s->smax;
          continue;
        }
        if (in_backCull && s->facing==1)
        {
          lastTop = (int)s->smax;
          continue;
        }

        if (m_colors)
        {
          m_colors[idx] = s->info & 0xffff; // it is assumed that we have HSV packed in 16bits
          m_areas[idx] = 1;   // the area bits are filled later.
        }
        else
          m_areas[idx] = s->info & 0x7f;  // area information

        if (s->facing == 0)
        {
          const int bot = (int)s->smax;
          const int top = s->m_next ? (int)s->m_next->smin : MAX_HEIGHT;
          m_spans[idx].y = (unsigned short)Clamp(bot, 0, 0xffff);
          m_spans[idx].h = (unsigned char)Clamp(top - bot, 0, 0xff);
          m_spans[idx].hv = (unsigned char)Clamp(s->smax - s->smin, 0, 0xff);
        }
        else
        {
          const int bot = (int)s->smin;
          const int top = lastTop;
          m_spans[idx].y = (unsigned short)Clamp(bot, 0, 0xffff);
          m_spans[idx].h = (unsigned char)Clamp(bot - lastTop, 0, 0xff);
          m_spans[idx].hv = (unsigned char)Clamp(s->smax - s->smin, 0, 0xff);
          if (m_areas[idx] != 0) // 0 area is kept 0 regarless of facing
            m_areas[idx] |= BACK_FACE;
        }
        idx++;
        c.count++;
        lastTop = (int)s->smax;
      }
    }
  }

  ComputeSurfaces(in_walkableHeight, in_walkableClimb);
  Log("Created compact height map size %d x %d.\n", m_width, m_height);
}


void CCompactHeightMap::ComputeSurfaces(float in_walkableHeight, float in_walkableClimb)
{
  const int w = m_width;
  const int h = m_height;
  const int NOT_CONNECTED = CCompactSpan::NOT_CONNECTED;
  CCompactCell* cells = (CCompactCell*)m_data;

  // TODO: move this outside
  int walkableClimb = (int)(in_walkableClimb * m_invResH);
  int walkableHeight = (int)(in_walkableHeight * m_invResH);


  // build surfaces by finding neighbour connections.
  const int MAX_LAYERS = NOT_CONNECTED-1;
  int tooHighNeighbour = 0;
  for (int y = 0; y < h; ++y)
  {
    for (int x = 0; x < w; ++x)
    {
      const CCompactCell& c = cells[x+y*w];
      for (int i = (int)c.index, ni = (int)(c.index+c.count); i < ni; ++i)
      {
        CCompactSpan& s = m_spans[i];
        int facing = m_areas[i] & BACK_FACE;

        for (int dir = 0; dir < 4; ++dir)
        {
          s.SetCon(dir, NOT_CONNECTED);
          const int nx = x + GetDirOffsetX(dir);
          const int ny = y + GetDirOffsetY(dir);
          // First check that the neighbour cell is in bounds.
          if (nx < 0 || ny < 0 || nx >= w || ny >= h)
            continue;
        
          // Iterate over all neighbour spans and check if any of them is
          // accessible from current cell.
          const CCompactCell& nc = cells[nx+ny*w];
          for (int k = (int)nc.index, nk = (int)(nc.index+nc.count); k < nk; ++k)
          {
            if ( (m_areas[k] & BACK_FACE) != facing)
              continue;
            const CCompactSpan& ns = m_spans[k];

            const int bot = MaxT(s.y, ns.y);
            const int top = MinT(s.y+s.h, ns.y+ns.h);

            // Check that the gap between the spans is walkable,
            // and that the climb height between the gaps is not too high.
            if ((top - bot) >= walkableHeight && Abs((int)ns.y - (int)s.y) <= walkableClimb)
            {
              // Mark direction as walkable.
              const int lidx = k - (int)nc.index;
              if (lidx < 0 || lidx > MAX_LAYERS)
              {
                tooHighNeighbour = MaxT(tooHighNeighbour, lidx);
                continue;
              }
              s.SetCon(dir, lidx);
              break;
            }
          }
        	
        }
      }
    }
  }
	
  if (tooHighNeighbour > MAX_LAYERS)
  {
    Log("Building of Compact height map failed: too many layers %d (max: %d)", tooHighNeighbour, MAX_LAYERS);
  }
}



CCompactHeightMap::~CCompactHeightMap()
{
  if (m_spans)
    delete[] m_spans;
  if (m_areas)
    delete[] m_areas;
  if (m_colors)
    delete[] m_colors;
}



int CCompactHeightMap::GetHeightAt(const CVec3& in_pos, int in_bufSize, float* out_heights) const
{
  CCompactCell l_cell = Get(in_pos);
  int l_n = 0;
  for (int i=l_cell.index+l_cell.count-1, ni=l_cell.index; i>=ni; --i)
  {
    if (l_n >= in_bufSize)
      break;
    CCompactSpan& l_span = m_spans[i];
    float z = l_span.y * m_resH + m_bbMin.z;
    out_heights[l_n++] = z;
  }
  return l_n;
}



int CCompactHeightMap::GetTopSpanId(const CVec3& in_pos) const
{
  CCompactCell l_cell = Get(in_pos);
  if (l_cell.count == 0)
    return -1;
  return l_cell.index + l_cell.count - 1;
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
}
