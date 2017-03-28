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


#include "SpatialHash.h"
#include "common.h"
#include <unordered_map>
#include <vector>


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

namespace tpcl{

/******************************************************************************
*                             INTERNAL CONSTANTS                              *
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

typedef TVec3<int> CInt3;


// hash function for 2D as needed by hash_map
class HashCompare2D
{
public:
  static const size_t bucket_size = 4;
  static const size_t min_buckets = 8;

  size_t hash( const CInt3& n) const  {size_t l_Val = (2166136261U * n.x) ^ (16777619U * n.y); return l_Val;}
  size_t operator() (const CInt3& v) const  {return hash(v);}
  bool operator() ( const CInt3& lhs, const CInt3& rhs) const 
  {
    return lhs.x==rhs.x && lhs.y == rhs.y;
  }
};


// internal node for 2D hashing containing both position and object pointer
struct Node2D {void* obj; CVec3 pt;};

// internal node for 1D hashimh containing both position and object pointer
struct Node1D {void* obj; float pt;};


/** Each cell is hashed to an array of object+position pairs */
typedef std::unordered_map<CInt3, std::vector<Node2D> , HashCompare2D, HashCompare2D> MapInt3;


/******************************************************************************
*                           EXPORTED CLASS METHODS                            *
******************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//
//                           GENGMTRX_SPAT_CSpatialHash2D
//
///////////////////////////////////////////////////////////////////////////////
/******************************************************************************
*                               Public methods                                *
******************************************************************************/
/******************************************************************************
*
*: Method name: CSpatialHash2D
*
******************************************************************************/
CSpatialHash2D::CSpatialHash2D (float res)
{
  m_data = new MapInt3();
  m_res = res;
  m_resInv = (float)(1.0 / m_res);
}

/******************************************************************************
*
*: Method name: ~CSpatialHash2D
*
******************************************************************************/
CSpatialHash2D::~CSpatialHash2D ()
{
  delete ((MapInt3*)m_data);
}


/******************************************************************************
*
*: Method name: Add
*
******************************************************************************/
void CSpatialHash2D::Add(const CVec3& Xi_pos, void* Xi_obj)
{
  MapInt3* l_data = (MapInt3*)m_data;
  if (l_data->size() == 0)
    m_pivot = Xi_pos;

  // convert coordinates to cell coordinates
  CVec3 l_v = (Xi_pos - m_pivot) * m_resInv;
  CInt3 l_cell = CInt3(l_v.x, l_v.y, l_v.z);

  // add to cell
  MapInt3::iterator l_it = l_data->find(l_cell);
  if (l_it == l_data->end()) // empty cell? create list
    (*l_data)[l_cell] = std::vector<Node2D>();
  Node2D n;
  n.obj = Xi_obj;
  n.pt = Xi_pos;
  (*l_data)[l_cell].push_back(n);
}


/******************************************************************************
*
*: Method name: FindNearest
*
* Find nearest object (using its associated point)
* @param Xo_pMinPt (optional) closest point
* @param max2DRadius maximum 2D radius to search in
******************************************************************************/

void* CSpatialHash2D::FindNearest(const CVec3& Xi_pos, CVec3* Xo_pMinPt, float Xi_max2DRadius) const
{
  MapInt3* l_data = (MapInt3*)m_data;
  const void* l_minObj = 0;
  const CVec3* l_minPt = 0;
  float l_minDist = 1E10;
  float s_epsilon = 0.01f * m_res;
  if (Xi_max2DRadius < s_epsilon)
    Xi_max2DRadius  = s_epsilon;
  float max2dRadSqr = Xi_max2DRadius * Xi_max2DRadius;
  // convert coordinates to cell coordinates
  CVec3 l_v = (Xi_pos - m_pivot) * m_resInv;
  CInt3 l_cell = CInt3(l_v.x, l_v.y, l_v.z);
  int rad = int(ceil(Xi_max2DRadius * m_resInv));

  // go over all cells and objects in them to find the minimal distance
  for (int x=-rad; x<=rad; x++)
  {
    for (int y=-rad; y<=rad; y++)
    {
      CInt3 l_c = CInt3(l_cell.x + x, l_cell.y + y, 0);
      MapInt3::const_iterator l_it = l_data->find(l_c);
      if (l_it == l_data->end()) 
        continue; // empty cell
      const std::vector<Node2D>& nodes = l_it->second;
      for (unsigned int i=0; i<nodes.size(); i++)
      {
        float l_dist2DSqr = DistSqr2D(nodes[i].pt, Xi_pos);
        if (l_dist2DSqr > max2dRadSqr)
          continue;
        float l_dist = DistSqr(nodes[i].pt, Xi_pos);
        if (l_dist >= l_minDist)
          continue;
        l_minDist = l_dist;
        l_minObj  = nodes[i].obj;
        l_minPt = &(nodes[i].pt);
      }
    }
  }

  if (Xo_pMinPt != 0)
    *Xo_pMinPt = l_minPt == 0 ? CVec3(0,0,0) : *l_minPt;
  return const_cast<void*>(l_minObj);
}



/******************************************************************************
*
*: Method name: GetNear
*
******************************************************************************/
int CSpatialHash2D::GetNear(const CVec3& Xi_pos, int Xi_bufSize, void** Xo_buf, CVec3* Xo_pos, float Xi_max2DRadius) const
{
  MapInt3* l_data = (MapInt3*)m_data;
  int l_n = 0;    // number of objects
  float l_minDist = 1E10;
  float s_epsilon = 0.01f * m_res;
  if (Xi_max2DRadius < s_epsilon)
    Xi_max2DRadius  = s_epsilon;
  float max2dRadSqr = Xi_max2DRadius * Xi_max2DRadius;
  // convert coordinates to cell coordinates
  CVec3 l_v = (Xi_pos - m_pivot) * m_resInv;
  CInt3 l_cell = CInt3(l_v.x, l_v.y, l_v.z);
  int rad = int(ceil(Xi_max2DRadius * m_resInv));

  // go over all cells and objects in them to find the minimal distance
  for (int x=-rad; x<=rad; x++)
  {
    for (int y=-rad; y<=rad; y++)
    {
      CInt3 l_c = CInt3(l_cell.x + x, l_cell.y + y, 0);
      MapInt3::const_iterator l_it = l_data->find(l_c);
      if (l_it == l_data->end()) 
        continue;   // empty cell
      const std::vector<Node2D>& nodes = l_it->second;
      for (unsigned int i=0; i<nodes.size(); i++)
      {
        float l_dist2DSqr = DistSqr2D(nodes[i].pt, Xi_pos);
        if (l_dist2DSqr > max2dRadSqr)
          continue;
        Xo_buf[l_n] = nodes[i].obj;
        Xo_pos[l_n++] = nodes[i].pt;
        if (l_n >= Xi_bufSize)
          return Xi_bufSize;
      }
    }
  }
  return l_n;
}


/******************************************************************************
*
*: Method name: Clear data
*
******************************************************************************/
void CSpatialHash2D::Clear()
{
  MapInt3* l_data = (MapInt3*)m_data;
  l_data->clear();
}



} // nsmaespace GenGmtrx
