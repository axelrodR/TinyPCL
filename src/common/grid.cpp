// File Location: S:\gen\gengmtrx\gengmtrx_grid.h

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
*: Package Name: tpcl_grid
*
******************************************************************************/
#include "vec.h"
//#include <ifr\ifrgen\ifrgen_stnd.h>
#include "grid.h"
#include "mesh.h"
#include "common.h"
//#include <ifr\ifrlog\ifrlog.h>
#include <stdio.h>
#include <string.h>

//#ifdef _DEBUG
//#define new_file DEBUG_NEW
//#endif


namespace tpcl
{
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
template<class T> inline T Clamp(T v, T mn, T mx) { return v < mn ? mn : (v > mx ? mx : v); }

void Log(...) {}

/******************************************************************************
*                           EXPORTED CLASS METHODS                            *
******************************************************************************/
///////////////////////////////////////////////////////////////////////////////
//
//                           TPCL_GRID_CGrid
//
///////////////////////////////////////////////////////////////////////////////
/******************************************************************************
*                               Public methods                                *
******************************************************************************/
/******************************************************************************
*
*: Method name: TPCL_GRID_CGrid
*
******************************************************************************/
CGrid2dBase::CGrid2dBase(const CVec3& in_bbMin, const CVec3& in_bbMax, float in_res, int in_stride)
{
  m_bbMin = in_bbMin;
  m_res = in_res;
  m_invRes = 1.0f / in_res;
  m_data = 0;

  if ( !IsValid(in_bbMin) || !IsValid(in_bbMin) )
  {
    Log("Error: bounding box is invalid (NaN or inf)\n");
    return;
  }
  m_width = (int)ceil((in_bbMax.x - in_bbMin.x) * m_invRes);
  m_height = (int)ceil((in_bbMax.y - in_bbMin.y) * m_invRes);
  if (m_width <= 0 || m_height <= 0)
  {
    Log("Error: box and resolution define a grid with no cells.\n");
    return;
  }
  m_stride = in_stride;
  m_lineSize = m_width * in_stride;
  m_data = new char[m_width*m_height*in_stride];
  Log("Trying to create grid from box (%4.2f x %4.2f x %4.2f)  <-->  (%4.2f x %4.2f x %4.2f)\n",
      in_bbMin.x, in_bbMin.y, in_bbMin.z, in_bbMax.x, in_bbMax.y, in_bbMax.z);
  if (!m_data)
    Log("Error: Could not create grid size %d x %d. Out of memory?\n", m_width, m_height);
  else
    Log("Created grid size %d x %d.\n", m_width, m_height);
  m_bbMax = in_bbMin + CVec3(m_width*in_res, m_height*in_res, 1000.0);
}


CGrid2dBase::CGrid2dBase(int in_width, int in_height, const CVec3& in_bbMin, float in_res, int in_stride)
{
  m_data = 0;
  m_bbMin = in_bbMin;
  m_res = in_res;
  m_invRes = 1.0f / in_res;

  m_width = in_width;
  m_height = in_height;
  m_stride = in_stride;
  m_lineSize = in_width * in_stride;

  // check that the bounding box is valid
  if ( !IsValid(in_bbMin) )
  {
    Log("Error: bounding box is invalid (NaN or inf)\n");
    return;
  }
  if (m_width <= 0 || m_height <= 0)
  {
    Log("Error: a grid with no cells.\n");
    return;
  }

  m_data = new char[in_width*in_height*in_stride];
  if (!m_data)
    Log("Could not create grid size %d x %d. Out of memory?\n", m_width, m_height);
  else
    Log("Created grid size %d x %d.\n", m_width, m_height);
  m_bbMax = in_bbMin + CVec3(in_width*in_res, in_height*in_res, 1000.0);
}


/******************************************************************************
*
*: Method name: ~TPCL_GRID_CGrid
*
******************************************************************************/
CGrid2dBase::~CGrid2dBase()
{
  if (m_data != 0)
    delete[] m_data;
  m_data = 0;
}


void CGrid2dBase::Clear()
{
  memset(m_data,0, m_width*m_height*m_stride);
}



/** get cell index using world coordinates */
int CGrid2dBase::GetIndex(const CVec3& in_pos) const
{
  CVec3 l_v = (in_pos - m_bbMin) * m_invRes;
  return GetIndex((int)l_v.x, (int)l_v.y);
}


/** get cell using world coordinates */
void* CGrid2dBase::Get(const CVec3& in_pos)
{
  CVec3 l_v = (in_pos - m_bbMin) * m_invRes;
  int l_x = (int)l_v.x;
  int l_y = (int)l_v.y;
  return Get(l_x, l_y);
}

/** get cell using world coordinates */
const void* CGrid2dBase::Get(const CVec3& in_pos) const
{
  CVec3 l_v = (in_pos - m_bbMin) * m_invRes;
  int l_x = (int)l_v.x;
  int l_y = (int)l_v.y;
  return Get(l_x, l_y);
}


void CGrid2dBase::Convert(const CVec3& in_pos, int& in_cellX, int& in_cellY) const
{
  CVec3 l_v = (in_pos - m_bbMin) * m_invRes;
  in_cellX = (int)l_v.x;
  in_cellY = (int)l_v.y;
}


void CGrid2dBase::ConvertSafe(const CVec3& in_pos, int& in_cellX, int& in_cellY) const
{
  CVec3 l_v = (in_pos - m_bbMin) * m_invRes;
  in_cellX = (int)l_v.x;
  in_cellY = (int)l_v.y;
  in_cellX = Clamp(in_cellX, 0, m_width);
  in_cellY = Clamp(in_cellY, 0, m_height);
}


void CGrid2dBase::Transform(const CVec3& in_shift, float scale)
{
  CVec3 l_extent = (m_bbMax - m_bbMin) * scale;
  m_bbMin += in_shift;
  m_bbMax = m_bbMin + l_extent;
}


bool CGrid2dBase::U_RasterizeTri(const CVec3& in_v0, const CVec3& in_v1, const CVec3& in_v2, int in_info)
{
  CVertexUV l_v0(in_v0, CVec3(0,0,0));
  CVertexUV l_v1(in_v1, CVec3(1,0,0));
  CVertexUV l_v2(in_v2, CVec3(1,1,0));
  return U_RasterizeTri(l_v0, l_v1, l_v2, 1, 1, (unsigned int*)&in_info);
}


bool CGrid2dBase::U_RasterizeTri(const CVertexUV& /*in_v0*/, const CVertexUV& /*in_v1*/, const CVertexUV& /*in_v2*/,
                               int /*in_imgWidth*/, int /*in_imgHeight*/, const unsigned int* /*in_img*/)
{
  return false;
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
} // namespace tpcl
