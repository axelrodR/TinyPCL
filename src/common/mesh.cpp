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
*: Package Name: tpcl_mesh
*
******************************************************************************/
//#include <d3dx9core.h>
//#include <ifr\ifrgen\ifrgen_stnd.h>
#include "vec.h"
#include "mesh.h"


#ifdef _DEBUG
#define new_file DEBUG_NEW
#endif

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

  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           TPCL_MESH_CMesh
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: TPCL_MESH_CMesh
  *
  ******************************************************************************/
  CMesh::CMesh ()
  {
    m_vertsUV = 0;
    m_verts = 0;
    m_faces = 0;
    m_numVerts = 0;
    m_numFaces = 0;
    m_XYZUV = 0;
    m_managed = 0;
  }

  /******************************************************************************
  *
  *: Method name: TPCL_MESH_CMesh
  *
  ******************************************************************************/
  CMesh::CMesh (int Xi_numVtx, int Xi_numFaces, bool Xi_XYZUV)
  {
    if (Xi_XYZUV)
    {
      m_vertsUV = new CVertexUV[Xi_numVtx];
      m_XYZUV = 1;
    }
    else
    {
      m_verts = new CVec3[Xi_numVtx];
      m_XYZUV = 0;
    }
    m_faces = new int[Xi_numFaces*3];
    m_numVerts = Xi_numVtx;
    m_numFaces = Xi_numFaces;
    m_managed = 1;
  }


  /******************************************************************************
  *
  *: Method name: ~TPCL_MESH_CMesh
  *
  ******************************************************************************/
  CMesh::~CMesh ()
  {
    if (!m_managed)
      return;

    if (m_vertsUV != 0)
      delete m_vertsUV;
    m_vertsUV = 0;

    if (m_faces != 0)
      delete m_faces;
    m_faces = 0;
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
