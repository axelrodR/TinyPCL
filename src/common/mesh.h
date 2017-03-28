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
*: Package Name: tpcl_mesh
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_mesh_H
#define __tpcl_mesh_H

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

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/

  ///////////////////////////////////////////////////////////////////////////////
  //
  //                            tpcl::CVertexUV
  //
  // a vertex with UV coordinates
  ///////////////////////////////////////////////////////////////////////////////
  struct CVertexUV
  {
    CVec3 m_vtx;      ///< vertex
    CVec3 m_uv;       ///< uv
    CVertexUV(){}
    CVertexUV(const CVec3& vtx, const CVec3& uv){m_vtx=vtx,m_uv=uv;}
  };


  ///////////////////////////////////////////////////////////////////////////////
  //
  //                            tpcl::CMesh
  //
  // Simple mesh container used mainly for interfacing between modules
  ///////////////////////////////////////////////////////////////////////////////
  struct CMesh
  {
    union
    {
      CVertexUV* m_vertsUV;       ///< vertices in XYZUV format
      CVec3* m_verts;       ///< vertices in XYZ format
    };
    int* m_faces;                 ///< indices (triangles)
    int m_numVerts;               ///< number of vertices
    int m_numFaces;               ///< number of faces

    //flags
    unsigned char m_XYZUV: 1;     ///< 1=vertices contain UVs, 0=vertices do not contain UVs
    unsigned char m_managed: 1;   ///< 1=memory is managed internally (i.e. destructor will delete arrays)

    /** constructor */
    CMesh ();
    CMesh (int Xi_numVtx, int Xi_numFaces, bool m_XYZUV);

    /** destructor: DOES NOT delete the data when deleted. */
    virtual ~CMesh ();
  };


  ///////////////////////////////////////////////////////////////////////////////
  //
  //                            tpcl::IRasterizable
  //
  // Provides interface for getting mesh 
  ///////////////////////////////////////////////////////////////////////////////

  class IRasterizable
  {
  public:
    /** Get vertices with texture */
    virtual const CMesh* U_GetMesh() = 0;

    /** Get bounding box */
    virtual void U_GetBBox(CVec3& Xo_min, CVec3& Xo_max) const = 0;

    /** get the texture / associated image */
    virtual const unsigned int* U_GetTexture(int& Xo_texWidth, int& Xo_texHeight) = 0;
  };


  /******************************************************************************
  *                            EXPORTED FUNCTIONS                               *
  ******************************************************************************/

} // namespace tpcl


#endif