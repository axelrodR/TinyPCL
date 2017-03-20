// File Location: S:\gen\gengmtrx\gengmtrx_grid.h
/******************************************************************************
*
*: Package Name: gengmtrx_grid
*
*: Title:
*
******************************************************************************/

#ifndef __gengmtrx_grid_H
#define __gengmtrx_grid_H

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

#pragma warning (disable : 4251)


namespace GenGmtrx
{
/******************************************************************************
*                             EXPORTED CONSTANTS                              *
******************************************************************************/

/******************************************************************************
*                        INCOMPLETE CLASS DECLARATIONS                        *
******************************************************************************/
class CInt3;
struct CVertexUV;

/******************************************************************************
*                              EXPORTED CLASSES                               *
******************************************************************************/

  /******************************************************************************
  *
  *: Class name: CGrid2dBase
  *
  *: Abstract: baseline of all grids. Use CGrid2D template instead
  ******************************************************************************/
  /** base for 2D grids. Use CGrid2D template instead */
class DLL_Entry CGrid2dBase : public IFRGEN_CBaseObject
  {
  public:
    /******************************************************************************
    *                               Public definitions                            *
    ******************************************************************************/

    /* internal data access */
    int GetWidth() const                            {return m_width;} ///< cells along x direction
    int GetHeight() const                           {return m_height;}///< cells along y direction
    const float GetRes() const                      {return m_res;}   ///< resolution (cell size)
    const int GetStride() const                     {return m_stride;}///< stride: bytes per cell
    const D3DXVECTOR3& GetBBoxMin() const           {return m_bbMin;} ///< min corner of bounding box
    const D3DXVECTOR3& GetBBoxMax() const           {return m_bbMax;} ///< max corner of bounding box

    /** get cell using 2D cell coordiantes */
    void* Get(int Xi_x, int Xi_y)                   {return (void*)Get(GetIndex(Xi_x,Xi_y));}
    const void* Get(int Xi_x, int Xi_y) const       {return (void*)Get(GetIndex(Xi_x,Xi_y));}

    /** get cell using world coordinates */
    void* Get(const D3DXVECTOR3& Xi_pos);
    const void* Get(const D3DXVECTOR3& Xi_pos) const;

    /** get the world position associated with a cell coordinates (corner of a cell) */
    D3DXVECTOR3 GetPos(int Xi_x, int Xi_y) const    {return m_bbMin + D3DXVECTOR3((float)Xi_x,(float)Xi_y,0)*m_res;}

    /** get 1D array index using cell coordinates */
    int GetIndex(int Xi_x, int Xi_y) const          {return Xi_x + Xi_y*m_width;}
    int GetIndex(const D3DXVECTOR3& Xi_pos) const; ///< get index using world position

    /** get cell using 1D array index */
    void* Get(int Xi_index)                         {return (void*)(m_data + (Xi_index*m_stride));}
    const void* Get(int Xi_index) const             {return (void*)(m_data + (Xi_index*m_stride));}
    
    /** convert world position to cell coordinates (changes Xi_cellCoord and returns it) */
    CInt3& Convert(const D3DXVECTOR3& Xi_pos, CInt3& Xi_cellCoord) const;
    CInt3& ConvertSafe(const D3DXVECTOR3& Xi_pos, CInt3& Xi_cellCoord) const; ///< convert position to cell coords with safety checks


    /** translate and scale the bounding box (changes resolution along the way) 
     * @param Xi_shift shift of the min box corner.
     */
    void Transform(const D3DXVECTOR3& Xi_shift, float scale=1.0f);

    /** clear contents to 0 */
    void Clear();


    /** Paint/rasterize a triangle to the grid (for height map the height is written) */
    virtual bool U_RasterizeTri(const D3DXVECTOR3& Xi_v0, const D3DXVECTOR3& Xi_v1, const D3DXVECTOR3& Xi_v2, int Xi_info);

    /** Paint/rasterize a triangle to the grid with info in "texture" image
     *  (for height map the height is written) */
    virtual bool U_RasterizeTri(const CVertexUV& Xi_v0, const CVertexUV& Xi_v1, const CVertexUV& Xi_v2,
                                int Xi_imgWidth, int Xi_imgHeight, const unsigned int* Xi_img);


    /** constructor */
    CGrid2dBase(int Xi_width, int Xi_height, const D3DXVECTOR3& Xi_bbmin, float Xi_res=1.0f, int Xi_stride=4);
    
    /** @brief constrctor 
     *
     * The width and height of the bounding box are computed from the bounding box corners.
     * The max corner is changed so it is exactly at the end of a grid cell
     */
    CGrid2dBase(const D3DXVECTOR3& Xi_bbMin, const D3DXVECTOR3& Xi_bbMax, float Xi_res=1.0f, int Xi_stride=4);
    
    /** destructor */
    virtual ~CGrid2dBase();

    /** @brief check if the array is allocated 
     * Inheriting classes may choose to differ creation or allocation might have failed.
     */
    bool IsDataAllocated()    { return m_data != 0;}

  protected:
    /******************************************************************************
    *                               proteced definitions                          *
    ******************************************************************************/

    char* m_data;           ///< data array
    int m_width;			      ///< The width of the heightfield. (Along the x-axis in cell units.)
	  int m_height;			      ///< The height of the heightfield. (Along the z-axis in cell units.)
    int m_stride;           ///< bytes per pixel
    int m_lineSize;         ///< line size in bytes
    float m_res, m_invRes;  ///< resolution: the size of each cell (m_invRes=1.0/m_res)
    D3DXVECTOR3 m_bbMin, m_bbMax; ///< bounding box of grid
  };




  /******************************************************************************
  *: Abstract: wrapper template for 2D grids
  ******************************************************************************/
  /** wrapper template for 2D grids */
  template <typename T> class CGrid2D : public CGrid2dBase
  {
  public:
    /** constructor */
    CGrid2D(int width, int height, const D3DXVECTOR3& bbmin, float res=1.0f)
      : CGrid2dBase(width, height, bbmin, res, sizeof(T)){}
    
    /** constructor */
    CGrid2D(const D3DXVECTOR3& bbMin, const D3DXVECTOR3& bbMax, float res=1.0f)
      : CGrid2dBase(bbMin, bbMax, res, sizeof(T)){}

    // wrapper to access cells 
    T& Get(int Xi_x, int Xi_y)                  {return *(T*)CGrid2dBase::Get(Xi_x,Xi_y);}    ///< the cell using cell coordinates
    T& Get(int Xi_ind)                          {return *(T*)CGrid2dBase::Get(Xi_ind);}       ///< the cell using an index
    T& Get(const D3DXVECTOR3& Xi_pos)           {return *(T*)CGrid2dBase::Get(Xi_pos);}       ///< get the cell at a world pos
    const T& Get(int Xi_x, int Xi_y) const      {return *(T*)CGrid2dBase::Get(Xi_x,Xi_y);}    ///< the cell using cell coordinates
    const T& Get(int Xi_ind) const              {return *(T*)CGrid2dBase::Get(Xi_ind);}       ///< the cell using an index
    const T& Get(const D3DXVECTOR3& Xi_pos) const {return *(T*)CGrid2dBase::Get(Xi_pos);}       ///< get the cell at a world pos
  };

} // namespace GenGmtrx

/******************************************************************************
*                            EXPORTED FUNCTIONS                               *
******************************************************************************/

  /// Gets the standard width (x-axis) offset for the specified direction.
  inline int GetDirOffsetX(int dir)
  {
	  static const int offset[4] = { -1, 0, 1, 0, };
	  return offset[dir&0x03];
  }

  /// Gets the standard height (y-axis) offset for the specified direction.
  inline int GetDirOffsetY(int dir)
  {
	  static const int offset[4] = { 0, 1, 0, -1 };
	  return offset[dir&0x03];
  }



#undef DLL_Entry
#endif

#pragma warning (default : 4251)

/******************************************************************************
*                    inline function implementation                           *
******************************************************************************/


/******************************************************************************
*                    OLD STYLE TYPEDEFS                                       *
******************************************************************************/


typedef GenGmtrx::CGrid2dBase GENGMTRX_RAST_CGrid2dBase;
