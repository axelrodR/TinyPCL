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
*: Package Name: sldrcr_dbg
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_dbg_H
#define __tpcl_dbg_H

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/


  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

namespace tpcl
{
  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: SLDR_DBG_CRegDebug
  *
  *: Abstract: debug tools for coarse registeration based on dictionary
  *
  ******************************************************************************/

  class CRegDebug
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor */
    CRegDebug();
    CRegDebug(char* in_Fpath);

    /** destructor */
    virtual ~CRegDebug(); 


    /** set output files path.
    * @param in_Fpath        output files path */
    void setPath(char* in_Fpath);

    /** save float image as bmp.
    * @param in_Fpath       output files name 
    * @param in_img         the float image */
    void SaveAsBmp(char* in_Fname, float* in_img, int in_Width, int in_Height, float minVal, float maxVal);
 

  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/


    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/


  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/
    char* m_Fpath;             ///< output files path.

   
    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/

  };



} // namespace SLDR

#endif

