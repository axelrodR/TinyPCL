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
*: Package Name: dbg
*
******************************************************************************/
#include "dbg.h"
#include "common.h"
#include <string>


//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


namespace tpcl
{

  /******************************************************************************
  *                             INTERNAL CONSTANTS  / Functions                 *
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

  // 24 bit blue, green, red used in BMPs
  struct BGR { unsigned char B, G, R; };

  /** bitmap (.BMP) file header */
#pragma pack(push,2)
  struct BITMAPFILEHEADER
  {
    unsigned short  bfType;
    unsigned int    bfSize;
    unsigned short  bfReserved1;
    unsigned short  bfReserved2;
    unsigned int    bfOffBits;
  };

  /** bitmap inforrmation (comes after BITMAPFILEHEADER in the BMP file */
  struct BITMAPINFOHEADER
  {
    unsigned int    biSize;
    int             biWidth;
    int             biHeight;
    unsigned short  biPlanes;
    unsigned short  biBitCount;
    unsigned int  biCompression;
    unsigned int  biSizeImage;
    int           biXPelsPerMeter;
    int           biYPelsPerMeter;
    unsigned int  biClrUsed;
    unsigned int  biClrImportant;
  };
#pragma pack(pop)



  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CRegDebug
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: SLDRCR_DBG_CRegDebug
  *
  ******************************************************************************/
  CRegDebug::CRegDebug()
  {
    m_Fpath = NULL;
  }

  CRegDebug::CRegDebug(char* Xi_Fpath)
  {
    m_Fpath = NULL;
    setPath(Xi_Fpath);
  }

  /******************************************************************************
  *
  *: Method name: ~SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CRegDebug::~CRegDebug()
  {
    delete[] m_Fpath;
  }


  /******************************************************************************
  *
  *: Method name: ~SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  void CRegDebug::setPath(char* Xi_Fpath)
  {
    delete[] m_Fpath;

    m_Fpath = new char[strlen(Xi_Fpath) + 3];
    strcpy(m_Fpath, Xi_Fpath);
    strcat(m_Fpath, "\\");
  }


  /******************************************************************************
  *
  *: Method name: SaveAsBmp
  *
  ******************************************************************************/
  /** save to disk as a BMP */
  void CRegDebug::SaveAsBmp(char* Xi_Fname, float* Xi_img, int Xi_Width, int Xi_Height, float minVal, float maxVal)
  {
    if (m_Fpath == NULL)
      return;

    char *l_FileFullPath = new char[strlen(m_Fpath) + strlen(Xi_Fname) + 1];
    strcpy(l_FileFullPath, m_Fpath);
    strcat(l_FileFullPath, Xi_Fname);


    ////IFRLOG_MSG_Assert(l_FileFullPath != NULL, "CDebugImg::SaveAsBMP: missing file name");
    FILE* l_pFile = fopen(l_FileFullPath, "wb");
    if (!l_pFile)
    {
      ////Log("CDebugImg::SaveAsBMP: file %s not created", l_FileFullPath);
      return;
    }

    BITMAPFILEHEADER l_BmpFileHdr;
    l_BmpFileHdr.bfType = 0x4d42;
    //int a = sizeof(BITMAPFILEHEADER);
    //int b = sizeof(BITMAPINFOHEADER);
    l_BmpFileHdr.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
    //l_BmpFileHdr.bfOffBits += sizeof(CGeoSimExtraInfo);

    BITMAPINFOHEADER l_BmpInfoHdr;
    l_BmpInfoHdr.biSize = sizeof(BITMAPINFOHEADER);
    l_BmpInfoHdr.biWidth = Xi_Width;
    l_BmpInfoHdr.biHeight = Xi_Height;
    l_BmpInfoHdr.biPlanes = 1;
    l_BmpInfoHdr.biBitCount = 24;
    l_BmpInfoHdr.biCompression = 0L;
    l_BmpInfoHdr.biSizeImage = 0;//m_Width*m_Height;
    l_BmpInfoHdr.biClrUsed = 0;
    l_BmpInfoHdr.biClrImportant = 0;
    l_BmpInfoHdr.biXPelsPerMeter = 1;
    l_BmpInfoHdr.biYPelsPerMeter = 1;

    // write headers
    fwrite(&l_BmpFileHdr, sizeof(BITMAPFILEHEADER), 1, l_pFile);
    fwrite(&l_BmpInfoHdr, sizeof(BITMAPINFOHEADER), 1, l_pFile);
    //fwrite(&l_geoSimExt, sizeof(CGeoSimExtraInfo), 1, l_pFile);

    //find parameters for float to unsigned char conversion:
    float bRange[2] = { minVal , maxVal };
    float multVal = 1.f / (bRange[1] - bRange[0]);

    //float bRange[2] = { 0 , 0 };
    //float minVal = 0;
    //float multVal = 0;

    //int totalSize = Xi_Width * Xi_Height;
    //int indexRange = 0;
    //for (indexRange; indexRange < totalSize; indexRange++)
    //{
    //  if (Xi_img[indexRange] != 0)
    //    break;
    //}

    //if (indexRange != totalSize)
    //{
    //  bRange[0] = bRange[1] = Xi_img[indexRange];
    //  for (indexRange; indexRange < totalSize; indexRange++)
    //  {
    //    if (Xi_img[indexRange] != 0)
    //    {
    //      bRange[0] = min(bRange[0], Xi_img[indexRange]);
    //      bRange[1] = max(bRange[1], Xi_img[indexRange]);
    //    }
    //  }
    //minVal = bRange[0];
    //multVal = float(UCHAR_MAX) / (bRange[1] - bRange[0]);
    //}

    // write data
    int l_bmpRowWidth = (l_BmpInfoHdr.biBitCount*Xi_Width + 31) / 32 * 4;
    unsigned char *l_row = new unsigned char[l_bmpRowWidth];
    BGR* l_RGBrow = new BGR[Xi_Width];

    for (int y = Xi_Height - 1; y >= 0; --y)
    {
      //float to RGB:
      for (int x = 0; x < Xi_Width; x++)
      {
        int index = y*Xi_Width + x;
        //unsigned char val = unsigned char(round((Xi_img[index] - minVal) * multVal));
        
        float clampedValNormalized = (MaxT(MinT(Xi_img[index], bRange[1]), bRange[0]) - bRange[0]) * multVal;
        //clampedValNormalized = 1.f - pow(1.f - clampedValNormalized, 2);
        if (clampedValNormalized <= 0.33)
        {
          l_RGBrow[x].B = l_RGBrow[x].G = 0;
          l_RGBrow[x].R = unsigned char(round((clampedValNormalized/ 0.33)*255));
        }
        else if (clampedValNormalized <= 0.66)
        {
          l_RGBrow[x].B = 0;
          l_RGBrow[x].G = 180;
          l_RGBrow[x].R = 255 - unsigned char(round(((clampedValNormalized - 0.33) / 0.33) * 255));
        }
        else
        {
          l_RGBrow[x].B = unsigned char(round(((clampedValNormalized - 0.66) / 0.34) * 255));
          l_RGBrow[x].G = 180;
          l_RGBrow[x].R = 0;
        }


        ////////////////////////////////////////
        //unsigned char val = unsigned char(round(clampedValNormalized * 255));
        //if (val < 90)
        //{
        //  l_RGBrow[x].B = l_RGBrow[x].G = 0;
        //  l_RGBrow[x].R = val;
        //}
        //else if (val < 90)
        //{
        //  l_RGBrow[x].B = 0;
        //  l_RGBrow[x].G = 90;
        //  l_RGBrow[x].R = val;
        //}
        //else
        //{
        //  l_RGBrow[x].B = 0;
        //  l_RGBrow[x].G = 180;
        //  l_RGBrow[x].R = val;
        //}
        //l_RGBrow[x].B = l_RGBrow[x].G = l_RGBrow[x].R = val;
        ////////////////////////////////////////
      }

      //write:
      memcpy(l_row, l_RGBrow, Xi_Width * sizeof(BGR));
      fwrite(l_row, l_bmpRowWidth, 1, l_pFile);
    }
    delete[] l_row;
    delete[] l_RGBrow;
    fclose(l_pFile);
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


 

} //namespace SLDR