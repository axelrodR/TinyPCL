// File Location: 

//
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


#ifndef __tpcl_sp_H
#define __tpcl_sp_H

#include "vec.h"


namespace tpcl
{
  class  CRegDictionary;         // registration Dictionary
  struct CPtCloud;

  /******************************************************************************
  *
  *: Class name: IRegister
  *
  ******************************************************************************/

  /* Interface class for registration classes */
  class IRegister
  {
  public:
    /** destructor */
    virtual ~IRegister() = 0;


    /** Set main cloud point.
     * Registration of secondary cloud points are done against this cloud using RegisterCloud()
     * @param in_pcl           point cloud.
     * @param in_append        if true, append points to the existing cloud
     */
    virtual void SetMainPtCloud(const CPtCloud& in_pcl, bool in_append = false) = 0;


    /** Get registration for a secondary point cloud against the main cloud
     * The second cloud is not stored
     * @param Xo_registration      best registration found.
     * @param Xi_pcl               secondary point cloud.
     * @param Xi_estimatedOrient   estimation of registration, if 0 then estimation is identity.
     * @return                     registration's grade/error - the lower the better. 
     */
    virtual float RegisterCloud(const CPtCloud& in_pcl, CMat4& out_registration, CMat4* in_estimatedOrient = 0) = 0;


  protected:

  };


} // namespace tpcl

#endif // __tpcl_sp_H