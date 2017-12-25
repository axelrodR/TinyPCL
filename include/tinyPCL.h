// File Location: 

//
// Copyright (c) 2016-2017 Geosim Ltd.
// 
// Written by Amit Henig
//            Ramon Axelrod
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


#ifndef  __tpcl_h
#define __tpcl_h

/** include all services */
#include "ptCloud.h"
#include "registration.h"
#include "iFeatures.h"


namespace tpcl
{
  // forward declaration
  enum ERegistrationType : char;


  /** factory for all tiny-PCL classes */
  class CTinyPCL
  {
  public:
    /** create a registration service 
     * @param in_type     type of registarion algorithm to use
     */
    IRegister* GenerateRegistration(ERegistrationType in_type);

    /** Create a feature manipulation (denoise, compute normals, etc.) servce */
    IFeatures* GenerateFeatures();

    // free created services
    void Free(IRegister* in_register);
    void Free(IFeatures* in_features);

    /** constructor */
    CTinyPCL();
  };


  /** defines the types of registration classes defined by the system */
  enum ERegistrationType : char
  {
    REGISTRATION_TYPE_ICP = 1,      ///< standatd ICP method
    REGISTRATION_TYPE_POV = 2,      ///< registration of single POV cloud against genral multi-pov clouds
  };




} // namespace tpcl


#endif // ! __tpcl_h
