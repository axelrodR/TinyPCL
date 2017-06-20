

#include "../include/tinyPCL.h"
#include "../registration/RegICP.h"
#include "../registration/RegistDict.h"
#include "features.h"

namespace tpcl
{


  CTinyPCL::CTinyPCL()
  {

  }



  IRegister* CTinyPCL::GenerateRegistration(ERegistrationType in_type)
  {
    switch (in_type)
    {
    case REGISTRATION_TYPE_ICP: return new ICP();
    case REGISTRATION_TYPE_POV: return new CCoarseRegister();
    }
  }


  IFeatures* CTinyPCL::GenerateFeatures()
  {
    return new Features();
  }



  void CTinyPCL::Free(IRegister* in_register)
  {
    delete in_register;
  }

  void CTinyPCL::Free(IFeatures* in_features)
  {
    delete in_features;
  }


} //namespace tpcl