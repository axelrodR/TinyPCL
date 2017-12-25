#include "../include/tinyPCL.h"

using namespace tpcl;


// Fill in a point cloud from file.
void LoadPointCloud(char* fileName, CPtCloud& cloud, bool globalCloud);
//Get estimated orientation
CMat4 GetEstimatedOrientation();



void main()
{
  CTinyPCL TPCL;

  // Get a Dictionary Registration registrator:
  IRegister* dicReg = TPCL.GenerateRegistration(ERegistrationType::REGISTRATION_TYPE_POV);

  // create global point cloud:
  CPtCloud globalPcl;
  LoadPointCloud("MyGlobalCloud.txt", globalPcl, true);

  // create dictionary from global point cloud:
  dicReg->SetMainPtCloud(globalPcl);

  // create local point cloud:
  CPtCloud localPcl;
  LoadPointCloud("MyLocalCloud.txt", localPcl, false);

  // Get registration estimation: 
  CMat4 estimatedRegisration = GetEstimatedOrientation();

  // get registration for local point cloud:
  CMat4 registrationResult;
  dicReg->RegisterCloud(localPcl, registrationResult, &estimatedRegisration);
  // registration result in registrationResult.

  //delete allocated arrays.
  delete[] globalPcl.m_pos;
  delete[] localPcl.m_pos;
}



// Fill in a point cloud from file.
void LoadPointCloud(char* fileName, CPtCloud& cloud, bool globalCloud)
{
  // this "header" needs to be read from file
  cloud.m_type = globalCloud ? PCL_TYPE_FUSED : PCL_TYPE_SINGLE_ORIGIN_SCAN;
  cloud.m_lineWidth = globalCloud ? 0 : 2;
  cloud.m_numPts = 10; //read number of points from file.

  // the points need to be read from file as well. For now - just fill 5 lines.
  cloud.m_pos = new CVec3[cloud.m_numPts];
  for (int i = 0; i < cloud.m_numPts; ++i)
  {
    CVec3 point(i%2, 0, i/2);
    cloud.m_pos[i] = point;
  }
}



//Get estimated orientation
CMat4 GetEstimatedOrientation()
{
  //Note: dictionary registration doesn't use rotation estimation.
  CVec3 laserPos(100, 100, 0);    // in example laser position assumed to be at (100, 100, 0)
  
  // set the point cloud oreintation matrix
  CMat4 estimatedRegisration;
  MatrixIdentity(&estimatedRegisration);
  estimatedRegisration.m[3][0] = laserPos[0];
  estimatedRegisration.m[3][1] = laserPos[1];
  estimatedRegisration.m[3][2] = laserPos[2];

  return estimatedRegisration;
}