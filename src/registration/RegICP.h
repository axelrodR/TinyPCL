// File Location: l:\SLDR\sldrcr\sldrcr_icp.h
/******************************************************************************
*
*: Package Name: sldrcr_icp
*
*: Title:
*
******************************************************************************/

#ifndef __tpcl_register_icp_H
#define __tpcl_register_icp_H

#include "../include/registration.h"


/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/


/******************************************************************************
*                        INCOMPLETE CLASS DECLARATIONS                        *
******************************************************************************/
struct CVec3;
struct CMat4;

namespace tpcl
{
  class CSpatialHash2D;
}

namespace tpcl
{

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/


  /******************************************************************************
  *
  *: Class name: SLDR_ICP_ICP
  *
  *: Abstract: basic ICP calculation on point cloud.
  *            see: <FILL PAPER REF HERE>
  *
  ******************************************************************************/

  class ICP : public IRegister
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor
    * @param Xi_regThresh         set registration threshold. */
    ICP();
    ICP(double Xi_regThresh);

    /** destructor */
    virtual ~ICP();

    /** update the main point cloud in which registration is searched for, for a secondary point cloud.
    *   if previous main point cloud was given from outside, that pointer it is forgotten.
    * @param Xi_pts           point cloud to add to existing main point cloud.
    * @param Xi_clean         if true deletes all previous information of main point cloud (if previous main point cloud wasn't given from outside) before adding input point cloud.
    * @param Xi_mainHashed    pointer to an already hashed point cloud. deletes any local main point cloud. expects data to be available whenever registration is called. */
    void MainPointCloudUpdate(int Xi_numPts, const CVec3* Xi_pts, bool Xi_clean = false);

    void MainPointCloudUpdate(void* Xi_mainHashed);

    /** Get hashed main point cloud.
    * return         pointer to hashed main point cloud. */
    void* getMainHashedPtr();

    /** Set registration threshold.
    * @param Xi_regThresh         registration threshold. */
    void setRegistrationThresh(double Xi_regThresh);

    /** Get best ICP registration for a secondary point cloud.
    * @param Xo_registration      best registration found.
    * @param Xi_pts               secondary point cloud.
    * @param Xi_numPts            total points in the secondary point cloud.
    * @param Xi_lineWidth         NOT USED! if it's an order point cloud then this is the width of the secondary point cloud.
    * @param Xi_estimatedOrient   estimation of registration, if NULL then estimation is identity.
    * return                      registration's grade/error - the lower the better.
    * !!!!currently no utilization of ordered point cloud!!!! */
    float SecondaryPointCloudRegistration(CMat4& Xo_registration, CVec3* Xi_pts, int Xi_numPts, int Xi_lineWidth = -1, CMat4* Xi_estimatedOrient = NULL);


  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/
    float m_voxelSize;                        ///< the voxel size parameter of the local hashed main point cloud.
    CSpatialHash2D* m_mainHashed;   ///< a hashed copy of the main point cloud.
    bool m_outsourceMainPC;                   ///< if true then hashed main point cloud used if given from outside (and will not be changed).
    double m_minDelta;                        ///< threshold of registration change between ICP iterations.
                                              /******************************************************************************
                                              *                             Protected methods                               *
                                              ******************************************************************************/
                                              /** Set default values to members. */
    void initMembers();
  };




} // namespace tpcl

#endif


