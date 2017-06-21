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
*                        INCOMPLETE CLASS DECLARATIONS                        *
******************************************************************************/
struct CVec3;
struct CMat4;


namespace tpcl
{
  class CSpatialHash2D;

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: ICP
  *
  *: Abstract: basic ICP calculation on point cloud.
  *            see: <FILL PAPER REF HERE>
  *
  ******************************************************************************/

  class ICP : public IRegister
  {
  public:
    /** Constructor 
    * @param in_regThresh         set registration threshold. */
    ICP();
    ICP(float in_regRes);

    /** destructor */
    virtual ~ICP();

    /** Set main cloud point.
    * Registration of secondary cloud points are done against this cloud using RegisterCloud()
    * @param in_pcl           point cloud.
    * @param in_append        if true, append points to the existing cloud
    */
    void SetMainPtCloud(const CPtCloud& in_pcl, bool in_append = false);

    /** Set main cloud point.
    * Registration of secondary cloud points are done against this cloud using RegisterCloud()
    * expects data to be available whenever registration is called!
    * @param in_mainHashed    pointer to an already hashed point cloud. deletes any local main point cloud. expects data to be available whenever registration is called.
    */
    void SetMainPtCloud(CSpatialHash2D* in_mainHashed);

    /** Get hashed main point cloud.
    * return         pointer to hashed main point cloud. */
    void* getMainHashedPtr();

    /** Set registration threshold.
    * @param in_regThresh         registration threshold. */
    void setRegistrationResolution(float in_regRes);

    /** Get registration for a secondary point cloud against the main cloud
    * The second cloud is not stored
    * @param out_registration      best registration found.
    * @param in_pcl               secondary point cloud.
    * @param in_estimatedOrient   estimation of registration, if 0 then estimation is identity.
    * @return                     registration's grade/error - the lower the better.
    */
    float RegisterCloud(const CPtCloud& in_pcl, CMat4& out_registration, CMat4* in_estimatedOrient = 0);


  protected:
    CSpatialHash2D* m_mainHashed;   ///< a hashed copy of the main point cloud.
    bool m_outsourceMainPC;         ///< if true then hashed main point cloud used if given from outside (and will not be changed).
    float m_regRes;                 ///< resolution of registration wanted.

    /** Set default values to members. */
    void initMembers();
  };

} // namespace tpcl

#endif


