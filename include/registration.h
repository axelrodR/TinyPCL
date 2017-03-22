// File Location: l:\SLDR\sldrcr\sldrcr_sp.h
/******************************************************************************
*
*: Package Name: sldrcr_sp
*
*: Title:
*
******************************************************************************/

#ifndef __sldrcr_sp_H
#define __sldrcr_sp_H

#include "vec.h"


/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#undef DLL_Entry 
#ifdef SLDRCR_EXPORTS
#define DLL_Entry __declspec(dllexport)
#else
#define DLL_Entry __declspec(dllimport)
#endif


  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/

namespace tpcl
{
  class  CRegDictionary;         // registration Dictionary

  /******************************************************************************
  *                              EXPORTED CLASSES                               *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Class name: SLDRCR_SP_CCoarseRegister
  *
  *: Abstract: Coarse registeration vased on dictionary
  *            Derived from the thesis of David Avidar (Malaach's group)
  *            see: <FILL PAPER REF HERE> 
  *
  ******************************************************************************/

  class DLL_Entry CCoarseRegister
  {
  public:
    /******************************************************************************
    *                               Public methods                                *
    ******************************************************************************/
    /** Constructor */
    CCoarseRegister();

    /** destructor */
    virtual ~CCoarseRegister();


    /** Build dictionary from global point cloud. Xi_pts IS being altered.
    * @param Xi_pts           point cloud */
    void BuildDictionary(int Xi_numPts, const CVec3* Xi_pts);


    /** Create and get list of registration matches for a local point cloud
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_pts               local point cloud.
    * @param Xi_originApprox      measured/approximate origin.
    * @param Xo_candidates        indices of the candidates in the dictionary.
    * @return                     number of candidates.*/
    int GetLocalRegistrationCandidates(int Xi_maxCandidates, int Xi_numPts, const CVec3* Xi_pts, CVec3 Xi_originApprox, int* Xo_candidates);


    /** Get best registration match from created candidates list.
    * @return grade of best match  */
    void GetLocalRegistration(D3DXMATRIX& Xo_best);



  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/

    int m_numPtsGlobal;
    CVec3* m_ptsGlobal;       ///< a copyu of the original point cloud
    CRegDictionary* m_dictionary;   ///< internal data used to store features of global cloud
    void* m_opts;                   ///< implementation specific options

    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/

    /** Create a descriptor from (laser) polar depth map
     * @param Xi_polarDistMap     polar distance map from which the descriptor is built
     * @param Xo_descriptor       Output descriptor */
    virtual void U_CreateDescriptor(float* Xi_polarDistMap, int* Xo_descriptor);
  };



  ///** another class with the 2nd algorithm proposed by David Avidar */
  //class DLL_Entry CCoarseRegister2 : public CCoarseRegister
  //{
  //public:
  //  CCoarseRegister2();
  //  virtual void U_CreateDescriptor(float* Xi_polarDistMap, int* Xo_descriptor);
  //};


#undef DLL_Entry
#endif


} // namespace tpcl

