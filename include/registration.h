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

/******************************************************************************
*                                   IMPORTED                                  *
******************************************************************************/

#undef DLL_Entry 
#ifdef SLDRCR_EXPORTS
#define DLL_Entry __declspec(dllexport)
#else
#define DLL_Entry __declspec(dllimport)
#endif

#include <vector>

namespace SLDR
{

  /******************************************************************************
  *                        INCOMPLETE CLASS DECLARATIONS                        *
  ******************************************************************************/
  struct CRegDictionaryEntry;    // registration Dictionary entry
  struct CRegDictionary;         // registration Dictionary


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


    /** Build dictionary from global point cloud 
    * @param Xi_pts           point cloud 
    * @param Xi_lineWidth     width of the polar depth map
    * @param Xi_numlines      height of the polar depth map */
    void BuildDictionary(int Xi_numPts, const D3DXVECTOR3* Xi_pts);


    /** Create and get list of registration matches for a local point cloud
    * @param Xi_maxCandidates     maximum number of matches to return.
    * @param Xi_pts               local point cloud
    * @param Xi_origin            measured/approximate origin.
    * @param Xi_origin            measured/approximate origin. */
    void GetLocalRegistrationCandidates(int Xi_maxCandidates, int Xi_numPts, const D3DXVECTOR3* Xi_pts, D3DXVECTOR3& Xi_origin, CRegDictionaryEntry* Xo_candidates, float* Xo_grades);


    /** Get best registration match from created candidates list.
     * @return grade of best match  */
    float GetLocalRegistration(CRegDictionaryEntry *best);



  protected:
    /******************************************************************************
    *                             Protected members                               *
    ******************************************************************************/

    CRegDictionary *m_dictionary;

    CRegDictionaryEntry* m_candidates;
    float* m_grades;

    int m_lineWidth; // width of the polar depth map
    int m_numlines;  // height of the polar depth map 

    /******************************************************************************
    *                             Protected methods                               *
    ******************************************************************************/

    /** Create a descriptor from (laser) polar depth map
     * @param Xi_polarDistMap     polar distance map from which the descriptor is built
     * @param Xo_descriptor       Output descriptor */
    virtual void U_CreateDescriptor(float* Xi_polarDistMap, int* Xo_descriptor);
  };



  /******************************************************************************
  *
  *: Class name: CRegDictionaryEntry
  *
  *: Abstract:
  *
  ******************************************************************************/

  struct DLL_Entry CRegDictionaryEntry
  {
    D3DXVECTOR3 m_pos;
    int* m_descriptor;
  };



  /** another class with the 2nd algorithm proposed by David Avidar */
  class DLL_Entry CCoarseRegister2 : public CCoarseRegister
  {
  public:
    CCoarseRegister2();
    virtual void U_CreateDescriptor(float* Xi_polarDistMap, int* Xo_descriptor);
  };


#undef DLL_Entry
#endif


} // namespace SLDR


/******************************************************************************
*                            old style typedefs                               *
******************************************************************************/

typedef SLDR::CCoarseRegister SLDRCR_SP_CCoarseRegister;
