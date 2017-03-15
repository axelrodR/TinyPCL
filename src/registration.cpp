/******************************************************************************
*
*: Package Name: sldrcr_sp
*
******************************************************************************/
#include <D3dx9core.h> // uncommenet if using DirectX
#include <ifr/ifrgen/ifrgen_stnd.h>
#include "sldrcr_sp.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

namespace SLDR
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


  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CRegDictionary
  //
  ///////////////////////////////////////////////////////////////////////////////
  class CRegDictionary
  {
  public:
    void AddDescriptor(const CRegDictionaryEntry&);

    /** return candidates suitable for the input descriptor.
    * @param Xi_descriptor                 input local descriptor.
    * @param Xo_Candidates                 list of the candidates.
    * @param Xo_Candidates                 list of the candidates grades. */
    void SearchDictionary(const int* Xi_descriptor, CRegDictionaryEntry* Xo_candidates, float* Xo_grades)
    {

    }

    /** Create a descriptor from (laser) polar depth map
    * @param Xi_polarDistMap     polar distance map from which the descriptor is built
    * @param Xo_descriptor       Output descriptor */
    virtual float U_MatchDescriptor(int* Xi_descriptor0, int* Xi_descriptor1)
    {

    }

    CRegDictionaryEntry* m_enteries;
  };



  /** denose a point cloud.
  * @param Xio_PointCloud    in: point cloud,  out: denoised point cloud  */
  void DenoisePointCloud(int Xi_numPts, const D3DXVECTOR3* Xi_pts)
  {

  }

  /** downsample a point cloud.
  * @param Xio_PointCloud    in: point cloud,  out: downsampled point cloud  */
  void DownSamplePointCloud(int Xi_numPts, const D3DXVECTOR3* Xi_pts)
  {

  }

  /** remove out liers of a point cloud.
  * @param Xio_PointCloud    in: point cloud,  out: point cloud without out liers  */
  void RemoveOutliersPointCloud(int Xi_numPts, const D3DXVECTOR3* Xi_pts)
  {

  }


  /** Convert a point cloud into a polar map of distance around a point
  * @param Xi_lineWidth      width of the polar depth map.
  * @param Xi_numlines       height of the polar depth map.
  * @param Xi_pts            list of points in the cloud.
  * @param Xi_origin         the center of the polar map.
  * @param Xo_polarDistMap   polar depth map.  */
  void ConvertPCL2polarMap(int Xi_lineWidth, int Xi_numlines, int Xi_numPts, const D3DXVECTOR3* Xi_pts, const D3DXVECTOR3& Xi_origin, float* Xo_polarDistMap)
  {

  }



  /** return a dictionary filled with only the grid - sets of location & normal per entry.
  * @param Xi_pts          point cloud.  */
  CRegDictionary* ViewpointGridCreation(int Xi_numPts, const D3DXVECTOR3* Xi_pts)
  {

  }



  /******************************************************************************
  *                           EXPORTED CLASS METHODS                            *
  ******************************************************************************/
  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CCoarseRegister
  //
  ///////////////////////////////////////////////////////////////////////////////
  /******************************************************************************
  *                               Public methods                                *
  ******************************************************************************/
  /******************************************************************************
  *
  *: Method name: SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::CCoarseRegister()
  {
    m_dictionary = new CRegDictionary();
  }

  /******************************************************************************
  *
  *: Method name: ~SLDRCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::~CCoarseRegister()
  {
  }


  /******************************************************************************
  *
  *: Method name: SetGlobalCloud
  *
  ******************************************************************************/
  void CCoarseRegister::BuildDictionary(int Xi_numPts, const D3DXVECTOR3* Xi_pts)
  {
    DenoisePointCloud(Xi_numPts, Xi_pts);
    DownSamplePointCloud(Xi_numPts, Xi_pts);
    RemoveOutliersPointCloud(Xi_numPts, Xi_pts);

    m_dictionary = ViewpointGridCreation(Xi_numPts, Xi_pts);

    float* PolarDistMap = new float[m_lineWidth*m_numlines];
    //go over grid points and create signature for each point
    {
      D3DXVECTOR3 l_pos; // fill me: position in grid
      ConvertPCL2polarMap(m_lineWidth, m_numlines, Xi_numPts, Xi_pts, l_pos, PolarDistMap);

      //put descriptor in dictionary.
      CRegDictionaryEntry l_entry;
      U_CreateDescriptor(PolarDistMap, l_entry.m_descriptor);
      l_entry.pos = l_pos;
      m_dictionary->AddDescriptor(l_entry);
    }

    delete[] PolarDistMap;
  }


  /******************************************************************************
  *
  *: Method name: GetLocalRegistrationCandidates
  *
  ******************************************************************************/
  void CCoarseRegister::GetLocalRegistrationCandidates(int Xi_maxCandidates, int Xi_numPts, const D3DXVECTOR3* Xi_pts, D3DXVECTOR3& Xi_origin, CRegDictionaryEntry* Xo_candidates, float* Xo_grades)
  {

    DenoisePointCloud(Xi_numPts, Xi_pts);

    float* PolarDistMap = new float[m_lineWidth*m_numlines];

    D3DXVECTOR3 localOrigin(0, 0, 0);
    ConvertPCL2polarMap(m_lineWidth, m_numlines, Xi_numPts, Xi_pts, localOrigin, PolarDistMap);

    int* localDescriptor;

    U_CreateDescriptor(PolarDistMap, localDescriptor);
    
    std::vector<CRegDictionaryEntry> candidates;
    //add use of GPS location filtering.
    m_dictionary->SearchDictionary(localDescriptor, Xo_candidates, Xo_grades);

    delete[] PolarDistMap;
  }



  /******************************************************************************
  *
  *: Method name: GetLocalRegistration
  *
  ******************************************************************************/
  float CCoarseRegister::GetLocalRegistration(CRegDictionaryEntry *best))
  {
    
  }

 
  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/


  /** creates the Height Profile descriptor.
  * @param Xi_pts        point cloud
  * @param Xo_descriptor   Height Profile descriptor */
  void CCoarseRegister::U_CreateDescriptor(float* Xi_polarDistMap, int* Xo_descriptor)
  {

  }


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