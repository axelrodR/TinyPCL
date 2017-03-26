/******************************************************************************
*
*: Package Name: sldrcr_sp
*
******************************************************************************/
#include "registration.h"
#include <vector>


//#ifdef _DEBUG
//#define new_file DEBUG_NEW
//#endif

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

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

  struct CRegOptions
  {
    //parameters:
    int m_k_normalEstimation; //the number of points used for local plane fitting.
    CVec3 m_viewNormalFlip; //normal vectors direction flipped towards the viewpoint viewNormalFlip.
    float m_voxelSize; //if >0, use voxel grid downsampling.
    CVec3 m_floodSource_xy; //x y coordinantes from which to start flood for ground segmentaion in ViewpointGridCreation.
    float m_heightMapRes; //resolution for height map used in ViewpointGridCreation.
    int m_lineWidth; // width of the polar depth map
    int m_numlines;  // height of the polar depth map 
    int m_searchRange; //range of grid locations to check. if not using GPS this needs to be set to inf.
    int m_medFiltSize; //deniseing median filter size.
    float m_distFromMedianThresh; //max distance between point and median filter's result.
  };



  class CRegDictionary
  {
  public:
    CMat4* m_Orient;      // 4x4 matrix, camera location & orientations.
    int m_size;
  };



  ///////////////////////////////////////////////////////////////////////////////
  //
  //                           CRegDictionary
  //
  ///////////////////////////////////////////////////////////////////////////////
  class CRegDictionaryE : public CRegDictionary
  {
  public:
    float** m_descriptors;          // descriptors array     (per orientation).
    float** m_descriptorsDFT;       // descriptors DFT array (per orientation).


                                    /** destructor */
    CRegDictionaryE::CRegDictionaryE()
    {
      m_Orient = NULL;
      m_descriptors = NULL;
      m_descriptorsDFT = NULL;
      m_size = -1;
    }

    /** destructor */
    CRegDictionaryE::~CRegDictionaryE()
    {
      DeleteDictionary();
    }

    /** free allocated memory of the dictionary. */
    void DeleteDictionary()
    {
      for (int dicIndex = 0; dicIndex < m_size; dicIndex++)
      {
        delete[] m_descriptors[dicIndex];
        delete[] m_descriptorsDFT[dicIndex];

      }
      delete[] m_descriptors;
      delete[] m_descriptorsDFT;
      m_descriptors = NULL;
      m_descriptorsDFT = NULL;

      delete[] m_Orient;
      m_Orient = NULL;

      m_size = -1;
    }


    /** deletes previous dictionary and fills the dictionary with the grid - sets of location & rotation per entry.
    * @param floodSource_xy   source for flood based ground segmentation.
    * @param heightMapRes     Resolution of the height map.
    * @param Xi_pts           global point cloud */
    void ViewpointGridCreation(CVec3 floodSource_xy, float heightMapRes, int Xi_numPts, const CVec3* Xi_pts)
    {
      DeleteDictionary();
      float invRes = 1.0f / heightMapRes;

      //DO: find points only on ground, create grod from that.
      float xMin = Xi_pts[0].x; float xMax = xMin;
      float yMin = Xi_pts[0].y; float yMax = yMin;

      for (int ptrIndex = 1; ptrIndex < Xi_numPts; ptrIndex++)
      {
        float x = Xi_pts[ptrIndex].x;
        float y = Xi_pts[ptrIndex].y;
        xMin = min(xMin, x);    yMin = min(yMin, y);
        xMax = max(xMax, x);    yMax = max(yMax, y);
      }

      int N = int(ceil((xMax - xMin) * invRes)); //height map width.
      int M = int(ceil((yMax - yMin) * invRes)); //height map height.
      float* heightMap = new float[M*N];
      for (int hIndex = 0; hIndex < M*N; hIndex++) //initialize height map with inf.
        heightMap[hIndex] = FLT_MAX;

      for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
      {
        int x = int(floor((Xi_pts[ptrIndex].x - xMin) * invRes));
        int y = int(floor((Xi_pts[ptrIndex].y - yMin) * invRes));
        int index = y*N + x;
        if (heightMap[index] == FLT_MAX)
          heightMap[index] = Xi_pts[ptrIndex].z;
        else if (Xi_pts[ptrIndex].z > heightMap[index])
          heightMap[index] = Xi_pts[ptrIndex].z;
      }


      //source in pixels.
      int floodSourceX = int(ceil((floodSource_xy.x - xMin) * invRes));
      int floodSourceY = int(ceil((floodSource_xy.y - yMin) * invRes));
      //DO: flood for grid.

      delete[] heightMap;


      ////////
      //CVec3 normal;
      ////DO: calculate normal (pcnormals - in matlab). using parameter m_k_normalEstimation.

      ////normal alignment
      //for (int ptsIndex = 0; ptsIndex < Xi_numPts; ptsIndex++)
      //{
      //  CVec3 view2pts = Xi_pts[ptsIndex] - m_viewNormalFlip;
      //  float dotProduct = D3DXVec3Dot(&view2pts, &normal);
      //  if (dotProduct > 0)
      //    normal = -normal;
      //}

      ////DO: find xyz: go to prid point position and take only points from global cloud that are IDX_minRange < at < IDX_maxRange.
      ////DO: Translate to grid: for each xyz, xyz = xyz - xyzPoV.
      //actually need to find whole rotation.

      ////////



      //DO:
      //m_size = number of points we got in grid.
      m_descriptors = new float*[m_size];
      m_descriptorsDFT = new float*[m_size];
      m_Orient = new CMat4[m_size];
      //put all positions and rotations in m_Orient.

    }



    /** Convert a point cloud into the descriptor - range image, a polar map of distance around a point
    * @param Xi_lineWidth      width of the polar depth map.
    * @param Xi_numlines       height of the polar depth map.
    * @param Xi_pts            list of points in the cloud.
    * @param Xi_origin         the center of the polar map.
    * @param Xo_RangeImage   polar depth map.  */
    void ConvertPCL2descriptor(int Xi_lineWidth, int Xi_numlines, int Xi_numPts, const CVec3* Xi_pts, float* Xo_RangeImage)
    {
      const double M_PI = 3.14159265358979323846264338327950288;

      float azimuthRes = 2 * float(M_PI) / Xi_lineWidth;
      float elevationRes = float(M_PI) / Xi_numlines;
      //fill range image with default value - 0.
      for (int elevation = 0; elevation < Xi_numlines; elevation++)
      {
        for (int azimuth = 0; azimuth < Xi_lineWidth; azimuth++)
        {
          Xo_RangeImage[elevation*Xi_lineWidth + azimuth] = 0;
        }
      }

      //put closest range per azimuth & elevation in the range image.
      for (int i = 0; i < Xi_numPts; i++)
      {
        float x = Xi_pts[i].x;
        float y = Xi_pts[i].y;
        float z = Xi_pts[i].z;
        float azimuth = atan2(y, x);
        float elevation = atan2(z, sqrt(x*x + y*y));
        float r = sqrt(x*x + y*y + z*z);


        int azimuthInds = int(floor((azimuth + M_PI) / azimuthRes));
        int elevationInds = int(floor((elevation + M_PI / 2) / elevationRes));

        int index = elevationInds*Xi_lineWidth + azimuthInds;
        if (r < Xo_RangeImage[index])
          Xo_RangeImage[index] = r;
      }

    }


    /** find FFT of a 2D descriptor.
    * @param Xi_desc1Dsize      size of the descriptor's first dimension (width of the range image - lineWidth)
    * @param Xi_desc2Dsize      size of the descriptor's second dimension (height of the range image - numlines).
    * @param Xi_Descriptor      2D descriptor - range image.
    * @param Xo_DescriptorFFT   descriptor's 2D FFT.  */
    void ConvertDescriptor2FFT(int Xi_desc1Dsize, int Xi_desc2Dsize, float* Xi_Descriptor, float* Xo_DescriptorFFT)
    {
      //DO:
    }


    /** fills the dictionary entries with their descriptors.
    * @param Xi_pts           global point cloud
    * @param Xi_desc1Dsize   size of the descriptor's first dimension (width of the range image - lineWidth)
    * @param Xi_desc2Dsize   size of the descriptor's second dimension (height of the range image - numlines).  */
    void DescriptorsCreation(int Xi_numPts, const CVec3* Xi_pts, int Xi_desc1Dsize, int Xi_desc2Dsize = 0)
    {
      CVec3* ptsTransformed = new CVec3[Xi_numPts];

      //go over grid points and create descriptors for each point:
      for (int gridIndex = 0; gridIndex < m_size; gridIndex++)
      {
        CMat4 orients = m_Orient[gridIndex];

        //DO: transform points around position so z axis is as normal and origin (0,0,0) be at position (point-position).
        //  from Xi_pts to ptsTransformed.

        //create descriptor - range image, FFT
        delete[] m_descriptors[gridIndex];
        m_descriptors[gridIndex] = new float[Xi_desc1Dsize*Xi_desc2Dsize];
        ConvertPCL2descriptor(Xi_desc1Dsize, Xi_desc2Dsize, Xi_numPts, ptsTransformed, m_descriptors[gridIndex]);

        delete[] m_descriptorsDFT[gridIndex];
        m_descriptorsDFT[gridIndex] = new float[Xi_desc1Dsize*Xi_desc2Dsize];
        ConvertDescriptor2FFT(Xi_desc1Dsize, Xi_desc2Dsize, m_descriptors[gridIndex], m_descriptorsDFT[gridIndex]);
      }

      delete[] ptsTransformed;
    }


    /** return candidates suitable for the input descriptor.
    * @param Xi_descriptor                 input local descriptor.
    * @param Xo_Candidates                 list of the candidates.
    * @param Xo_Candidates                 list of the candidates grades. */
    void SearchDictionary(const int* Xi_descriptor, int* Xo_candidates, float* Xo_grades)
    {

    }

    /** Create a descriptor from (laser) polar depth map
    * @param Xi_descriptor0     polar distance map from which the descriptor is built
    * @param Xi_descriptor1       Output descriptor */
    virtual float U_MatchDescriptor(int* Xi_descriptor0, int* Xi_descriptor1)
    {
      float tempForRunning = 0;
      return tempForRunning;
    }

  };





  /** denose a point cloud.
  * @param Xi_lineWidth                width of the polar depth map.
  * @param Xi_numlines                 height of the polar depth map
  * @param Xi_medFiltSize              median filter size.
  * @param Xi_distFromMedianThresh     max distance between point and median filter's result.
  * @param Xio_pts                     in: point cloud,  out: denoised point cloud.
  * @param Xo_mask                     denoised point cloud mask on orig input.
  * @param removeIsolatedPixels        flag if to remove isolated pixels.  */
  void DenoisePointCloud(int Xi_lineWidth, int Xi_numlines, int Xi_medFiltSize, float Xi_distFromMedianThresh, int& Xio_numPts, const CVec3* Xio_pts, bool* Xo_mask, bool removeIsolatedPixels = false)
  {
    //TODO: waiting for final algorithm version.

    //float* RangeImage;
    ////calc range image.

    ////TODO: will it be better NOT to allocate every time this function called but rewrite over the same memory over and over?
    ////      or use only one matrix of each type?
    //float* Median_rangeImage = new float[Xi_lineWidth*Xi_numlines];
    //bool* distFromMedianMask = new bool[Xi_lineWidth*Xi_numlines];
    //bool* median_distFromMedianMask = new bool[Xi_lineWidth*Xi_numlines];
    //bool* denoisingMask = new bool[Xi_lineWidth*Xi_numlines]; //newMedianMask
    //bool* denoisingMask_median = new bool[Xi_lineWidth*Xi_numlines];
    ////median filter on RangeImage to Median_rangeImage, filter size medFiltSize.

    //for (int row = 0; row < Xi_numlines; row++)
    //{
    //  for (int col = 0; col < Xi_lineWidth; col++)
    //  {
    //    int index = row * Xi_lineWidth + col;
    //    float distFromMedian = abs(RangeImage[index] - Median_rangeImage[index]);
    //    distFromMedianMask[index] = distFromMedian < Xi_distFromMedianThresh;
    //  }
    //}

    ////median filter on distFromMedianMask to median_distFromMedianMask, filter size medFiltSize.
    //for (int row = 0; row < Xi_numlines; row++)
    //{
    //  for (int col = 0; col < Xi_lineWidth; col++)
    //  {
    //    int index = row * Xi_lineWidth + col;
    //    bool changed1to0Mask = (distFromMedianMask[index] - median_distFromMedianMask[index]) == 1;
    //    denoisingMask[index] = distFromMedianMask[index] && !changed1to0Mask;
    //  }
    //}

    //if (removeIsolatedPixels)
    //{
    //  //median filter on denoisingMask to denoisingMask_median, filter size medFiltSize.
    //  for (int row = 0; row < Xi_numlines; row++)
    //  {
    //    for (int col = 0; col < Xi_lineWidth; col++)
    //    {
    //      int index = row * Xi_lineWidth + col;
    //      bool isolatedPixels = (denoisingMask[index] - denoisingMask_median[index] == 1);
    //      if (isolatedPixels)
    //        denoisingMask[index] = false;
    //    }
    //  }

    //}




    //delete[] Median_rangeImage;
    //delete[] distFromMedianMask;
  }

  /** downsample a point cloud. Divides to grid from minXYZ (of pts) to max XYZ, of size m_voxelSize.
  *  If more than one point in same grid index, takes only first one.
  * @param Xi_voxelSize            size of a voxel in grid. if equal zero - NO downsampling.
  * @param Xio_pts                 in: point cloud,  out: downsampled point cloud. */
  void DownSamplePointCloud(float Xi_voxelSize, int& Xio_numPts, CVec3* Xio_pts)
  {
    float InvVoxelSize = 1.0f / Xi_voxelSize;

    //find min/max x/y/z:
    float xMin = Xio_pts[0].x; float xMax = xMin;
    float yMin = Xio_pts[0].y; float yMax = yMin;
    float zMin = Xio_pts[0].z; float zMax = zMin;

    for (int ptrIndex = 1; ptrIndex < Xio_numPts; ptrIndex++)
    {
      float x = Xio_pts[ptrIndex].x;
      float y = Xio_pts[ptrIndex].y;
      float z = Xio_pts[ptrIndex].z;

      xMin = min(xMin, x);    yMin = min(yMin, y);    zMin = min(zMin, z);
      xMax = max(xMax, x);    yMax = max(yMax, y);    zMax = max(zMax, z);
    }

    //create downsampling vector:
    int Mx = int(ceil((xMax - xMin) * InvVoxelSize));
    int My = int(ceil((yMax - yMin) * InvVoxelSize));
    int Mz = int(ceil((zMax - zMin) * InvVoxelSize));
    int Mxy = Mx*My;

    bool* VisitedVoxel = new bool[Mxy*Mz];
    memset(VisitedVoxel, false, Mxy*Mz * sizeof(bool));
    std::vector<CVec3> DownSampledPts;

    for (int ptrIndex = 0; ptrIndex < Xio_numPts; ptrIndex++)
    {
      float x = Xio_pts[ptrIndex].x;
      float y = Xio_pts[ptrIndex].y;
      float z = Xio_pts[ptrIndex].z;

      int xInd = int(floor((x - xMin) * InvVoxelSize));
      int yInd = int(floor((y - yMin) * InvVoxelSize));
      int zInd = int(floor((z - zMin) * InvVoxelSize));

      if (xInd == Mx) xInd--;
      if (yInd == My) yInd--;
      if (zInd == Mz) zInd--;

      int index = zInd*Mxy + yInd*Mx + xInd;
      if (!VisitedVoxel[index])
      {
        DownSampledPts.push_back(Xio_pts[ptrIndex]);
        VisitedVoxel[index] = true;
      }
    }
    delete[] VisitedVoxel;


    //save selected points:
    for (int VecIndex = 0; VecIndex < DownSampledPts.size(); VecIndex++)
    {
      Xio_pts[VecIndex] = DownSampledPts[VecIndex];
    }
    Xio_numPts = int(DownSampledPts.size());
  }


  /** remove out liers of a point cloud.
  * @param Xio_pts                 in: point cloud,  out: point cloud without outliers.
  * @param Xio_mask                in: mask of downsampled orig PC. out: mask on downsampled input without outliers. */
  void RemoveOutliersPointCloud(int& Xio_numPts, const CVec3* Xio_pts, bool* Xio_mask)
  {
    //TODO: might be removed.
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
  *: Method name: tpclCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::CCoarseRegister()
  {
    m_opts = new CRegOptions;
    CRegOptions* optsP = (CRegOptions*)m_opts;
    optsP->m_k_normalEstimation = 0; //the number of points used for local plane fitting.
                                     //    optsP->m_viewNormalFlip       = ; //normal vectors direction flipped towards the viewpoint viewNormalFlip.
    optsP->m_voxelSize = 0; //if >0, use voxel grid downsampling.
                            //    optsP->m_floodSource_xy       = ; //x y coordinantes from which to start flood for ground segmentaion in ViewpointGridCreation.
    optsP->m_heightMapRes = 0; //resolution for height map used in ViewpointGridCreation.
    optsP->m_lineWidth = 0; // width of the polar depth map
    optsP->m_numlines = 0; // height of the polar depth map 
    optsP->m_searchRange = 0; //range of grid locations to check. if not using GPS this needs to be set to inf.
    optsP->m_medFiltSize = 0; //deniseing median filter size.
    optsP->m_distFromMedianThresh = 0; //max distance between point and median filter's result.

    m_dictionary = new CRegDictionaryE;
  }

  /******************************************************************************
  *
  *: Method name: ~tpclCR_SP_CCoarseRegister
  *
  ******************************************************************************/
  CCoarseRegister::~CCoarseRegister()
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    delete optsP;

    CRegDictionaryE* dictionaryP = (CRegDictionaryE*)m_dictionary;
    delete dictionaryP;

    delete[] m_ptsGlobal;
  }


  /******************************************************************************
  *
  *: Method name: SetGlobalCloud
  *
  ******************************************************************************/
  void CCoarseRegister::BuildDictionary(int Xi_numPts, const CVec3* Xi_pts)
  {
    CRegOptions* optsP = (CRegOptions*)m_opts;
    CRegDictionaryE* dictionaryP = (CRegDictionaryE*)m_dictionary;

    //copy global cloud to local variable:
    delete[] m_ptsGlobal;
    m_ptsGlobal = new CVec3[Xi_numPts];
    for (int ptrIndex = 0; ptrIndex < Xi_numPts; ptrIndex++)
    {
      m_ptsGlobal[ptrIndex] = Xi_pts[ptrIndex];
    }
    m_numPtsGlobal = Xi_numPts;

    //preprocess global cloud:
    DownSamplePointCloud(optsP->m_voxelSize, m_numPtsGlobal, m_ptsGlobal);



    //create grid's position and orientation: 
    dictionaryP->ViewpointGridCreation(optsP->m_floodSource_xy, optsP->m_heightMapRes, Xi_numPts, Xi_pts); //TODO: complete function ViewpointGridCreation.

                                                                                                           //fill the dictionary entries with their descriptors:
    dictionaryP->DescriptorsCreation(Xi_numPts, Xi_pts, optsP->m_lineWidth, optsP->m_numlines);
  }


  /******************************************************************************
  *
  *: Method name: GetLocalRegistrationCandidates
  *
  ******************************************************************************/
  int CCoarseRegister::GetLocalRegistrationCandidates(int Xi_maxCandidates, int Xi_numPts, const CVec3* Xi_pts, CVec3 Xi_originApprox, int* Xo_candidates)
  {
    int NumOfCandidates = 0;

    ////Preprocess local cloud:
    //bool* local_mask = NULL; //TODO
    //DenoisePointCloud(m_lineWidth, m_numlines, m_medFiltSize, m_distFromMedianThresh, Xi_numPts, Xi_pts, local_mask);
    //DownSamplePointCloud(m_voxelSize, Xi_numPts, Xi_pts);


    ////create descriptor - range image:
    //float* LocalRangeImage = new float[m_lineWidth*m_numlines];
    //ConvertPCL2polarMap(m_lineWidth, m_numlines, Xi_numPts, Xi_pts, LocalRangeImage);

    ////create range image's 2D FFT:
    //float* LocalRangeImageFFT = new float[m_lineWidth*m_numlines * 2]; // *2 - unreal values.
    ////DO: create 2D fft. FFT2D(m_lineWidth, m_numlines, LocalRangeImage, LocalRangeImageFFT);


    //int* localDescriptor;

    //U_CreateDescriptor(RangeImage, localDescriptor);
    //
    //std::vector<CRegDictionaryEntry> candidates;
    ////add use of GPS location filtering.
    //m_dictionary->SearchDictionary(localDescriptor, Xo_candidates, Xo_grades); //DO: should return candidates so need Xi_maxCandidates.

    //delete[] RangeImage;

    return NumOfCandidates;
  }



  /******************************************************************************
  *
  *: Method name: GetLocalRegistration
  *
  ******************************************************************************/
  void CCoarseRegister::GetLocalRegistration(D3DXMATRIX& Xo_best)
  {

  }


  /******************************************************************************
  *                             Protected methods                               *
  ******************************************************************************/


  /** creates the Height Profile descriptor.
  * @param Xi_pts        point cloud
  * @param Xo_descriptor   Height Profile descriptor */
  void CCoarseRegister::U_CreateDescriptor(float* Xi_RangeImage, int* Xo_descriptor)
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


} //namespace tpcl