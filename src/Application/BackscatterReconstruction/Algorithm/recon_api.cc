
#include <Application/BackscatterReconstruction/Algorithm/recon_api.h>


namespace BackscatterReconstruction
{
  
  
// calibration - segment image stack of calibration pattern to find disks
void CalibrationSegment(const ReconImageVolumeType::Pointer images,
                        ReconMaterialIdVolumeType::Pointer diskIds)
{
}


// calibration - given images and disk id masks (after mapping to match calibration pattern), fit detector geometry
void CalibrationFitGeometry(const ReconImageVolumeType::Pointer images,
                            const ReconMaterialIdVolumeType::Pointer diskIds,
                            const char *workDirectory,
                            const char *initialGeometryConfigFile,
                            const char *outputGeometryConfigFile)
{
}


// reconstruction - start a thread internally so this returns immediately
void ReconstructionStart(const ReconImageVolumeType::Pointer images,
                         const ReconMaterialIdVolumeType::Pointer initialGuess, // possibly NULL
                         const char *geometryConfigFile,
                         int iterations)
{
}

// reconstruction - stop it if it's currently running
void ReconstructionAbort()
{
}


// reconstruction - get the latest reconstructed volume material ids, can be called any time during reconstruction
void ReconstructionGetMaterialVolume(ReconMaterialIdVolumeType::Pointer reconVolume)
{
}

}