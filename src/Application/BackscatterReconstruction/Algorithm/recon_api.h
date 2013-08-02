#ifndef APPLICATION_BACKSCATTERRECONSTRUCTION_ALGORITHM_RECON_API_H
#define APPLICATION_BACKSCATTERRECONSTRUCTION_ALGORITHM_RECON_API_H

#include <itkImage.h>

typedef itk::Image<float, 3> ReconImageVolumeType;
typedef itk::Image<unsigned char, 3> ReconMaterialIdVolumeType;

// TODO: helper classes would be better here, especially to manage status
// variables and to do cleanup on error

// TODO: add Seg3D log output

// calibration - segment image stack of calibration pattern to find disks
void CalibrationSegment(const ReconImageVolumeType::Pointer images,
                        std::vector<unsigned char>& diskIds);

// calibration - given images and disk id masks (after mapping to match calibration pattern), fit detector geometry
void CalibrationFitGeometry(const ReconImageVolumeType::Pointer images,
                            const std::vector<unsigned char>& diskIds,
                            const char *workDirectory,
                            const char *initialGeometryConfigFile,
                            const char *outputGeometryConfigFile);


// reconstruction - don't return until reconstruction is completed
void ReconstructionStart(const ReconImageVolumeType::Pointer images,
                         const ReconMaterialIdVolumeType::Pointer initialGuess, // possibly NULL
                         const char *geometryConfigFile,
                         const char *sourceFalloffMapFile,
                         const float voxelSizeCM[3],
                         int iterations,
                         ReconMaterialIdVolumeType::Pointer finalReconVolume);

// reconstruction - stop it if it's currently running
void ReconstructionAbort();


// reconstruction - get the latest reconstructed volume material ids, can be called any time during reconstruction
// returns a normalized progress value in the range [0,1]
double ReconstructionGetMaterialVolume(ReconMaterialIdVolumeType::Pointer reconVolume);

#endif
