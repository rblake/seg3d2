
#ifdef WIN32
#pragma warning(disable: 4996)
#endif

#include <string.h>
#ifndef WIN32
#define sscanf_s sscanf
#define _stricmp strcasecmp
#define sprintf_s snprintf
#define fopen_s(f,n,m) (*f) = fopen(n,m)
#endif

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "geometry.h"
#include "sart.h"
#include "markov.h"
#include "calibration.h"
#include "material.h"

// include this stuff for poisson distributed variables!
#include "nr3.h"
#include "gamma.h"
#include "ran.h"
#include "deviates.h"

#include <tiff.h>
#include <tiffio.h>

#define XRAY_ENERGY_LOW  45
#define XRAY_ENERGY_HIGH 45


typedef itk::RGBPixel<float> RGBPixelType;

typedef itk::Image<RGBPixelType, 3> RGBVolumeType;
typedef itk::ImageFileWriter< RGBVolumeType > RGBWriterType;

typedef itk::Image<float, 4> Float4VolumeType;
typedef itk::ImageFileWriter< Float4VolumeType > Float4WriterType;

typedef itk::ImageFileReader< FloatVolumeType > FloatReaderType;
typedef itk::ImageFileWriter< FloatVolumeType > FloatWriterType;

typedef itk::ImageFileReader< ByteVolumeType > ByteReaderType;
typedef itk::ImageFileWriter< ByteVolumeType > ByteWriterType;



int RunMarkovForwardProject(int argc, char **argv) {
  std::cerr<<"Running Markov Forward Project"<<std::endl;
  if (argc != 11) {
    std::cerr<<"usage: backscatter markov_forward_project <geometry.txt> <air_material.txt> <foam_material.txt> <aluminum_material.txt> <volume.nrrd> <samples_per_pixel> <voxel_step_size> <sourceatten_output.nrrd> <forward_output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *airName = argv[3];
  char *foamName = argv[4];
  char *aluminumName = argv[5];
  char *volumeName = argv[6];
  int samplesPerPixel = atoi(argv[7]);
  float voxelStepSize = atof(argv[8]);
  char *sourceAttenName = argv[9];
  char *outputName = argv[10];


  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read materials
  vector<Material> materials(3);
  materials[0].InitFromFile(airName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[1].InitFromFile(foamName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[2].InitFromFile(aluminumName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);

  // read volume to project from
  ByteReaderType::Pointer reader = ByteReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  ByteVolumeType::Pointer volume = reader->GetOutput();
  ByteVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != geometry.GetVolumeNodeSamplesX() ||
      volumeSize[1] != geometry.GetVolumeNodeSamplesY() ||
      volumeSize[2] != geometry.GetVolumeNodeSamplesZ()) {
    std::cerr<<"input volume dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  MarkovContext markovContext(geometry, materials, samplesPerPixel, voxelStepSize, 0);
  markovContext.SetCurrentVolume(volume);

  // cast rays from source into volume
  markovContext.ComputeSourceAttenuation();

  FloatVolumeType::Pointer sourceAttenVol = FloatVolumeType::New();
  markovContext.GetCurrentVolumeSourceAttenuation(sourceAttenVol);

  FloatWriterType::Pointer sourceAttenWriter = FloatWriterType::New();
  sourceAttenWriter->SetInput(sourceAttenVol);
  sourceAttenWriter->SetFileName(sourceAttenName);
  sourceAttenWriter->Update();


  // cast rays from detector
  markovContext.ComputeForwardProjection();

  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  markovContext.GetCurrentProjection(outVol);

  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}


int RunMarkovInitialGuess(int argc, char **argv) {
  std::cerr<<"Running Markov Initial Guess"<<std::endl;
  if (argc != 10) {
    std::cerr<<"usage: backscatter markov_initial_guess <geometry.txt> <air_material.txt> <foam_material.txt> <aluminum_material.txt> <forward_proj.nrrd> <samples_per_pixel> <voxel_step_size> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *airName = argv[3];
  char *foamName = argv[4];
  char *aluminumName = argv[5];
  char *forwardName = argv[6];
  int samplesPerPixel = atoi(argv[7]);
  float voxelStepSize = atof(argv[8]);
  char *outputName = argv[9];


  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read materials
  vector<Material> materials(3);
  materials[0].InitFromFile(airName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[1].InitFromFile(foamName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[2].InitFromFile(aluminumName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);

  // read baseline projections
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(forwardName);
  reader->Update();
  FloatVolumeType::Pointer forward = reader->GetOutput();
  FloatVolumeType::SizeType projSize = forward->GetLargestPossibleRegion().GetSize();

  if (projSize[0] != geometry.GetDetectorSamplesWidth() ||
      projSize[1] != geometry.GetDetectorSamplesHeight() ||
      projSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input forward projection dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }


  MarkovContext markovContext(geometry, materials, samplesPerPixel, voxelStepSize, 0);
  markovContext.SetBaselineProjection(forward);
  markovContext.FindInitialGuess();

  ByteVolumeType::Pointer outVol = ByteVolumeType::New();
  markovContext.GetCurrentVolume(outVol);

  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}


int RunMarkovScaleVolume(int argc, char **argv) {
  std::cerr<<"Running Markov Scale Volume"<<std::endl;
  
  if (argc != 6) {
    std::cerr<<"usage: backscatter markov_scale_volume <geometry_from.txt> <geometry_to.txt> <vol_from.nrrd> <vol_to.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryFromName = argv[2];
  char *geometryToName = argv[3];
  char *volFromName = argv[4];
  char *volToName = argv[5];

  Geometry fromGeometry, toGeometry;
  fromGeometry.LoadFromFile(geometryFromName);
  toGeometry.LoadFromFile(geometryToName);

  // read volume to sample from
  ByteReaderType::Pointer reader = ByteReaderType::New();
  reader->SetFileName(volFromName);
  reader->Update();
  ByteVolumeType::Pointer volume = reader->GetOutput();
  ByteVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != fromGeometry.GetVolumeNodeSamplesX() ||
      volumeSize[1] != fromGeometry.GetVolumeNodeSamplesY() ||
      volumeSize[2] != fromGeometry.GetVolumeNodeSamplesZ()) {
    std::cerr<<"input volume dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }


  ByteVolumeType::Pointer outVol = ByteVolumeType::New();
  ByteVolumeType::SizeType size;
  size.SetElement(0, toGeometry.GetVolumeNodeSamplesX());
  size.SetElement(1, toGeometry.GetVolumeNodeSamplesY());
  size.SetElement(2, toGeometry.GetVolumeNodeSamplesZ());
  outVol->SetRegions(ByteVolumeType::RegionType(size));
  outVol->Allocate();


  for (int z=0; z<toGeometry.GetVolumeNodeSamplesZ(); z++) {
    int fz = (int)(((z+0.5)/toGeometry.GetVolumeNodeSamplesZ()) * fromGeometry.GetVolumeNodeSamplesZ());
    for (int y=0; y<toGeometry.GetVolumeNodeSamplesY(); y++) {
      int fy = (int)(((y+0.5)/toGeometry.GetVolumeNodeSamplesY()) * fromGeometry.GetVolumeNodeSamplesY());
      for (int x=0; x<toGeometry.GetVolumeNodeSamplesX(); x++) {
        int fx = (int)(((x+0.5)/toGeometry.GetVolumeNodeSamplesX()) * fromGeometry.GetVolumeNodeSamplesX());

        int fromVi = fromGeometry.VolumeNodeCoordToIndex(fx, fy, fz);
        int toVi = toGeometry.VolumeNodeCoordToIndex(x, y, z);

        outVol->GetBufferPointer()[toVi] = volume->GetBufferPointer()[fromVi];
      }
    }
  }
  
  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(volToName);
  outWriter->Update();

  return 0;
}


int RunMarkovScaleProjections(int argc, char **argv) {
  std::cerr<<"Running Markov Scale Projections"<<std::endl;
  
  if (argc != 6) {
    std::cerr<<"usage: backscatter markov_scale_projections <geometry_from.txt> <geometry_to.txt> <vol_from.nrrd> <vol_to.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryFromName = argv[2];
  char *geometryToName = argv[3];
  char *volFromName = argv[4];
  char *volToName = argv[5];

  Geometry fromGeometry, toGeometry;
  fromGeometry.LoadFromFile(geometryFromName);
  toGeometry.LoadFromFile(geometryToName);

  // read volume to sample from
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volFromName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != fromGeometry.GetDetectorSamplesWidth() ||
      volumeSize[1] != fromGeometry.GetDetectorSamplesHeight() ||
      volumeSize[2] != fromGeometry.GetNumProjectionAngles()) {
    std::cerr<<"input volume dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }


  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  FloatVolumeType::SizeType size;
  size.SetElement(0, toGeometry.GetDetectorSamplesWidth());
  size.SetElement(1, toGeometry.GetDetectorSamplesHeight());
  size.SetElement(2, toGeometry.GetNumProjectionAngles());
  outVol->SetRegions(FloatVolumeType::RegionType(size));
  outVol->Allocate();


  if (2*toGeometry.GetDetectorSamplesHeight() == fromGeometry.GetDetectorSamplesHeight() &&
      2*toGeometry.GetDetectorSamplesWidth() == fromGeometry.GetDetectorSamplesWidth()) {
    
    for (int z=0; z<toGeometry.GetNumProjectionAngles(); z++) {
      for (int y=0; y<toGeometry.GetDetectorSamplesHeight(); y++) {
        for (int x=0; x<toGeometry.GetDetectorSamplesWidth(); x++) {

          int fromPi = fromGeometry.ProjectionCoordToIndex(z, 2*x, 2*y);
          int toPi = toGeometry.ProjectionCoordToIndex(z, x, y);

          outVol->GetBufferPointer()[toPi] = (volume->GetBufferPointer()[fromGeometry.ProjectionCoordToIndex(z, 2*x, 2*y)] +
                                              volume->GetBufferPointer()[fromGeometry.ProjectionCoordToIndex(z, 2*x+1, 2*y)] +
                                              volume->GetBufferPointer()[fromGeometry.ProjectionCoordToIndex(z, 2*x+1, 2*y+1)] +
                                              volume->GetBufferPointer()[fromGeometry.ProjectionCoordToIndex(z, 2*x, 2*y+1)]) * 0.25;
        }
      }
    }
  }
  else {
    for (int z=0; z<toGeometry.GetNumProjectionAngles(); z++) {
      for (int y=0; y<toGeometry.GetDetectorSamplesHeight(); y++) {
        int fy = (int)(((y+0.5)/toGeometry.GetDetectorSamplesHeight()) * fromGeometry.GetDetectorSamplesHeight());
        for (int x=0; x<toGeometry.GetDetectorSamplesWidth(); x++) {
          int fx = (int)(((x+0.5)/toGeometry.GetDetectorSamplesWidth()) * fromGeometry.GetDetectorSamplesWidth());

          int fromPi = fromGeometry.ProjectionCoordToIndex(z, fx, fy);
          int toPi = toGeometry.ProjectionCoordToIndex(z, x, y);

          outVol->GetBufferPointer()[toPi] = volume->GetBufferPointer()[fromPi];
        }
      }
    }
  }


  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(volToName);
  outWriter->Update();

  return 0;
}



int RunGibbs(int argc, char **argv) {
  std::cerr<<"Running Gibbs"<<std::endl;

  if (argc != 19) {
    std::cerr<<"usage: backscatter gibbs <geometry.txt> <air_material.txt> <foam_material.txt> <aluminum_material.txt> <forward.nrrd> <initial_vol.nrrd> <samples_per_pixel> <voxel_step_size> <num_iter> <regularization_weight> <start_temp> <end_temp> <numThreads>  <numSubThreads> <outputErrors.nrrd> <outputProbs.nrrd>  <outputRecon.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *airName = argv[3];
  char *foamName = argv[4];
  char *aluminumName = argv[5];
  char *forwardName = argv[6];
  char *initialName = argv[7];
  int samplesPerPixel = atoi(argv[8]);
  float voxelStepSize = atof(argv[9]);
  int niter = atoi(argv[10]);
  float regWeight = atof(argv[11]);
  float startTemp = atof(argv[12]);
  float endTemp = atof(argv[13]);
  int numThreads = atoi(argv[14]);
  int numSubThreads = atoi(argv[15]);
  char *outputErrorsName = argv[16];
  char *outputProbsName = argv[17];
  char *outputReconName = argv[18];


  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read materials
  vector<Material> materials(3);
  materials[0].InitFromFile(airName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[1].InitFromFile(foamName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);
  materials[2].InitFromFile(aluminumName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);

  // read baseline projections
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(forwardName);
  reader->Update();
  FloatVolumeType::Pointer forward = reader->GetOutput();
  FloatVolumeType::SizeType projSize = forward->GetLargestPossibleRegion().GetSize();

  if (projSize[0] != geometry.GetDetectorSamplesWidth() ||
      projSize[1] != geometry.GetDetectorSamplesHeight() ||
      projSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input forward projection dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }


  MarkovContext markovContext(geometry, materials, samplesPerPixel, voxelStepSize, regWeight);
  markovContext.SetBaselineProjection(forward);

  if (strcmp(initialName, "init")==0) {
    markovContext.SetCurrentVolume(1);
  }

  else {
    ByteReaderType::Pointer initreader = ByteReaderType::New();
    initreader->SetFileName(initialName);
    initreader->Update();
    ByteVolumeType::Pointer initvolume = initreader->GetOutput();
    ByteVolumeType::SizeType initvolumeSize = initvolume->GetLargestPossibleRegion().GetSize();

    if (initvolumeSize[0] != geometry.GetVolumeNodeSamplesX() ||
        initvolumeSize[1] != geometry.GetVolumeNodeSamplesY() ||
        initvolumeSize[2] != geometry.GetVolumeNodeSamplesZ()) {
      std::cerr<<"initial volume dimensions do not match geometry configuration file!"<<std::endl;
      return 1;
    }

    markovContext.SetCurrentVolume(initvolume);
  }

  markovContext.Gibbs(numThreads, numSubThreads, niter, startTemp, endTemp);

  // save reconstruction
  ByteVolumeType::Pointer outVol = ByteVolumeType::New();
  markovContext.GetCurrentVolume(outVol);

  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputReconName);
  outWriter->Update();


  // save material probabilities
  vector< vector<float> > matErrors;
  vector< vector<float> > matProbs;

  double T = startTemp/niter + endTemp * (niter-1)/niter;
  T = startTemp / log(endTemp*(niter-1)+2.0);

  markovContext.GetGibbsConditionalProbabilities(numThreads*numSubThreads, T, matErrors, matProbs);

  RGBVolumeType::Pointer rgbProbVol = RGBVolumeType::New();
  RGBVolumeType::SizeType rgbProbSize;
  rgbProbSize.SetElement(0, geometry.GetVolumeNodeSamplesX());
  rgbProbSize.SetElement(1, geometry.GetVolumeNodeSamplesY());
  rgbProbSize.SetElement(2, geometry.GetVolumeNodeSamplesZ());
  rgbProbVol->SetRegions(RGBVolumeType::RegionType(rgbProbSize));
  rgbProbVol->Allocate();

  RGBVolumeType::Pointer rgbErrorVol = RGBVolumeType::New();
  RGBVolumeType::SizeType rgbErrorSize;
  rgbErrorSize.SetElement(0, geometry.GetVolumeNodeSamplesX());
  rgbErrorSize.SetElement(1, geometry.GetVolumeNodeSamplesY());
  rgbErrorSize.SetElement(2, geometry.GetVolumeNodeSamplesZ());
  rgbErrorVol->SetRegions(RGBVolumeType::RegionType(rgbErrorSize));
  rgbErrorVol->Allocate();

  for (int m=0; m<materials.size(); m++) {

    FloatVolumeType::SizeType outProbSize;
    outProbSize.SetElement(0, geometry.GetVolumeNodeSamplesX());
    outProbSize.SetElement(1, geometry.GetVolumeNodeSamplesY());
    outProbSize.SetElement(2, geometry.GetVolumeNodeSamplesZ());

    FloatVolumeType::Pointer outProbVol = FloatVolumeType::New();
    outProbVol->SetRegions(FloatVolumeType::RegionType(outProbSize));
    outProbVol->Allocate();

    FloatVolumeType::Pointer outErrorVol = FloatVolumeType::New();
    outErrorVol->SetRegions(FloatVolumeType::RegionType(outProbSize));
    outErrorVol->Allocate();

    for (int vi=0; vi<geometry.GetTotalVolumeNodeSamples(); vi++) {
      outErrorVol->GetBufferPointer()[vi] = matErrors[m][vi];
      outProbVol->GetBufferPointer()[vi] = matProbs[m][vi];
      rgbErrorVol->GetBufferPointer()[vi].SetNthComponent(m, matErrors[m][vi]);
      rgbProbVol->GetBufferPointer()[vi].SetNthComponent(m, matProbs[m][vi]*255);
    }

    FloatWriterType::Pointer outProbWriter = FloatWriterType::New();
    outProbWriter->SetInput(outProbVol);
    char outputProbsName2[1024];
    sprintf_s(outputProbsName2, 1024, "%s_%d.nrrd", outputProbsName, m);
    outProbWriter->SetFileName(outputProbsName2);
    outProbWriter->Update();

    FloatWriterType::Pointer outErrorWriter = FloatWriterType::New();
    outErrorWriter->SetInput(outErrorVol);
    char outputErrorsName2[1024];
    sprintf_s(outputErrorsName2, 1024, "%s_%d.nrrd", outputErrorsName, m);
    outErrorWriter->SetFileName(outputErrorsName2);
    outErrorWriter->Update();
  }


  RGBWriterType::Pointer rgbProbWriter = RGBWriterType::New();
  rgbProbWriter->SetInput(rgbProbVol);
  char outputProbsName2[1024];
  sprintf_s(outputProbsName2, 1024, "%s.nrrd", outputProbsName);
  rgbProbWriter->SetFileName(outputProbsName2);
  rgbProbWriter->Update();

  RGBWriterType::Pointer rgbErrorWriter = RGBWriterType::New();
  rgbErrorWriter->SetInput(rgbErrorVol);
  char outputErrorsName2[1024];
  sprintf_s(outputErrorsName2, 1024, "%s.nrrd", outputErrorsName);
  rgbErrorWriter->SetFileName(outputErrorsName2);
  rgbErrorWriter->Update();

  return 0;
}




int RunSARTForwardProject(int argc, char **argv) {
  std::cerr<<"Running SART Forward Project"<<std::endl;

  if (argc != 9) {
    std::cerr<<"usage: backscatter forward_project <geometry.txt> <foam_material.txt> <volume.nrrd> <sample_step> <samples_per_pixel> <numThreads> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *foamName = argv[3];
  char *volumeName = argv[4];
  float sampleSize = atof(argv[5]);
  int samplesPerPixel = atoi(argv[6]);
  int numThreads = atoi(argv[7]);
  char *outputName = argv[8];


  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read foam material configuration
  Material foamMaterial;
  foamMaterial.InitFromFile(foamName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);

  // read volume to project from
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != geometry.GetVolumeSamplesX() ||
      volumeSize[1] != geometry.GetVolumeSamplesY() ||
      volumeSize[2] != geometry.GetVolumeSamplesZ()) {
    std::cerr<<"input volume dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  SARTContext sartContext(geometry, foamMaterial);
  sartContext.SetCurrentVolume(volume);
  sartContext.CreateWeights(sampleSize, samplesPerPixel, numThreads);
  sartContext.ForwardProject();

  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  sartContext.GetCurrentProjection(outVol);

  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();


  return 0;
}



int RunSARTReconstruct(int argc, char **argv) {
  std::cerr<<"Running SART Reconstruction"<<std::endl;

  if (argc != 13) {
    std::cerr<<"usage: backscatter reconstruct <geometry.txt> <foam_material.txt> <projections.nrrd> <initial_value> <sample_step> <samples_per_pixel> <lambda> <total_iter> <weight_iter> <num_threads> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *foamName = argv[3];
  char *volumeName = argv[4];
  float initialValue = atof(argv[5]);
  float sampleSize = atof(argv[6]);
  int samplesPerPixel = atoi(argv[7]);
  float lambda = atof(argv[8]);
  int totalIter = atoi(argv[9]);
  int weightIter = atoi(argv[10]);
  int numThreads = atoi(argv[11]);
  char *outputName = argv[12];


  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read foam material configuration
  Material foamMaterial;
  foamMaterial.InitFromFile(foamName, XRAY_ENERGY_LOW, XRAY_ENERGY_HIGH);


  // read volume to project from
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != geometry.GetDetectorSamplesWidth() ||
      volumeSize[1] != geometry.GetDetectorSamplesHeight() ||
      volumeSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input projection dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  SARTContext sartContext(geometry, foamMaterial);
  sartContext.SetBaselineProjection(volume);
  sartContext.SetCurrentVolume(initialValue);

  for (int i=0; i<totalIter; i++) {
    if (i%weightIter == 0)
      sartContext.CreateWeights(sampleSize, samplesPerPixel, numThreads);

    sartContext.ForwardProject();
    sartContext.BackwardProject(lambda, samplesPerPixel);


    // save forward projection
    FloatVolumeType::Pointer outVolForwardProj = FloatVolumeType::New();
    sartContext.GetCurrentProjection(outVolForwardProj);

    char forwardName[1024];
    sprintf_s(forwardName, 1024, "forward_%d.nrrd", i);
    FloatWriterType::Pointer outWriterForwardProj = FloatWriterType::New();
    outWriterForwardProj->SetInput(outVolForwardProj);
    outWriterForwardProj->SetFileName(forwardName);
    outWriterForwardProj->Update();


    // save delta r
    FloatVolumeType::Pointer outVolForwardProjDR = FloatVolumeType::New();
    sartContext.GetCurrentDeltaR(outVolForwardProjDR);

    char forwardNameDR[1024];
    sprintf_s(forwardNameDR, 1024, "deltar_%d.nrrd", i);
    FloatWriterType::Pointer outWriterForwardProjDR = FloatWriterType::New();
    outWriterForwardProjDR->SetInput(outVolForwardProjDR);
    outWriterForwardProjDR->SetFileName(forwardNameDR);
    outWriterForwardProjDR->Update();


    // save back projection
    FloatVolumeType::Pointer outVolBackProj = FloatVolumeType::New();
    sartContext.GetCurrentVolume(outVolBackProj);

    char backName[1024];
    sprintf_s(backName, 1024, "back_%d.nrrd", i);
    FloatWriterType::Pointer outWriterBackProj = FloatWriterType::New();
    outWriterBackProj->SetInput(outVolBackProj);
    outWriterBackProj->SetFileName(backName);
    outWriterBackProj->Update();
  }

  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  sartContext.GetCurrentVolume(outVol);

  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}


int ScaleVolume(int argc, char **argv) {

  std::cerr<<"Scaling Volume"<<std::endl;

  if (argc != 6) {
    std::cerr<<"usage: backscatter scale_vol <in.nrrd> <scale> <offset> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *volumeName = argv[2];
  float scale = atof(argv[3]);
  float offset = atof(argv[4]);
  char *outputName = argv[5];


  // read volume
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  itk::ImageRegionIterator<FloatVolumeType> it(volume,
                                               volume->GetLargestPossibleRegion());
  for (it=it.Begin(); !it.IsAtEnd(); ++it)
    it.Set(it.Value() * scale + offset);


  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(volume);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}



int PoissonCorrupt(int argc, char **argv) {

  std::cerr<<"Poisson"<<std::endl;

  if (argc != 5) {
    std::cerr<<"usage: backscatter poisson <in.nrrd> <scale> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *volumeName = argv[2];
  float scale = atof(argv[3]);
  char *outputName = argv[4];


  // read volume
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  itk::ImageRegionIterator<FloatVolumeType> it(volume,
                                               volume->GetLargestPossibleRegion());
  for (it=it.Begin(); !it.IsAtEnd(); ++it) {
    Poissondev pd(it.Value() * scale, rand());
    it.Set((float)pd.dev() / scale);
  }


  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(volume);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}


int GenerateSynthetic(int argc, char **argv) {
  std::cerr<<"Generating synthetic volume"<<std::endl;
  if (argc != 5) {
    std::cerr<<"usage: backscatter gen_synth <geometry.txt> <outdensity.nrrd>  <outmaterial.nrrd>" << std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *outputDensityName = argv[3];
  char *outputMaterialName = argv[4];

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // init output density volume
  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  FloatVolumeType::SizeType size;
  size.SetElement(0, geometry.GetVolumeNodeSamplesX());
  size.SetElement(1, geometry.GetVolumeNodeSamplesY());
  size.SetElement(2, geometry.GetVolumeNodeSamplesZ());
  outVol->SetRegions(FloatVolumeType::RegionType(size));
  outVol->Allocate();

  ByteVolumeType::Pointer matVol = ByteVolumeType::New();
  ByteVolumeType::SizeType matsize;
  matsize.SetElement(0, geometry.GetVolumeNodeSamplesX());
  matsize.SetElement(1, geometry.GetVolumeNodeSamplesY());
  matsize.SetElement(2, geometry.GetVolumeNodeSamplesZ());
  matVol->SetRegions(ByteVolumeType::RegionType(matsize));
  matVol->Allocate();


  std::vector<Vec3f> voids;
  voids.push_back(Vec3f(-0.0508, -0.0508, -0.0408));
  voids.push_back(Vec3f(-0.0508,  0.0508, -0.0408));
  voids.push_back(Vec3f(0.0508, -0.0508, 0));
  voids.push_back(Vec3f(0.0508,  0.0508, 0));

  int highestAlZ = 0;
  //for (int iz=0; iz<geometry.GetVolumeNodeSamplesZ(); iz++) {
  for (int iz=geometry.GetVolumeNodeSamplesZ()-1; iz>=0; iz--) {
    for (int iy=0; iy<geometry.GetVolumeNodeSamplesY(); iy++) {
      for (int ix=0; ix<geometry.GetVolumeNodeSamplesX(); ix++) {

        Vec3f pixelcoord(ix, iy, iz);
        Vec3f worldcoord;
        geometry.VolumeToWorld(pixelcoord, worldcoord);

        // determine which material type this is
        int mat = 1; // default to foam

#if 0
        // check for air in voids
        for (size_t v=0; v<voids.size(); v++) {
          //if ((voids[v] - worldcoord).Length() < 0.02) // diameter 40mm
          if ((voids[v] - worldcoord).Length() < 0.025)
            mat = 0;
        }

        // the top half inch is air
        if (worldcoord[2] > 0.0508)
          mat = 0;

        // alluminum on bottom
        if (iz<=1 || worldcoord[2] < -0.0448)
          mat = 2;

        // aluminum bar
        if (((worldcoord[0] - -0.03)*(worldcoord[0] - -0.03) + (worldcoord[2] - -0.03)*(worldcoord[2] - -0.03)) < (0.01*0.01))
          mat = 2;

        if (iz < 2)
          mat = 2;

#else

        worldcoord[2] -= geometry.GetVolumeCenter()[2];

        // cross void
        if (worldcoord[0] > -1.5*2.54/100 && 
            worldcoord[0] <  1.5*2.54/100 && 
            worldcoord[1] > -0.5*2.54/100 && 
            worldcoord[1] <  0.5*2.54/100 &&
            worldcoord[2] > -0.5*2.69/100 && 
            worldcoord[2] <  0.5*2.69/100)
          mat = 0;

        if (worldcoord[0] > -0.5*2.54/100 && 
            worldcoord[0] <  0.5*2.54/100 && 
            worldcoord[1] > -1.5*2.54/100 && 
            worldcoord[1] <  1.5*2.54/100 &&
            worldcoord[2] > -0.5*2.69/100 && 
            worldcoord[2] <  0.5*2.69/100 )
          mat = 0;

        // square block of foam on aluminum
        if (worldcoord[2] < -0.5*5.39/100)
          mat = 2;

        if ((worldcoord[0] < -0.1525 ||
             worldcoord[0] >  0.1525 ||
             worldcoord[1] < -0.1525 ||
             worldcoord[1] >  0.1525 ||
             worldcoord[2] >  0.5*5.39/100 ||
             iz < highestAlZ))
          mat = 0;

        if (mat == 2)
          highestAlZ = std::max(highestAlZ, iz);


#endif


        float density = 0;
        if (mat == 1) // foam
          density = 0.0384443121 * 100*100*100;
        else if (mat == 0) // air
          density = 0.00122521 * 100*100*100;
        else // aluminum
          density = 2.7 * 100*100*100;

        outVol->GetBufferPointer()[geometry.VolumeNodeCoordToIndex(ix,iy,iz)] = density;
        matVol->GetBufferPointer()[geometry.VolumeNodeCoordToIndex(ix,iy,iz)] = mat;

      }
    }
  }

  // save it
  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputDensityName);
  outWriter->Update();

  ByteWriterType::Pointer matWriter = ByteWriterType::New();
  matWriter->SetInput(matVol);
  matWriter->SetFileName(outputMaterialName);
  matWriter->Update();
  
  return 0;
}



int SegmentCalibPoints(int argc, char **argv) {
  std::cerr<<"Segment Calibration Points"<<std::endl;

  if (argc < 10) {
    std::cerr<<"usage: backscatter segment_calib_points <geometry.txt> <forward.nrrd> <threshold> <mask_size/2> <id_size/2> <border> <light/dark> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *volumeName = argv[3];
  int numPoints = atof(argv[4]);
  int maskSize = atof(argv[5]);
  int idSize = atof(argv[6]);
  int border = atof(argv[7]);
  bool dark = (stricmp(argv[8], "dark")==0);
  char *outputName = argv[9];

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read input volume
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != geometry.GetDetectorSamplesWidth() ||
      volumeSize[1] != geometry.GetDetectorSamplesHeight() ||
      volumeSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input projection dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  // create output volume
  ByteVolumeType::Pointer idVol = ByteVolumeType::New();
  ByteVolumeType::SizeType size;
  size.SetElement(0, geometry.GetDetectorSamplesWidth());
  size.SetElement(1, geometry.GetDetectorSamplesHeight());
  size.SetElement(2, geometry.GetNumProjectionAngles());
  idVol->SetRegions(ByteVolumeType::RegionType(size));
  idVol->Allocate();


  // copy input into vector
  vector<float> forward(geometry.GetTotalProjectionSamples());
  for (size_t i=0; i<forward.size(); i++) {
    forward[i] = volume->GetBufferPointer()[i];
  }

  vector<unsigned char> ids;
  SegmentCalibPoints(geometry.GetDetectorSamplesWidth(),
                     geometry.GetDetectorSamplesHeight(),
                     geometry.GetNumProjectionAngles(),
                     forward, ids, numPoints, maskSize, idSize, border, dark);

  for (size_t i=0; i<ids.size(); i++) {
    idVol->GetBufferPointer()[i] = ids[i];
  }


  // save it
  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(idVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}




int RunLocalizeCalibPoints(int argc, char **argv) {
  std::cerr<<"Localizing Calibration Points"<<std::endl;

  if (argc < 7) {
    std::cerr<<"usage: backscatter localize_calib_points <geometry.txt> <forward.nrrd> <mask.nrrd> <light/dark> <output.txt>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *volumeName = argv[3];
  char *maskName = argv[4];
  bool dark = (stricmp(argv[5], "dark")==0);
  char *outputName = argv[6];

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);

  // read input volume
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  if (volumeSize[0] != geometry.GetDetectorSamplesWidth() ||
      volumeSize[1] != geometry.GetDetectorSamplesHeight() ||
      volumeSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input projection dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  // read mask volume
  ByteReaderType::Pointer maskreader = ByteReaderType::New();
  maskreader->SetFileName(maskName);
  maskreader->Update();
  ByteVolumeType::Pointer mask = maskreader->GetOutput();
  ByteVolumeType::SizeType maskVolumeSize = mask->GetLargestPossibleRegion().GetSize();

  if (maskVolumeSize[0] != geometry.GetDetectorSamplesWidth() ||
      maskVolumeSize[1] != geometry.GetDetectorSamplesHeight() ||
      maskVolumeSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input mask dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }

  // copy input into vectors
  vector<float> forward(geometry.GetTotalProjectionSamples());
  for (size_t i=0; i<forward.size(); i++) {
    forward[i] = volume->GetBufferPointer()[i];
  }

  vector<unsigned char> maskv(geometry.GetTotalProjectionSamples());
  for (size_t i=0; i<maskv.size(); i++) {
    maskv[i] = mask->GetBufferPointer()[i];
  }


  // run localization
  vector<int> proj_i, point_i;
  vector<float> point_x, point_y;
  vector<float> residuals;
  LocalizeCalibPoints(geometry.GetDetectorSamplesWidth(),
                      geometry.GetDetectorSamplesHeight(),
                      geometry.GetNumProjectionAngles(),
                      forward, maskv,
                      proj_i, point_i, point_x, point_y,
                      residuals,
                      dark);

  std::cerr<<"found "<<proj_i.size()<<"points"<<std::endl;


  // write output
  FILE *fout=NULL;
  fopen_s(&fout, outputName, "w");
  if (!fout) {
    printf("error opening file %s for output\n", outputName);
    return 1;
  }
  for (size_t i=0; i<proj_i.size(); i++) {
    fprintf(fout, "%d, %d, %g, %g\n", proj_i[i], point_i[i]-1, point_x[i], point_y[i]);
  }
  fclose(fout);


  // create a new mask volume for debugging with point locations marked
  ByteVolumeType::Pointer outVol = ByteVolumeType::New();
  ByteVolumeType::SizeType size;
  size.SetElement(0, geometry.GetDetectorSamplesWidth());
  size.SetElement(1, geometry.GetDetectorSamplesHeight());
  size.SetElement(2, geometry.GetNumProjectionAngles());
  outVol->SetRegions(ByteVolumeType::RegionType(size));
  outVol->Allocate();

  for (size_t i=0; i<maskv.size(); i++) {
    outVol->GetBufferPointer()[i] = 0;
  }

  for (size_t p=0; p<proj_i.size(); p++) {
    int x = std::max(0, std::min(geometry.GetDetectorSamplesWidth()-1, (int)point_x[p]));
    int y = std::max(0, std::min(geometry.GetDetectorSamplesHeight()-1, (int)point_y[p]));

    outVol->GetBufferPointer()[proj_i[p] * (geometry.GetDetectorSamplesWidth() * geometry.GetDetectorSamplesHeight()) +
                               y * geometry.GetDetectorSamplesWidth() + x] = point_i[p];
  }

  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(outVol);

  char outvolname[1024];
  sprintf_s(outvolname, 1024, "%s.id.nrrd", outputName);
  outWriter->SetFileName(outvolname);
  outWriter->Update();


  // create a volume for the residual of each particle fit
  FloatVolumeType::Pointer routVol = FloatVolumeType::New();
  FloatVolumeType::SizeType rsize;
  rsize.SetElement(0, geometry.GetDetectorSamplesWidth());
  rsize.SetElement(1, geometry.GetDetectorSamplesHeight());
  rsize.SetElement(2, geometry.GetNumProjectionAngles());
  routVol->SetRegions(FloatVolumeType::RegionType(rsize));
  routVol->Allocate();

  for (size_t i=0; i<maskv.size(); i++) {
    routVol->GetBufferPointer()[i] = 0;
  }

  for (size_t p=0; p<proj_i.size(); p++) {
    int x = std::max(0, std::min(geometry.GetDetectorSamplesWidth()-1, (int)point_x[p]));
    int y = std::max(0, std::min(geometry.GetDetectorSamplesHeight()-1, (int)point_y[p]));

    routVol->GetBufferPointer()[proj_i[p] * (geometry.GetDetectorSamplesWidth() * geometry.GetDetectorSamplesHeight()) +
                                y * geometry.GetDetectorSamplesWidth() + x] = residuals[p];
  }

  FloatWriterType::Pointer routWriter = FloatWriterType::New();
  routWriter->SetInput(routVol);

  char routvolname[1024];
  sprintf_s(routvolname, 1024, "%s.res.nrrd", outputName);
  routWriter->SetFileName(routvolname);
  routWriter->Update();

  return 0;
}


int RunCalibrateDetector(int argc, char **argv) {
  std::cerr<<"Calibrating Detector"<<std::endl;

  if (argc < 12 || (argc-9)%3 != 0) {
    std::cerr<<"usage: backscatter calibrate_detector <geometry.txt> <points.txt> <foreground> <background> <radius> <px py pz ...> <out_geometry.txt> <output_vol.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *pointsName = argv[3];
  int background = atoi(argv[4]);
  int foreground = atoi(argv[5]);
  float radius = atof(argv[6]);
  char *outputName = argv[argc-2];
  char *outputVolName = argv[argc-1];



  vector<Vec3f> cpoints;
  for (int arg=7; arg<argc-2; arg+=3) {
    float px = atof(argv[arg+0]);
    float py = atof(argv[arg+1]);
    float pz = atof(argv[arg+2]);
    Vec3f p(px,py,pz);

    cpoints.push_back(p);
  }

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);


  // read the localized point positions
  std::ifstream file;
  file.open(pointsName);
  if (!file.is_open()) {
    printf("error opening points file %s\n", pointsName);
    return 1;
  }

  
  vector<int> proj_i, point_i;
  vector<float> point_x, point_y;
  char line[1024];
  while (file.getline(line, 1024)) {
    int tproj_i, tpoint_i;
    float tpoint_x, tpoint_y;

    if (sscanf_s(line, "%d, %d, %g, %g", &tproj_i, &tpoint_i, &tpoint_x, &tpoint_y) == 4) {
      proj_i.push_back(tproj_i);
      point_i.push_back(tpoint_i);
      point_x.push_back(tpoint_x);
      point_y.push_back(tpoint_y);
    }
  }
  file.close();

  // do the linear algebra
  Vec3f baseCenter, baseX, baseY, baseZ;
  CalibrateDetector(geometry, proj_i, point_i, point_x, point_y, cpoints,
                    baseCenter, baseX, baseY, baseZ);

  // save the new geometry file
  geometry.SaveToFile(outputName);


  //
  // generate a matching synthetic volume
  //

  // init output volume
  ByteVolumeType::Pointer matVol = ByteVolumeType::New();
  ByteVolumeType::SizeType bsize;
  bsize.SetElement(0, geometry.GetVolumeNodeSamplesX());
  bsize.SetElement(1, geometry.GetVolumeNodeSamplesY());
  bsize.SetElement(2, geometry.GetVolumeNodeSamplesZ());
  matVol->SetRegions(ByteVolumeType::RegionType(bsize));
  matVol->Allocate();

  // set background
  for (int v=0; v<geometry.GetTotalVolumeNodeSamples(); v++) {
    matVol->GetBufferPointer()[v] = background;
  }

  // create beads
  for (int c=0; c<cpoints.size(); c++) {
    Vec3f p = (baseCenter + 
               Vec3f(baseX.Dot(cpoints[c]),
                     baseY.Dot(cpoints[c]),
                     baseZ.Dot(cpoints[c])));


    for (int v=0; v<geometry.GetTotalVolumeNodeSamples(); v++) {
      int vx, vy, vz;
      geometry.VolumeIndexToNodeCoord(v, vx, vy, vz);

      Vec3f world;
      geometry.VolumeToWorld(Vec3f(vx, vy, vz), world);

      Vec3f diff = world-p;
      diff[2] *= 3;
      if (diff.Length() < radius) {
        matVol->GetBufferPointer()[v] = foreground;
      }
    }
  }

  // save it
  ByteWriterType::Pointer matWriter = ByteWriterType::New();
  matWriter->SetInput(matVol);
  matWriter->SetFileName(outputVolName);
  matWriter->Update();


  /*
  // project each of the beads at each angle
  FILE *fideal=NULL;
  fopen_s(&fideal, idealName, "w");
  if (!fideal) {
    printf("error opening file %s for output\n", idealName);
    return 1;
  }
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int c=0; c<(int)beads.size(); c++) {

      Vec3f ir = geometry.GetInverseRotatedPoint(beads[c], p);
      Vec3f irrel = ir - geometry.GetDetectorCenter();
      
      float mx = irrel.Dot(geometry.GetDetectorRight());
      float my = irrel.Dot(geometry.GetDetectorUp());

      float px, py;
      geometry.ProjectionPhysicalOffsetToPixel(mx, my, px, py);

      fprintf(fideal, "%d, %d, %g, %g\n", p, c, px, py);
    }
  }
  fclose(fideal);
  */


  return 0;
}


int RunCalibrateIntensity(int argc, char **argv) {
  std::cerr<<"Calibrating Detector Gain/Offset"<<std::endl;

  if (argc != 6) {
    std::cerr<<"usage: backscatter calibrate_intensity <geometry.txt> <detector_images.nrrd> <forward.nrrd> <out_geometry.txt>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *imagesName = argv[3];
  char *forwardName = argv[4];
  char *outputName = argv[5];

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);


  // read acquired images
  FloatReaderType::Pointer imreader = FloatReaderType::New();
  imreader->SetFileName(imagesName);
  imreader->Update();
  FloatVolumeType::Pointer images = imreader->GetOutput();
  FloatVolumeType::SizeType improjSize = images->GetLargestPossibleRegion().GetSize();

  if (improjSize[0] != geometry.GetDetectorSamplesWidth() ||
      improjSize[1] != geometry.GetDetectorSamplesHeight() ||
      improjSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input image dimensions do not match geometry configuration file!"<<std::endl;
    return 1;
  }


  // read our forward projections
  FloatReaderType::Pointer freader = FloatReaderType::New();
  freader->SetFileName(forwardName);
  freader->Update();
  FloatVolumeType::Pointer forward = freader->GetOutput();
  FloatVolumeType::SizeType fprojSize = forward->GetLargestPossibleRegion().GetSize();

  if (fprojSize[0] != geometry.GetDetectorSamplesWidth() ||
      fprojSize[1] != geometry.GetDetectorSamplesHeight() ||
      fprojSize[2] != geometry.GetNumProjectionAngles()) {
    std::cerr<<"input forward projection dimensions do not match geometry configuration file!"<<std::endl;
    std::cerr<<fprojSize[0]<<std::endl;
    std::cerr<<fprojSize[1]<<std::endl;
    std::cerr<<fprojSize[2]<<std::endl;
    return 1;
  }


  // do the linear algebra
  CalibrateIntensity(geometry, images->GetBufferPointer(), forward->GetBufferPointer());

  // save the new geometry file
  geometry.SaveToFile(outputName);

  return 0;
}


int MapVolume(int argc, char **argv) {
  std::cerr<<"Mapping Volume Values"<<std::endl;

  if (argc < 5) {
    std::cerr<<"usage: backscatter map_volume <input.nrrd> <map_file> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *volumeName = argv[2];
  char *mapName = argv[3];
  char *outputName = argv[4];

  // read input volume
  ByteReaderType::Pointer reader = ByteReaderType::New();
  reader->SetFileName(volumeName);
  reader->Update();
  ByteVolumeType::Pointer volume = reader->GetOutput();
  ByteVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  // read value map
  vector<int> map(256);
  for (int i=0; i<256; i++)
    map[i] = i;

  FILE *mapf=NULL;
  fopen_s(&mapf, mapName, "r");
  if (!mapf) {
    printf("error opening map file %s\n", mapName);
    return 1;
  }

  char line[1024];
  while (fgets(line, 1024, mapf)) {
    if (strchr(line, '\n'))
      *strchr(line, '\n') = 0;
    if (strchr(line, '\r'))
      *strchr(line, '\r') = 0;

    if (strchr(line, ' ')) {
      char *space = strchr(line, ' ');
      *space = 0;
      space++;

      int map_from = atoi(line);
      int map_to = atoi(space);

      map[map_from] = map_to;
    }
  }
  fclose(mapf);


  // map input values into output
  ByteVolumeType::Pointer outVol = ByteVolumeType::New();
  outVol->SetRegions(FloatVolumeType::RegionType(volumeSize));
  outVol->Allocate();

  for (int i=0; i<volumeSize[0]*volumeSize[1]*volumeSize[2]; i++) {
    outVol->GetBufferPointer()[i] = map[volume->GetBufferPointer()[i]];
  }


  // save it
  ByteWriterType::Pointer outWriter = ByteWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}



int ProjectionsToNrrd(int argc, char **argv) {

  std::cerr<<"Converting projection .tiff files to .nrrd"<<std::endl;

  if (argc < 6) {
    std::cerr<<"usage: backscatter projections_to_nrrd <geometry.txt> <input_prefix> <index_style> <output.nrrd>"<<std::endl;
    return 1;
  }

  char *geometryName = argv[2];
  char *inputPrefix = argv[3];
  bool indexAngle = (stricmp(argv[4], "angle")==0);
  char *outputName = argv[5];

  // read geometry configuration
  Geometry geometry;
  geometry.LoadFromFile(geometryName);


  uint32 width=0;
  uint32 height=0;
  vector< vector<float> > images;
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    images.resize(p+1);
    vector<float> &pimage = images[p];

    char projName[1024];
    if (!indexAngle)
      sprintf(projName, "%s%d.tiff", inputPrefix, (int)(0.5 + geometry.GetProjectionAngle(p) * 180 / M_PI)/10);
    else
      sprintf(projName, "%s%03d.tiff", inputPrefix, (int)(0.5 + geometry.GetProjectionAngle(p) * 180 / M_PI)+10);

    TIFF *tiff = TIFFOpen(projName, "r");
    if (!tiff) {
      std::cerr<<"error opening file "<<projName<<std::endl;
      return 1;
    }
    TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &height);

    if (width!=1024 || height!=1024) {
      std::cerr<<"expected input image dimensions to be 1024x1024, found "<<width<<"x"<<height<<std::endl;
    }

    tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tiff));
    for (int row=0; row<height; row++) {
	    TIFFReadScanline(tiff, buf, row);

      for (int col=0; col<width; col++) {
        pimage.push_back(((float*)buf)[col]);
      }
    }
    _TIFFfree(buf);
    TIFFClose(tiff);
  }

  // remove stuck pixels
  vector<int> mask(1024*1024, 0);
  mask[1024*12+49] = 1;
  mask[1024*32+24] = 1;
  mask[1024*12+82] = 1;
  mask[1024*29+211] = 1;
  mask[1024*19+243] = 1;
  mask[1024*20+335] = 1;
  mask[1024*50+460] = 1;
  mask[1024*41+712] = 1;
  mask[1024*22+803] = 1;
  mask[1024*52+57] = 1;
  mask[1024*72+49] = 1;
  mask[1024*76+111] = 1;
  mask[1024*122+60] = 1;
  mask[1024*101+120] = 1;
  mask[1024*93+145] = 1;
  mask[1024*120+133] = 1;
  mask[1024*131+147] = 1;
  mask[1024*77+265] = 1;
  mask[1024*79+458] = 1;
  mask[1024*164+574] = 1;
  mask[1024*182+611] = 1;
  mask[1024*191+727] = 1;
  mask[1024*151+908] = 1;
  mask[1024*168+284] = 1;
  mask[1024*189+289] = 1;
  mask[1024*217+325] = 1;
  mask[1024*268+24] = 1;
  mask[1024*260+46] = 1;
  mask[1024*257+130] = 1;
  mask[1024*292+192] = 1;
  mask[1024*292+193] = 1;
  mask[1024*310+59] = 1;
  mask[1024*286+249] = 1;
  mask[1024*253+333] = 1;
  mask[1024*217+325] = 1;
  mask[1024*229+310] = 1;
  mask[1024*212+377] = 1;
  mask[1024*296+571] = 1;
  mask[1024*330+593] = 1;
  mask[1024*422+686] = 1;
  mask[1024*409+716] = 1;
  mask[1024*431+727] = 1;
  mask[1024*432+727] = 1;
  mask[1024*401+780] = 1;
  mask[1024*379+998] = 1;
  mask[1024*485+579] = 1;
  mask[1024*475+613] = 1;
  mask[1024*500+597] = 1;
  mask[1024*532+658] = 1;
  mask[1024*544+681] = 1;
  mask[1024*603+816] = 1;
  mask[1024*593+982] = 1;
  mask[1024*772+29] = 1;
  mask[1024*763+140] = 1;
  mask[1024*781+113] = 1;
  mask[1024*802+18] = 1;
  mask[1024*804+146] = 1;
  mask[1024*846+11] = 1;
  mask[1024*858+188] = 1;
  mask[1024*890+58] = 1;
  mask[1024*936+61] = 1;
  mask[1024*953+59] = 1;
  mask[1024*970+56] = 1;
  mask[1024*985+17] = 1;
  mask[1024*994+17] = 1;
  mask[1024*1002+29] = 1;
  mask[1024*1011+157] = 1;
  mask[1024*994+202] = 1;
  mask[1024*964+204] = 1;
  mask[1024*583+104] = 1;
  mask[1024*1023+93] = 1;
  mask[1024*661+78] = 1;
  mask[1024*611+174] = 1;
  mask[1024*642+199] = 1;
  mask[1024*642+205] = 1;
  mask[1024*656+260] = 1;
  mask[1024*707+257] = 1;
  mask[1024*571+412] = 1;
  mask[1024*584+405] = 1;
  mask[1024*756+546] = 1;
  mask[1024*737+639] = 1;
  mask[1024*724+801] = 1;
  mask[1024*724+802] = 1;
  mask[1024*723+802] = 1;
  mask[1024*728+844] = 1;
  mask[1024*729+844] = 1;
  mask[1024*730+844] = 1;
  mask[1024*741+972] = 1;
  mask[1024*825+979] = 1;
  mask[1024*775+949] = 1;
  mask[1024*776+949] = 1;
  mask[1024*777+949] = 1;
  mask[1024*871+802] = 1;
  mask[1024*872+802] = 1;
  mask[1024*871+936] = 1;
  mask[1024*872+936] = 1;
  mask[1024*871+937] = 1;
  mask[1024*872+937] = 1;
  mask[1024*882+836] = 1;
  mask[1024*883+836] = 1;
  mask[1024*990+823] = 1;
  mask[1024*884+732] = 1;
  mask[1024*819+788] = 1;
  mask[1024*940+995] = 1;
  mask[1024*511+195] = 1;
  mask[1024*512+195] = 1;
  mask[1024*337+18] = 1;
  mask[1024*352+27] = 1;
  mask[1024*332+117] = 1;
  mask[1024*404+177] = 1;
  mask[1024*144+13] = 1;
  mask[1024*84+85] = 1;
  mask[1024*151+48] = 1;
  mask[1024*151+78] = 1;
  mask[1024*171+95] = 1;
  mask[1024*172+99] = 1;
  mask[1024*217+41] = 1;
  mask[1024*223+119] = 1;
  mask[1024*241+101] = 1;
  mask[1024*258+115] = 1;
  mask[1024*280+136] = 1;
  mask[1024*330+166] = 1;
  mask[1024*229+212] = 1;
  mask[1024*177+159] = 1;
  mask[1024*19+167] = 1;
  mask[1024*15+197] = 1;
  mask[1024*46+244] = 1;
  mask[1024*348+307] = 1;
  mask[1024*325+375] = 1;
  mask[1024*426+385] = 1;
  mask[1024*436+402] = 1;
  mask[1024*532+311] = 1;
  mask[1024*537+352] = 1;
  mask[1024*622+93] = 1;
  mask[1024*617+106] = 1;
  mask[1024*687+220] = 1;
  mask[1024*705+250] = 1;
  mask[1024*707+256] = 1;
  for (int y=707; y<=712; y++)
    mask[1024*y+257] = 1;
  mask[1024*622+93] = 1;
  mask[1024*617+106] = 1;
  mask[1024*311+715] = 1;
  mask[1024*229+212] = 1;
  mask[1024*223+119] = 1;
  mask[1024*217+41] = 1;
  mask[1024*177+159] = 1;
  mask[1024*392+99] = 1;
  mask[1024*403+71] = 1;
  mask[1024*282+133] = 1;
  mask[1024*299+140] = 1;
  mask[1024*242+281] = 1;
  mask[1024*241+301] = 1;
  mask[1024*206+348] = 1;
  mask[1024*181+611] = 1;
  mask[1024*266+750] = 1;
  mask[1024*375+838] = 1;
  mask[1024*376+838] = 1;
  mask[1024*922+545] = 1;
  mask[1024*932+536] = 1;
  mask[1024*909+612] = 1;
  mask[1024*913+871] = 1;
  mask[1024*817+838] = 1;
  mask[1024*380+588] = 1;
  mask[1024*437+499] = 1;
  mask[1024*459+504] = 1;
  mask[1024*181+611] = 1;
  mask[1024*125+735] = 1;
  mask[1024*116+778] = 1;
  mask[1024*49+460] = 1;
  mask[1024*206+348] = 1;
  mask[1024*170+318] = 1;
  mask[1024*140+322] = 1;
  mask[1024*136+281] = 1;
  mask[1024*142+201] = 1;
  mask[1024*180+166] = 1;
  mask[1024*132+135] = 1;
  mask[1024*163+85] = 1;
  mask[1024*82+83] = 1;
  mask[1024*84+36] = 1;
  mask[1024*2+148] = 1;
  mask[1024*843+62] = 1;
  mask[1024*794+42] = 1;
  mask[1024*769+22] = 1;
  mask[1024*747+29] = 1;
  mask[1024*680+28] = 1;
  mask[1024*691+58] = 1;
  mask[1024*550+85] = 1;
  mask[1024*526+114] = 1;
  mask[1024*507+188] = 1;
  mask[1024*833+373] = 1;
  mask[1024*941+992] = 1;
  mask[1024*691+540] = 1;
  mask[1024*679+524] = 1;
  mask[1024*422+687] = 1;
  mask[1024*1023+473] = 1;

  for (int r=728; r<1024; r++)
    mask[1024*r+845] = 1;
  for (int r=775; r<1024; r++)
    mask[1024*r+950] = 1;


  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int xy=0; xy<images[p].size(); xy++) {
      if (mask[xy]) {
        int x = xy%width;
        int y = xy/width;

        float nv = 0;
        int ni = 0;
        if (x>0 && mask[xy-1]==0) {
          nv += images[p][xy-1];
          ni++;
        }

        if (x<width-1 && mask[xy+1]==0) {
          nv += images[p][xy+1];
          ni++;
        }

        if (y>0 && mask[xy-width]==0) {
          nv += images[p][xy-width];
          ni++;
        }

        if (y<height-1 && mask[xy+width]==0) {
          nv += images[p][xy+width];
          ni++;
        }

        if (ni>0)
          images[p][xy] = nv / ni;
      }
    }
  }


  // rotate 90 degrees
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    int d = width;
    for (int y=0; y<d/2; y++) {
      for (int x=0; x<d/2; x++) {

        float vals[4];
        vals[0] = images[p][y*width+x];
        vals[1] = images[p][x*width+(d-y-1)];
        vals[2] = images[p][(d-y-1)*width+(d-x-1)];
        vals[3] = images[p][(d-x-1)*width+y];

        images[p][y*width+x] = vals[1];
        images[p][x*width+(d-y-1)] = vals[2];
        images[p][(d-y-1)*width+(d-x-1)] = vals[3];
        images[p][(d-x-1)*width+y] = vals[0];
      }
    }
  }


  // remove 64 pixel border - do it in place
  int border = 64;
  int nwidth = width - border*2;
  int nheight = height - border*2;
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int r=0; r<nheight; r++) {
      for (int c=0; c<nwidth; c++) {
        images[p][r*nwidth+c] = images[p][(r+border)*width + (c+border)];
      }
    }
  }
  width = nwidth;
  height = nheight;

  // resample to the geometry size
  if (width % geometry.GetDetectorSamplesWidth() ||
      height % geometry.GetDetectorSamplesHeight()) {
    std::cerr<<"cannot evenly downsample projections!"<<std::endl;
    return 1;
  }

  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    int wfactor = width / geometry.GetDetectorSamplesWidth();
    int hfactor = height / geometry.GetDetectorSamplesHeight();

    for (int y=0; y<geometry.GetDetectorSamplesHeight(); y++) {
      for (int x=0; x<geometry.GetDetectorSamplesWidth(); x++) {

        float v=0;
        for (int y2=0; y2<hfactor; y2++) {
          for (int x2=0; x2<wfactor; x2++) {
            v += images[p][(y*hfactor+y2)*width + (x*wfactor+x2)];
          }
        }

        images[p][y*geometry.GetDetectorSamplesWidth() + x] = v / (wfactor*hfactor);
      }
    }
    images[p].resize(geometry.GetDetectorSamplesHeight()*geometry.GetDetectorSamplesWidth());
  }


  // compute median intensities
  vector<float> medians;
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    vector<float> sim;// = images[p];
    for (int y=0; y<geometry.GetDetectorSamplesHeight(); y++) {
      float fy = 2*((y+0.5) / geometry.GetDetectorSamplesHeight() - 0.5);
      for (int x=0; x<geometry.GetDetectorSamplesWidth(); x++) {
        float fx = 2*((x+0.5) / geometry.GetDetectorSamplesWidth() - 0.5);

        if (fx*fx + fy*fy < 1) {
          sim.push_back(images[p][y*geometry.GetDetectorSamplesWidth() + x]);
        }

      }      
    }
    std::sort(sim.begin(), sim.end());
    medians.push_back(sim[sim.size()/2]);

    std::cerr<<p<<", "<<medians.back()<<std::endl;
  }

#if 0
  // normalize images to have same median value
  float avemedian = 0;
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    avemedian += medians[p];
  }  
  avemedian /= geometry.GetNumProjectionAngles();

  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int y=0; y<geometry.GetDetectorSamplesHeight(); y++) {
      for (int x=0; x<geometry.GetDetectorSamplesWidth(); x++) {
        images[p][y*geometry.GetDetectorSamplesWidth() + x] *= avemedian / medians[p];
      }      
    }
  }


  // try to adjust for vertical background gradient
  float etop = 540;
  float ebottom = 480;
  float emid = (etop+ebottom)/2;
  etop /= emid;
  ebottom /= emid;
  
  float atop = 1740;
  float abottom = 950;
  float amid = (atop+abottom)/2;
  atop /= amid;
  abottom /= amid;

  float topscale = etop/atop;
  float bottomscale = ebottom/abottom;


  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int y=0; y<geometry.GetDetectorSamplesHeight(); y++) {
      float fy = (y+0.5) / geometry.GetDetectorSamplesHeight();
      float efrac = ebottom + fy*(etop-ebottom);
      float afrac = abottom + fy*(atop-abottom);
      float scale = efrac / afrac;

      scale = bottomscale + fy*(topscale-bottomscale);


      for (int x=0; x<geometry.GetDetectorSamplesWidth(); x++) {
        images[p][y*geometry.GetDetectorSamplesWidth() + x] *= scale;
      }      
    }
  }



  // try to adjust for vertical background gradient
  float eright = 518;
  float eleft = 517;
  emid = (eright+eleft)/2;
  eright /= emid;
  eleft /= emid;
  
  float aright = 1540;
  float aleft = 1750;
  amid = (aright+aleft)/2;
  aright /= amid;
  aleft /= amid;

  float rightscale = eright/aright;
  float leftscale = eleft/aleft;


  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int x=0; x<geometry.GetDetectorSamplesWidth(); x++) {
      float fx = (x+0.5) / geometry.GetDetectorSamplesHeight();
      float efrac = eleft + fx*(eright-eleft);
      float afrac = aleft + fx*(aright-aleft);
      float scale = efrac / afrac;
      scale = leftscale + fx*(rightscale-leftscale);

      for (int y=0; y<geometry.GetDetectorSamplesHeight(); y++) {
        images[p][y*geometry.GetDetectorSamplesWidth() + x] *= scale;
      }      
    }
  }
#endif  



  // init output volume
  FloatVolumeType::Pointer outVol = FloatVolumeType::New();
  FloatVolumeType::SizeType size;
  size.SetElement(0, geometry.GetDetectorSamplesWidth());
  size.SetElement(1, geometry.GetDetectorSamplesHeight());
  size.SetElement(2, geometry.GetNumProjectionAngles());
  outVol->SetRegions(FloatVolumeType::RegionType(size));
  outVol->Allocate();

  int i=0;
  for (int p=0; p<geometry.GetNumProjectionAngles(); p++) {
    for (int xy=0; xy<images[p].size(); xy++) {
      outVol->GetBufferPointer()[p*geometry.GetDetectorSamplesWidth()*geometry.GetDetectorSamplesHeight() + xy] = images[p][xy];
    }
  }

  // save it
  FloatWriterType::Pointer outWriter = FloatWriterType::New();
  outWriter->SetInput(outVol);
  outWriter->SetFileName(outputName);
  outWriter->Update();

  return 0;
}


int Nrrd2Tiff(int argc, char **argv) {

  std::cerr<<"Converting .nrrd to .tiff series"<<std::endl;

  if (argc < 4) {
    std::cerr<<"usage: backscatter nrrd2tiff <input.nrrd> <output_prefix>"<<std::endl;
    return 1;
  }

  char *inputName = argv[2];
  char *outputPrefix = argv[3];

  // read input volume
  FloatReaderType::Pointer reader = FloatReaderType::New();
  reader->SetFileName(inputName);
  reader->Update();
  FloatVolumeType::Pointer volume = reader->GetOutput();
  FloatVolumeType::SizeType volumeSize = volume->GetLargestPossibleRegion().GetSize();

  int width = volumeSize[0];
  int height = volumeSize[1];
  for (int p=0; p<volumeSize[2]; p++) {
    char outputName[1024];
    sprintf(outputName, "%s%d.tiff", outputPrefix, p);

    TIFF *tif = TIFFOpen(outputName, "w");
    if ( !tif ) {
      std::cerr<<"error writing file!"<<std::endl;
      return 1;
    }

    int rowLength = sizeof(float) * width;
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_IEEEFP);
    char *outPtr = (char*)&volume->GetBufferPointer()[width*height*p];

    uint32 w = width;
    uint32 h = height;

    TIFFSetDirectory(tif, 0);
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(tif, -1));

    int row = 0;
    for (unsigned int idx2 = 0; idx2 < height; idx2++) {

      for (int i=0; i<width; i++) {
        ((float*)outPtr)[i] = idx2;
      }


      if ( TIFFWriteScanline(tif, const_cast<char*>(outPtr), row, 0) < 0) {
        break;
      }
      outPtr += rowLength;
      row++;
    }
    TIFFClose(tif);
  }


  return 0;
}



int Test(int argc, char **argv) {

  vector<double> errors;
  for (int i=0; i<5; i++) {
    errors.push_back((double)rand()*100/RAND_MAX);
  }

  for (int j=0; j<10; j++) {
    double adjustment = (double)rand()*100/RAND_MAX;

    double probsum = 0;
    vector<double> probs;
    for (int i=0; i<(int)errors.size(); i++) {
      probs.push_back(exp(-errors[i] + adjustment));
      probsum += probs.back();
    }

    for (int i=0; i<(int)errors.size(); i++) {
      std::cerr<<(probs[i]/probsum)<<" ";
    }
    std::cerr<<std::endl;

  }

  return 0;
}


int main(int argc, char **argv) {

  if (argc>=2 && _stricmp(argv[1], "markov_forward_project")==0) {
    return RunMarkovForwardProject(argc, argv);
  }

  if (argc>=2 && _stricmp(argv[1], "markov_scale_volume")==0) {
    return RunMarkovScaleVolume(argc, argv);
  }

  if (argc>=2 && _stricmp(argv[1], "markov_scale_projections")==0) {
    return RunMarkovScaleProjections(argc, argv);
  }

  if (argc>=2 && _stricmp(argv[1], "markov_initial_guess")==0) {
    return RunMarkovInitialGuess(argc, argv);
  }

  if (argc>=2 && _stricmp(argv[1], "gibbs")==0) {
    return RunGibbs(argc, argv);
  }

  if (argc>=2 && _stricmp(argv[1], "sart_forward_project")==0) {
    return RunSARTForwardProject(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "sart_reconstruct")==0) {
    return RunSARTReconstruct(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "scale")==0) {
    return ScaleVolume(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "poisson")==0) {
    return PoissonCorrupt(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "gen_synth")==0) {
    return GenerateSynthetic(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "segment_calib_points")==0) {
    return SegmentCalibPoints(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "map_volume")==0) {
    return MapVolume(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "localize_calib_points")==0) {
    return RunLocalizeCalibPoints(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "calibrate_detector")==0) {
    return RunCalibrateDetector(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "calibrate_intensity")==0) {
    return RunCalibrateIntensity(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "projections_to_nrrd")==0) {
    return ProjectionsToNrrd(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "nrrd2tiff")==0) {
    return Nrrd2Tiff(argc, argv);
  }

  else if (argc>=2 && _stricmp(argv[1], "test")==0) {
    return Test(argc, argv);
  }


  else {
    std::cerr<<"usage: backscatter markov_forward_project ..."<<std::endl;
    std::cerr<<"usage: backscatter markov_scale_volume ..."<<std::endl;
    std::cerr<<"usage: backscatter markov_scale_projections ..."<<std::endl;
    std::cerr<<"usage: backscatter markov_initial_guess ..."<<std::endl;
    std::cerr<<"usage: backscatter gibbs ..."<<std::endl;
    std::cerr<<"usage: backscatter sart_forward_project ..."<<std::endl;
    std::cerr<<"usage: backscatter sart_reconstruct ..."<<std::endl;
    std::cerr<<"usage: backscatter scale ..."<<std::endl;
    std::cerr<<"usage: backscatter poisson ..."<<std::endl;
    std::cerr<<"usage: backscatter gen_synth ..."<<std::endl;
    std::cerr<<"usage: backscatter segment_calib_points ..."<<std::endl;
    std::cerr<<"usage: backscatter map_volume ..."<<std::endl;
    std::cerr<<"usage: backscatter localize_calib_points ..."<<std::endl;
    std::cerr<<"usage: backscatter calibrate_detector ..."<<std::endl;
    std::cerr<<"usage: backscatter calibrate_intensity ..."<<std::endl;
    std::cerr<<"usage: backscatter projections_to_nrrd ..."<<std::endl;
    std::cerr<<"usage: backscatter nrrd2tiff ..."<<std::endl;
    return 1;
  }
}
