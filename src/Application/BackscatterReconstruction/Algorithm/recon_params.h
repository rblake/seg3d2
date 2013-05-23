
// the number of steps to take in voxel coordinates when finding material interfaces
#define RECON_VOXEL_STEP_SIZE  0.33

// how much of the regularization term to include in the final forward projection error
#define RECON_REGULARIZATION_WEIGHT 2e-5

// start and end temperature for annealing in Gibbs optimization
#define RECON_TEMP_START 8.4e-12
#define RECON_TEMP_END 8.4e-13

// high and low xray energies (KeV) - scattering factors are averaged across the range
#define XRAY_ENERGY_LOW  10
#define XRAY_ENERGY_HIGH 10

// the positions of the disks in the calibration pattern
#define RECON_CALIB_PATTERN_DEF(name)             \
  vector<Vec3f> name;                             \
  name.push_back(Vec3f( 0.03175,-0.03175, 0));    \
  name.push_back(Vec3f(-0.03175, 0.03175, 0));    \
  name.push_back(Vec3f(-0.03175,-0.03175, 0));
