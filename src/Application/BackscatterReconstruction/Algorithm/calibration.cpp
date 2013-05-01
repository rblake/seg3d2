
#include <map>
#include "sart.h"
#include "calibration.h"
#include "lmmin.h"
#include "gradient.h"

//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_linalg.h>

void SegmentCalibPoints(int detWidth, int detHeight, int detAngles,
                        const vector<float> &volume,
                        vector<unsigned char> &ids,
                        int numPoints,
                        int maskSize,
                        int idSize,
                        int border,
                        bool dark) {

  int totalSamples = detWidth * detHeight * detAngles;

  vector<unsigned char> mask(volume.size(), 0);
  ids.resize(volume.size());
  for (size_t i=0; i<ids.size(); i++)
    ids[i] = 0;

  // mask the points
  for (int p=0; p<detAngles; p++) {

    for (int point=0; point<numPoints; point++) {

      // find darkest or lightest pixel
      int px=-1;
      int py=-1;
      float pv=dark?1e9:-1e9;
      for (int y=border; y<detHeight-border; y++) {
        for (int x=border; x<detWidth-border; x++) {
          int i = p*detWidth*detHeight + y*detWidth + x;
          if (mask[i] == 0) {
            if ((px<0 || py<0) || 
                (dark && (volume[i] < pv)) ||
                (!dark && (volume[i] > pv))) {
              px = x;
              py = y;
              pv = volume[i];
            }
          }
        }
      }


      // mask the area around it out
      for (int y2=std::max(0,py-idSize); y2<=std::min(detHeight-1, py+idSize); y2++) {
        for (int x2=std::max(0,px-idSize); x2<=std::min(detWidth-1, px+idSize); x2++) {
          if ((px-x2)*(px-x2) + (py-y2)*(py-y2) < idSize*idSize) {
            ids[p*detWidth*detHeight + y2*detWidth + x2] = 255;
          }
        }
      }

      // mask the area around it out
      for (int y2=std::max(0,py-maskSize); y2<=std::min(detHeight-1, py+maskSize); y2++) {
        for (int x2=std::max(0,px-maskSize); x2<=std::min(detWidth-1, px+maskSize); x2++) {
          if ((px-x2)*(px-x2) + (py-y2)*(py-y2) < maskSize*maskSize) {
            mask[p*detWidth*detHeight + y2*detWidth + x2] = 255;
          }
        }
      }
    }

  }


  // flood all of the points
  int np = 1;
  for (int start_i=0; start_i<totalSamples; start_i++) {

    if (ids[start_i] != 255)
      continue;

    vector<int> stack;
    ids[start_i] = np;
    stack.push_back(start_i);

    while (!stack.empty()) {

      int next_i = stack.back();
      stack.pop_back();

      int p, x, y;
      p = next_i / (detWidth*detHeight);
      int temp = next_i % (detWidth*detHeight);
      y = temp / detWidth;
      x = temp % detWidth;

      if (p<detAngles-1) {
        int i = (p+1)*detWidth*detHeight + y*detWidth + x;
        if (ids[i] == 255) {
          ids[i] = np;
          stack.push_back(i);
        }
      }

      if (x>0) {
        int i = p*detWidth*detHeight + y*detWidth + (x-1);
        if (ids[i] == 255) {
          ids[i] = np;
          stack.push_back(i);
        }
      }
      if (x<detWidth-1) {
        int i = p*detWidth*detHeight + y*detWidth + (x+1);
        if (ids[i] == 255) {
          ids[i] = np;
          stack.push_back(i);
        }
      }

      if (y>0) {
        int i = p*detWidth*detHeight + (y-1)*detWidth + x;
        if (ids[i] == 255) {
          ids[i] = np;
          stack.push_back(i);
        }
      }
      if (y<detHeight-1) {
        int i = p*detWidth*detHeight + (y+1)*detWidth + x;
        if (ids[i] == 255) {
          ids[i] = np;
          stack.push_back(i);
        }
      }
    }

    np++;
  }


}




void FitGaussian1D(const vector<int> &xs, vector<float> &ys,
                   float &mean, float &std, float &amp, float &back) {
  mean = std = amp = back = 0;
 
  // Initial guess for amplitude, center, std deviation, and offset
  double params[4] = {1.0, 0.0, xs.size()*0.5, 0.0 };

  for (int i = 0; i < xs.size(); i++) {
    //std::cerr<<xs[i]<<" "<<ys[i]<<std::endl;
  }

  // Find min/max 
  double min = ys[0], max = ys[0];
  int minyi = 0;
  int maxyi = 0;
  for (int i = 0; i < xs.size(); i++) {
    if (ys[i] < min) {
      min = ys[i];
      minyi = i;
    }
    if (ys[i] > max) {
      max = ys[i];
      maxyi = i;
    }
  }

  int xi0, xi1;
  if (1 || -min > max) { // flipped
    params[0] = min-max;
    params[1] = xs[minyi];
    params[3] = max;
    xi0 = std::max(0, minyi-3);
    xi1 = std::min((int)xs.size()-1, minyi+3);
  }
  else {
    params[0] = max-min;
    params[1] = xs[maxyi];
    params[3] = min;
    xi0 = std::max(0, maxyi-3);
    xi1 = std::min((int)xs.size()-1, maxyi+3);
  }


  // copy out data we're actually going to use
  double *X = new double[xs.size()];
  double *Y = new double[ys.size()];
  for (int i=xi0; i<=xi1; i++) {
    X[i-xi0] = (double)xs[i];
    Y[i-xi0] = (double)ys[i];
  }


  // Setup minimization
  lm_control_type control;
  lm_data_1d data;
  lm_initialize_control(&control);
  data.f = lm_gaussian_fit_1d;
  data.xvec = X;
  data.fvec = Y;


  mean = params[1];
  std = params[2];
  amp = params[0];
  back = params[3];
  std::cerr<<mean<<" "<<std<<" "<<amp<<" "<<back<<std::endl;
  
  // Perform least squares minimization
  lm_minimize(xi1-xi0+1, 4, params, lm_evaluation_1d, lm_print_1d, &data, &control);
  std::cerr << "status: " << lm_shortmsg[control.info] << " after " << control.nfev << " evaluations" << std::endl;

  delete [] X;
  delete [] Y;

  mean = params[1];
  std = params[2];
  amp = params[0];
  back = params[3];
  std::cerr<<mean<<" "<<std<<" "<<amp<<" "<<back<<std::endl;
}



float FitGaussian2D(const vector<float> &xs, vector<float> &ys, vector<float> &fs,
                    float &meanx, float &meany,
                    float &stdx, float &stdy, float &cross,
                    float &amp, float &back,
                    bool dark) {
  meanx = meany = stdx = stdy = amp = back = 0;
 
  // Initial guess for amplitude, center, std deviation, and offset
  // 0: amp
  // 1: meanx
  // 2: stdx
  // 3: meany
  // 4: stdy
  // 5: back
  double params[12] = { 1.0, 
                        0.0, 4.5,
                        0.0, 4.5,
                        0.0,
                        0.0,
                        0.0, 0.0,
                        0.0, 0.0, 0.0 };

  for (int i = 0; i < xs.size(); i++) {
    //std::cerr<<xs[i]<<" "<<ys[i]<<" "<<fs[i]<<std::endl;
  }

  // Find min/max 
  double min = fs[0], max = fs[0];
  int mini = 0;
  int maxi = 0;
  for (int i = 0; i < xs.size(); i++) {
    if (fs[i] < min) {
      min = fs[i];
      mini = i;
    }
    if (fs[i] > max) {
      max = fs[i];
      maxi = i;
    }
  }

  if (dark) { // flipped
    params[0] = min-max;
    params[1] = xs[mini];
    params[3] = ys[mini];
    params[5] = max;
  }
  else {
    params[0] = max-min;
    params[1] = xs[maxi];
    params[3] = ys[maxi];
    params[5] = min;
  }


  // copy out data we're actually going to use
  double *X = new double[xs.size()];
  double *Y = new double[ys.size()];
  double *F = new double[ys.size()];
  for (int i=0; i<xs.size(); i++) {
    X[i] = (double)xs[i];
    Y[i] = (double)ys[i];
    F[i] = (double)fs[i];
  }

  // Setup minimization
  lm_control_type control;
  lm_initialize_control(&control);
  control.epsilon = 1e-5;
  control.stepbound = 10;
  lm_data_2d data;
  data.f = lm_gaussian_fit_2dg;
  data.xvec = X;
  data.yvec = Y;
  data.fvec = F;


  meanx = params[1];
  stdx = params[2];
  meany = params[3];
  stdy = params[4];
  amp = params[0];
  back = params[5];
  cross = params[6];
  double dx = params[7];
  double dy = params[8];
  double dxx = params[9];
  double dxy = params[10];
  double dyy = params[11];
  std::cerr<<meanx<<" "<<meany<<" "<<stdx<<" "<<stdy<<" "<<cross<<" "<<amp<<" "<<back<<" "<<std::endl;
  std::cerr<<dx<<" "<<dy<<std::endl;
  std::cerr<<dxx<<" "<<dxy<<" "<<dyy<<std::endl;
  
  // Perform least squares minimization
  lm_minimize(fs.size(), 12, params, lm_evaluation_2d, lm_print_2d, &data, &control);
  std::cerr << "status: " << lm_shortmsg[control.info] << " after " << control.nfev << " evaluations" << std::endl;
  std::cerr << "norm: " <<control.fnorm << std::endl;

  delete [] X;
  delete [] Y;
  delete [] F;

  meanx = params[1];
  stdx = params[2];
  meany = params[3];
  stdy = params[4];
  amp = params[0];
  back = params[5];
  cross = params[6];
  dx = params[7];
  dy = params[8];
  dxx = params[9];
  dxy = params[10];
  dyy = params[11];
  std::cerr<<meanx<<" "<<meany<<" "<<stdx<<" "<<stdy<<" "<<cross<<" "<<amp<<" "<<back<<" "<<std::endl;
  std::cerr<<dx<<" "<<dy<<std::endl;
  std::cerr<<dxx<<" "<<dxy<<" "<<dyy<<std::endl;

  return control.fnorm;
}


void LocalizeCalibPoints(int detWidth, int detHeight, int detAngles,
                         const vector<float> &volume,
                         const vector<unsigned char> &mask,
                         vector<int> &proj_i,
                         vector<int> &point_i,
                         vector<float> &point_x,
                         vector<float> &point_y,
                         vector<float> &residuals,
                         bool dark) {

  for (int p=0; p<detAngles; p++) {
    const float *volp = &volume[p * detWidth * detHeight];
    const unsigned char *maskp = &mask[p * detWidth * detHeight];

    // find the max point id in this projection
    unsigned char maxmask = 0;
    for (int i=0; i<detWidth * detHeight; i++)
      maxmask = std::max(maxmask, maskp[i]);
    std::cerr<<(int)maxmask<<std::endl;

    // for each mask id/point - id 0 is background
    for (int m=1; m<=maxmask; m++) {

      std::cerr<<"projection "<<p<<std::endl;
      std::cerr<<"mask "<<m<<std::endl;

      vector<float> xs, ys, fs;

      // find the center of the mask for this point
      float centerx=0;
      float centery=0;
      int n = 0;
      int minx = 99999;
      int miny = 99999;
      int maxx = 0;
      int maxy = 0;

      for (int y=0; y<detHeight; y++) {
        for (int x=0; x<detWidth; x++) {

          if (maskp[y*detWidth + x] == m) {
            centerx += x;
            centery += y;
            minx = std::min(minx, x);
            miny = std::min(miny, y);
            maxx = std::max(maxx, x);
            maxy = std::max(maxy, y);

            xs.push_back(x);
            ys.push_back(y);
            fs.push_back(volp[y*detWidth + x]);
            
            n++;
          }
        }
      }

      // this point doesn't appear in this projecion
      if (n<5)
        continue;
      centerx /= n;
      centery /= n;


      // cutout a box around the particle - directly project in x and y
      const int cutout_size = (int)(2.2*sqrt(n/M_PI));
      std::cerr<<"cutout_size: "<<cutout_size<<std::endl;
      int cminx = std::max(0, (int)(centerx - cutout_size/2));
      int cmaxx = std::min(detWidth, (int)(centerx + cutout_size/2)); // exclusive
      int cminy = std::max(0, (int)(centery - cutout_size/2));
      int cmaxy = std::min(detHeight, (int)(centery + cutout_size/2)); // exclusive
      std::cerr<<"cutout: "<<cminx<<" "<<cmaxx<<" "<<cminy<<" "<<cmaxy<<std::endl;

      vector<int> xi(cmaxx-cminx,0), yi(cmaxy-cminy,0);
      vector<float> xv(cmaxx-cminx,0), yv(cmaxy-cminy,0);

      for (int y=cminy; y<cmaxy; y++) {
        yi[y-cminy] = y;
        for (int x=cminx; x<cmaxx; x++) {
          xi[x-cminx] = x;

          float v = volp[y*detWidth + x];
          yv[y-cminy] += v;
          xv[x-cminx] += v;
        }
      }


      // fit a gaussian to both x and y
      float meanx, meany;
      float stdx, stdy, cross;
      float amp, back;

      float res = FitGaussian2D(xs, ys, fs,
                                meanx, meany, 
                                stdx, stdy, cross,
                                amp, back,
                                dark);

      proj_i.push_back(p);
      point_i.push_back(m);
      point_x.push_back((meanx+0.5) / detWidth); // measurements are normalized [0,1]
      point_y.push_back((meany+0.5) / detHeight);
      residuals.push_back(res);
    }
  }
}


void OrthoNormalize(double M[3][3]) {
  Vec3d v1(M[0][0], M[1][0], M[2][0]);
  Vec3d v2(M[0][1], M[1][1], M[2][1]);
  Vec3d v3(M[0][2], M[1][2], M[2][2]);
  v1.Normalize();

  // remove v1 component from v2 and v3
  v2 = v2 - v1.Dot(v2)*v1;
  v3 = v3 - v1.Dot(v3)*v1;

  // remove v2 component from v3
  v2.Normalize();
  v3 = v3 - v2.Dot(v3)*v2;

  v3.Normalize();

  M[0][0] = v1[0];  M[1][0] = v1[1];  M[2][0] = v1[2];
  M[0][1] = v2[0];  M[1][1] = v2[1];  M[2][1] = v2[2];
  M[0][2] = v3[0];  M[1][2] = v3[1];  M[2][2] = v3[2];
}

void OrthoNormalizeCross(double M[3][4]) {
  Vec3d v1(M[0][0], M[1][0], M[2][0]);
  Vec3d v2(M[0][1], M[1][1], M[2][1]);
  Vec3d v3(M[0][2], M[1][2], M[2][2]);
  v1.Normalize();

  // remove v1 component from v2
  v2 = v2 - v1.Dot(v2)*v1;
  v2.Normalize();

  // v3 is cross product of v1 and v2
  v3 = v1.Cross(v2);
  v3.Normalize();

  M[0][0] = v1[0];  M[1][0] = v1[1];  M[2][0] = v1[2];
  M[0][1] = v2[0];  M[1][1] = v2[1];  M[2][1] = v2[2];
  M[0][2] = v3[0];  M[1][2] = v3[1];  M[2][2] = v3[2];
}


void OrthoNormalizeCross(double M[3][3]) {
  Vec3d v1(M[0][0], M[1][0], M[2][0]);
  Vec3d v2(M[0][1], M[1][1], M[2][1]);
  Vec3d v3(M[0][2], M[1][2], M[2][2]);
  v1.Normalize();

  // remove v1 component from v2
  v2 = v2 - v1.Dot(v2)*v1;
  v2.Normalize();

  // v3 is cross product of v1 and v2
  v3 = v1.Cross(v2);
  v3.Normalize();

  M[0][0] = v1[0];  M[1][0] = v1[1];  M[2][0] = v1[2];
  M[0][1] = v2[0];  M[1][1] = v2[1];  M[2][1] = v2[2];
  M[0][2] = v3[0];  M[1][2] = v3[1];  M[2][2] = v3[2];
}



template <typename T>
void MatrixToQuaternion(const T mat[3][3], T quat[4]) {

  int u,v,w;
  if (mat[0][0]>=mat[1][1] && mat[0][0]>=mat[2][2]) {
    u = 0;
    v = 1;
    w = 2;
  }

  else if (mat[1][1]>=mat[0][0] && mat[1][1]>=mat[2][2]) {
    u = 1;
    v = 2;
    w = 0;
  }

  else {
    u = 2;
    v = 0;
    w = 1;
  }

  T r = sqrt(1 + mat[u][u] - mat[v][v] - mat[w][w]);
  quat[0] = (mat[w][v] - mat[v][w]) / (2*r);
  quat[u+1] = r/2;
  quat[v+1] = (mat[u][v] + mat[v][u]) / (2*r);
  quat[w+1] = (mat[w][u] + mat[u][w]) / (2*r);

}


template <typename T>
void QuaternionToMatrix(const T quat[4], T mat[3][3]) {

  const T &a = quat[0];
  const T &b = quat[1];
  const T &c = quat[2];
  const T &d = quat[3];

  mat[0][0] = a*a + b*b - c*c - d*d;
  mat[0][1] = 2*b*c - 2*a*d;
  mat[0][2] = 2*b*d + 2*a*c;

  mat[1][0] = 2*b*c + 2*a*d;
  mat[1][1] = a*a - b*b + c*c - d*d;
  mat[1][2] = 2*c*d - 2*a*b;

  mat[2][0] = 2*b*d - 2*a*c;
  mat[2][1] = 2*c*d + 2*a*b;
  mat[2][2] = a*a - b*b - c*c + d*d;
}


template <typename T>
T EvalResidual(Geometry &geometry,
               const T detQuat[4], const Vec3<T> &_detPos,
               const T baseQuat[4], const Vec3<T> &_baseOffset,
               const vector<int> &proj_i,
               const vector<int> &point_i,
               const vector<float> &point_mx,
               const vector<float> &point_my,
               const vector<Vec3f> &cpoints,
               float detectorDistance,
               T *residuals=NULL) {

  T residual = 0;

  T detMat[3][3];
  QuaternionToMatrix(detQuat, detMat);

  T baseMat[3][3];
  QuaternionToMatrix(baseQuat, baseMat);


  // get right and up vectors from the quaternion
  Vec3<T> right(detMat[0][0], detMat[1][0], detMat[2][0]);
  Vec3<T> up(detMat[0][1], detMat[1][1], detMat[2][1]);
  Vec3<T> forward(-detMat[0][2], -detMat[1][2], -detMat[2][2]);


  // keep detector at fixed distance from origin
  Vec3<T> detPos = _detPos;
  Vec3<T> baseOffset = _baseOffset;
  detPos[2] -= baseOffset[2];
  baseOffset[2] = 0;
  if (1 || fabs(geometry.GetColimatorFocalLength()) > 10) {
    // find closest point to origin along forward direction
    Vec3<T> p0 = detPos - forward.Dot(detPos)*forward;
    T osq = p0.LengthSquared();
    T discrim = detectorDistance*detectorDistance - osq;
    if (discrim < 0) {
      std::cerr<<"calibration error - original detector distance not attainable!"<<std::endl;
      std::cerr<<_baseOffset[2]<<std::endl;
      std::cerr<<_detPos[2]<<std::endl;
    }
    else {
      T y = sqrt(discrim);
      detPos = p0 - y*forward;
    }
  }


  for (size_t ij=0; ij<proj_i.size(); ij++) {
    int i = proj_i[ij];
    int j = point_i[ij];

    // bead location in it's own coordinate space
    Vec3f c_j  = cpoints[j];

    // account for the rotation / translation of the calibration grid
    Vec3<T> tc_j(c_j[0]*baseMat[0][0] + c_j[1]*baseMat[1][0] + c_j[2]*baseMat[2][0] + baseOffset[0],
                 c_j[0]*baseMat[0][1] + c_j[1]*baseMat[1][1] + c_j[2]*baseMat[2][1] + baseOffset[1],
                 c_j[0]*baseMat[0][2] + c_j[1]*baseMat[1][2] + c_j[2]*baseMat[2][2] + baseOffset[2]);

    // rotate the bead into the current projection angle
    Vec3<T> tc = geometry.GetInverseRotatedPoint(tc_j, i);

    // physical position of this projection on the detector
    Vec3<T> tp(point_mx[ij], point_my[ij], 0);

    // relative position of the transformed bead to the detector origin
    Vec3<T> tcrel = tc - detPos;

    // l2 distance of transformed bead and projection localization, after projection of the bead
    T focalScale = geometry.GetColimatorFocalLength() / (forward.Dot(tcrel) + geometry.GetColimatorFocalLength());
    T rightComp = focalScale * right.Dot(tcrel) - tp[0];
    T upComp = focalScale * up.Dot(tcrel) - tp[1];

    residual += rightComp*rightComp + upComp*upComp;
    if (residuals) {
      residuals[ij*2+0] = rightComp;
      residuals[ij*2+1] = upComp;
    }
        

    //std::cerr<<i<<" "<<j<<" "<<sqrt((rightComp*rightComp + upComp*upComp).value)<<std::endl;
  }
  
  return residual;
}


typedef struct {
  Geometry *geometry;
  const vector<int> *proj_i;
  const vector<int> *point_i;
  const vector<float> *point_mx;
  const vector<float> *point_my;
  const vector<Vec3f> *cpoints;
  float detectorDistance;
} lm_calib_data;

void lm_calib_eval(double *par, int m_dat, double *fvec, void *data, int *info) {
  lm_calib_data *mydata = (lm_calib_data*)data;

  double tmp;

  double detQuat[4] = { par[0], par[1], par[2], par[3] };
  tmp = sqrt(detQuat[0]*detQuat[0] + detQuat[1]*detQuat[1] + detQuat[2]*detQuat[2] + detQuat[3]*detQuat[3]);
  detQuat[0]/=tmp; detQuat[1]/=tmp; detQuat[2]/=tmp; detQuat[3]/=tmp;

  Vec3<double> detPos(par[4], par[5], par[6]);

  double baseQuat[4] = { par[7], par[8], par[9], par[10] };
  tmp = sqrt(baseQuat[0]*baseQuat[0] + baseQuat[1]*baseQuat[1] + baseQuat[2]*baseQuat[2] + baseQuat[3]*baseQuat[3]);
  baseQuat[0]/=tmp; baseQuat[1]/=tmp; baseQuat[2]/=tmp; baseQuat[3]/=tmp;
  
  Vec3<double> basePos(par[11], par[12], par[13]);


  EvalResidual(*mydata->geometry,
               detQuat, detPos,
               baseQuat, basePos,
               *mydata->proj_i,
               *mydata->point_i,
               *mydata->point_mx,
               *mydata->point_my,
               *mydata->cpoints,
               mydata->detectorDistance,
               fvec);
}




void CalibrateDetector(Geometry &geometry, 
                       const vector<int> &proj_i,
                       const vector<int> &point_i,
                       const vector<float> &point_x,
                       const vector<float> &point_y,
                       const vector<Vec3f> &cpoints,
                       Vec3f &baseCenter,
                       Vec3f &baseX,
                       Vec3f &baseY,
                       Vec3f &baseZ) {

  std::cerr<<"before: "<<std::endl;
  std::cerr<<geometry.GetDetectorCenter()[0]<<" "
           <<geometry.GetDetectorCenter()[1]<<" "
           <<geometry.GetDetectorCenter()[2]<<std::endl;
  std::cerr<<geometry.GetDetectorForward()[0]<<" "
           <<geometry.GetDetectorForward()[1]<<" "
           <<geometry.GetDetectorForward()[2]<<std::endl;
  std::cerr<<geometry.GetDetectorRight()[0]<<" "
           <<geometry.GetDetectorRight()[1]<<" "
           <<geometry.GetDetectorRight()[2]<<std::endl;
  std::cerr<<geometry.GetDetectorUp()[0]<<" "
           <<geometry.GetDetectorUp()[1]<<" "
           <<geometry.GetDetectorUp()[2]<<std::endl;


  // first convert point positions into physical distances relative to the center of the detector
  vector<float> point_mx, point_my;
  for (size_t i=0; i<point_x.size(); i++) {
    float mx, my;
    geometry.ProjectionPixelToPhysicalOffset(point_x[i] * geometry.GetDetectorSamplesWidth(),
                                             point_y[i] * geometry.GetDetectorSamplesHeight(),
                                             mx, my);
    point_mx.push_back(mx);
    point_my.push_back(my);

    //std::cerr<<mx<<" "<<my<<std::endl;
  }


  Vec3f forwardOrig = geometry.GetDetectorForward();
  Vec3f rightOrig = geometry.GetDetectorRight();
  Vec3f upOrig = geometry.GetDetectorUp();
  Vec3f posOrig = geometry.GetDetectorCenter();
  float matOrig[3][3] = { { rightOrig[0], upOrig[0], -forwardOrig[0] },
                          { rightOrig[1], upOrig[1], -forwardOrig[1] },
                          { rightOrig[2], upOrig[2], -forwardOrig[2] } };
  float origDetDist = posOrig.Length();


  // get our current detector orientation as a quaternion
  Vec3f forward = geometry.GetDetectorForward();
  Vec3f right = geometry.GetDetectorRight();
  Vec3f up = geometry.GetDetectorUp();
  Vec3f detPos = geometry.GetDetectorCenter();
  double detMat[3][3] = { { right[0], up[0], -forward[0] },
                          { right[1], up[1], -forward[1] },
                          { right[2], up[2], -forward[2] } };

  double detQuat[4];
  MatrixToQuaternion(detMat, detQuat);
  QuaternionToMatrix(detQuat, detMat);

  double baseQuat[4] = { 1,0,0,0 };
  Vec3f baseOffset(0,0,0.01*(float)rand()/RAND_MAX);


  lm_calib_data lmdata;
  lmdata.geometry = &geometry;
  lmdata.proj_i = &proj_i;
  lmdata.point_i = &point_i;
  lmdata.point_mx = &point_mx;
  lmdata.point_my = &point_my;
  lmdata.cpoints = &cpoints;
  lmdata.detectorDistance = origDetDist;

  lm_control_type control;
  lm_initialize_control(&control);
  //control.epsilon = 1e-5;

  double params[14];
  params[0] = detQuat[0];
  params[1] = detQuat[1];
  params[2] = detQuat[2];
  params[3] = detQuat[3];
  params[4] = detPos[0];
  params[5] = detPos[1];
  params[6] = detPos[2];
  params[7] = baseQuat[0];
  params[8] = baseQuat[1];
  params[9] = baseQuat[2];
  params[10] = baseQuat[3];
  params[11] = baseOffset[0];
  params[12] = baseOffset[1];
  params[13] = baseOffset[2];

  lm_minimize(proj_i.size()*2, 14, params, lm_calib_eval, NULL, &lmdata, &control);

  std::cerr << "status: " << lm_shortmsg[control.info] << " after " << control.nfev << " evaluations" << std::endl;
  std::cerr << "norm: " <<control.fnorm << std::endl;
  std::cerr << "sqnorm: " <<control.fnorm*control.fnorm << std::endl;
  std::cerr << "rms: " <<sqrt(control.fnorm*control.fnorm / (proj_i.size()*2)) << std::endl;
  

  detQuat[0] = params[0];
  detQuat[1] = params[1];
  detQuat[2] = params[2];
  detQuat[3] = params[3];
  detPos[0] = params[4];
  detPos[1] = params[5];
  detPos[2] = params[6];
  baseQuat[0] = params[7];
  baseQuat[1] = params[8];
  baseQuat[2] = params[9];
  baseQuat[3] = params[10];
  baseOffset[0] = params[11];
  baseOffset[1] = params[12];
  baseOffset[2] = params[13];

  double qDetLen = sqrt(detQuat[0]*detQuat[0] + detQuat[1]*detQuat[1] + detQuat[2]*detQuat[2] + detQuat[3]*detQuat[3]);
  detQuat[0] /= qDetLen;
  detQuat[1] /= qDetLen;
  detQuat[2] /= qDetLen;
  detQuat[3] /= qDetLen;

  double qBaseLen = sqrt(baseQuat[0]*baseQuat[0] + baseQuat[1]*baseQuat[1] + baseQuat[2]*baseQuat[2] + baseQuat[3]*baseQuat[3]);
  baseQuat[0] /= qBaseLen;
  baseQuat[1] /= qBaseLen;
  baseQuat[2] /= qBaseLen;
  baseQuat[3] /= qBaseLen;

  
  // convert the quat back into the detector geometry
  QuaternionToMatrix(detQuat, detMat);

  Vec3f forwardNew(-detMat[0][2], -detMat[1][2], -detMat[2][2]);
  Vec3f rightNew(detMat[0][0], detMat[1][0], detMat[2][0]);
  Vec3f upNew(detMat[0][1], detMat[1][1], detMat[2][1]);

  // take any z offset from the base transform and add it to the detector transform
  // sets the origin to be at the height of the calibration pattern
  detPos[2] -= baseOffset[2];
  baseOffset[2] = 0;

  // if close to ortho projection,
  // move the detector back to the same distance from the origin 
  // (this must be measured and setup correctly as the initial guess)
  // do this before rotation!
  if (1 || fabs(geometry.GetColimatorFocalLength()) > 10) {
    // find closest point to origin along forward direction
    Vec3f p0 = detPos - forwardNew.Dot(detPos)*forwardNew;
    float osq = p0.LengthSquared();
    float discrim = origDetDist*origDetDist - osq;
    if (discrim < 0) {
      std::cerr<<"calibration error - original detector distance not attainable!"<<std::endl;
    }
    else {
      float y = sqrt(discrim);
      detPos = p0 - y*forwardNew;
    }
  }

  // rotate everything so the y component of the detector position is 0
  double theta = atan2(detPos[1], detPos[0]);// + M_PI;
  detPos = Geometry::RotateVector(detPos, -theta);
  forwardNew = Geometry::RotateVector(forwardNew, -theta);
  rightNew = Geometry::RotateVector(rightNew, -theta);
  upNew = Geometry::RotateVector(upNew, -theta);


  geometry.SetDetectorForward(forwardNew);
  geometry.SetDetectorRight(rightNew);
  geometry.SetDetectorUp(upNew);
  geometry.SetDetectorCenter(detPos);

  std::cerr<<"after: "<<std::endl;
  std::cerr<<"center: "
           <<geometry.GetDetectorCenter()[0]<<" "
           <<geometry.GetDetectorCenter()[1]<<" "
           <<geometry.GetDetectorCenter()[2]<<std::endl;
  std::cerr<<"forward: "
           <<geometry.GetDetectorForward()[0]<<" "
           <<geometry.GetDetectorForward()[1]<<" "
           <<geometry.GetDetectorForward()[2]<<std::endl;
  std::cerr<<"right: "
           <<geometry.GetDetectorRight()[0]<<" "
           <<geometry.GetDetectorRight()[1]<<" "
           <<geometry.GetDetectorRight()[2]<<std::endl;
  std::cerr<<"up: "
           <<geometry.GetDetectorUp()[0]<<" "
           <<geometry.GetDetectorUp()[1]<<" "
           <<geometry.GetDetectorUp()[2]<<std::endl;


  // print the base rotation for debugging
  // convert the quat back into the detector geometry
  double baseMat[3][3];
  QuaternionToMatrix(baseQuat, baseMat);

  Vec3f bX(baseMat[0][0], baseMat[1][0], baseMat[2][0]);
  Vec3f bY(baseMat[0][1], baseMat[1][1], baseMat[2][1]);
  Vec3f bZ(baseMat[0][2], baseMat[1][2], baseMat[2][2]);

  bX = Geometry::RotateVector(bX, theta);
  bY = Geometry::RotateVector(bY, theta);
  bZ = Geometry::RotateVector(bZ, theta);
  baseOffset = Geometry::RotateVector(baseOffset, theta);

  std::cerr<<"theta: "<<theta<<std::endl;

  // return the base transform so we can synthesize a volume for intensity calibration
  baseCenter = baseOffset;
  baseX = bX;
  baseY = bY;
  baseZ = bZ;
}


void CalibrateIntensity(Geometry &geometry, 
                        const float *images,
                        const float *forward) {

  double A[2][2] = { { 0,0 }, { 0,0 } };
  double b[2] = { 0,0 };

  for (int i=0; i<geometry.GetTotalProjectionSamples(); i++) {
    A[0][0] += 2*forward[i]*forward[i];
    A[0][1] += 2*forward[i];
    b[0] += 2*forward[i]*images[i];

    A[1][0] += 2*forward[i];
    A[1][1] += 2;
    b[1] += 2*images[i];
  }

  double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
  double AInv[2][2];
  AInv[0][0] =  A[1][1] / det;
  AInv[0][1] = -A[0][1] / det;
  AInv[1][0] = -A[1][0] / det;
  AInv[1][1] =  A[0][0] / det;

  geometry.SetDetectorGain(AInv[0][0]*b[0] + AInv[0][1]*b[1]);
  geometry.SetDetectorOffset(AInv[1][0]*b[0] + AInv[1][1]*b[1]);
}
