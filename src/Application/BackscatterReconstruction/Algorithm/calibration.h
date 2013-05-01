
void SegmentCalibPoints(int detWidth, int detHeight, int detAngles,
                        const vector<float> &volume,
                        vector<unsigned char> &ids,
                        int numPoints,
                        int maskSize,
                        int idSize,
                        int border,
                        bool dark);

void LocalizeCalibPoints(int detWidth, int detHeight, int detAngles,
                         const vector<float> &volume,
                         const vector<unsigned char> &mask,
                         vector<int> &proj_i,
                         vector<int> &point_i,
                         vector<float> &point_x,
                         vector<float> &point_y,
                         vector<float> &residuals,
                         bool dark);

void CalibrateDetector(Geometry &geometry, 
                       const vector<int> &proj_i,
                       const vector<int> &point_i,
                       const vector<float> &point_x,
                       const vector<float> &point_y,
                       const vector<Vec3f> &cpoints,
                       Vec3f &baseCenter,
                       Vec3f &baseX,
                       Vec3f &baseY,
                       Vec3f &baseZ);

void CalibrateIntensity(Geometry &geometry, 
                        const float *images,
                        const float *forward);
