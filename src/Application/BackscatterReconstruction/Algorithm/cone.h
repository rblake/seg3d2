
#ifndef CONE_H
#define CONE_H

#include <iostream>
#include <Application/BackscatterReconstruction/Algorithm/vec3.h>

#define COSTHETA_EPS 0.01f
#define MINDIST_EPS 0.001f

class Cone {
public:
  Cone()  {  }


  Cone(const Vec3f &origin, const Vec3f &dir, const float cosTheta)
    : mOrigin(origin),
      mDir(dir),
      mCosTheta(cosTheta),
      mMinDist(0) { 

    if (mCosTheta < 0) {
      std::cerr<<"cone with greater than 90 degrees!"<<std::endl;
      exit(0);
    }
  }


  // initialize a cone at the origin that contains the box - might not be the minimum angle cone, but should be close
  Cone(const Vec3f &origin, const Vec3f &min, const Vec3f &max) {
    mOrigin = origin;
    mDir = ((min+max)/2 - origin).Normalized();
    mCosTheta = 1;
    mMinDist = 1e10;

    for (int i=0; i<2; i++) {
      float x = (i==0) ? min[0] : max[0];
      for (int j=0; j<2; j++) {
        float y = (j==0) ? min[1] : max[1];
        for (int k=0; k<2; k++) {
          float z = (k==0) ? min[2] : max[2];

          Vec3f p(x,y,z);
          mCosTheta = std::min(mCosTheta, mDir.Dot((p-origin).Normalized()));
          mMinDist = std::min(mMinDist,mDir.Dot(p-mOrigin));
        }
      }
    }

    if (mCosTheta < 0) {
      std::cerr<<"cone with greater than 90 degrees!"<<std::endl;
      exit(0);
    }

    if (mMinDist<0)
      std::cerr<<"min dist 0!!!!"<<std::endl;
  }


  // initialize a cone at the origin that contains the sphere
  Cone(const Vec3f &origin, float sphereRad, const Vec3f &sphereOrigin) {
    float dist = (sphereOrigin-origin).Length();
    if (dist < sphereRad) {
      std::cerr<<"cone with greater than 90 degrees!"<<std::endl;
      exit(0);
    }

    mOrigin = origin;
    mDir = (sphereOrigin - origin) / dist;
    mCosTheta = cos(atan(sphereRad / dist));
    mMinDist = mDir.Dot(sphereOrigin-mOrigin) - sphereRad;
    if (mMinDist<0)
      std::cerr<<"min dist 0!!!!"<<std::endl;
  }

  
  bool ContainsPoint(const Vec3f &p) const {
    Vec3f offset = p - mOrigin;
    return ((mDir.Dot(offset) >= mMinDist) &&
            (mDir.Dot(offset.Normalized()) >= mCosTheta));
  }

  void MergeSameAxisOrigin(const Cone &rhs) {
    mCosTheta = std::min(mCosTheta, rhs.mCosTheta);
    mMinDist = std::min(mMinDist, rhs.mMinDist);
  }

  void Expand(const Vec3f &min, const Vec3f &max) {

    for (int i=0; i<2; i++) {
      float x = (i==0) ? min[0] : max[0];
      for (int j=0; j<2; j++) {
        float y = (j==0) ? min[1] : max[1];
        for (int k=0; k<2; k++) {
          float z = (k==0) ? min[2] : max[2];

          Vec3f p(x,y,z);
          mCosTheta = std::min(mCosTheta, mDir.Dot((p-mOrigin).Normalized()));
          mMinDist = std::min(mMinDist,mDir.Dot(p-mOrigin));
        }
      }
    }

    if (mCosTheta < 0) {
      std::cerr<<"cone with greater than 90 degrees!"<<std::endl;
      exit(0);
    }

    if (mMinDist<0)
      std::cerr<<"min dist 0!!!!"<<std::endl;
  }
  

  bool Intersect(const Vec3f &_lineOrigin, const Vec3f &_lineDir,
                 float &t0, float &t1) const {

    return FastIntersect(_lineOrigin, _lineDir, t0, t1);
#if 0

    Vec3d coneOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);
    Vec3d coneDir(mDir[0], mDir[1], mDir[2]);
    Vec3d lineDir(_lineDir[0], _lineDir[1], _lineDir[2]);
    Vec3d lineOrigin(_lineOrigin[0], _lineOrigin[1], _lineOrigin[2]);

    double AdD = coneDir.Dot(lineDir);
    double cosSqr = (mCosTheta-COSTHETA_EPS)*(mCosTheta-COSTHETA_EPS);
    Vec3d E = lineOrigin - coneOrigin;
    double AdE = coneDir.Dot(E);
    double DdE = lineDir.Dot(E);
    double EdE = E.Dot(E);
    double c2 = AdD*AdD - cosSqr;
    double c1 = AdD*AdE - cosSqr*DdE;
    double c0 = AdE*AdE - cosSqr*EdE;
    double dot;

    Vec3d point;

    // Solve the quadratic.  Keep only those X for which Dot(A,X-V) >= 0.
    if (fabs(c2) >= 1e-4) {
      // c2 != 0
      double discr = c1*c1 - c0*c2;
      if (discr < 0) {
        // Q(t) = 0 has no real-valued roots.  The line does not
        // intersect the double-sided cone.
        return false;
      }
      
      else if (discr > 1e-4) {
        // Q(t) = 0 has two distinct real-valued roots.  However, one or
        // both of them might intersect the portion of the double-sided
        // cone "behind" the vertex.  We are interested only in those
        // intersections "in front" of the vertex.
        double root = sqrt(discr);
        double invC2 = ((double)1)/c2;

        double t = (-c1 - root)*invC2;
        point = lineOrigin + t*lineDir;
        E = point - coneOrigin;
        dot = E.Dot(coneDir);
        if (dot > mMinDist-MINDIST_EPS) {
          return true;
        }

        t = (-c1 + root)*invC2;
        point = lineOrigin + t*lineDir;
        E = point - coneOrigin;
        dot = E.Dot(coneDir);
        if (dot > mMinDist-MINDIST_EPS) {
          return true;
        }

        return false;
      }
      else {
        // One repeated real root (line is tangent to the cone).
        point = lineOrigin - (c1/c2)*lineDir;
        E = point - coneOrigin;
        if (E.Dot(coneDir) > mMinDist-MINDIST_EPS) {
          return true;
        }
        else {
          return false;
        }
      }
    }

    else if (fabs(c1) >= 1e-4) {
      // c2 = 0, c1 != 0 (D is a direction vector on the cone boundary)
      point = lineOrigin - (0.5*c0/c1)*lineDir;
      E = point - coneDir;
      dot = E.Dot(coneOrigin);
      if (dot > mMinDist-MINDIST_EPS) {
        return true;
      }
      else {
        return false;
      }
    }

    else if (fabs(c0) >= 1e-4) {
      // c2 = c1 = 0, c0 != 0
      return false;
    }
    else {
      // c2 = c1 = c0 = 0, cone contains ray V+t*D where V is cone vertex
      // and D is the line direction.
      return true;
    }
#endif
  }  




  // quick intersection test - only checks for double intersection
  bool FastIntersect(const Vec3f &_lineOrigin, const Vec3f &_lineDir,
                     float &t0, float &t1) const {
    Vec3d coneOrigin(mOrigin[0], mOrigin[1], mOrigin[2]);
    Vec3d coneDir(mDir[0], mDir[1], mDir[2]);
    Vec3d lineDir(_lineDir[0], _lineDir[1], _lineDir[2]);
    Vec3d lineOrigin(_lineOrigin[0], _lineOrigin[1], _lineOrigin[2]);

    double AdD = coneDir.Dot(lineDir);
    double cosSqr = (mCosTheta-COSTHETA_EPS)*(mCosTheta-COSTHETA_EPS);
    Vec3d E = lineOrigin - coneOrigin;
    double AdE = coneDir.Dot(E);
    double DdE = lineDir.Dot(E);
    double EdE = E.Dot(E);
    double c2 = AdD*AdD - cosSqr;
    double c1 = AdD*AdE - cosSqr*DdE;
    double c0 = AdE*AdE - cosSqr*EdE;
    double dot;

    Vec3d point;

    // Solve the quadratic.  Keep only those X for which Dot(A,X-V) >= 0.
    if (fabs(c2) >= 1e-4) {
      // c2 != 0
      double discr = c1*c1 - c0*c2;
      if (discr > 1e-4) {
        // Q(t) = 0 has two distinct real-valued roots.  However, one or
        // both of them might intersect the portion of the double-sided
        // cone "behind" the vertex.  We are interested only in those
        // intersections "in front" of the vertex.
        double root = sqrt(discr);
        double invC2 = ((double)1)/c2;
        bool hit = false;

        double t = (-c1 - root)*invC2;
        point = lineOrigin + t*lineDir;
        E = point - coneOrigin;
        dot = E.Dot(coneDir);
        if (dot > mMinDist-MINDIST_EPS) {
          hit = true;
        }
        t0 = (float)t;

        t = (-c1 + root)*invC2;
        point = lineOrigin + t*lineDir;
        E = point - coneOrigin;
        dot = E.Dot(coneDir);
        if (dot > mMinDist-MINDIST_EPS) {
          hit = true;
        }
        t1 = (float)t;

        if (t1<t0)
          std::swap(t0, t1);

        return hit;
      }
    }

    return false;
  }  



  //private:
  Vec3f mOrigin;
  Vec3f mDir;
  float mCosTheta;
  float mMinDist;
};


#endif
