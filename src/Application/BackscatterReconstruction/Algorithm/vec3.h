
#ifndef VEC3_H
#define VEC3_H

#include <math.h>

template <typename T>
class Vec3 {
public:
  template <typename T2>
  Vec3(const Vec3<T2> &rhs) { v[0]=(T)rhs[0]; v[1]=(T)rhs[1]; v[2]=(T)rhs[2]; }
  Vec3() { v[0]=0; v[1]=0; v[2]=0; }
  Vec3(T x, T y, T z) { v[0]=x; v[1]=y; v[2]=z; }
  Vec3(const T *b) { v[0]=b[0]; v[1]=b[1]; v[2]=b[2]; }
  
  T Dot(const Vec3 &b) const { return v[0]*b.v[0] + v[1]*b.v[1] + v[2]*b.v[2]; }
  Vec3 Cross(const Vec3 &b) const { return Vec3(v[1]*b.v[2] - v[2]*b.v[1],
                                                v[2]*b.v[0] - v[0]*b.v[2],
                                                v[0]*b.v[1] - v[1]*b.v[0]); }
  T Length() const { return sqrtf(Dot(*this)); }
  T LengthSquared() const { return Dot(*this); }
  void Normalize() { T l=Length(); v[0]/=l; v[1]/=l; v[2]/=l; }
  Vec3 Normalized() const { T l=Length(); return Vec3(v[0]/l, v[1]/l, v[2]/l); }

  T& operator[](int i) { return v[i]; }
  const T& operator[](int i) const { return v[i]; }

  Vec3 operator+(const Vec3 &b) const { return Vec3(v[0]+b.v[0], v[1]+b.v[1], v[2]+b.v[2]); }
  void operator+=(const Vec3 &b) { v[0]+=b.v[0]; v[1]+=b.v[1]; v[2]+=b.v[2]; }
  Vec3 operator-(const Vec3 &b) const { return Vec3(v[0]-b.v[0], v[1]-b.v[1], v[2]-b.v[2]); }
  void operator-=(const Vec3 &b) { v[0]-=b.v[0]; v[1]-=b.v[1]; v[2]-=b.v[2]; }
  void operator*=(const T x) { v[0]*=x; v[1]*=x; v[2]*=x; }
  void operator/=(const T x) { v[0]/=x; v[1]/=x; v[2]/=x; }

  friend Vec3 operator*(const Vec3 &lhs, T rhs) { return Vec3(lhs.v[0]*rhs, lhs.v[1]*rhs, lhs.v[2]*rhs); }
  friend Vec3 operator*(T lhs, const Vec3 &rhs) { return Vec3(lhs*rhs.v[0], lhs*rhs.v[1], lhs*rhs.v[2]); }
  friend Vec3 operator/(const Vec3 &lhs, T rhs) { return Vec3(lhs.v[0]/rhs, lhs.v[1]/rhs, lhs.v[2]/rhs); }

  friend Vec3 operator-(const Vec3 &lhs) { return Vec3(-lhs.v[0], -lhs.v[1], -lhs.v[2]); }

private:
  T v[3];
};


typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3<int> Vec3i;

#endif
