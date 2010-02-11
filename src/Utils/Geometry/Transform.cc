/*
For more information, please see: http://software.sci.utah.edu

The MIT License

Copyright (c) 2009 Scientific Computing and Imaging Institute,
University of Utah.


Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/

#include <Utils/Geometry/Transform.h>
#include <Utils/Geometry/Point.h>

namespace Utils {

Transform::Transform()
{
	load_identity();
}

Transform::Transform(const Transform& copy) 
{
	this->mat_ = copy.mat_;
	this->inverse_mat_ = copy.inverse_mat_;
	this->inverse_valid_ = copy.inverse_valid_;
}

Transform& Transform::operator=(const Transform& copy) 
{
	this->mat_ = copy.mat_;
	this->inverse_mat_ = copy.inverse_mat_;
	this->inverse_valid_ = copy.inverse_valid_;

	return (*this);
}

Transform::Transform(const Point& p, const Vector& i, 
									const Vector& j, const Vector& k)
{
	load_basis(p, i, j, k);
}

void Transform::load_basis(const Point &p,
											const Vector &x,
											const Vector &y,
											const Vector &z)
{
	load_frame(x,y,z);
	post_translate(Vector(-p));
}

void	Transform::load_frame(const Vector& x, 
											const Vector& y, 
											const Vector& z)
{
	mat_(3, 3) = 1.0;
	mat_(0, 3) = mat_(1, 3) = mat_(2, 3) = 0.0;
	mat_(3, 0) = mat_(3, 1) = mat_(3, 2) = 0.0;

	mat_(0, 0) = x.x();
	mat_(0, 1) = x.y();
	mat_(0, 2) = x.z();

	mat_(1, 0) = y.x();
	mat_(1, 1) = y.y();
	mat_(1, 2) = y.z();

	mat_(2, 0) = z.x();
	mat_(2, 1) = z.y();
	mat_(2, 2) = z.z();

	this->inverse_valid_ = false;
}

void Transform::post_transform(const Transform& trans)
{
	this->mat_ *= trans.mat_;
	this->inverse_valid_ = false;
}

void Transform::pre_transform(const Transform& trans)
{
	this->mat_ = trans.mat_ * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::post_mult_matrix(const Matrix& m)
{
	this->mat_ *= m;
	this->inverse_valid_ = false;
}

void Transform::pre_mult_matrix(const Matrix& m)
{
	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::pre_scale(const Vector& v)
{
	Matrix m;
	Transform::BuildScale(m, v);

	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::post_scale(const Vector& v)
{
	Matrix m;
	Transform::BuildScale(m, v);

	this->mat_ *= m;
	this->inverse_valid_ = false;
}

void Transform::pre_shear(const Vector& s, const Plane& p)
{
	Matrix m;
	Transform::BuildShear(m, s, p);

	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::post_shear(const Vector& s, const Plane& p)
{
	Matrix m;
	Transform::BuildShear(m, s, p);

	this->mat_ *= m;
	this->inverse_valid_ = false;
}

void Transform::pre_translate(const Vector& v)
{
	Matrix m;
	Transform::BuildTranslate(m, v);

	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::post_translate(const Vector& v)
{
	Matrix m;
	Transform::BuildTranslate(m, v);

	this->mat_ *= m;
	this->inverse_valid_ = false;
}

void Transform::pre_rotate(double angle, const Vector& axis)
{
	Matrix m;
	Transform::BuildRotate(m, angle, axis);

	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}       

void Transform::post_rotate(double angle, const Vector& axis)
{
	Matrix m;
	Transform::BuildRotate(m, angle, axis);

	this->mat_ *= m;
	this->inverse_valid_ = false;
}

bool Transform::rotate(const Vector& from, const Vector& to)
{
	Vector t(to); 
	t.normalize();
	Vector f(from); 
	f.normalize();
	Vector axis(Cross(f, t));

	if (axis.length2() < 1.0e-8) 
	{
		// Vectors are too close to each other to get a stable axis of
		// rotation, so return.
		return (false);
	}

	double sinth=axis.length();
	double costh=Dot(f,t);
	if(Abs(sinth) < 1.0e-9)
	{
		if(costh > 0.0)
			return (false); // no rotate;
		else 
		{
			// from and to are in opposite directions, find an axis of rotation
			// Try the Z axis first.  This will fail if from is along Z, so try
			// Y next.  Then rotate 180 degrees.
			axis = Cross(from, Vector(0.0,0.0,1.0));
			if(axis.length2() < 1.0e-9)
				axis = Cross(from, Vector(0.0,1.0,0.0));
			axis.normalize();
			post_rotate(Pi(), axis);
		}
	} 
	else 
	{
		post_rotate(Atan2(sinth, costh), axis.normal());
	}
	return (true);
}

void Transform::pre_permute(int xmap, int ymap, int zmap)
{
	Matrix m;
	Transform::BuildPermute(m, xmap, ymap, zmap, true);

	this->mat_ = m * this->mat_;
	this->inverse_valid_ = false;
}

void Transform::post_permute(int xmap, int ymap, int zmap)
{
	Matrix m;
	Transform::BuildPermute(m, xmap, ymap, zmap, false);

	this->mat_ *= m;
	this->inverse_valid_ = false;
}

Point Transform::project(const Point& p) const
{
	return this->mat_ * p;
}

PointF Transform::project(const PointF& p) const
{
	return this->mat_ * p;
}

Vector Transform::project(const Vector& v) const
{
	return this->mat_ * v;
}

VectorF Transform::project(const VectorF& v) const
{
	return this->mat_ * v;
}

Point Transform::unproject(const Point& p) const
{
	compute_inverse();
	return this->inverse_mat_ * p;
}

Vector Transform::unproject(const Vector& v) const
{
	compute_inverse();
	return this->inverse_mat_ * v;
}

PointF Transform::unproject(const PointF& p) const
{
	compute_inverse();
	return this->inverse_mat_ * p;
}

VectorF Transform::unproject(const VectorF& v) const
{
	compute_inverse();
	return this->inverse_mat_ * v;
}

const Matrix& Transform::get_matrix() const
{
	return this->mat_;
}

const Matrix& Transform::get_inverse_matrix() const
{
	compute_inverse();
	return this->inverse_mat_;
}

void Transform::get(double* data) const
{
	std::memcpy(data, this->mat_.data(), sizeof(double) * 16);
}

void Transform::set(const double* data)
{
	std::memcpy(this->mat_.data(), data, sizeof(double) * 16);
	this->inverse_valid_ = false;
}

void Transform::load_identity()
{
	mat_ = Matrix::IDENTITY_C;
	inverse_mat_ = Matrix::IDENTITY_C;
	this->inverse_valid_ = true;
}

void Transform::invert()
{
	compute_inverse();

	Matrix tmp = mat_;
	mat_ = inverse_mat_;
	inverse_mat_ = tmp;
}

void Transform::compute_inverse() const
{
	if ( !this->inverse_valid_ )
	{
		Invert(this->mat_, this->inverse_mat_);
		this->inverse_valid_ = true;
	}
}

void Transform::perspective( const Point& eyep, const Point& lookat, const Vector& up, 
											 double fovy, double znear, double zfar, double aspect )
{
 	Vector z(eyep - lookat); 
	z.normalize();

 	Vector x(Cross(up, z)); 
	x.normalize();

 	Vector y(Cross(z, x));

	// View transformation
	Transform tf(eyep, x, y, z);
	pre_transform(tf);

	// Perspective projection
	double f = Cot(DegreeToRadian(fovy * 0.5));

	Matrix proj = Matrix::ZERO_C;
	proj(0, 0) = f / aspect;
	proj(1, 1) = f;
	proj(2, 2) = (zfar + znear) / (znear - zfar);
	proj(2, 3) = 2 * zfar * znear / (znear - zfar);
	proj(3, 2) = -1;

	this->mat_ = proj * this->mat_;
	this->inverse_valid_ = false;
 }

void Transform::BuildTranslate( Matrix& m, const Vector& v )
{
	m = Matrix::IDENTITY_C;
	m(0, 3) = v.x();
	m(1, 3) = v.y();
	m(2, 3) = v.z();
}

// rotate into a new frame (z=shear-fixed-plane, y=projected shear vector),
// shear in y (based on value of z), rotate back to original frame
void Transform::BuildShear( Matrix& m, const Vector& s, const Plane& p )
{    
	m = Matrix::IDENTITY_C;

	Vector sv(p.project(s));      // s projected onto p
	Vector dn(s-sv);				// difference (in normal direction) between s and sv
	double d = Dot(dn,p.normal());

	// shear vector lies in shear fixed plane, return identity.
	if (Abs(d)<1.0e-8) 
	{ 
		return;
	}

	double yshear = sv.length()/d; // compute the length of the shear vector,
	// after the normal-to-shear-plane component
	// has been made unit-length.
	Vector svn(sv);
	svn.normalize();      // normalized vector for building orthonormal basis

	Vector su(Cross(svn, p.normal()));

	// the rotation to take the z-axis to the shear normal
	// and the y-axis to the projected shear vector
	Transform r;  
	r.load_frame(su, svn, p.normal());

	// the shear matrix in the new frame
	Matrix shear = Matrix::IDENTITY_C;
	shear(1, 2) = yshear;
	shear(1, 3) = -yshear * p.distance();

	m = r.inverse_mat_*shear*r.mat_;
}

void Transform::BuildScale( Matrix& m, const Vector& v )
{
	m = Matrix::IDENTITY_C;
	m(0, 0) = v.x();
	m(1, 1) = v.y();
	m(2, 2) = v.z();
}

void Transform::BuildRotate( Matrix& m, double angle, const Vector& axis )
{
	double sintheta = Sin(angle);
	double costheta = Cos(angle);
	double ux=axis.x();
	double uy=axis.y();
	double uz=axis.z();

	m(0, 0) = ux*ux + costheta*(1.0-ux*ux);
	m(0, 1) = ux*uy*(1.0-costheta) - uz*sintheta;
	m(0, 2) = uz*ux*(1.0-costheta) + uy*sintheta;
	m(0, 3) = 0.0;

	m(1, 0) = ux*uy*(1.0-costheta) + uz*sintheta;
	m(1, 1) = uy*uy + costheta*(1-uy*uy);
	m(1, 2) = uy*uz*(1.0-costheta)-ux*sintheta;
	m(1, 3) = 0.0;

	m(2, 0) = uz*ux*(1.0-costheta) - uy*sintheta;
	m(2, 1) = uy*uz*(1.0-costheta) + ux*sintheta;
	m(2, 2) = uz*uz + costheta*(1-uz*uz);
	m(2, 3) = 0.0;

	m(3, 0) = 0.0;
	m(3, 1) = 0.0;
	m(3, 2) = 0.0;
	m(3, 3) = 1.0;
}

void Transform::BuildPermute( Matrix& m, int xmap, int ymap, int zmap, bool pre )
{
	m = Matrix::ZERO_C;

	m(3, 3)=1.0;

	int x = xmap < 0 ? (-1 - xmap) : (xmap - 1);
	int y = ymap < 0 ? (-1 - ymap) : (ymap - 1);
	int z = zmap < 0 ? (-1 - zmap) : (zmap - 1);

	if (pre) 
	{    
		// for each row, set the mapped row
		m(0, x) = Sign(xmap) * 1.0;
		m(1, y) = Sign(ymap) * 1.0;
		m(2, z) = Sign(zmap) * 1.0;
	} 
	else 
	{      
		// for each column, set the mapped column
		m(x, 0) = Sign(xmap) * 1.0;
		m(y, 1) = Sign(ymap) * 1.0;
		m(z, 2) = Sign(zmap) * 1.0;
	}
}

Point operator*( const Transform& t, const Point& d )
{
	return t.project(d);
}

Vector operator*( const Transform& t, const Vector& d )
{
	return t.project(d);
}

PointF operator*( const Transform& t, const PointF& d )
{
	return t.project(d);
}

VectorF operator*( const Transform& t, const VectorF& d )
{
	return t.project(d);
}



} // namespace Utils
