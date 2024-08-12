#include "symx_quaternion.h"

using namespace stark;

using Scalar = symx::Scalar;
using Vector = symx::Vector;
using Matrix = symx::Matrix;

Quaternion Quaternion::from_vec4(const Vector& q) {
    assert(q.size() == 4);
    return Quaternion { q };
}

symx::Vector Quaternion::to_vec3(symx::Energy& e) const {
	symx::Vector v = e.make_zero_vector(3);
	auto& q = this->coeffs;

	v[0] = q[1];
	v[1] = q[2];
	v[2] = q[3];
	return v;
}

// Normalized quaternion to rotation matrix as in Eigen
symx::Matrix Quaternion::to_rotation(symx::Energy& e) const {
	symx::Matrix R = e.make_zero_matrix({ 3, 3 });
	auto& q = this->coeffs;

	const Scalar& qw = q[0];
	const Scalar& qx = q[1];
	const Scalar& qy = q[2];
	const Scalar& qz = q[3];

	Scalar tx = 2.0 * qx;
	Scalar ty = 2.0 * qy;
	Scalar tz = 2.0 * qz;
	Scalar twx = tx * qw;
	Scalar twy = ty * qw;
	Scalar twz = tz * qw;
	Scalar txx = tx * qx;
	Scalar txy = ty * qx;
	Scalar txz = tz * qx;
	Scalar tyy = ty * qy;
	Scalar tyz = tz * qy;
	Scalar tzz = tz * qz;

	R(0, 0) = 1.0 - (tyy + tzz);
	R(0, 1) = txy - twz;
	R(0, 2) = txz + twy;
	R(1, 0) = txy + twz;
	R(1, 1) = 1.0 - (txx + tzz);
	R(1, 2) = tyz - twx;
	R(2, 0) = txz - twy;
	R(2, 1) = tyz + twx;
	R(2, 2) = 1.0 - (txx + tyy);

	return R;
}

Quaternion Quaternion::conjugated() const {
	auto& q = this->coeffs;
	// w, x, y, z = q
	const Scalar& a = q[0];
	const Scalar& b = q[1];
	const Scalar& c = q[2];
	const Scalar& d = q[3];
	return Quaternion::from_vec4(Vector({ a, -b, -c, -d }));
}

Quaternion Quaternion::normalized() const {
	auto& q = this->coeffs;
	return Quaternion::from_vec4(q / q.norm());
}

Quaternion Quaternion::inverse() const {
	// TODO: Is this correct?
    return this->conjugated();
}

// Quaternion multiplication
Quaternion Quaternion::prod(const Quaternion& other) const {
	auto& q1 = this->coeffs;
	auto& q2 = other.coeffs;

	// w, x, y, z = q
	const Scalar& a = q1[0];
	const Scalar& b = q1[1];
	const Scalar& c = q1[2];
	const Scalar& d = q1[3];

	const Scalar& e = q2[0];
	const Scalar& f = q2[1];
	const Scalar& g = q2[2];
	const Scalar& h = q2[3];

	return Quaternion::from_vec4(Vector({
		a*e - b*f - c*g - d*h,
		b*e + a*f + c*h - d*g,
		a*g - b*h + c*e + d*f,
		a*h + b*g - c*f + d*e
	}));
}

Quaternion Quaternion::time_integration_with_global_w(const Vector& w_glob, const Scalar& dt) const {
	const Quaternion& q_start = *this;
	const Quaternion w = Quaternion::from_vec4(Vector({ dt.get_zero(), w_glob[0], w_glob[1], w_glob[2] }));
	return q_start + 0.5 * dt * w.prod(*this);
}
/*
Quaternion Quaternion::time_integration_exp_with_global_w(const Vector& w_glob, const Scalar& dt) const {
	const Quaternion& q_start = *this;

	Vector w_({ dt.get_zero(), w_glob[0], w_glob[1], w_glob[2] });
	//w_[1] += 1e-7;
	//w_[2] += 1e-7;
	//w_[3] += 1e-7;

	Vector arg = 0.5 * dt * w_;
	Scalar arg_norm = arg.squared_norm();

	//Scalar sinc_arg_norm = symx::branch(arg_norm - 1e-14, sin(arg_norm) / arg_norm, dt.get_one());
	//Scalar sinc_arg_norm = sin(arg_norm) / arg_norm;
	Scalar sinc_arg_norm = arg_norm.sqrt_sinc();

	Quaternion q_exp = Quaternion::from_vec4(Vector({
		arg_norm.sqrt_cos(),
		sinc_arg_norm * arg[1],
		sinc_arg_norm * arg[2],
		sinc_arg_norm * arg[3],
	}));

	return q_exp.prod(q_start);
}
*/
Quaternion Quaternion::operator*(const Quaternion& other) const {
    return this->prod(other);
}

Quaternion Quaternion::operator+(const Quaternion& other) const {
    return Quaternion::from_vec4(this->coeffs + other.coeffs);
}

Quaternion Quaternion::operator*(const symx::Scalar& other) const {
    return Quaternion::from_vec4(this->coeffs * other);
}

Quaternion stark::operator*(const symx::Scalar& scalar, const Quaternion& quat) {
    return Quaternion::from_vec4(scalar * quat.coeffs);
}
