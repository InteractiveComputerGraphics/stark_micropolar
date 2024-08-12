#pragma once
#include <symx>

namespace stark
{
    class Quaternion
    {
    public:
        symx::Vector coeffs;

        static Quaternion from_vec4(const symx::Vector& q);

        /// Returns the vector (or imaginary) part of the quaternion
        symx::Vector to_vec3(symx::Energy& e) const;
        /// Converts the quaternion to a rotation matrix (does not normalize)
        symx::Matrix to_rotation(symx::Energy& e) const;

        Quaternion conjugated() const;
        Quaternion normalized() const;
        Quaternion inverse() const;

        /// Computes the quaternion product with the other quaternion
        Quaternion prod(const Quaternion& other) const;
        Quaternion dot(const Quaternion& other) const;

        Quaternion time_integration_with_global_w(const symx::Vector& w_glob, const symx::Scalar& dt) const;
        //Quaternion time_integration_exp_with_global_w(const symx::Vector& w_glob, const symx::Scalar& dt) const;

        Quaternion operator+(const Quaternion& other) const;
        Quaternion operator*(const Quaternion& other) const;
        Quaternion operator*(const symx::Scalar& other) const;
    };

    Quaternion operator*(const symx::Scalar& scalar, const Quaternion& quat);

    template<typename STATIC_VECTOR>
    Quaternion make_quaternion(symx::Energy& energy, std::vector<STATIC_VECTOR>& data, symx::Index conn) {
        symx::Vector q_vec = energy.make_vector(data, conn);
        return Quaternion::from_vec4(q_vec);
    }

    template<typename STATIC_VECTOR>
    std::vector<Quaternion> make_quaternions(symx::Energy& energy, std::vector<STATIC_VECTOR>& data, const std::vector<symx::Index>& conn) {
        std::vector<Quaternion> q;
        std::vector<symx::Vector> q_vec = energy.make_vectors(data, conn);
        for (const auto& q_vec_i : q_vec) {
            q.push_back(Quaternion::from_vec4(q_vec_i));
        }
        return q;
    }
}
