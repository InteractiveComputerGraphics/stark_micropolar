#pragma once

#include "../core/Stark.h"

namespace stark
{
    /// Stores a matrix in row-major ordering for use with symx and easy conversion to Eigen matrices
    template<typename T, std::size_t Rows, std::size_t Cols>
    struct MatrixWrapper : std::array<T, Rows*Cols> {
        MatrixWrapper() = default;

        MatrixWrapper(std::array<T, Rows*Cols> array) : std::array<T, Rows * Cols>(array) {}

        MatrixWrapper(Eigen::Matrix<T, Rows, Cols> matrix) {
            const Eigen::Matrix<T,Rows,Cols> mat_t = matrix.transpose();
            for (int i = 0; i < mat_t.size(); ++i) {
                this->data()[i] = mat_t.data()[i];
            }
        }

        static MatrixWrapper zero() {
            MatrixWrapper matrix;
            matrix.fill(T(0.0));
            return matrix;
        }

        static MatrixWrapper identity() {
            const auto I = Eigen::Matrix<T, Rows, Cols>::Identity();
            return MatrixWrapper::from_matrix(I);
        }

        static MatrixWrapper<T, Rows, Cols> from_array(const std::array<T, Rows*Cols>& array) {
            return MatrixWrapper<T, Rows, Cols>(array);
        }

        static MatrixWrapper<T, Rows, Cols> from_matrix(const Eigen::Matrix<T, Rows, Cols> matrix) {
            return MatrixWrapper<T, Rows, Cols>(matrix);
        }

        Eigen::Matrix<T, Rows, Cols> to_matrix() const {
            return Eigen::Matrix<T, Cols, Rows>(this->data()).transpose();
        }

        constexpr std::size_t rows() const {
            return Rows;
        }

        constexpr std::size_t cols() const {
            return Cols;
        }

        T& operator()(const int row, const int col) {
            return this->at(row * Cols + col);
        }

        const T& operator()(const int row, const int col) const {
            return this->at(row * Cols + col);
        }
    };

    /// Stores a quaternion as an array [w,x,y,z] for use with symx and easy conversion to Eigen quaternions
    template <typename T>
    struct QuaternionWrapper : std::array<T, 4> {
        QuaternionWrapper() : std::array<T, 4>({T(1.0), T(0.0), T(0.0), T(0.0)}) {}
        QuaternionWrapper(std::array<T, 4> array) : std::array<T, 4>(array) {}
        QuaternionWrapper(Eigen::Quaternion<T> q) : std::array<T, 4>({q.w(), q.x(), q.y(), q.z()}) {}

        operator Eigen::Quaternion<T>() const {
            return this->to_quaternion();
        }

        static QuaternionWrapper identity() {
            QuaternionWrapper q;
            q.fill(T(0.0));
            q[0] = T(1.0);
            return q;
        }

        static QuaternionWrapper from_array(const std::array<T, 4>& array) {
            return QuaternionWrapper(array);
        }

        static QuaternionWrapper from_quaternion(const Eigen::Quaternion<T>& q) {
            return QuaternionWrapper(q);
        }

        Eigen::Quaternion<T> to_quaternion() const {
            const auto& arr = this->data();
            return Eigen::Quaterniond(arr[0], arr[1], arr[2], arr[3]);
        }

        Eigen::Matrix3<T> to_rotation_matrix() const {
            return this->to_quaternion().toRotationMatrix();
        }
    };

    using MatrixWrapper3d = MatrixWrapper<double, 3, 3>;
    using MatrixWrapper2d = MatrixWrapper<double, 2, 2>;
    using QuaternionWrapperd = QuaternionWrapper<double>;
}
