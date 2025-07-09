#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include <Eigen/Core>
#include <g2o/core/base_vertex.h>

#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_multi_edge.h>
#include <g2o/core/base_unary_edge.h>

#include <g2o/core/block_solver.h>
#include <g2o/core/sparse_block_matrix.h>
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>
#include <g2o/solvers/pcg/linear_solver_pcg.h>

#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>

#include <g2o/core/robust_kernel_impl.h>

class PositionVertex : public g2o::BaseVertex<3, Eigen::Vector3d> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    PositionVertex() = default;
    ~PositionVertex() = default;

    void setToOriginImpl() override
    {
        _estimate << 0.5, 0.5, 0.5;
    }

    void oplusImpl(const double* update) override
    {
        _estimate += Eigen::Vector3d(update);
    }

    bool read(std::istream& in) override { return false; }
    bool write(std::ostream& out) const override { return false; }
};

class OptimLocEdge : public g2o::BaseUnaryEdge<1, double, PositionVertex> {
private:
    Eigen::Vector3d p_known_;
    double bias_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    explicit OptimLocEdge(const Eigen::Vector3d& p_known, double bias = 0)
        : p_known_(p_known)
        , bias_(bias)
    {
    }

    void computeError() override
    {
        const auto* v = dynamic_cast<const PositionVertex*>(_vertices[0]);
        const Eigen::Vector3d p = v->estimate();
        _error(0, 0) = _measurement - (1 + bias_) * (p - p_known_).norm();
    }

    void linearizeOplus() override
    {
        const auto* v = dynamic_cast<const PositionVertex*>(_vertices[0]);
        const Eigen::Vector3d p = v->estimate();

        const Eigen::Vector3d difference = p - p_known_;
        const double distance = difference.norm();

        for (int i = 0; i < 3; i++) {
            _jacobianOplusXi(0, i) = -(1 + bias_) * difference[i] / distance;
        }
    }

    bool read(std::istream& in) override { return false; }
    bool write(std::ostream& out) const override { return false; }
};

#endif