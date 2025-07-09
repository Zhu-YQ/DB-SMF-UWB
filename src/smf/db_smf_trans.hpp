#include <Eigen/Core>
#include <Eigen/StdVector>
#include <vector>

#include "optimization.hpp"

// # define BINARY_SEARCH

using Vec9d = Eigen::Matrix<double, 9, 1>;
using Mat9d = Eigen::Matrix<double, 9, 9>;
using Vec3d = Eigen::Vector3d;
using Mat3d = Eigen::Matrix3d;
using Mat93d = Eigen::Matrix<double, 9, 3>;

class DBSMFtrans {
private:
    static const int n_ = 9;
    int m_;

    double NLOS_THRESH_;

    double calMinWeightedTraceBetaNumeric(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2)
    {
        auto ObjFunc = [P1, P2, this](double beta) {
            const double w_x = P_(0, 0) + P_(1, 1) + P_(2, 2);
            const double tr_P1_x = P1(0, 0) + P1(1, 1) + P1(2, 2);
            const double tr_P2_x = P2(0, 0) + P2(1, 1) + P2(2, 2);

            const double w_y = P_(3, 3) + P_(4, 4) + P_(5, 5);
            const double tr_P1_y = P1(3, 3) + P1(4, 4) + P1(5, 5);
            const double tr_P2_y = P2(3, 3) + P2(4, 4) + P2(5, 5);

            const double w_z = P_(6, 6) + P_(7, 7) + P_(8, 8);
            const double tr_P1_z = P1(6, 6) + P1(7, 7) + P1(8, 8);
            const double tr_P2_z = P2(6, 6) + P2(7, 7) + P2(8, 8);

            const double val_x = w_x * ((1.0 / (1 - beta)) * tr_P1_x + (1.0 / beta) * tr_P2_x);
            const double val_y = w_y * ((1.0 / (1 - beta)) * tr_P1_y + (1.0 / beta) * tr_P2_y);
            const double val_z = w_z * ((1.0 / (1 - beta)) * tr_P1_z + (1.0 / beta) * tr_P2_z);

            return val_x + val_y + val_z;
        };

        double a = 1e-30, b = 1.0 - 1e-30;
#ifdef BINARY_SEARCH

        while (1) {
            double x1 = a + (b - a) / 3;
            double x2 = a + 2.0 * (b - a) / 3;

            if (std::abs(x1 - x2) < 1e-5) {
                return (x1 + x2) / 2;
            }

            const double fx1 = ObjFunc(x1);
            const double fx2 = ObjFunc(x2);

            if (fx1 < fx2) {
                b = x2;
            } else {
                a = x1;
            }
        }
#else
        const double tau = (std::sqrt(5) + 1) / 2;
        double y = a + (b - a) / std::pow(tau, 2);
        double z = a + (b - a) / tau;
        while (b - a > 1e-5) {
            const double fy = ObjFunc(y);
            const double fz = ObjFunc(z);

            if (fy <= fz) {
                b = z;
                z = y;
                y = a + (b - a) / std::pow(tau, 2);
            }

            else {
                a = y;
                y = z;
                z = a + (b - a) / tau;
            }
        }
#endif
        return (a + b) / 2;
    }

    Eigen::MatrixXd calHi(const Eigen::Vector3d& q_i) const
    {
        const Eigen::Vector3d p = getTagPosition();

        Eigen::Vector3d dp = p - q_i;
        double distance = dp.norm();

        Eigen::MatrixXd H_i;
        H_i.setZero(1, n_);
        H_i(0, 0) = dp(0) / distance;
        H_i(0, 3) = dp(1) / distance;
        H_i(0, 6) = dp(2) / distance;

        return H_i;
    }

    double calMinWeightedTraceRho(const Eigen::MatrixXd& P, const Eigen::MatrixXd& H, const Eigen::MatrixXd& R, const Eigen::MatrixXd& e)
    {
        auto ObjFunc = [P, H, R, e](double rho) {
            const auto linearized_P = H * P * H.transpose();

            const auto W = (1.0 / (1.0 - rho)) * linearized_P + (1.0 / rho) * R;
            const auto K = (1.0 / (1.0 - rho)) * P * H.transpose() * W.inverse();

            Mat9d P_rho = (1.0 / (1.0 - rho)) * P - K * H * (1.0 / (1.0 - rho)) * P;

            double delta = 1.0 - (e.transpose() * W.inverse() * e)(0, 0);
            if (delta > 0) {
                P_rho = delta * P_rho;
            }

            double w_x = P(0, 0) + P(1, 1) + P(2, 2);
            double w_y = P(3, 3) + P(4, 4) + P(5, 5);
            double w_z = P(6, 6) + P(7, 7) + P(8, 8);

            double tr_P_x = P_rho(0, 0) + P_rho(1, 1) + P_rho(2, 2);
            double tr_P_y = P_rho(3, 3) + P_rho(4, 4) + P_rho(5, 5);
            double tr_P_z = P_rho(6, 6) + P_rho(7, 7) + P_rho(8, 8);
            return w_x * tr_P_x + w_y * tr_P_y + w_z * tr_P_z;
        };

        double a = 1e-30, b = 1.0 - 1e-30;
#ifdef BINARY_SEARCH
        while (1) {
            double x1 = a + (b - a) / 3;
            double x2 = a + 2.0 * (b - a) / 3;

            if (std::abs(x1 - x2) < 1e-5) {
                return (x1 + x2) / 2;
            }

            double fx1 = ObjFunc(x1);
            double fx2 = ObjFunc(x2);

            if (fx1 < fx2) {
                b = x2;
            } else {
                a = x1;
            }
        }
#else
        const double tau = (std::sqrt(5) + 1) / 2;
        double y = a + (b - a) / std::pow(tau, 2);
        double z = a + (b - a) / tau;
        while (b - a > 1e-5) {
            const double fy = ObjFunc(y);
            const double fz = ObjFunc(z);

            if (fy <= fz) {
                b = z;
                z = y;
                y = a + (b - a) / std::pow(tau, 2);
            }

            else {
                a = y;
                y = z;
                z = a + (b - a) / tau;
            }
        }
#endif
        return (a + b) / 2;
    }

public:
    Vec9d x_;
    Mat9d P_;

    DBSMFtrans(int num_anchors, double nlos_thresh, double initial_robot_px_bound,
        double initial_robot_py_bound, double initial_robot_pz_bound,
        double initial_robot_v_bound, double initial_robot_a_bound)
        : m_(num_anchors)
        , x_(Vec9d::Zero())
        , P_(Mat9d::Zero())
        , NLOS_THRESH_(nlos_thresh)
    {
        P_(0, 0) = n_ * std::pow(initial_robot_px_bound, 2);
        P_(1, 1) = n_ * std::pow(initial_robot_v_bound, 2);
        P_(2, 2) = n_ * std::pow(initial_robot_a_bound, 2);

        P_(3, 3) = n_ * std::pow(initial_robot_py_bound, 2);
        P_(4, 4) = n_ * std::pow(initial_robot_v_bound, 2);
        P_(5, 5) = n_ * std::pow(initial_robot_a_bound, 2);

        P_(6, 6) = n_ * std::pow(initial_robot_pz_bound, 2);
        P_(7, 7) = n_ * std::pow(initial_robot_v_bound, 2);
        P_(8, 8) = n_ * std::pow(initial_robot_a_bound, 2);
    }

    Vec3d getTagPosition() const { return { x_(0), x_(3), x_(6) }; }

    std::vector<double> getTagPositionBound() const
    {
        return { std::sqrt(P_(0, 0)), std::sqrt(P_(3, 3)), std::sqrt(P_(6, 6)) };
    }

    static Eigen::Vector3d
    LocateOnce(const Eigen::VectorXd& y,
        const std::vector<Eigen::Vector3d>& aps_vec,
        const Eigen::Vector3d& p0 = Eigen::Vector3d::Ones(),
        const std::vector<double>& bias_vec = {})
    {

        auto* v = new PositionVertex();
        v->setId(0);
        v->setEstimate(p0);

        using BlockSolver = g2o::BlockSolver<g2o::BlockSolverTraits<3, 1>>;
        using LinearSolverType = g2o::LinearSolverEigen<BlockSolver::PoseMatrixType>;
        auto solver = new g2o::OptimizationAlgorithmDogleg(g2o::make_unique<BlockSolver>(g2o::make_unique<LinearSolverType>()));
        g2o::SparseOptimizer optimizer;
        optimizer.setAlgorithm(solver);

        optimizer.addVertex(v);

        for (int i = 0; i < aps_vec.size(); i++) {
            double bias = 0;
            if (!bias_vec.empty()) {
                bias = bias_vec.at(i);
            }

            auto edge = new OptimLocEdge(aps_vec.at(i), bias);
            edge->setId(i);
            edge->setVertex(0, v);
            edge->setMeasurement(y(i));
            edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

            auto robust_kernel = new g2o::RobustKernelGemanMcClure();
            robust_kernel->setDelta(0.5);
            edge->setRobustKernel(robust_kernel);

            optimizer.addEdge(edge);
        }

        optimizer.initializeOptimization();
        optimizer.optimize(20);

        return v->estimate();
    }

    void initTagPosition(const Eigen::VectorXd& y, double tag_height,
        const std::vector<Eigen::Vector3d>& aps_vec)
    {
        Eigen::Vector3d p0 = LocateOnce(y, aps_vec);
        x_(0) = p0(0);
        x_(3) = p0(1);
        x_(6) = tag_height;
        // x_(6) = p0(2);
    }

    void predict(double dt, double motion_robot_jx_bound, double motion_robot_jy_bound, double motion_robot_jz_bound)
    {
        Mat9d F = Mat9d::Identity();
        F(0, 1) = F(3, 4) = F(6, 7) = dt;
        F(0, 2) = F(3, 5) = F(6, 8) = 0.5 * dt * dt;
        F(1, 2) = F(4, 5) = F(7, 8) = dt;

        x_ = F * x_;

        Mat93d G = Mat93d::Zero();
        G(0, 0) = G(3, 1) = G(6, 2) = (1.0 / 6) * dt * dt * dt;
        G(1, 0) = G(4, 1) = G(7, 2) = (1.0 / 2) * dt * dt;
        G(2, 0) = G(5, 1) = G(8, 2) = dt;

        Mat3d Q = Mat3d::Zero();
        Q(0, 0) = 3 * std::pow(motion_robot_jx_bound, 2);
        Q(1, 1) = 3 * std::pow(motion_robot_jy_bound, 2);
        Q(2, 2) = 3 * std::pow(motion_robot_jz_bound, 2);

        const auto P1 = F * P_ * F.transpose();
        const auto P2 = G * Q * G.transpose();

        const double beta = calMinWeightedTraceBetaNumeric(P1, P2);

        P_ = (1.0 / (1 - beta)) * P1 + (1.0 / beta) * P2;
    }

    void update_seq(const Eigen::VectorXd& y, double noise_bound, const std::vector<Eigen::Vector3d>& aps_vec)
    {
        for (int i = 0; i < m_; i++) {
            const double y_hat_i = (getTagPosition() - aps_vec[i]).norm();

            Eigen::Matrix<double, 1, 1> e_i;
            e_i << y(i) - y_hat_i;

            if (std::abs(e_i(0, 0)) > NLOS_THRESH_) {
                // std::cout << "NLOS!\n";
                continue;
            }

            const Eigen::MatrixXd H_i = calHi(aps_vec.at(i));

            const auto linearized_P = H_i * P_ * H_i.transpose();

            Eigen::Matrix<double, 1, 1> R_i;

            // const double REMAINDER_RADIUS = 0.04;
            const double REMAINDER_RADIUS = 0.03;
            R_i(0, 0) = std::pow(noise_bound + REMAINDER_RADIUS, 2);

            const double rho = calMinWeightedTraceRho(P_, H_i, R_i, e_i);

            const auto W_i = ((1.0 / (1.0 - rho))) * linearized_P + (1.0 / rho) * R_i;
            const auto K_i = (1.0 / (1.0 - rho)) * P_ * H_i.transpose() * W_i.inverse();

            double delta = 1.0 - e_i.transpose() * W_i.inverse()(0, 0) * e_i;
            if (delta <= 0) {
                continue;
            }

            x_ = x_ + K_i * e_i;
            P_ = (1.0 / (1.0 - rho)) * P_ - K_i * H_i * (1.0 / (1.0 - rho)) * P_;

            if (delta > 0) {
                P_ = delta * P_;
            }
        }
    }
};
