#include <fstream>
#include <iostream>
#include <queue>

#include <Eigen/Core>

#include <geometry_msgs/PoseWithCovarianceStamped.h>
#include <nav_msgs/Path.h>
#include <ros/ros.h>

#include "db_smf_uwb/MyUWB.h"
#include "smf/db_smf_trans.hpp"

std::queue<db_smf_uwb::MyUWB> uwb_msg_queue;
void UWBMsgCallback(const db_smf_uwb::MyUWB& msg)
{
    uwb_msg_queue.push(msg);
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "db_smf_node");
    ros::NodeHandle nh;

    // subscribe sensor data
    std::string uwb_topic;
    nh.param<std::string>("uwb_topic", uwb_topic, "/");
    ros::Subscriber uwb_sub = nh.subscribe(uwb_topic, 2000, UWBMsgCallback);

    // publish
    ros::Publisher ellipse_pub = nh.advertise<geometry_msgs::PoseWithCovarianceStamped>("/ellipse_smf", 100);

    ros::Publisher path_pub = nh.advertise<nav_msgs::Path>("/est_path", 20);
    nav_msgs::Path path;
    path.header.frame_id = "world";

    // save settings
    std::string save_path;
    nh.param<std::string>("save_path", save_path, "/");

    const std::string center_save_path = save_path + "/est_center.txt";
    std::ofstream center_save_fs;
    center_save_fs.open(center_save_path, std::ios::out);
    center_save_fs.setf(std::ios::fixed, std::ios::floatfield);
    center_save_fs.precision(16);

    const std::string extrema_save_path = save_path + "/est_extrema.txt";
    std::ofstream extrema_save_fs;
    extrema_save_fs.open(extrema_save_path, std::ios::out);
    extrema_save_fs.setf(std::ios::fixed, std::ios::floatfield);
    extrema_save_fs.precision(16);

    /* parameters */
    // anchor positions
    int num_anchors = 0;
    nh.param<int>("num_anchors", num_anchors, 8);

    std::vector<double> aps_val_vec;
    nh.param<std::vector<double>>("anchor_positions", aps_val_vec, std::vector<double>());

    ROS_ASSERT_MSG(aps_val_vec.size() == (num_anchors * 3), "Please check the anchor configuration!");

    std::vector<Eigen::Vector3d> aps_vec;
    for (size_t i = 0; i < num_anchors; i++) {
        aps_vec.emplace_back(aps_val_vec[i * 3], aps_val_vec[i * 3 + 1], aps_val_vec[i * 3 + 2]);
    }

    // bias
    std::vector<double> bias_vec;
    nh.param<std::vector<double>>("bias", bias_vec, std::vector<double>());

    // NLOS
    double nlos_thresh = 0;
    nh.param<double>("nlos_thresh", nlos_thresh, 10000);

    // shape matrix
    double initial_robot_px_bound = 0;
    nh.param<double>("initial_robot_px_bound", initial_robot_px_bound, 1);
    double initial_robot_py_bound = 0;
    nh.param<double>("initial_robot_py_bound", initial_robot_py_bound, 1);
    double initial_robot_pz_bound = 0;
    nh.param<double>("initial_robot_pz_bound", initial_robot_pz_bound, 1);

    double initial_robot_v_bound = 0;
    nh.param<double>("initial_robot_v_bound", initial_robot_v_bound, 1);
    double initial_robot_a_bound = 0;
    nh.param<double>("initial_robot_a_bound", initial_robot_a_bound, 1);

    double motion_robot_jx_bound = 0;
    nh.param<double>("motion_robot_jx_bound", motion_robot_jx_bound, 1);
    double motion_robot_jy_bound = 0;
    nh.param<double>("motion_robot_jy_bound", motion_robot_jy_bound, 1);
    double motion_robot_jz_bound = 0;
    nh.param<double>("motion_robot_jz_bound", motion_robot_jz_bound, 1);

    // other
    double tag_height = 0;
    nh.param<double>("tag_height", tag_height, 0);

    double noise_bound = 0;
    nh.param<double>("noise_bound", noise_bound, 0);

    /// core
    DBSMFtrans smf(num_anchors, nlos_thresh, initial_robot_px_bound, initial_robot_py_bound, initial_robot_pz_bound, initial_robot_v_bound, initial_robot_a_bound);

    double last_stamp = 0;

    // const double DT_THRESH = 0.05; // 20Hz
    // const double DT_THRESH = 0.1; // 10Hz
    const double DT_THRESH = 0.2; // 5Hz
    // const double DT_THRESH = 0.5; // 2Hz

    bool init_done = false;
    const int INIT_Y_THRESH = 10;
    std::vector<Eigen::VectorXd> y_for_init_vec;

    while (ros::ok()) {
        ros::spinOnce();

        if (uwb_msg_queue.empty() || aps_vec.empty()) {
            continue;
        }

        // get data
        auto uwb_msg = uwb_msg_queue.front();
        uwb_msg_queue.pop();

        auto dis_arr = uwb_msg.dis_arr;
        Eigen::VectorXd y;
        y.setZero(num_anchors, 1);
        for (size_t i = 0; i < num_anchors; i++) {
            y(i) = dis_arr.at(i);
        }

        // initialization
        if (!init_done) {
            y_for_init_vec.push_back(y);
            if (y_for_init_vec.size() < INIT_Y_THRESH) {
                continue;
            }

            Eigen::VectorXd y_mean;
            y_mean.setZero(y.rows(), y.cols());
            for (auto& y_i : y_for_init_vec) {
                y_mean += y_i;
            }
            y_mean /= INIT_Y_THRESH;

            smf.initTagPosition(y_mean, tag_height, aps_vec);

            init_done = true;
        }

        // sampling period constraints
        double stamp = uwb_msg.stamp;

        if (last_stamp == 0) {
            last_stamp = stamp;
        }
        if (stamp - last_stamp < DT_THRESH) {
            continue;
        }

        double dt = stamp - last_stamp;
        last_stamp = stamp;

        // bias
        for (int i = 0; i < num_anchors; i++) {
            y(i) /= (1 + bias_vec.at(i));
        }

        smf.predict(dt, motion_robot_jx_bound, motion_robot_jy_bound, motion_robot_jz_bound);
        smf.update_seq(y, noise_bound, aps_vec);

        /* publish */
        // bound
        geometry_msgs::PoseWithCovarianceStamped ellipse_msg;
        ellipse_msg.header.stamp = ros::Time::now();
        ellipse_msg.header.frame_id = "world";
        ellipse_msg.pose.pose.position.x = smf.x_(0);
        ellipse_msg.pose.pose.position.y = smf.x_(3);
        ellipse_msg.pose.pose.position.z = smf.x_(6);
        ellipse_msg.pose.pose.orientation.x = 0;
        ellipse_msg.pose.pose.orientation.y = 0;
        ellipse_msg.pose.pose.orientation.z = 0;
        ellipse_msg.pose.pose.orientation.w = 1;
        ellipse_msg.pose.covariance[0] = smf.P_(0, 0);
        ellipse_msg.pose.covariance[1] = ellipse_msg.pose.covariance[6] = smf.P_(0, 3);
        ellipse_msg.pose.covariance[2] = ellipse_msg.pose.covariance[12] = smf.P_(0, 6);
        ellipse_msg.pose.covariance[7] = smf.P_(3, 3);
        ellipse_msg.pose.covariance[8] = ellipse_msg.pose.covariance[13] = smf.P_(3, 6);
        ellipse_msg.pose.covariance[14] = smf.P_(6, 6);
        ellipse_pub.publish(ellipse_msg);

        // path
        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header.stamp = ros::Time::now();
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose.position.x = smf.x_(0);
        pose_stamped.pose.position.y = smf.x_(3);
        pose_stamped.pose.position.z = smf.x_(6);

        path.header.stamp = ros::Time(ros::Time::now());
        path.poses.push_back(pose_stamped);
        path_pub.publish(path);

        /* record */
        // center
        std::stringstream ss;
        ss.setf(std::ios::fixed, std::ios::floatfield);
        ss.precision(9);
        ss << uwb_msg.stamp << " ";
        ss.precision(16);
        ss << smf.x_(0) << " " << smf.x_(3) << " " << smf.x_(6) << " 0 0 0 1\n";

        center_save_fs << ss.str();

        // shape matrix
        ss.str("");
        ss.precision(9);
        ss << uwb_msg.stamp << " ";
        ss.precision(16);
        ss << smf.P_(0, 0) << " " << smf.P_(0, 3) << " " << smf.P_(0, 6) << " ";
        ss << smf.P_(3, 0) << " " << smf.P_(3, 3) << " " << smf.P_(3, 6) << " ";
        ss << smf.P_(6, 0) << " " << smf.P_(6, 3) << " " << smf.P_(6, 6) << "\n";

        extrema_save_fs << ss.str();

        // time
        ss.str("");
        ss.precision(9);
        ss << uwb_msg.stamp << " ";
        ss.precision(16);
    }

    center_save_fs.close();
    extrema_save_fs.close();

    return 0;
}