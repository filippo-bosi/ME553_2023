#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include "Eigen/Geometry"

// Define a struct to represent a joint
struct Joint {

    enum class Type {
        FIXED,
        FLOATING,
        REVOLUTE,
        PRISMATIC
    }type;

    // robot definition
    Eigen::Vector3d jointPosition_B;
    Eigen::Vector3d jointAxis_B;
    double jointGC;
    double jointGV;
    Eigen::Vector3d jointPositionCOM_B;
    double mass;
    Eigen::Matrix3d inertiaMatrix_B;
    Eigen::Matrix3d jointRotation_B;

    std::vector<Joint*> children; // Pointers to child joints
    std::vector<Joint*> parent; // Pointers to parent joints

    // variables in the World Frame
    Eigen::Vector3d jointPosition_W;
    Eigen::Vector3d jointPositionCOM_W;
    Eigen::Vector3d jointLinearVelocity;
    Eigen::Vector3d jointAngularVelocity;
    Eigen::Vector3d jointLinearAcceleration;
    Eigen::Vector3d jointAngularAcceleration;
    Eigen::Matrix3d inertiaMatrix_W;
    Eigen::Vector3d jointPositionCOM_comp;
    double mass_comp;
    Eigen::Matrix3d inertiaMatrix_comp;
    Eigen::Vector3d S_matrix_temp;
    Eigen::VectorXd S_matrix_W = Eigen::VectorXd::Zero(6,1);
    Eigen::Vector3d jointPosToEnd;
    Eigen::VectorXd generalizedForce = Eigen::VectorXd::Zero(6,1);
};

// Define a class to represent the tree structure
class JointTree {
public:
    JointTree(Joint* root) : root(root) {}
    bool trunk_reached = false;
    int i = 0;
    Eigen::Vector3d jointPosition_W;
    Eigen::Vector3d jointPositionCOM_W;
    Eigen::Vector3d jointLinearVelocity;
    Eigen::Vector3d jointAngularVelocity;
    Eigen::Matrix3d inertiaMatrix_W;

    // Compute the position of a joint in Cartesian coordinates
    void computeKinematicsDownTheTree(Joint* joint, std::vector<Eigen::Matrix3d>& jointsRotationMatrices, int& counter) {
        if (joint->type == Joint::Type::FLOATING)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED) {
                joint->jointPositionCOM_B = joint->jointPosition_B;
                joint->jointRotation_B = Eigen::Matrix3d::Identity(3,3);
                Eigen::Matrix3d multipliedRotations;
                multipliedRotations = Eigen::Matrix3d::Identity();
                for (const auto &e: jointsRotationMatrices) {
                    multipliedRotations *= e;
                }
                jointPosition_W = jointPosition_W + multipliedRotations *  joint->jointPosition_B;
                joint->jointPosition_W = jointPosition_W;
                joint->jointPositionCOM_W = jointPosition_W;
                joint->inertiaMatrix_W = multipliedRotations*joint->inertiaMatrix_B*multipliedRotations.transpose();
                // Initialization for composite body
                joint->jointPositionCOM_comp = joint->jointPositionCOM_W;
                joint->mass_comp = joint->mass;
                joint->inertiaMatrix_comp = joint->inertiaMatrix_W;
                joint->jointAngularVelocity = jointAngularVelocity;
                joint->S_matrix_temp << Eigen::VectorXd::Zero(3, 1);
                joint->S_matrix_W << Eigen::VectorXd::Zero(6, 1);
                joint->jointPosToEnd = multipliedRotations*joint->jointPosition_B;
                joint->generalizedForce << 0,0,0,0,0,0;
            } else
            if (joint->type == Joint::Type::FLOATING) {
                joint->jointPosition_W = joint->jointPosition_B;
                //joint->S_matrix = Eigen::MatrixXd::Identity(6,6);
                jointPosition_W = joint->jointPosition_W;
                joint->jointPositionCOM_W = joint->jointPosition_W + joint->jointRotation_B*joint->jointPositionCOM_B;
                joint->inertiaMatrix_W = joint->jointRotation_B*joint->inertiaMatrix_B*joint->jointRotation_B.transpose();
                jointAngularVelocity = joint->jointAngularVelocity;
                joint->jointPosToEnd = joint->jointPosition_B;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                Eigen::AngleAxisd angAx;
                angAx.angle() = joint->jointGC; angAx.axis() = joint->jointAxis_B;
                joint->jointRotation_B = angAx.toRotationMatrix();
                Eigen::Matrix3d multipliedRotations;
                multipliedRotations = Eigen::Matrix3d::Identity();
                for (const auto &e: jointsRotationMatrices) {
                    multipliedRotations *= e;
                }
                jointPosition_W = jointPosition_W + multipliedRotations * joint->jointPosition_B;
                joint->jointPosition_W = jointPosition_W;
                joint->jointPositionCOM_W = jointPosition_W + multipliedRotations * joint->jointRotation_B * joint->jointPositionCOM_B;
                inertiaMatrix_W = multipliedRotations * joint->jointRotation_B * joint->inertiaMatrix_B *(multipliedRotations * joint->jointRotation_B).transpose();
                joint->inertiaMatrix_W = inertiaMatrix_W;
                joint->S_matrix_temp = multipliedRotations * joint->jointAxis_B;
                joint->S_matrix_W << 0, 0, 0, joint->S_matrix_temp;
                jointAngularVelocity = jointAngularVelocity + multipliedRotations*joint->jointRotation_B*joint->jointGV*joint->jointAxis_B; ;
                joint->jointAngularVelocity = jointAngularVelocity;
                joint->jointPosToEnd = multipliedRotations*joint->jointPosition_B;
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                joint->jointRotation_B = Eigen::Matrix3d::Identity();
                Eigen::Matrix3d multipliedRotations;
                multipliedRotations = Eigen::Matrix3d::Identity();
                for (const auto &e: jointsRotationMatrices) {
                    multipliedRotations *= e;
                }
                jointPosition_W = jointPosition_W + multipliedRotations*(joint->jointGC * joint->jointAxis_B + joint->jointPosition_B);
                joint->jointPosition_W = jointPosition_W;
                joint->jointPositionCOM_W = jointPosition_W + multipliedRotations * joint->jointPositionCOM_B;
                inertiaMatrix_W = multipliedRotations * joint->inertiaMatrix_B * multipliedRotations.transpose();
                joint->inertiaMatrix_W = inertiaMatrix_W;
                joint->S_matrix_temp = multipliedRotations * joint->jointAxis_B;
                joint->S_matrix_W << joint->S_matrix_temp, 0, 0, 0;
                jointAngularVelocity = jointAngularVelocity;
                joint->jointAngularVelocity = jointAngularVelocity;
                joint->jointPosToEnd = multipliedRotations*joint->jointPosition_B;
            }
            for (Joint *child: joint->children) {
                // In case we are not interested in the position of the foot -> stop before
                if (i == counter)
                    break;
                jointsRotationMatrices.push_back(joint->jointRotation_B);
                i++;
                computeKinematicsDownTheTree(child, jointsRotationMatrices, counter);
            }
        }
        if (!trunk_reached) {
            for (Joint *parent: joint->parent) {
                counter++;
                computeKinematicsDownTheTree(parent, jointsRotationMatrices, counter);
            }
        }
    }

    Eigen::Vector3d parentAngVelocity;
    //std::vector<Eigen::Vector3d> childJointPosToEnd;
    std::vector<Eigen::Matrix3d> jointsRotationMatrices2;


    void computeJointLinearVelocity(Joint* joint, int& counter) {
        if (joint->type == Joint::Type::FLOATING)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED) {
                jointLinearVelocity = jointLinearVelocity + parentAngVelocity.cross(joint->jointPosToEnd);
                joint->jointLinearVelocity = jointLinearVelocity;
            } else
            if (joint->type == Joint::Type::FLOATING) {
                jointLinearVelocity = joint->jointLinearVelocity;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                jointLinearVelocity = jointLinearVelocity + parentAngVelocity.cross(joint->jointPosToEnd);
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                Eigen::Matrix3d multipliedRotations;
                multipliedRotations = Eigen::Matrix3d::Identity();
                for (const auto &e: jointsRotationMatrices2) {
                    multipliedRotations *= e;
                }
                jointLinearVelocity = jointLinearVelocity + multipliedRotations*joint->jointGV*joint->jointAxis_B;
            }

            for (Joint *child: joint->children) {
                // In case we are not interested in the position of the foot -> stop before
                if (i == counter) {
                    joint->jointLinearVelocity = jointLinearVelocity;
                    break;
                }
                parentAngVelocity = joint->jointAngularVelocity;
                jointsRotationMatrices2.push_back(joint->jointRotation_B);
                i++;
                computeJointLinearVelocity(child,counter);
            }
        }
        if (!trunk_reached) {
            //childJointPosToEnd[counter] = joint->jointPosToEnd;
            for (Joint *parent: joint->parent) {
                counter++;
                parent->jointPosToEnd = parent->jointPosToEnd + joint->jointPosToEnd;
                computeJointLinearVelocity(parent, counter);
            }
        }
    }
    Eigen::Vector3d parentAngularAcceleration;
    Eigen::Vector3d parentLinearAcceleration;
    Eigen::Vector3d jointLinearAcceleration;
    Eigen::Vector3d parentJointPosition;


    void computeAcceleration(Joint* joint, int& counter) {
        if (joint->type == Joint::Type::FLOATING)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED) {
                Eigen::Matrix3d alfa_x, omegaParent_x, omega_x;
                alfa_x << 0, -parentAngularAcceleration(2), parentAngularAcceleration(1),
                        parentAngularAcceleration(2), 0, -parentAngularAcceleration(0),
                        -parentAngularAcceleration(1), parentAngularAcceleration(0), 0;
                omegaParent_x << 0, -parentAngVelocity(2), parentAngVelocity(1),
                        parentAngVelocity(2), 0, -parentAngVelocity(0),
                        -parentAngVelocity(1), parentAngVelocity(0), 0;
                omega_x << 0, -joint->jointAngularVelocity(2), joint->jointAngularVelocity(1),
                        joint->jointAngularVelocity(2), 0, -joint->jointAngularVelocity(0),
                        -joint->jointAngularVelocity(1), joint->jointAngularVelocity(0), 0;
                Eigen::Vector3d temp_lin_acc;
                Eigen::VectorXd temp_acc, acc;
                Eigen::MatrixXd S_dot;
                temp_acc.setZero(6,1);
                acc.setZero(6,1);
                S_dot.setZero(6,1);
                temp_lin_acc = parentLinearAcceleration + alfa_x*(joint->jointPosition_W-parentJointPosition) +
                               omegaParent_x*omegaParent_x*(joint->jointPosition_W-parentJointPosition);
                temp_acc << temp_lin_acc, parentAngularAcceleration;
                S_dot << 0, 0, 0, omega_x*joint->S_matrix_temp;
                acc << temp_acc+S_dot*joint->jointGV;
                joint->jointLinearAcceleration << acc(0), acc(1), acc(2);
                joint->jointAngularAcceleration << acc(3), acc(4), acc(5);
            } else
            if (joint->type == Joint::Type::FLOATING) {
                jointLinearAcceleration = joint->jointLinearAcceleration;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                Eigen::Matrix3d alfa_x, omegaParent_x, omega_x;
                alfa_x << 0, -parentAngularAcceleration(2), parentAngularAcceleration(1),
                        parentAngularAcceleration(2), 0, -parentAngularAcceleration(0),
                        -parentAngularAcceleration(1), parentAngularAcceleration(0), 0;
                omegaParent_x << 0, -parentAngVelocity(2), parentAngVelocity(1),
                        parentAngVelocity(2), 0, -parentAngVelocity(0),
                        -parentAngVelocity(1), parentAngVelocity(0), 0;
                omega_x << 0, -joint->jointAngularVelocity(2), joint->jointAngularVelocity(1),
                        joint->jointAngularVelocity(2), 0, -joint->jointAngularVelocity(0),
                        -joint->jointAngularVelocity(1), joint->jointAngularVelocity(0), 0;
                Eigen::Vector3d temp_lin_acc;
                Eigen::VectorXd temp_acc, acc;
                Eigen::MatrixXd S_dot;
                temp_acc.setZero(6,1);
                acc.setZero(6,1);
                S_dot.setZero(6,1);
                temp_lin_acc = parentLinearAcceleration + alfa_x*(joint->jointPosition_W-parentJointPosition) +
                        omegaParent_x*omegaParent_x*(joint->jointPosition_W-parentJointPosition);
                temp_acc << temp_lin_acc, parentAngularAcceleration;
                S_dot << 0, 0, 0, omegaParent_x*joint->S_matrix_temp;
                acc << temp_acc+S_dot*joint->jointGV;
                joint->jointLinearAcceleration << acc(0), acc(1), acc(2);
                joint->jointAngularAcceleration << acc(3), acc(4), acc(5);
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                Eigen::Matrix3d alfa_x, omegaParent_x, omega_x;
                alfa_x << 0, -parentAngularAcceleration(2), parentAngularAcceleration(1),
                        parentAngularAcceleration(2), 0, -parentAngularAcceleration(0),
                        -parentAngularAcceleration(1), parentAngularAcceleration(0), 0;
                omegaParent_x << 0, -parentAngVelocity(2), parentAngVelocity(1),
                        parentAngVelocity(2), 0, -parentAngVelocity(0),
                        -parentAngVelocity(1), parentAngVelocity(0), 0;
                omega_x << 0, -joint->jointAngularVelocity(2), joint->jointAngularVelocity(1),
                        joint->jointAngularVelocity(2), 0, -joint->jointAngularVelocity(0),
                        -joint->jointAngularVelocity(1), joint->jointAngularVelocity(0), 0;
                Eigen::Vector3d temp_lin_acc;
                Eigen::VectorXd temp_acc, acc;
                Eigen::MatrixXd S_dot;
                temp_acc.setZero(6,1);
                acc.setZero(6,1);
                S_dot.setZero(6,1);
                temp_lin_acc = parentLinearAcceleration + alfa_x*(joint->jointPosition_W-parentJointPosition) +
                               omegaParent_x*omegaParent_x*(joint->jointPosition_W-parentJointPosition);
                temp_acc << temp_lin_acc, parentAngularAcceleration;
                S_dot << omegaParent_x*joint->S_matrix_temp, 0, 0, 0;
                acc << temp_acc+2*S_dot*joint->jointGV;
                joint->jointLinearAcceleration << acc.head(3);
                joint->jointAngularAcceleration << acc(3), acc(4), acc(5);
            }

            for (Joint *child: joint->children) {
                // In case we are not interested in the position of the foot -> stop before
                if (i == counter) {
                    break;
                }
                parentLinearAcceleration = joint->jointLinearAcceleration;
                parentAngularAcceleration = joint->jointAngularAcceleration;
                parentAngVelocity = joint->jointAngularVelocity;
                parentJointPosition = joint->jointPosition_W;
                i++;
                computeAcceleration(child,counter);
            }
        }
        if (!trunk_reached) {
            for (Joint *parent: joint->parent) {
                counter++;
                computeAcceleration(parent, counter);
            }
        }
    }

    void getCompositeBody(Joint* joint1, Joint* joint2) {
        Eigen::Matrix3d Inertia_combined;
        Eigen::Vector3d r_COM, r1, r2;
        Eigen::Matrix3d r1_x, r2_x;
        double m;

        m = joint1->mass + joint2->mass;
        r_COM = 1/m*(joint1->mass*joint1->jointPositionCOM_W + joint2->mass*joint2->jointPositionCOM_W);
        r1 = joint1->jointPositionCOM_W - r_COM;
        r2 = joint2->jointPositionCOM_W - r_COM;
        r1_x << 0, -r1(2), r1(1),
                r1(2), 0, -r1(0),
                -r1(1), r1(0), 0;
        r2_x << 0, -r2(2), r2(1),
                r2(2), 0, -r2(0),
                -r2(1), r2(0), 0;
        Inertia_combined = joint1->inertiaMatrix_W+joint2->inertiaMatrix_W-joint1->mass*r1_x*r1_x-joint2->mass*r2_x*r2_x;
//        joint1->mass_comp = m;
//        joint1->jointPositionCOM_comp = r_COM;
//        joint1->inertiaMatrix_comp = Inertia_combined;
        joint1->mass = m;
        joint1->jointPositionCOM_W = r_COM;
        joint1->inertiaMatrix_W = Inertia_combined;
    }

    inline Eigen::MatrixXd getSpatialInertiaMatrix(Joint* joint_j) {
        Eigen::MatrixXd spatialInertiaMatrix(6,6);
        Eigen::Vector3d r_COM;
        Eigen::Matrix3d rCOM_x;
        double m;
        m = joint_j->mass;
        r_COM = joint_j->jointPositionCOM_W - joint_j->jointPosition_W;
        rCOM_x << 0, -r_COM(2), r_COM(1),
                r_COM(2), 0, -r_COM(0),
                -r_COM(1), r_COM(0), 0;
        spatialInertiaMatrix << m*Eigen::Matrix3d::Identity(3,3), -m*rCOM_x,
                m*rCOM_x,                 joint_j->inertiaMatrix_W-m*rCOM_x*rCOM_x;
        return spatialInertiaMatrix;
    }

    inline Eigen::MatrixXd getBias(Joint* joint_j) {
        Eigen::VectorXd bias(6,1);
        Eigen::Vector3d omega, r_COM;
        Eigen::Matrix3d omega_x, rCOM_x;
        double m;
        m = joint_j->mass;
        r_COM = joint_j->jointPositionCOM_W - joint_j->jointPosition_W;
        rCOM_x << 0, -r_COM(2), r_COM(1),
                r_COM(2), 0, -r_COM(0),
                -r_COM(1), r_COM(0), 0;
        omega = joint_j->jointAngularVelocity;
        omega_x << 0, -omega(2), omega(1),
                omega(2), 0, -omega(0),
                -omega(1), omega(0), 0;

        bias << m*omega_x*omega_x*r_COM, omega_x*(joint_j->inertiaMatrix_W-m*rCOM_x*rCOM_x)*omega;
        return bias;
    }

    void getGeneralizedForces(Joint* joint) {
        Eigen::VectorXd acceleration(6,1), bias(6,1), childMatrix(6,1);
        Eigen::MatrixXd spatialInertia(6,6);
        Eigen::Vector3d force_child, torque_child;
        force_child << joint->generalizedForce.head(3);
        torque_child << joint->generalizedForce[3], joint->generalizedForce[4], joint->generalizedForce[5];

        for (Joint *parent: joint->parent) {
            if (parent->type == Joint::Type::FLOATING)
                break;
            childMatrix << force_child, torque_child+(joint->jointPosition_W-parent->jointPosition_W).cross(force_child);
            acceleration << parent->jointLinearAcceleration, parent->jointAngularAcceleration;
            spatialInertia = getSpatialInertiaMatrix(parent);
            bias = getBias(parent);
            parent->generalizedForce = spatialInertia * acceleration + bias + childMatrix;
            getGeneralizedForces(parent);
        }
    }

private:
    Joint* root;
};

inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    double m_trunk, m_imu, m_hip, m_thigh, m_calf, m_foot;
    m_trunk = 9.041;
    m_imu = 0.001;
    m_hip = 1.993;
    m_thigh = 0.639;
    m_calf = 0.207;
    m_foot = 0.06;

    double ixx, ixy, ixz, iyy, iyz, izz;
    Eigen::Matrix3d I_trunk, I_imu, I_FR_hip, I_FL_hip, I_RR_hip, I_RL_hip,
            I_FR_thigh, I_FL_thigh, I_RR_thigh, I_RL_thigh,
            I_FR_calf, I_FL_calf, I_RR_calf, I_RL_calf,
            I_FR_foot, I_FL_foot, I_RR_foot, I_RL_foot;
    // inertia matrices of bodies from URDF
    ixx=0.033260231; ixy=-0.000451628; ixz=0.000487603; iyy=0.16117211; iyz=4.8356e-05; izz=0.17460442;
    I_trunk << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    ixx=0.0001; ixy=0; ixz=0; iyy=0.000001; iyz=0; izz=0.0001;
    I_imu << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    ixx=0.002903894; ixy=7.185e-05; ixz=-1.262e-06; iyy=0.004907517; iyz=1.75e-06; izz=0.005586944;
    I_FR_hip << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    I_FL_hip << ixx, -ixy, ixz, -ixy, iyy, -iyz, ixz, -iyz, izz;
    I_RR_hip << ixx, -ixy, -ixz, -ixy, iyy, iyz, -ixz, iyz, izz;
    I_RL_hip << ixx, ixy, -ixz, ixy, iyy, -iyz, -ixz, -iyz, izz;
    ixx=0.005666803; ixy=-3.597e-06; ixz=0.000491446; iyy=0.005847229; iyz=-1.0086e-05; izz=0.000369811;
    I_FR_thigh << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    I_FL_thigh << ixx, -ixy, ixz, -ixy, iyy, -iyz, ixz, -iyz, izz;
    I_RR_thigh << I_FR_thigh;
    I_RL_thigh << I_FL_thigh;
    ixx=0.006341369; ixy=-3e-09; ixz=-8.7951e-05; iyy=0.006355157; iyz=-1.336e-06; izz=3.9188e-05;
    I_FR_calf << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    I_FL_calf << I_FR_calf;
    I_RR_calf << I_FR_calf;
    I_RL_calf << I_FR_calf;
    ixx=1.6854e-05; ixy=0.0; ixz=0.0; iyy=1.6854e-05; iyz=0.0; izz=1.6854e-05;
    I_FR_foot << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    I_FL_foot << I_FR_foot;
    I_RR_foot << I_FR_foot;
    I_RL_foot << I_FR_foot;

    // Create the joints for the quadruped robot

    // Front-Right Leg
    Joint* FR_hip = new Joint{Joint::Type::REVOLUTE, {0.2399, -0.051, 0},{1,0,0},
                              gc[7], gv[6], {-0.022191, -0.015144, -1.5e-05}, m_hip, I_FR_hip};
    Joint* FR_thigh = new Joint{Joint::Type::REVOLUTE, {0, -0.083, 0},{0,1,0},
                                gc[8], gv[7], {-0.005607, 0.003877, -0.048199}, m_thigh, I_FR_thigh};
    Joint* FR_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[9], gv[8], {0.002781, 6.3e-05, -0.142518}, m_calf, I_FR_calf};
    Joint* FR_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {},{},{}, m_foot, I_FR_foot, Eigen::Matrix3d::Identity()};

    FR_hip->children.push_back(FR_thigh);
    FR_thigh->parent.push_back(FR_hip);
    FR_thigh->children.push_back(FR_calf);
    FR_calf->parent.push_back(FR_thigh);
    FR_calf->children.push_back(FR_foot);
    FR_foot->parent.push_back(FR_calf);

    // Front-Left Leg
    Joint* FL_hip = new Joint{Joint::Type::REVOLUTE, {0.2399, 0.051, 0},{1,0,0},
                              gc[10],gv[9],{-0.022191, 0.015144, -1.5e-05}, m_hip, I_FL_hip};
    Joint* FL_thigh = new Joint{Joint::Type::PRISMATIC, {0, 0.083, 0},{0,1,0},
                                gc[11], gv[10],{-0.005607, -0.003877, -0.048199}, m_thigh, I_FL_thigh};
    Joint* FL_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[12],gv[11], {0.002781, 6.3e-05, -0.142518}, m_calf, I_FL_calf};
    Joint* FL_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {},{}, m_foot, I_FL_foot, Eigen::Matrix3d::Identity()};
    FL_hip->children.push_back(FL_thigh);
    FL_thigh->parent.push_back(FL_hip);
    FL_thigh->children.push_back(FL_calf);
    FL_calf->parent.push_back(FL_thigh);
    FL_calf->children.push_back(FL_foot);
    FL_foot->parent.push_back(FL_calf);

    // Rear-Right Leg
    Joint* RR_hip = new Joint{Joint::Type::REVOLUTE, {-0.2399, -0.051, 0},{1,0,0},
                              gc[13], gv[12], {0.022191, -0.015144, -1.5e-05}, m_hip, I_RR_hip};
    Joint* RR_thigh = new Joint{Joint::Type::REVOLUTE, {0, -0.083, 0},{0,1,0},
                                gc[14], gv[13],{-0.005607, 0.003877, -0.048199}, m_thigh, I_RR_thigh};
    Joint* RR_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[15], gv[14],{0.002781, 6.3e-05, -0.142518}, m_calf, I_RR_calf};
    Joint* RR_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {},{}, m_foot, I_RR_foot, Eigen::Matrix3d::Identity()};
    RR_hip->children.push_back(RR_thigh);
    RR_thigh->parent.push_back(RR_hip);
    RR_thigh->children.push_back(RR_calf);
    RR_calf->parent.push_back(RR_thigh);
    RR_calf->children.push_back(RR_foot);
    RR_foot->parent.push_back(RR_calf);

    // Rear-Left Leg
    Joint* RL_hip = new Joint{Joint::Type::PRISMATIC, {-0.2399, 0.051, 0},{1,0,0},
                              gc[16], gv[15],{0.022191, 0.015144, -1.5e-05}, m_hip, I_RL_hip};
    Joint* RL_thigh = new Joint{Joint::Type::PRISMATIC, {0, 0.083, 0},{0,1,0},
                                gc[17], gv[16],{-0.005607, -0.003877, -0.048199}, m_thigh, I_RL_thigh};
    Joint* RL_calf = new Joint{Joint::Type::PRISMATIC, {0, 0, -0.25},{0,1,0},
                               gc[18], gv[17],{0.002781, 6.3e-05, -0.142518}, m_calf, I_RL_calf};
    Joint* RL_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {}, {},m_foot, I_RL_foot, Eigen::Matrix3d::Identity()};
    RL_hip->children.push_back(RL_thigh);
    RL_thigh->parent.push_back(RL_hip);
    RL_thigh->children.push_back(RL_calf);
    RL_calf->parent.push_back(RL_thigh);
    RL_calf->children.push_back(RL_foot);
    RL_foot->parent.push_back(RL_calf);

    // Create the tree structure for the quadruped robot
    Joint* imu = new Joint{Joint::Type::FIXED, {0,0,0},{0,0,0},
                           0, 0, {0,0,0}, m_imu, I_imu, Eigen::Matrix3d::Identity()};

    // convert quaternion to rotation matrix
    Eigen::Quaterniond q;
    q.w() = gc[3];
    q.x() = gc[4];
    q.y() = gc[5];
    q.z() = gc[6];
    //q.toRotationMatrix();

    Joint* trunk = new Joint{Joint::Type::FLOATING, gc.head(3),{0,0,0},
                             0, 0, {0.008465, 0.004045, -0.000763}, m_trunk, I_trunk,
                             q.toRotationMatrix(),{}};
    trunk->jointLinearVelocity << gv[0], gv[1], gv[2];
    trunk->jointAngularVelocity << gv[3], gv[4], gv[5];
    trunk->jointLinearAcceleration << 0, 0, 9.81;
    trunk->jointAngularAcceleration << 0, 0, 0;

    ///// Compute FWD Kinematics, COM and Inertia for each joint //////////////
    // Idea: attach first leg, compute kinematics - COM - Inertia expressed in
    // the World frame, detach leg, attach next leg and repeat until the direct
    // kinematics has been computed for every joint

    // Variables declaration
    std::vector<Eigen::Matrix3d> jointsRotationMatrices;
    int counter = 0;

    // attach FR leg
    trunk->children.push_back(FR_hip);
    FR_hip->parent.push_back(trunk);
    JointTree tree_FR{trunk};
    tree_FR.computeKinematicsDownTheTree(FR_foot,jointsRotationMatrices, counter);
    counter = 0;
    tree_FR.trunk_reached = false;
    tree_FR.i = 0;
    tree_FR.computeAcceleration(FR_foot, counter);
    counter = 0;
    tree_FR.trunk_reached = false;
    tree_FR.i = 0;
    // detach FR leg
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    // attach FL leg
    trunk->children.push_back(FL_hip);
    FL_hip->parent.push_back(trunk);
    JointTree tree_FL{trunk};
    tree_FL.computeKinematicsDownTheTree(FL_foot, jointsRotationMatrices,counter);
    counter = 0;
    tree_FL.trunk_reached = false;
    tree_FL.i = 0;
    tree_FL.computeAcceleration(FL_foot, counter);
    counter = 0;
    tree_FL.trunk_reached = false;
    tree_FL.i = 0;
    // detach FL leg
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    // attach RR leg
    trunk->children.push_back(RR_hip);
    RR_hip->parent.push_back(trunk);
    JointTree tree_RR{trunk};
    tree_RR.computeKinematicsDownTheTree(RR_foot,jointsRotationMatrices, counter);
    counter = 0;
    tree_RR.trunk_reached = false;
    tree_RR.i = 0;
    tree_RR.computeAcceleration(RR_foot, counter);
    counter = 0;
    tree_RR.trunk_reached = false;
    tree_RR.i = 0;
    // detach RR leg
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    // attach RL leg
    trunk->children.push_back(RL_hip);
    RL_hip->parent.push_back(trunk);
    JointTree tree_RL{trunk};
    tree_RL.computeKinematicsDownTheTree(RL_foot,jointsRotationMatrices, counter);
    counter = 0;
    tree_RL.trunk_reached = false;
    tree_RL.i = 0;
    tree_RL.computeAcceleration(RL_foot, counter);
    counter = 0;
    tree_RL.trunk_reached = false;
    tree_RL.i = 0;

    // detach RL leg
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    // Combine Foot and Calf to get the right mass, inertia and COM
    tree_FR.getCompositeBody(FR_calf, FR_foot);
    tree_FL.getCompositeBody(FL_calf, FL_foot);
    tree_RR.getCompositeBody(RR_calf, RR_foot);
    tree_RL.getCompositeBody(RL_calf, RL_foot);

    trunk->children.push_back(imu);
    imu->parent.push_back(trunk);
    JointTree tree_imu{trunk};
    tree_imu.computeKinematicsDownTheTree(imu, jointsRotationMatrices, counter);
    counter = 0;
    tree_FR.trunk_reached = false;
    tree_FR.i = 0;
    tree_imu.getCompositeBody(trunk, imu);
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    Eigen::VectorXd generalizedForces;
    tree_FR.getGeneralizedForces(FR_foot);
    tree_FL.getGeneralizedForces(FL_foot);
    tree_RR.getGeneralizedForces(RR_foot);
    tree_RL.getGeneralizedForces(RL_foot);

    Eigen::VectorXd acceleration(6,1), bias(6,1), childMatrix_FR(6,1), childMatrix_FL(6,1), childMatrix_RR(6,1), childMatrix_RL(6,1);
    Eigen::MatrixXd spatialInertia(6,6);
    Eigen::Vector3d force_child, torque_child;
    force_child << FR_hip->generalizedForce.head(3);
    torque_child << FR_hip->generalizedForce[3], FR_hip->generalizedForce[4], FR_hip->generalizedForce[5];
    childMatrix_FR << force_child, torque_child+(FR_hip->jointPosition_W-trunk->jointPosition_W).cross(force_child);

    force_child << FL_hip->generalizedForce.head(3);
    torque_child << FL_hip->generalizedForce[3], FL_hip->generalizedForce[4], FL_hip->generalizedForce[5];
    childMatrix_FL << force_child, torque_child+(FL_hip->jointPosition_W-trunk->jointPosition_W).cross(force_child);

    force_child << RR_hip->generalizedForce.head(3);
    torque_child << RR_hip->generalizedForce[3], RR_hip->generalizedForce[4], RR_hip->generalizedForce[5];
    childMatrix_RR << force_child, torque_child+(RR_hip->jointPosition_W-trunk->jointPosition_W).cross(force_child);

    force_child << RL_hip->generalizedForce.head(3);
    torque_child << RL_hip->generalizedForce[3], RL_hip->generalizedForce[4], RL_hip->generalizedForce[5];
    childMatrix_RL << force_child, torque_child+(RL_hip->jointPosition_W-trunk->jointPosition_W).cross(force_child);

    acceleration << trunk->jointLinearAcceleration, trunk->jointAngularAcceleration;
    spatialInertia = tree_imu.getSpatialInertiaMatrix(trunk);
    bias = tree_imu.getBias(trunk);

    trunk->generalizedForce = spatialInertia * acceleration + bias
            + childMatrix_FR + childMatrix_FL + childMatrix_RR + childMatrix_RL;

    Eigen::VectorXd biasFinal;
    biasFinal.setZero(18,1);
    biasFinal << Eigen::MatrixXd::Identity(6,6)*trunk->generalizedForce, FR_hip->S_matrix_W.transpose()*FR_hip->generalizedForce,
            FR_thigh->S_matrix_W.transpose()*FR_thigh->generalizedForce, FR_calf->S_matrix_W.transpose()*FR_calf->generalizedForce,
            FL_hip->S_matrix_W.transpose()*FL_hip->generalizedForce, FL_thigh->S_matrix_W.transpose()*FL_thigh->generalizedForce,
            FL_calf->S_matrix_W.transpose()*FL_calf->generalizedForce, RR_hip->S_matrix_W.transpose()*RR_hip->generalizedForce,
            RR_thigh->S_matrix_W.transpose()*RR_thigh->generalizedForce, RR_calf->S_matrix_W.transpose()*RR_calf->generalizedForce,
            RL_hip->S_matrix_W.transpose()*RL_hip->generalizedForce, RL_thigh->S_matrix_W.transpose()*RL_thigh->generalizedForce,
            RL_calf->S_matrix_W.transpose()*RL_calf->generalizedForce;

//    std::cout << "My nonlinearities: " << std::endl;
//    std::cout << biasFinal.transpose() << std::endl;

//    std::cout << trunk->generalizedForce << std::endl;
//    std::cout << FR_hip->generalizedForce << std::endl;
//    std::cout << FR_thigh->generalizedForce << std::endl;
//    std::cout << FR_calf->generalizedForce << std::endl<<std::endl;

    return biasFinal.transpose();
}