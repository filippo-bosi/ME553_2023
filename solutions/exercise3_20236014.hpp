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
    double jointGF;
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
    Eigen::Vector3d jointAngularVelocity_W;
    Eigen::Vector3d jointLinearAcceleration;
    Eigen::Vector3d jointAngularAcceleration;
    Eigen::Matrix3d inertiaMatrix_W;
    Eigen::Vector3d jointPositionCOM_comp;
    double mass_comp;
    Eigen::Matrix3d inertiaMatrix_comp;
    Eigen::Vector3d S_matrix_temp;
    Eigen::VectorXd S_matrix_W = Eigen::VectorXd::Zero(6,1);
    Eigen::VectorXd generalizedForce = Eigen::VectorXd::Zero(6,1);
    Eigen::VectorXd S_matrix_dot;
    Eigen::MatrixXd articulatedMassMatrix;
    Eigen::VectorXd articulatedBias;
    Eigen::MatrixXd X_BP;
    Eigen::MatrixXd X_BP_dot;
    Eigen::VectorXd jointTwist = Eigen::VectorXd::Zero(6,1);
    Eigen::VectorXd jointTwist_dot = Eigen::VectorXd::Zero(6,1);
    double generalizedAcceleration;
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
                joint->jointAngularVelocity_W = jointAngularVelocity;
                joint->S_matrix_temp << Eigen::VectorXd::Zero(3, 1);
                joint->S_matrix_W << Eigen::VectorXd::Zero(6, 1);
                joint->generalizedForce << 0,0,0,0,0,0;
            } else
            if (joint->type == Joint::Type::FLOATING) {
                joint->jointPosition_W = joint->jointPosition_B;
                //joint->S_matrix = Eigen::MatrixXd::Identity(6,6);
                jointPosition_W = joint->jointPosition_W;
                joint->jointPositionCOM_W = joint->jointPosition_W + joint->jointRotation_B*joint->jointPositionCOM_B;
                joint->inertiaMatrix_W = joint->jointRotation_B*joint->inertiaMatrix_B*joint->jointRotation_B.transpose();
                jointAngularVelocity = joint->jointAngularVelocity;
                joint->jointAngularVelocity_W = jointAngularVelocity;
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
                joint->inertiaMatrix_W = multipliedRotations * joint->jointRotation_B * joint->inertiaMatrix_B *(multipliedRotations * joint->jointRotation_B).transpose();
                joint->S_matrix_temp = multipliedRotations * joint->jointAxis_B;
                joint->S_matrix_W << 0, 0, 0, joint->S_matrix_temp;
                jointAngularVelocity = jointAngularVelocity + multipliedRotations*joint->jointRotation_B*joint->jointGV*joint->jointAxis_B; ;
                joint->jointAngularVelocity = jointAngularVelocity;
                joint->jointAngularVelocity_W = multipliedRotations*joint->jointRotation_B*joint->jointGV*joint->jointAxis_B;
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
                joint->inertiaMatrix_W = multipliedRotations * joint->inertiaMatrix_B * multipliedRotations.transpose();
                joint->S_matrix_temp = multipliedRotations * joint->jointAxis_B;
                joint->S_matrix_W << joint->S_matrix_temp, 0, 0, 0;
                joint->jointAngularVelocity = jointAngularVelocity;
                joint->jointAngularVelocity_W << 0,0,0;
            }
            for (Joint *child: joint->children) {
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

    Eigen::Vector3d parentAngularVelocity;
    inline Eigen::Matrix3d skew(Eigen::Vector3d& vector) {
        Eigen::Matrix3d vector_x;
        vector_x << 0, -vector(2), vector(1),
                vector(2), 0, -vector(0),
                -vector(1), vector(0), 0;
        return vector_x;
    }

    Eigen::Vector3d jointEEPosition_W;
    Eigen::Vector3d parentWorldPosition;
    Joint::Type parentType;

    void computeJointLinearVelocity(Joint* joint, int& counter) {
        if (joint->type == Joint::Type::FLOATING)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED && parentType == Joint::Type::PRISMATIC) {
                joint->jointLinearVelocity = jointLinearVelocity;
                joint->jointTwist << joint->jointLinearVelocity, joint->jointAngularVelocity;
            } else
            if (joint->type == Joint::Type::FIXED) {
                jointLinearVelocity = jointLinearVelocity + parentAngularVelocity.cross(jointEEPosition_W-parentWorldPosition);
                joint->jointLinearVelocity = jointLinearVelocity;
                joint->jointTwist << joint->jointLinearVelocity, joint->jointAngularVelocity;
            }
            if (joint->type == Joint::Type::FLOATING) {
                jointLinearVelocity = joint->jointLinearVelocity;
            } else
            if (joint->type == Joint::Type::REVOLUTE && parentType == Joint::Type::PRISMATIC) {
                jointLinearVelocity = jointLinearVelocity;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                jointLinearVelocity = jointLinearVelocity + parentAngularVelocity.cross(jointEEPosition_W-parentWorldPosition);
            } else
            if ((joint->type == Joint::Type::PRISMATIC) &&  (parentType == Joint::Type::FLOATING || parentType == Joint::Type::REVOLUTE)) {
                jointLinearVelocity = jointLinearVelocity + parentAngularVelocity.cross(jointEEPosition_W-parentWorldPosition) + joint->S_matrix_temp*joint->jointGV;
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                jointLinearVelocity = jointLinearVelocity + joint->S_matrix_temp*joint->jointGV;
            }

            for (Joint *child: joint->children) {
                if (i == counter) {
                    joint->jointLinearVelocity = jointLinearVelocity;
                    joint->jointTwist << joint->jointLinearVelocity, joint->jointAngularVelocity;
                    break;
                }
                parentType = joint->type;
                parentAngularVelocity = joint->jointAngularVelocity_W;
                parentWorldPosition = joint->jointPosition_W;
                i++;
                computeJointLinearVelocity(child,counter);
            }
        }
        if (!trunk_reached) {
            for (Joint *parent: joint->parent) {
                if (counter == 0) {
                    jointEEPosition_W = joint->jointPosition_W;
                }
                counter++;
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
                alfa_x = skew(parentAngularAcceleration);
                omegaParent_x = skew(parentAngularVelocity);
                omega_x = skew(joint->jointAngularVelocity);
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
                joint->S_matrix_dot = S_dot;
                acc << temp_acc+S_dot*joint->jointGV;
                joint->jointLinearAcceleration << acc(0), acc(1), acc(2);
                joint->jointAngularAcceleration << acc(3), acc(4), acc(5);
            } else
            if (joint->type == Joint::Type::FLOATING) {
                jointLinearAcceleration = joint->jointLinearAcceleration;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                Eigen::Matrix3d alfa_x, omegaParent_x, omega_x;
                alfa_x = skew(parentAngularAcceleration);
                omegaParent_x = skew(parentAngularVelocity);
                omega_x = skew(joint->jointAngularVelocity);
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
                joint->S_matrix_dot = S_dot;
                acc << temp_acc+S_dot*joint->jointGV;
                joint->jointLinearAcceleration << acc(0), acc(1), acc(2);
                joint->jointAngularAcceleration << acc(3), acc(4), acc(5);
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                Eigen::Matrix3d alfa_x, omegaParent_x, omega_x;
                alfa_x = skew(parentAngularAcceleration);
                omegaParent_x = skew(parentAngularVelocity);
                omega_x = skew(joint->jointAngularVelocity);
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
                joint->S_matrix_dot = S_dot;
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
                parentAngularVelocity = joint->jointAngularVelocity;
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
        r1_x << skew(r1);
        r2_x << skew(r2);
        Inertia_combined = joint1->inertiaMatrix_W+joint2->inertiaMatrix_W-joint1->mass*r1_x*r1_x-joint2->mass*r2_x*r2_x;

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
        rCOM_x << skew(r_COM);
        spatialInertiaMatrix << m*Eigen::Matrix3d::Identity(3,3), -m*rCOM_x,
                m*rCOM_x,                 joint_j->inertiaMatrix_W-m*rCOM_x*rCOM_x;
        return spatialInertiaMatrix;
    }

    void getCompositeBody_CRBA(Joint* joint1, Joint* joint2) {
        Eigen::Matrix3d Inertia_combined;
        Eigen::Vector3d r_COM, r1, r2;
        Eigen::Matrix3d r1_x, r2_x;
        double m;

        m = joint1->mass + joint2->mass_comp;
        r_COM = 1/m*(joint1->mass*joint1->jointPositionCOM_W + joint2->mass_comp*joint2->jointPositionCOM_comp);
        r1 = joint1->jointPositionCOM_W - r_COM;
        r2 = joint2->jointPositionCOM_comp - r_COM;
        r1_x << skew(r1);
        r2_x << skew(r2);
        Inertia_combined = joint1->inertiaMatrix_W+joint2->inertiaMatrix_comp-joint1->mass*r1_x*r1_x-joint2->mass_comp*r2_x*r2_x;
        joint1->mass_comp = m;
        joint1->jointPositionCOM_comp = r_COM;
        joint1->inertiaMatrix_comp = Inertia_combined;

        if (joint1->type == Joint::Type::FLOATING) {
            joint1->mass = joint1->mass_comp;
            joint1->jointPositionCOM_W = joint1->jointPositionCOM_comp;
            joint1->inertiaMatrix_W = joint1->inertiaMatrix_comp;
        }
    }

    inline Eigen::MatrixXd getSpatialInertiaMatrix_CRBA(Joint* joint_j) {
        Eigen::MatrixXd spatialInertiaMatrix(6,6);
        Eigen::Vector3d r_COM;
        Eigen::Matrix3d rCOM_x;
        double m;
        m = joint_j->mass_comp;
        r_COM = joint_j->jointPositionCOM_comp - joint_j->jointPosition_W;
        rCOM_x << skew(r_COM);
        spatialInertiaMatrix << m*Eigen::Matrix3d::Identity(3,3), -m*rCOM_x,
                m*rCOM_x,                 joint_j->inertiaMatrix_comp-m*rCOM_x*rCOM_x;
        return spatialInertiaMatrix;
    }

    inline Eigen::MatrixXd getBias(Joint* joint_j) {
        Eigen::VectorXd bias(6,1);
        Eigen::Vector3d omega, r_COM;
        Eigen::Matrix3d omega_x, rCOM_x;
        double m;
        m = joint_j->mass;
        r_COM = joint_j->jointPositionCOM_W - joint_j->jointPosition_W;
        rCOM_x << skew(r_COM);
        omega = joint_j->jointAngularVelocity;
        omega_x << skew(omega);

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

    void getX_BP(Joint* joint) {
        Eigen::Vector3d jointPositionRelative_W, jointVelocityRelative_W;
        for (Joint *parent: joint->parent) {
            jointPositionRelative_W = joint->jointPosition_W - parent->jointPosition_W;
            joint->X_BP.setZero(6, 6);
            joint->X_BP.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity(3, 3);
            joint->X_BP.block<3, 3>(3, 0) = skew(jointPositionRelative_W);
            joint->X_BP.block<3, 3>(3, 3) = Eigen::Matrix3d::Identity(3, 3);

            jointVelocityRelative_W = joint->jointLinearVelocity - parent->jointLinearVelocity;
            joint->X_BP_dot.setZero(6, 6);
            joint->X_BP_dot.block<3, 3>(3, 0) = skew(jointVelocityRelative_W);

            getX_BP(parent);
        }
    }

    void computeArticulatedBody(Joint* joint) {
        for (Joint *parent: joint->parent) {
            parent->articulatedMassMatrix = getSpatialInertiaMatrix(parent) + joint->X_BP*joint->articulatedMassMatrix*
                                                                              (-joint->S_matrix_W*(joint->S_matrix_W.transpose()*joint->articulatedMassMatrix*joint->S_matrix_W).inverse()*
                                                                               (joint->S_matrix_W.transpose()*joint->articulatedMassMatrix*joint->X_BP.transpose())+joint->X_BP.transpose());
            parent->articulatedBias = getBias(parent) + joint->X_BP*(joint->articulatedMassMatrix*(joint->S_matrix_W*
                                                                                                   (joint->S_matrix_W.transpose()*joint->articulatedMassMatrix*joint->S_matrix_W).inverse()*
                                                                                                   (joint->jointGF-joint->S_matrix_W.transpose()*joint->articulatedMassMatrix*(joint->S_matrix_dot*joint->jointGV
                                                                                                                                                                               +joint->X_BP_dot.transpose()*parent->jointTwist)-joint->S_matrix_W.transpose()*joint->articulatedBias)
                                                                                                   +joint->S_matrix_dot*joint->jointGV+joint->X_BP_dot.transpose()*parent->jointTwist)+joint->articulatedBias);
            computeArticulatedBody(parent);
        }
    }

    void getGeneralizedAcceleration(Joint* joint) {
        for (Joint *child: joint->children) {
            if (child->type == Joint::Type::FIXED)
                break;
            child->generalizedAcceleration = 1/(child->S_matrix_W.transpose()*child->articulatedMassMatrix*child->S_matrix_W)*
                                             (child->jointGF - child->S_matrix_W.transpose()*child->articulatedMassMatrix*(child->S_matrix_dot*child->jointGV
                                                                                                                           +child->X_BP_dot.transpose()*joint->jointTwist + child->X_BP.transpose()*joint->jointTwist_dot)-child->S_matrix_W.transpose()*child->articulatedBias);
            child->jointTwist_dot = child->S_matrix_W*child->generalizedAcceleration +
                                    child->S_matrix_dot*child->jointGV +
                                    child->X_BP.transpose()*joint->jointTwist_dot +
                                    child->X_BP_dot.transpose()*joint->jointTwist;
            getGeneralizedAcceleration(child);
        }
    }

private:
    Joint* root;
};
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
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
                              gc[7], 0,0,{-0.022191, -0.015144, -1.5e-05}, m_hip, I_FR_hip};
    Joint* FR_thigh = new Joint{Joint::Type::REVOLUTE, {0, -0.083, 0},{0,1,0},
                                gc[8], 0,0,{-0.005607, 0.003877, -0.048199}, m_thigh, I_FR_thigh};
    Joint* FR_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[9],0,0, {0.002781, 6.3e-05, -0.142518}, m_calf, I_FR_calf};
    Joint* FR_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {},{},{}, {}, m_foot, I_FR_foot, Eigen::Matrix3d::Identity()};

    FR_hip->children.push_back(FR_thigh);
    FR_thigh->parent.push_back(FR_hip);
    FR_thigh->children.push_back(FR_calf);
    FR_calf->parent.push_back(FR_thigh);
    FR_calf->children.push_back(FR_foot);
    FR_foot->parent.push_back(FR_calf);

    // Front-Left Leg
    Joint* FL_hip = new Joint{Joint::Type::REVOLUTE, {0.2399, 0.051, 0},{1,0,0},
                              gc[10],0,0,{-0.022191, 0.015144, -1.5e-05}, m_hip, I_FL_hip};
    Joint* FL_thigh = new Joint{Joint::Type::PRISMATIC, {0, 0.083, 0},{0,1,0},
                                gc[11], 0,0,{-0.005607, -0.003877, -0.048199}, m_thigh, I_FL_thigh};
    Joint* FL_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[12], 0, 0, {0.002781, 6.3e-05, -0.142518}, m_calf, I_FL_calf};
    Joint* FL_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {}, {}, {}, m_foot, I_FL_foot, Eigen::Matrix3d::Identity()};
    FL_hip->children.push_back(FL_thigh);
    FL_thigh->parent.push_back(FL_hip);
    FL_thigh->children.push_back(FL_calf);
    FL_calf->parent.push_back(FL_thigh);
    FL_calf->children.push_back(FL_foot);
    FL_foot->parent.push_back(FL_calf);

    // Rear-Right Leg
    Joint* RR_hip = new Joint{Joint::Type::REVOLUTE, {-0.2399, -0.051, 0},{1,0,0},
                              gc[13], 0,0,{0.022191, -0.015144, -1.5e-05}, m_hip, I_RR_hip};
    Joint* RR_thigh = new Joint{Joint::Type::REVOLUTE, {0, -0.083, 0},{0,1,0},
                                gc[14], 0,0,{-0.005607, 0.003877, -0.048199}, m_thigh, I_RR_thigh};
    Joint* RR_calf = new Joint{Joint::Type::REVOLUTE, {0, 0, -0.25},{0,1,0},
                               gc[15],0,0, {0.002781, 6.3e-05, -0.142518}, m_calf, I_RR_calf};
    Joint* RR_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {}, {}, {}, m_foot, I_RR_foot, Eigen::Matrix3d::Identity()};
    RR_hip->children.push_back(RR_thigh);
    RR_thigh->parent.push_back(RR_hip);
    RR_thigh->children.push_back(RR_calf);
    RR_calf->parent.push_back(RR_thigh);
    RR_calf->children.push_back(RR_foot);
    RR_foot->parent.push_back(RR_calf);

    // Rear-Left Leg
    Joint* RL_hip = new Joint{Joint::Type::PRISMATIC, {-0.2399, 0.051, 0},{1,0,0},
                              gc[16],0,0, {0.022191, 0.015144, -1.5e-05}, m_hip, I_RL_hip};
    Joint* RL_thigh = new Joint{Joint::Type::PRISMATIC, {0, 0.083, 0},{0,1,0},
                                gc[17], 0,0,{-0.005607, -0.003877, -0.048199}, m_thigh, I_RL_thigh};
    Joint* RL_calf = new Joint{Joint::Type::PRISMATIC, {0, 0, -0.25},{0,1,0},
                               gc[18], 0,0, {0.002781, 6.3e-05, -0.142518}, m_calf, I_RL_calf};
    Joint* RL_foot = new Joint{Joint::Type::FIXED, {0, 0, -0.25},{0,0,0},
                               {}, {}, {}, {}, m_foot, I_RL_foot, Eigen::Matrix3d::Identity()};
    RL_hip->children.push_back(RL_thigh);
    RL_thigh->parent.push_back(RL_hip);
    RL_thigh->children.push_back(RL_calf);
    RL_calf->parent.push_back(RL_thigh);
    RL_calf->children.push_back(RL_foot);
    RL_foot->parent.push_back(RL_calf);

    // Create the tree structure for the quadruped robot
    Joint* imu = new Joint{Joint::Type::FIXED, {0,0,0},{0,0,0},
                           0, 0, 0,{0,0,0}, m_imu, I_imu, Eigen::Matrix3d::Identity()};

    // convert quaternion to rotation matrix
    Eigen::Quaterniond q;
    q.w() = gc[3];
    q.x() = gc[4];
    q.y() = gc[5];
    q.z() = gc[6];
    //q.toRotationMatrix();

    Joint* trunk = new Joint{Joint::Type::FLOATING, gc.head(3),{0,0,0},
                             0, 0, 0,{0.008465, 0.004045, -0.000763}, m_trunk, I_trunk,
                             q.toRotationMatrix(),{}};


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

    // detach RL leg
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    trunk->children.push_back(imu);
    imu->parent.push_back(trunk);
    JointTree tree_imu{trunk};
    tree_imu.computeKinematicsDownTheTree(imu, jointsRotationMatrices, counter);
    counter = 0;
    tree_imu.trunk_reached = false;
    tree_imu.i = 0;
    tree_imu.getCompositeBody_CRBA(trunk, imu);
    trunk->children.pop_back();
    jointsRotationMatrices.clear();

    Eigen::MatrixXd spatialInertiaMatrix(6,6), conversionMatrix(6,6);
    Eigen::Matrix3d r_x;
    Eigen::Vector3d r;

    tree_FR.getCompositeBody_CRBA(FR_calf, FR_foot);
    tree_FR.getCompositeBody_CRBA(FR_thigh, FR_calf);
    tree_FR.getCompositeBody_CRBA(FR_hip, FR_thigh);
    tree_FR.getCompositeBody_CRBA(trunk, FR_hip);
    Eigen::Matrix3d M_FR = Eigen::Matrix3d::Zero(3,3);
    Eigen::MatrixXd coupling_FR = Eigen::MatrixXd::Zero(6,3);
    std::vector<Joint*> joint_listFR = {FR_hip, FR_thigh, FR_calf};
    for (int i=static_cast<int>(joint_listFR.size())-1; i>=0; i--) {
        for (int j=i; j>=0; j--) {
            spatialInertiaMatrix = tree_FR.getSpatialInertiaMatrix_CRBA(joint_listFR[i]);
            r = joint_listFR[i]->jointPosition_W-joint_listFR[j]->jointPosition_W;
            r_x << tree_FR.skew(r);
            conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                    Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);

            M_FR.coeffRef(j,i) = joint_listFR[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*joint_listFR[j]->S_matrix_W;
            M_FR.coeffRef(i,j) = M_FR.coeffRef(j,i);
        }
        // Coupling with trunk
        r = joint_listFR[i]->jointPosition_W - trunk->jointPosition_W;
        r_x << tree_FR.skew(r);
        conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
        coupling_FR.col(i) = (joint_listFR[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*Eigen::MatrixXd::Identity(6,6)).transpose();
    }

    tree_FL.getCompositeBody_CRBA(FL_calf, FL_foot);
    tree_FL.getCompositeBody_CRBA(FL_thigh, FL_calf);
    tree_FL.getCompositeBody_CRBA(FL_hip, FL_thigh);
    tree_FL.getCompositeBody_CRBA(trunk, FL_hip);
    Eigen::Matrix3d M_FL = Eigen::Matrix3d::Zero(3,3);
    Eigen::MatrixXd coupling_FL = Eigen::MatrixXd::Zero(6,3);
    std::vector<Joint*> joint_listFL = {FL_hip, FL_thigh, FL_calf};
    for (int i=static_cast<int>(joint_listFL.size())-1; i>=0; i--) {
        for (int j=i; j>=0; j--) {
            spatialInertiaMatrix = tree_FL.getSpatialInertiaMatrix_CRBA(joint_listFL[i]);
            r = joint_listFL[i]->jointPosition_W-joint_listFL[j]->jointPosition_W;
            r_x << tree_FL.skew(r);

            conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                    Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);

            M_FL.coeffRef(j,i) = joint_listFL[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*joint_listFL[j]->S_matrix_W;
            M_FL.coeffRef(i,j) = M_FL.coeffRef(j,i);
        }
        // Coupling with trunk
        r = joint_listFL[i]->jointPosition_W - trunk->jointPosition_W;
        r_x << tree_FL.skew(r);
        conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
        coupling_FL.col(i) = (joint_listFL[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*Eigen::MatrixXd::Identity(6,6)).transpose();
    }

    tree_RR.getCompositeBody_CRBA(RR_calf, RR_foot);
    tree_RR.getCompositeBody_CRBA(RR_thigh, RR_calf);
    tree_RR.getCompositeBody_CRBA(RR_hip, RR_thigh);
    tree_RR.getCompositeBody_CRBA(trunk, RR_hip);
    Eigen::Matrix3d M_RR = Eigen::Matrix3d::Zero(3,3);
    Eigen::MatrixXd coupling_RR = Eigen::MatrixXd::Zero(6,3);
    std::vector<Joint*> joint_listRR = {RR_hip, RR_thigh, RR_calf};
    for (int i=static_cast<int>(joint_listRR.size())-1; i>=0; i--) {
        for (int j=i; j>=0; j--) {
            spatialInertiaMatrix = tree_RR.getSpatialInertiaMatrix_CRBA(joint_listRR[i]);
            r = joint_listRR[i]->jointPosition_W-joint_listRR[j]->jointPosition_W;
            r_x << tree_RR.skew(r);
            conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                    Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
            M_RR.coeffRef(j,i) = joint_listRR[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*joint_listRR[j]->S_matrix_W;
            M_RR.coeffRef(i,j) = M_RR.coeff(j,i);
        }
        // Coupling with trunk
        r = joint_listRR[i]->jointPosition_W - trunk->jointPosition_W;
        r_x << tree_RR.skew(r);
        conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
        coupling_RR.col(i) = (joint_listRR[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*Eigen::MatrixXd::Identity(6,6)).transpose();
    }

    tree_RL.getCompositeBody_CRBA(RL_calf, RL_foot);
    tree_RL.getCompositeBody_CRBA(RL_thigh, RL_calf);
    tree_RL.getCompositeBody_CRBA(RL_hip, RL_thigh);
    tree_RL.getCompositeBody_CRBA(trunk, RL_hip);
    Eigen::Matrix3d M_RL = Eigen::Matrix3d::Zero(3,3);
    Eigen::MatrixXd coupling_RL = Eigen::MatrixXd::Zero(6,3);
    std::vector<Joint*> joint_listRL = {RL_hip, RL_thigh, RL_calf};
    for (int i=static_cast<int>(joint_listRL.size())-1; i>=0; i--) {
        for (int j=i; j>=0; j--) {
            spatialInertiaMatrix = tree_RL.getSpatialInertiaMatrix_CRBA(joint_listRL[i]);
            r = joint_listRL[i]->jointPosition_W-joint_listRL[j]->jointPosition_W;
            r_x << tree_RL.skew(r);

            conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                    Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
            M_RL.coeffRef(j,i) = joint_listRL[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*joint_listRL[j]->S_matrix_W;
            M_RL.coeffRef(i,j) = M_RL.coeff(j,i);
        }
        // Coupling with trunk
        r = joint_listRL[i]->jointPosition_W - trunk->jointPosition_W;
        r_x << tree_RL.skew(r);
        conversionMatrix << Eigen::Matrix3d::Identity(3,3), -r_x,
                Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Identity(3,3);
        coupling_RL.col(i) = (joint_listRL[i]->S_matrix_W.transpose()*spatialInertiaMatrix*conversionMatrix*Eigen::MatrixXd::Identity(6,6)).transpose();
    }

    Eigen::MatrixXd M_legs(12,12);
    M_legs = Eigen::MatrixXd::Zero(12,12);
    M_legs.block<3,3>(0,0) = M_FR;
    M_legs.block<3,3>(3,3) = M_FL;
    M_legs.block<3,3>(6,6) = M_RR;
    M_legs.block<3,3>(9,9) = M_RL;

    Eigen::MatrixXd M_trunk(6,6);
    M_trunk = tree_imu.getSpatialInertiaMatrix_CRBA(trunk);

    Eigen::MatrixXd MassMatrix(18,18);
    MassMatrix = Eigen::MatrixXd::Zero(18,18);
    MassMatrix.block<6,6>(0,0) = M_trunk;
    MassMatrix.block<12,12>(6,6) = M_legs;
    MassMatrix.block<6,3>(0,6) = coupling_FR;
    MassMatrix.block<6,3>(0,9) = coupling_FL;
    MassMatrix.block<6,3>(0,12) = coupling_RR;
    MassMatrix.block<6,3>(0,15) = coupling_RL;
    MassMatrix.block<3,6>(6,0) = coupling_FR.transpose();
    MassMatrix.block<3,6>(9,0) = coupling_FL.transpose();
    MassMatrix.block<3,6>(12,0) = coupling_RR.transpose();
    MassMatrix.block<3,6>(15,0) = coupling_RL.transpose();

    return MassMatrix;
}