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
        if (joint->type == Joint::Type::FIXED)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED) {
                joint->jointPosition_W = joint->jointPosition_B;
                joint->S_matrix_W = Eigen::MatrixXd::Zero(6,1);
                jointPosition_W = joint->jointPosition_W;
                joint->jointPositionCOM_W = joint->jointPosition_W + joint->jointRotation_B*joint->jointPositionCOM_B;
                joint->inertiaMatrix_W = joint->jointRotation_B*joint->inertiaMatrix_B*joint->jointRotation_B.transpose();
                jointAngularVelocity = joint->jointAngularVelocity;
                joint->jointAngularVelocity_W = jointAngularVelocity;
//                joint->jointPositionCOM_B = joint->jointPosition_B;
//                joint->jointRotation_B = Eigen::Matrix3d::Identity(3,3);
//                jointPosition_W = jointPosition_W + multipliedRotations *  joint->jointPosition_B;
//                joint->jointPosition_W = jointPosition_W;
//                joint->jointPositionCOM_W = jointPosition_W;
//                joint->inertiaMatrix_W = multipliedRotations*joint->inertiaMatrix_B*multipliedRotations.transpose();
//                // Initialization for composite body
//                joint->jointPositionCOM_comp = joint->jointPositionCOM_W;
//                joint->mass_comp = joint->mass;
//                joint->inertiaMatrix_comp = joint->inertiaMatrix_W;
//                joint->jointAngularVelocity = jointAngularVelocity;
//                joint->jointAngularVelocity_W = jointAngularVelocity;
//                joint->S_matrix_temp << Eigen::VectorXd::Zero(3, 1);
//                joint->S_matrix_W << Eigen::VectorXd::Zero(6, 1);
//                joint->generalizedForce << 0,0,0,0,0,0;
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
        if (joint->type == Joint::Type::FIXED)
            trunk_reached = true;
        if (trunk_reached) {
            if (joint->type == Joint::Type::FIXED) {
                jointLinearVelocity = joint->jointLinearVelocity;
                joint->jointTwist << joint->jointLinearVelocity, joint->jointAngularVelocity;
            } else
            if (joint->type == Joint::Type::REVOLUTE) {
                jointLinearVelocity = jointLinearVelocity + parentAngularVelocity.cross(jointEEPosition_W-parentWorldPosition);
            } else
            if (joint->type == Joint::Type::PRISMATIC) {
                jointLinearVelocity = jointLinearVelocity + parentAngularVelocity.cross(jointEEPosition_W-parentWorldPosition) + joint->S_matrix_temp*joint->jointGV;
                joint->jointLinearVelocity = jointLinearVelocity;
                joint->jointTwist << joint->jointLinearVelocity, joint->jointAngularVelocity;
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
        if (joint->type == Joint::Type::FIXED)
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
//            if (child->type == Joint::Type::FIXED)
//                break;
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

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
    double m1, m2, m3, m4;
    m1 = 0;
    m2 = 1;
    m3 = 1;
    m4 = 1;

    double ixx, ixy, ixz, iyy, iyz, izz;
    Eigen::Matrix3d I1, I2, I3, I4;
    I1 = Eigen::Matrix3d::Zero();
    ixx=0.001; ixy=0; ixz=0; iyy=0.001; iyz=0; izz=0.001;
    I2 << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    ixx=0.001; ixy=0; ixz=0; iyy=0.001; iyz=0; izz=0.001;
    I3 << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
    ixx=0.001; ixy=0; ixz=0; iyy=0.001; iyz=0; izz=0.001;
    I4 << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;

    Joint* base = new Joint{Joint::Type::FIXED, {0, 0, 0},{0,0,0},
                               {},{},{}, {}, m1, I1, Eigen::Matrix3d::Identity()};
    Joint* Body2 = new Joint{Joint::Type::REVOLUTE, {0, 0, 0.3},{1,0,0},
                              gc[0], gv[0], gf[0], {0, 0, 0.2}, m2, I2};
    Joint* Body3 = new Joint{Joint::Type::REVOLUTE, {0, 0, 0.3},{1,0,0},
                             gc[1], gv[1], gf[1], {0, 0, 0.2}, m3, I3};
    Joint* Body4 = new Joint{Joint::Type::PRISMATIC, {0, 0, 0.3},{0,1,0},
                             gc[2], gv[2], gf[2], {0, 0, 0.2}, m4, I4};
    base->children.push_back(Body2);
    Body2->parent.push_back(base);
    Body2->children.push_back(Body3);
    Body3->parent.push_back(Body2);
    Body3->children.push_back(Body4);
    Body4->parent.push_back(Body3);

    base->jointLinearVelocity << 0, 0, 0;
    base->jointAngularVelocity << 0, 0, 0;
    base->jointTwist << base->jointLinearVelocity, base->jointAngularVelocity;
    base->jointLinearAcceleration << 0, 0, 9.81;
    base->jointAngularAcceleration << 0, 0, 0;

    std::vector<Eigen::Matrix3d> jointsRotationMatrices;
    Eigen::MatrixXd trunkArticulatedMassMatrix = Eigen::MatrixXd::Zero(6,6);
    Eigen::VectorXd trunkArticulatedBias = Eigen::VectorXd::Zero(6,1);
    int counter = 0;

    JointTree tree_1{base};
    tree_1.computeKinematicsDownTheTree(Body4,jointsRotationMatrices, counter);
    counter = 0;
    tree_1.trunk_reached = false;
    tree_1.i = 0;
    tree_1.computeAcceleration(Body4, counter);
    counter = 0;
    tree_1.trunk_reached = false;
    tree_1.i = 0;
    tree_1.computeJointLinearVelocity(Body2, counter);
    counter = 0;
    tree_1.trunk_reached = false;
    tree_1.i = 0;
    tree_1.computeJointLinearVelocity(Body3, counter);
    counter = 0;
    tree_1.trunk_reached = false;
    tree_1.i = 0;
    tree_1.computeJointLinearVelocity(Body4, counter);
    tree_1.getX_BP(Body4);
    Body4->articulatedMassMatrix = tree_1.getSpatialInertiaMatrix(Body4);
    Body4->articulatedBias = tree_1.getBias(Body4);
    tree_1.computeArticulatedBody(Body4);

    Eigen::VectorXd generalizedAccelerations;
    base->jointTwist_dot << 0, 0, 9.81, 0, 0, 0;
    tree_1.getGeneralizedAcceleration(base);
    generalizedAccelerations.setZero(3);
    generalizedAccelerations << Body2->generalizedAcceleration, Body3->generalizedAcceleration, Body4->generalizedAcceleration;

    return generalizedAccelerations;
}