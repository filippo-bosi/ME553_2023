#ifndef ME553_2022_SOLUTIONS_EXERCISE1_20236014_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_20236014_HPP_

#include <Eigen/Core>
#include <iostream>
#include <raisim/math.hpp>
#include "Eigen/Geometry"

inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {

    // variables declaration
    raisim::Vec<3> rw_0, rw_1, rw_2, rw_3, rw_00p, r0p_0p1, r1p_1p2, r2p_2p3, r3p_e, rw_e;
    raisim::Vec<4> q = {gc[3], gc[4], gc[5], gc[6]};
    raisim::Mat<3, 3> R_00p, R_0p1, R_11p, R_1p2, R_22p, R_2p3, R_33p;
    raisim::Vec<3> a1, a2, a3;
    // translation vector from world to trunk
    rw_00p = gc.head(3);

    // quaternion to rotation matrix conversion
    raisim::quatToRotMat(q, R_00p);

    // identity matrices from URDF
    R_0p1.setIdentity();
    R_1p2.setIdentity();
    R_2p3.setIdentity();

    // translation vectors from aliengo URDF file
    r0p_0p1 = {0.2399, -0.051, 0};
    r1p_1p2 = {0, -0.083, 0};
    r2p_2p3 = {0, 0, -0.25};
    r3p_e = {0, 0, -0.25};

    // rotation axis from URDF
    a1 = {1,0,0};
    a2 = {0,1,0};
    a3 = {0,1,0};

    // rotation matrices dependent on robot's configuration
    raisim::angleAxisToRotMat(a1, gc[7], R_11p);
    raisim::angleAxisToRotMat(a2, gc[8], R_22p);
    raisim::angleAxisToRotMat(a3, gc[9], R_33p);

    // get the end-effector (FR_foot_fixed) position
    // in a floating-base system, frame 0 is the World frame
    rw_0 = rw_00p;
    rw_1 = rw_0+R_00p*r0p_0p1;
    rw_2 = rw_1+R_00p*R_0p1*R_11p*r1p_1p2;
    rw_3 = rw_2+R_00p*R_0p1*R_11p*R_1p2*R_22p*r2p_2p3;
    rw_e = rw_3+R_00p*R_0p1*R_11p*R_1p2*R_22p*R_2p3*R_33p*r3p_e;


    std::cout << rw_0 << std::endl;
    std::cout << rw_1 << std::endl;
    std::cout << rw_2 << std::endl;
    std::cout << rw_3 << std::endl;
    std::cout << rw_e << std::endl;
    return rw_e.e();
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_20236014_HPP_
