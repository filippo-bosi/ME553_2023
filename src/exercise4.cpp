//
// Created by Jemin Hwangbo on 2022/04/08.
//

#include "raisim/RaisimServer.hpp"
#include "exercise4_20236014.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world

  // kinova
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo_modified.urdf");

    // kinova configuration
  Eigen::VectorXd gc(aliengo->getGeneralizedCoordinateDim()), gv(aliengo->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  aliengo->setState(gc, gv);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  aliengo->getMassMatrix();
    std::cout.precision(8);
    std::cout<<"nonlinearities should be \n"<< aliengo->getNonlinearities({0,0,-9.81}).e().transpose()<<std::endl;
//    raisim::Vec<3> force, torque;
//    aliengo->setComputeInverseDynamics(true);
//    std::cout << aliengo->getGeneralizedForce() << std::endl;
//
//    aliengo->setComputeInverseDynamics(true);
//    force = aliengo->getForceAtJointInWorldFrame(2);
//    torque = aliengo->getTorqueAtJointInWorldFrame(2);
//    std::cout<< force <<std::endl;
//    std::cout<< torque <<std::endl;

  if((getNonlinearities(gc, gv) - aliengo->getNonlinearities({0,0,-9.81}).e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else {
      std::cout << "failed " << std::endl;
      std::cout << "Error:" << std::endl;
      std::cout << (getNonlinearities(gc, gv) - aliengo->getNonlinearities({0,0,-9.81}).e()).norm() << std::endl;
  }

//  raisim::Vec<3> rootVel, hipVel, thighVel, calfVel, footVel;
//    aliengo->getFrameAngularVelocity("ROOT", rootVel);
//    aliengo->getFrameAngularVelocity("FR_hip_joint", hipVel);
//    aliengo->getFrameAngularVelocity("FR_thigh_joint", thighVel);
//    aliengo->getFrameAngularVelocity("FR_calf_joint", calfVel);
//    aliengo->getFrameAngularVelocity("FR_foot_fixed", footVel);
//    std::cout << " Velocity " << std::endl;
//    std::cout << rootVel << std::endl;
//    std::cout << hipVel << std::endl;
//    std::cout << thighVel << std::endl;
//    std::cout << calfVel << std::endl;
//    std::cout << footVel << std::endl;
//    raisim::Vec<3> hipAcc, rootAcc, thighAcc, calfAcc, footAcc;
//    aliengo->getFrameAcceleration("ROOT", rootAcc);
//    aliengo->getFrameAcceleration("FR_hip_joint", hipAcc);
//    aliengo->getFrameAcceleration("FR_thigh_joint", thighAcc);
//    aliengo->getFrameAcceleration("FR_calf_joint", calfAcc);
//    aliengo->getFrameAcceleration("FR_foot_fixed", footAcc);
//    std::cout << " Acceleration " << std::endl;
//    std::cout << rootAcc << std::endl;
//    std::cout << hipAcc << std::endl;
//    std::cout << thighAcc << std::endl;
//    std::cout << calfAcc << std::endl;
//    std::cout << footAcc << std::endl;
    return 0;
}
