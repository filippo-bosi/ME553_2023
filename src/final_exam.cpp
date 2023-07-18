//
// Created by jemin on 23. 6. 15.
//

#include "raisim/RaisimServer.hpp"
#include "finalexam_20236014.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // kinova
  auto robotArm = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/2DrobotArm/robot_3D.urdf");

  // kinova configuration
  Eigen::VectorXd gc(robotArm->getGeneralizedCoordinateDim()), gv(robotArm->getDOF()), gf(robotArm->getDOF());
  gc << 0.1, 0.2, 0.3; /// Jemin: I'll randomize the gc, gv, gf when grading
  gv << 0.1, 0.2, 0.3;
  gf << 0.15, 0.21, 0.36;
  robotArm->setState(gc, gv);
  robotArm->setGeneralizedForce(gf);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  Eigen::VectorXd nonlinearity(robotArm->getDOF());
  Eigen::MatrixXd massMatrix(robotArm->getDOF(), robotArm->getDOF());
  massMatrix = robotArm->getMassMatrix().e();
  nonlinearity = robotArm->getNonlinearities({0,0,-9.81}).e();
    Eigen::VectorXd result(robotArm->getDOF());
    result = massMatrix.inverse() * (gf-nonlinearity);
    std::cout << result.transpose() << std::endl;
  if((computeGeneralizedAcceleration(gc, gv, gf) - massMatrix.inverse() * (gf-nonlinearity)).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  std::cout << computeGeneralizedAcceleration(gc, gv, gf).transpose() << std::endl;
  server.launchServer();

//  raisim::Vec<3> pos;
//  raisim::Vec<3> vel1, vel2, vel3, vel4;
//  raisim::Vec<3> ang1, ang2, ang3, ang4;

//  std::cout << "Position: " << std::endl;
//  robotArm->getFramePosition("worldTolink1", pos);
//  std::cout << pos << std::endl;
//  robotArm->getFramePosition("link1Tolink2", pos);
//  std::cout << pos << std::endl;
//  robotArm->getFramePosition("link2Tolink3", pos);
//  std::cout << pos << std::endl;
//  robotArm->getFramePosition("link3Tolink4", pos);
//  std::cout << pos << std::endl;

//  robotArm->getFrameVelocity("worldTolink1", vel1);
//  robotArm->getFrameVelocity("link1Tolink2", vel2);
//  robotArm->getFrameVelocity("link2Tolink3", vel3);
//  robotArm->getFrameVelocity("link3Tolink4", vel4);
//
//  std::cout << "Linear vel:" << std::endl;
//
//  std::cout << vel1 << std::endl;
//  std::cout << vel2 << std::endl;
//  std::cout << vel3 << std::endl;
//  std::cout << vel4 << std::endl;

//    robotArm->getFrameAngularVelocity("worldTolink1", ang1);
//    robotArm->getFrameAngularVelocity("link1Tolink2", ang2);
//    robotArm->getFrameAngularVelocity("link2Tolink3", ang3);
//    robotArm->getFrameAngularVelocity("link3Tolink4", ang4);

//    std::cout << "Angular vel:" << std::endl;
//    std::cout << ang1 << std::endl;
//    std::cout << ang2 << std::endl;
//    std::cout << ang3 << std::endl;
//    std::cout << ang4 << std::endl;


  /// this is for visualization only
  while (true) {
    RS_TIMED_LOOP(1)
  }

  server.closeConnection();
  return 0;
}
