//
// Created by jemin on 23. 6. 15.
//

#include "raisim/RaisimServer.hpp"
#include "finalexam_STUDENTID.hpp"

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

  if((computeGeneralizedAcceleration(gc, gv, gf) - massMatrix.inverse() * (gf-nonlinearity)).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  server.launchServer();

  /// this is for visualization only
  while (true) {
    RS_TIMED_LOOP(1)
  }

  server.closeConnection();
  return 0;
}
