//
// Created by Jemin Hwangbo on 2022/03/17.
//


#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "midterm_STUDENTID.hpp"

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  // a1
  // aliengo
  auto cartPole = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/cartPole/cartpole.urdf");
  cartPole->setName("cartpole");
  server.focusOn(cartPole);

  // a1 configuration
  Eigen::VectorXd gc(cartPole->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(cartPole->getDOF());

  gc << 0, 0.5;
  gv << 0, 0;
  cartPole->setState(gc, gv);

  // visualization
  server.launchServer();
  bool answerCorrect = true;

  for (int i=0; i<3; i++) {
    RS_TIMED_LOOP(world.getTimeStep()*1e6)

    cartPole->getState(gc, gv);
    if((cartPole->getMassMatrix().e() - getMassMatrix(gc)).norm() < 1e-8) {
      std::cout<<"the mass matrix is correct "<<std::endl;
    } else {
      std::cout<<"the mass matrix is not correct "<<std::endl;
      answerCorrect = false;
    }

    server.integrateWorldThreadSafe();
  }

  server.killServer();

  if(answerCorrect) {
    std::cout<<"\n\nThe solution is correct "<<std::endl;
  } else {
    std::cout<<"\n\nThe solution is not correct "<<std::endl;
  }

  return 0;
}
