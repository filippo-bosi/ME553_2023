//
// Created by Jemin Hwangbo on 2022/03/17.
//


#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "exercise2_20236014.hpp"

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  // a1
  // aliengo
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo_modified.urdf");
  aliengo->setName("aliengo");
  server.focusOn(aliengo);

  // a1 configuration
  Eigen::VectorXd gc(aliengo->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(aliengo->getDOF());
    Eigen::VectorXd gf(aliengo->getDOF());

  gc << 0, 0, 10.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  gv << 0.1, 0.2, 0.3, 0.1, 0.4, 0.3, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4;
    gf << 0.15, 0.21, 0.36, 0.24, 0.35, 0.46, 0.57, 0.18, 0.29, 1.0, 1.1, 1.5, 1.1, 1.2, 1.3, 1.6, 1.7, 1.8;

    aliengo->setState(gc, gv);

  // visualization
  server.launchServer();
  raisim::Vec<3> rootAngVel, footVel, footAngVel, calfAngVel, thighAngVel, hipAngVel;
    raisim::Vec<3> foot2Vel, foot2AngVel, calf2AngVel, thigh2AngVel, hip2AngVel;
    raisim::Vec<3> footAcc, calfAcc, thighAcc, hipAcc, rootAcc;

    bool answerCorrect = true;

  for (int i=0; i<200; i++) {
    RS_TIMED_LOOP(world.getTimeStep()*1e6);

//      aliengo->getFrameVelocity("ROOT", rootAngVel);
//
//      aliengo->getFrameVelocity("RL_foot_fixed", footAngVel);
//      aliengo->getFrameVelocity("RL_calf_joint", calfAngVel);
//      aliengo->getFrameVelocity("RL_thigh_joint", thighAngVel);
//      aliengo->getFrameVelocity("RL_hip_joint", hipAngVel);

//      aliengo->getFrameVelocity("RL_foot_fixed", foot2AngVel);
//      aliengo->getFrameVelocity("RL_calf_joint", calf2AngVel);
//      aliengo->getFrameVelocity("RL_thigh_joint", thigh2AngVel);
//      aliengo->getFrameVelocity("RL_hip_joint", hip2AngVel);

//    if((thighAngVel.e() - getFootLinearVelocity(gc, gv)).norm() < 1e-8) {
//      std::cout<<"the linear velocity is correct "<<std::endl;
//    } else {
//      std::cout<<"the linear velocity is not correct "<<std::endl;
//      answerCorrect = false;
//    }

//      std::cout << " Actual Linear velocity" << std::endl;
//
//      std::cout << footAngVel << std::endl;
//      std::cout << calfAngVel << std::endl;
//      std::cout << thighAngVel << std::endl;
//      std::cout << hipAngVel << std::endl << std::endl;
//      std::cout << rootAngVel << std::endl << std::endl;


    aliengo->getMassMatrix();
    aliengo->getNonlinearities({0,0,-9.81});
    aliengo->getFrameAcceleration("FR_hip_joint", hipAcc);
    aliengo->getFrameAcceleration("ROOT", rootAcc);

      std::cout << " Acceleration " << std::endl;
      std::cout << rootAcc << std::endl;
      std::cout << hipAcc << std::endl;

      std::cout << " My Acc" << std::endl;
      if((hipAcc.e() - getFootLinearVelocity(gc, gv)).norm() < 1e-8) {
      std::cout<<"the acc is correct "<<std::endl;
    } else {
      std::cout<<"the acc is not correct "<<std::endl;
      answerCorrect = false;
    }

      server.integrateWorldThreadSafe();
    aliengo->getState(gc, gv);
  }

  server.killServer();

//  if(answerCorrect) {
//    std::cout<<"The solution is correct "<<std::endl;
//  } else {
//    std::cout<<"The solution is not correct "<<std::endl;
//  }

//    raisim::Vec<3> hipAcc;
//    aliengo->getFrameAcceleration("FR_hip_joint", hipAcc);
//    std::cout << hipAcc << std::endl;

//    std::cout << foot2AngVel << std::endl;
//    std::cout << calf2AngVel << std::endl;
//    std::cout << thigh2AngVel << std::endl;
//    std::cout << hip2AngVel << std::endl << std::endl;
//

    aliengo->getMassMatrix();

    return 0;
}
