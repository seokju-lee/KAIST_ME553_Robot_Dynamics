//
// Created by Jemin Hwangbo on 2022/03/17.
//


#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

#include "raisim/RaisimServer.hpp"
#include "exercise2_20233536.hpp"

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world); // visualization server
  world.addGround();
  world.setTimeStep(0.001);

  // a1
  // aliengo
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo.urdf");
  aliengo->setName("aliengo");
  server.focusOn(aliengo);
  
  // a1 configuration
  Eigen::VectorXd gc(aliengo->getGeneralizedCoordinateDim());
  Eigen::VectorXd gv(aliengo->getDOF());

  gc << 0, 0, 10.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8;
  gv << 0.1, 0.2, 0.3, 0.1, 0.4, 0.3, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4;
  aliengo->setState(gc, gv);

  // visualization
  server.launchServer();
  raisim::Vec<3> footVel, footAngVel;
  bool answerCorrect = true;

  for (int i=0; i<2000; i++) {
    RS_TIMED_LOOP(world.getTimeStep()*1e6);

    aliengo->getFrameVelocity("FR_foot_fixed", footVel);
    aliengo->getFrameAngularVelocity("FR_foot_fixed", footAngVel);

    if((footVel.e() - getFootLinearVelocity(gc, gv)).norm() < 1e-8) {
      std::cout<<"the linear velocity is correct "<<std::endl;
    } else {
      std::cout<<"the linear velocity is not correct "<<std::endl;
      answerCorrect = false;
    }

    if((footAngVel.e() - getFootAngularVelocity(gc, gv)).norm() < 1e-8) {
      std::cout<<"the angular velocity is correct "<<std::endl;
    } else {
      std::cout<<"the angular velocity is not correct "<<std::endl;
      answerCorrect = false;
    }

    server.integrateWorldThreadSafe();
    aliengo->getState(gc, gv);
  }

  server.killServer();

  if(answerCorrect) {
    std::cout<<"The solution is correct "<<std::endl;
  } else {
    std::cout<<"The solution is not correct "<<std::endl;
  }

  return 0;
}
