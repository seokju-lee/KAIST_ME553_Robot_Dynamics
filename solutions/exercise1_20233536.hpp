//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE1_20233536_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_20233536_HPP_

#include <Eigen/Core>
#include <iostream>

namespace seokju {
class Body {
  public:
    struct Joint{
      /// variables
      Eigen::Vector3d jointPosition_W;
      Eigen::Vector3d jointAxis_W;
      Eigen::Matrix3d rotation_W;

      /// robot definition
      Eigen::Vector3d jointAxis_P;
      Eigen::Vector3d jointPosition_P;
      int i;
      enum class Type {
        FIXED,
        FLOATING,
        REVOLUTE,
        PRISMATC
      } type = Type::REVOLUTE;
      Eigen::Matrix3d jointRotation_P;

    };

    Body(const Joint::Type type,
         const Eigen::Vector3d& jointAxis_P, 
         const Eigen::Vector3d& jointPosition_P,
         const Eigen::Matrix3d& jointRotation_P,
         int i) {
      joint_.type = type;
      joint_.jointAxis_P = jointAxis_P;
      joint_.jointPosition_P = jointPosition_P;
      joint_.jointRotation_P = jointRotation_P;
      joint_.i = i;
    }

    void computeKinematicsDownTheTree(const Eigen::VectorXd& gc) {
      // Floating base type
      if (joint_.type == Joint::Type::FLOATING) {
        // set the joint postion
        joint_.jointPosition_W[0] = gc[0];
        joint_.jointPosition_W[1] = gc[1];
        joint_.jointPosition_W[2] = gc[2];

        // convert quaternion to rotation matrix
        joint_.rotation_W << 1-2*pow(gc[5],2)-2*pow(gc[6],2), 2*gc[4]*gc[5]-2*gc[3]*gc[6], 2*gc[4]*gc[6]+2*gc[3]*gc[5],
                                  2*gc[4]*gc[5]+2*gc[3]*gc[6], 1-2*pow(gc[4],2)-2*pow(gc[6],2), 2*gc[5]*gc[6]-2*gc[3]*gc[4],
                                  2*gc[4]*gc[6]-2*gc[3]*gc[5], 2*gc[5]*gc[6]+2*gc[3]*gc[4], 1-2*pow(gc[4],2)-2*pow(gc[5],2);
      }
      // Revolute joint type
      else if (joint_.type == Joint::Type::REVOLUTE) {
        // Rotation about x-axis
        if (joint_.jointAxis_P[0] == 1 && joint_.jointAxis_P[1] == 0 && joint_.jointAxis_P[2] == 0){
          joint_.jointRotation_P << 1,        0,         0,
                                    0, cos(gc[joint_.i]), -sin(gc[joint_.i]),
                                    0, sin(gc[joint_.i]),  cos(gc[joint_.i]);
        }
        // Rotation about y-axis
        else if (joint_.jointAxis_P[0] == 0 && joint_.jointAxis_P[1] == 1 && joint_.jointAxis_P[2] == 0){
          joint_.jointRotation_P << cos(gc[joint_.i]), 0, sin(gc[joint_.i]),
                                       0,       1,      0,
                                   -sin(gc[joint_.i]), 0, cos(gc[joint_.i]);
        }
      }

      for (auto& par: parent_){
        par->computeKinematicsDownTheTree(gc);
        joint_.rotation_W = par->joint_.rotation_W*joint_.jointRotation_P;
        joint_.jointPosition_W = par->joint_.rotation_W*joint_.jointPosition_P + par->joint_.jointPosition_W;
      }    
    }
    
    void setParent(std::vector<Body*> parent){
      parent_ = parent;
    }

    Joint joint_;
    std::vector<Body*> children_, parent_;
    Body* parent__ = nullptr;
};
}

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  // set the body
  Eigen::Vector3d tar_pos;
  Eigen::Matrix3d identity;
  
  identity << 1, 0, 0,
          0, 1, 0,
          0, 0, 1;
  
  seokju::Body trunk (seokju::Body::Joint::Type::FLOATING,
                      {0, 0, 0}, 
                      {gc[0], gc[1], gc[2]},
                      identity, 0);
  seokju::Body FR_hip (seokju::Body::Joint::Type::REVOLUTE,
                      {1, 0, 0}, 
                      {0.2399, -0.051, 0},
                      identity, 7);
  seokju::Body FR_thigh (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0},
                        {0, -0.083, 0},
                        identity, 8);
  seokju::Body FR_calf (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0},
                        {0, 0, -0.25},
                        identity, 9);
  seokju::Body FR_foot (seokju::Body::Joint::Type::FIXED,
                        {0, 0, 0},
                        {0, 0, -0.25},
                        identity, 0);
  
  // set the parent body
  trunk.setParent({});
  FR_hip.setParent({&trunk});
  FR_thigh.setParent({&FR_hip});
  FR_calf.setParent({&FR_thigh});
  FR_foot.setParent({&FR_calf});

  FR_foot.computeKinematicsDownTheTree(gc);
  
  tar_pos = FR_foot.joint_.jointPosition_W;
  return tar_pos; /// replace this

}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_20233536_HPP_