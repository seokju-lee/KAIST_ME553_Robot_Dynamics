//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_
#define ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_

#pragma once
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
      Eigen::Vector3d jointAngularvel_P, jointLinearvel_P, jointLinearvel_W, jointAngularvel_W;

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
    // set the body for saving the defined value
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

    void computeKinematicsDownTheTree(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
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

        // set the joint linear and angular velocity
        joint_.jointLinearvel_W[0] = gv[0];
        joint_.jointLinearvel_W[1] = gv[1];
        joint_.jointLinearvel_W[2] = gv[2];
        joint_.jointAngularvel_W << gv[3], gv[4], gv[5];
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
        // set the angular velocity in body frame
        joint_.jointAngularvel_P = gv[joint_.i-1]*joint_.jointAxis_P;
      }
      
      // using recursive method
      for (auto& par: parent_) {
        par->computeKinematicsDownTheTree(gc, gv);
        // compute the kinematics from the parent value
        joint_.rotation_W = par->joint_.rotation_W*joint_.jointRotation_P;
        joint_.jointPosition_W = par->joint_.rotation_W*joint_.jointPosition_P + par->joint_.jointPosition_W;
        // compute the velocity using rotation 
        joint_.jointLinearvel_W = par->joint_.jointLinearvel_W + par->joint_.jointAngularvel_W.cross(joint_.jointPosition_W - par->joint_.jointPosition_W);
        joint_.jointAngularvel_W = par->joint_.jointAngularvel_W + par->joint_.rotation_W*joint_.jointAngularvel_P;
      }
      
    }  
    // set the parent node
    void setParent(std::vector<Body*> parent){
      parent_ = parent;
    }

    Joint joint_;
    std::vector<Body*> children_, parent_;
    Body* parent__ = nullptr;
};
}


/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return Eigen::MatrixXd::Ones(18,18);
}

#endif // ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_