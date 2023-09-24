//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2023_SOLUTIONS_FINALEXAM_20233536_HPP_
#define ME553_2023_SOLUTIONS_FINALEXAM_20233536_HPP_

#pragma once
#include <Eigen/Core>
#include <iostream>

#include <vector>

namespace seokju {
class Body {
  public:
    struct Joint{
      /// variables
      Eigen::Vector3d jointPosition_W, jointAxis_W;
      Eigen::Matrix3d rotation_W;
      Eigen::Vector3d jointAngularvel_P, jointLinearvel_P, jointLinearvel_W, jointAngularvel_W, jointCom_P;
      Eigen::Vector3d jointAngularacc_P, jointLinearacc_P, jointLinearacc_W, jointAngularacc_W;
      std::vector<Eigen::MatrixXd> inertia_C, Si, Sd, X_BP, X_BP_dot, MA;
      std::vector<Eigen::VectorXd> sA, sF, sP, sV, Ab;
      Eigen::MatrixXd jointI, motion_W, X, Xd, motion_W_d, Arti_M;
      Eigen::VectorXd motion_P, motion_P_d, v, a, nonlinear, ff, u_dot, Arti_b, ga;

      /// robot definition
      Eigen::Vector3d jointAxis_P;
      Eigen::Vector3d jointPosition_P;
      int i;
      enum class Type {
        FIXED,
        FLOATING,
        REVOLUTE,
        PRISMATIC
      } type = Type::REVOLUTE;
      Eigen::Matrix3d jointRotation_P, jointMass, jointInertia;

    };
    // set the body for saving the defined value
    Body(const Joint::Type type,
         const Eigen::Vector3d& jointAxis_P, 
         const Eigen::Vector3d& jointPosition_P,
         const Eigen::Matrix3d& jointRotation_P,
         const Eigen::Matrix3d& jointMass,
         const Eigen::Matrix3d& jointInertia,
         const Eigen::Vector3d& jointCom_P,
         int i) {
      joint_.type = type;
      joint_.jointAxis_P = jointAxis_P;
      joint_.jointPosition_P = jointPosition_P;
      joint_.jointRotation_P = jointRotation_P;
      joint_.jointMass = jointMass;
      joint_.jointInertia = jointInertia;
      joint_.jointCom_P = jointCom_P;
      joint_.i = i;
    }

    // function of the skew-symmetric matrix
    Eigen::Matrix3d skew(const Eigen::VectorXd& c) {
      Eigen::Matrix3d skew_symmetric;
      skew_symmetric << 0, -c[2], c[1],
                        c[2], 0, -c[0],
                        -c[1], c[0], 0;
      return skew_symmetric;
    }

    void computeKinematicsDownTheTree(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, std::vector<Eigen::VectorXd> sA, std::vector<Eigen::VectorXd> sV, std::vector<Eigen::VectorXd> sP, std::vector<Eigen::VectorXd> sF, std::vector<Eigen::MatrixXd> si, std::vector<Eigen::MatrixXd> inertia_C, std::vector<Eigen::MatrixXd> MA, std::vector<Eigen::VectorXd> Ab, std::vector<Eigen::MatrixXd> sd, std::vector<Eigen::MatrixXd> X_BP, std::vector<Eigen::MatrixXd> X_BP_dot) {
      joint_.a = Eigen::VectorXd::Zero(6,1);
      joint_.v = Eigen::VectorXd::Zero(6,1);
      joint_.jointI = Eigen::MatrixXd::Zero(6,6);
      joint_.Arti_M = Eigen::MatrixXd::Zero(6,6);
      joint_.Arti_b = Eigen::VectorXd::Zero(6,1);
      joint_.ff = Eigen::VectorXd::Zero(6,1);
      joint_.motion_P = Eigen::MatrixXd::Zero(6,1);
      joint_.motion_W = Eigen::MatrixXd::Zero(6,1);
      joint_.motion_P_d = Eigen::MatrixXd::Zero(6,1);
      joint_.motion_W_d = Eigen::MatrixXd::Zero(6,1);
      joint_.X = Eigen::MatrixXd::Zero(6,6);
      joint_.Xd = Eigen::MatrixXd::Zero(6,6);
      
      // Revolute joint type
      if (joint_.type == Joint::Type::REVOLUTE) {
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
        
        // set the velocity and accleration in body frame
        joint_.jointAngularvel_P = gv[joint_.i]*joint_.jointAxis_P;
        joint_.jointLinearvel_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointAngularacc_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearacc_P << 0, 0, 9.81;
        
        // set the motion subspace
        joint_.motion_P << 0, 0, 0, joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2];
      }
      // Prismatic joint type
      else if(joint_.type == Joint::Type::PRISMATIC) {
        // set the velocity and accleration in body frame
        joint_.jointAngularvel_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearvel_P = gv[joint_.i]*joint_.jointAxis_P;
        joint_.jointAngularacc_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearacc_P = Eigen::Vector3d::Zero(3,1);
        // in the case of the prismatic, joint position will be changed
        joint_.jointPosition_P = joint_.jointPosition_P + gc[joint_.i]*joint_.jointAxis_P;
        // set the motion subspace
        joint_.motion_P << joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2], 0, 0, 0;
      }
      // Fixed joint type
      else if(joint_.type == Joint::Type::FIXED) {
        joint_.jointPosition_W = Eigen::VectorXd::Zero(3,1);
        joint_.rotation_W = Eigen::Matrix3d::Identity(3,3);
        joint_.jointAngularvel_W = Eigen::VectorXd::Zero(3,1);
        joint_.jointLinearvel_W = Eigen::VectorXd::Zero(3,1);
        joint_.jointAngularacc_W = Eigen::VectorXd::Zero(3,1);
        joint_.jointLinearacc_W << 0, 0, 9.81;
        joint_.v << joint_.jointLinearvel_W, joint_.jointAngularvel_W;
        joint_.a << joint_.jointLinearacc_W, joint_.jointAngularacc_W;
        joint_.motion_W = Eigen::VectorXd::Zero(6,1);
        joint_.motion_W_d = Eigen::VectorXd::Zero(6,1);
        joint_.X = Eigen::MatrixXd::Identity(6,6);
        joint_.Xd = Eigen::MatrixXd::Zero(6,6);
        joint_.ff << joint_.jointMass*skew(joint_.jointAngularvel_W)*skew(joint_.jointAngularvel_W)*joint_.jointCom_P, skew(joint_.jointAngularvel_W)*(joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P))*joint_.jointAngularvel_W;
        joint_.jointI << joint_.jointMass, -joint_.jointMass*skew(joint_.jointCom_P),
                          joint_.jointMass*skew(joint_.jointCom_P), joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P);

        joint_.sA.push_back(joint_.a);
        joint_.sV.push_back(joint_.v);
        joint_.sF.push_back(joint_.ff);
        joint_.sP.push_back(joint_.jointPosition_W);
        joint_.Si.push_back(joint_.motion_W);
        joint_.Sd.push_back(joint_.motion_W_d);
        joint_.X_BP.push_back(joint_.X);
        joint_.X_BP_dot.push_back(joint_.Xd);
        joint_.inertia_C.push_back(joint_.jointI);
        joint_.MA.push_back(joint_.Arti_M);
        joint_.Ab.push_back(joint_.Arti_b);
      }
      
      // using recursive method
      for (auto& par: parent_) {
        par->computeKinematicsDownTheTree(gc, gv, sA, sV, sP, sF, si, inertia_C, MA, Ab, sd, X_BP, X_BP_dot);
        // compute the kinematics from the parent value
        joint_.rotation_W = par->joint_.rotation_W*joint_.jointRotation_P;
        joint_.jointPosition_W = par->joint_.rotation_W*joint_.jointPosition_P + par->joint_.jointPosition_W;
        joint_.jointCom_P = joint_.rotation_W*joint_.jointCom_P;
        joint_.jointAngularvel_W = par->joint_.jointAngularvel_W + par->joint_.rotation_W*joint_.jointAngularvel_P;
        joint_.jointAngularacc_W = par->joint_.jointAngularacc_W + skew(par->joint_.jointAngularvel_W)*par->joint_.rotation_W*joint_.jointAngularvel_P;
        
        // compute the velocity, accleration, and motion subspace in world coordinate using rotation
        if(joint_.type == Joint::Type::REVOLUTE) {
          joint_.jointLinearvel_W = par->joint_.jointLinearvel_W + par->joint_.jointAngularvel_W.cross(joint_.jointPosition_W - par->joint_.jointPosition_W);
          joint_.jointLinearacc_W = par->joint_.jointLinearacc_W + par->joint_.jointAngularacc_W.cross(joint_.jointPosition_W - par->joint_.jointPosition_W)
                                  + skew(par->joint_.jointAngularvel_W)*skew(par->joint_.jointAngularvel_W)*(joint_.jointPosition_W - par->joint_.jointPosition_W);
          joint_.motion_P.tail(3) = par->joint_.rotation_W*joint_.motion_P.tail(3);
          joint_.motion_P_d.tail(3) = skew(par->joint_.jointAngularvel_W)*joint_.rotation_W*joint_.jointAxis_P;

        }
        else if(joint_.type == Joint::Type::PRISMATIC) {
          joint_.jointLinearvel_W = par->joint_.jointLinearvel_W + par->joint_.jointAngularvel_W.cross(par->joint_.rotation_W*joint_.jointPosition_P) + par->joint_.rotation_W*gv[joint_.i]*joint_.jointAxis_P;
          joint_.jointLinearacc_W = par->joint_.jointLinearacc_W + par->joint_.jointAngularacc_W.cross(par->joint_.rotation_W*joint_.jointPosition_P)
                                  + skew(par->joint_.jointAngularvel_W)*skew(par->joint_.jointAngularvel_W)*(par->joint_.rotation_W*joint_.jointPosition_P)
                                  + 2*par->joint_.jointAngularvel_W.cross(par->joint_.rotation_W*gv[joint_.i]*joint_.jointAxis_P);
          joint_.motion_P.head(3) = par->joint_.rotation_W*joint_.motion_P.head(3);
          joint_.motion_P_d.head(3) = skew(par->joint_.jointAngularvel_W)*joint_.rotation_W*joint_.jointAxis_P;
        }
        // set the spatial velocity, accleration, spatial inertia, fictitious force
        joint_.X << Eigen::Matrix3d::Identity(3,3), Eigen::Matrix3d::Zero(3,3),
                    skew(joint_.jointPosition_W - par->joint_.jointPosition_W), Eigen::Matrix3d::Identity(3,3);
        joint_.Xd << Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Zero(3,3),
                     skew(joint_.jointLinearvel_W - par->joint_.jointLinearvel_W), Eigen::Matrix3d::Zero(3,3);
        joint_.v << joint_.jointLinearvel_W, joint_.jointAngularvel_W;
        joint_.a << joint_.jointLinearacc_W, joint_.jointAngularacc_W;
        joint_.jointI << joint_.jointMass, -joint_.jointMass*skew(joint_.jointCom_P),
                        joint_.jointMass*skew(joint_.jointCom_P), joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P);
        joint_.ff << joint_.jointMass*skew(joint_.jointAngularvel_W)*skew(joint_.jointAngularvel_W)*joint_.jointCom_P, skew(joint_.jointAngularvel_W)*(joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P))*joint_.jointAngularvel_W;
        
        // save the set for computed value
        joint_.sA.insert(joint_.sA.begin(), par->joint_.sA.begin(), par->joint_.sA.end());
        joint_.sA.push_back(joint_.a);
        joint_.sV.insert(joint_.sV.begin(), par->joint_.sV.begin(), par->joint_.sV.end());
        joint_.sV.push_back(joint_.v);
        joint_.sF.insert(joint_.sF.begin(), par->joint_.sF.begin(), par->joint_.sF.end());
        joint_.sF.push_back(joint_.ff);
        joint_.sP.insert(joint_.sP.begin(), par->joint_.sP.begin(), par->joint_.sP.end());
        joint_.sP.push_back(joint_.jointPosition_W);
        joint_.Si.insert(joint_.Si.begin(), par->joint_.Si.begin(), par->joint_.Si.end());
        joint_.Si.push_back(joint_.motion_P);
        joint_.Sd.insert(joint_.Sd.begin(), par->joint_.Sd.begin(), par->joint_.Sd.end());
        joint_.Sd.push_back(joint_.motion_P_d);
        joint_.X_BP.insert(joint_.X_BP.begin(), par->joint_.X_BP.begin(), par->joint_.X_BP.end());
        joint_.X_BP.push_back(joint_.X);
        joint_.X_BP_dot.insert(joint_.X_BP_dot.begin(), par->joint_.X_BP_dot.begin(), par->joint_.X_BP_dot.end());
        joint_.X_BP_dot.push_back(joint_.Xd);
        joint_.inertia_C.insert(joint_.inertia_C.begin(), par->joint_.inertia_C.begin(), par->joint_.inertia_C.end());
        joint_.inertia_C.push_back(joint_.jointI);
        joint_.MA.insert(joint_.MA.begin(), par->joint_.MA.begin(), par->joint_.MA.end());
        joint_.MA.push_back(joint_.Arti_M);
        joint_.Ab.insert(joint_.Ab.begin(), par->joint_.Ab.begin(), par->joint_.Ab.end());
        joint_.Ab.push_back(joint_.Arti_b);
      } 
    }  

    // set the articulated body variable
    void computeArticulatedBody(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
      Eigen::MatrixXd M_A, M_A_;
      Eigen::VectorXd b_A, b_A_, W;
      std::vector<Eigen::MatrixXd> M_A_leg;
      std::vector<Eigen::VectorXd> b_A_leg;

      for (int i=3; i>1; i--) {
        // set the initial value for recursion
        joint_.MA[3] = joint_.inertia_C[2];
        joint_.Ab[3] = joint_.sF[2];

        M_A_ = joint_.MA[i];
        b_A_ = joint_.Ab[i];
        W = joint_.sV[i-1];
        // compute articulated inertia
        M_A = joint_.X_BP[i]*M_A_*(-joint_.Si[i]*(joint_.Si[i].transpose()*M_A_*joint_.Si[i]).inverse()*(joint_.Si[i].transpose()*M_A_*joint_.X_BP[i].transpose()) + joint_.X_BP[i].transpose());
        // compute articulated fictitious force
        b_A = joint_.X_BP[i]*(M_A_*(joint_.Si[i]*(joint_.Si[i].transpose()*M_A_*joint_.Si[i]).inverse()*(gf.segment(i-1,1)-joint_.Si[i].transpose()*M_A_*(gv[i-1]*joint_.Sd[i]+joint_.X_BP_dot[i].transpose()*W)
              -joint_.Si[i].transpose()*b_A_)+joint_.Sd[i]*gv[i-1]+joint_.X_BP_dot[i].transpose()*W)+b_A_);
        // update the articulated value from children node
        joint_.MA[i-1] = joint_.inertia_C[i-1] + M_A;
        joint_.Ab[i-1] = joint_.sF[i-1] + b_A;
      }
    }

    // compute generalized acceleration
    void computeGA(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf){
      Eigen::VectorXd W, Wd;
      W = Eigen::VectorXd::Zero(6,1);
      Wd = Eigen::VectorXd::Zero(6,1);
      // initialize the generalized accleration
      joint_.ga = Eigen::VectorXd::Zero(3,1);
      // compute generalized acceleration
      for(int j=1; j<=3; j++){
        // update the value from parent
        W = joint_.sV[j-1];
        Wd = joint_.sA[j-1];
        // compute generalized accleration of each joint
        joint_.ga.segment(j-1,1) = (joint_.Si[j].transpose()*joint_.MA[j]*joint_.Si[j]).inverse()*(gf.segment(j-1,1)-joint_.Si[j].transpose()*joint_.MA[j]*(joint_.Sd[j]*gv[j-1]+joint_.X_BP_dot[j].transpose()*W+joint_.X_BP[j].transpose()*Wd)-joint_.Si[j].transpose()*joint_.Ab[j]);
        joint_.sA[j] = joint_.Si[j]*joint_.ga.segment(j-1,1)+joint_.Sd[j]*gv[j-1]+joint_.X_BP_dot[j].transpose()*W+joint_.X_BP[j].transpose()*Wd;
      }
    }

    // set the parent node
    void setParent(std::vector<Body*> parent){
      parent_ = parent;
    }

    //set the children node
    void setChildren(Body* child){
      children_.push_back(child);
      child->parent_.push_back(this);
    }

    Joint joint_;
    std::vector<Body*> children_, parent_;
    Body* parent__ = nullptr;
};
}

inline Eigen::Matrix3d skew1(const Eigen::Vector3d& c){
  Eigen::Matrix3d skew;
  skew << 0, -c[2], c[1],
          c[2], 0, -c[0],
          -c[1], c[0], 0;
  return skew;
}

// compute composite body for fixed joint
inline seokju::Body computecompositebody(const seokju::Body& parent, const seokju::Body& children) {
  seokju::Body composite = parent;
  // compute the composite rigid body mass, center of mass, and inertia
  composite.joint_.jointMass = parent.joint_.jointMass + children.joint_.jointMass;
  composite.joint_.jointCom_P = (parent.joint_.jointMass*parent.joint_.jointCom_P + children.joint_.jointMass*(children.joint_.jointPosition_P+children.joint_.jointCom_P))/composite.joint_.jointMass(0,0);
  composite.joint_.jointInertia = parent.joint_.jointInertia + children.joint_.jointInertia - parent.joint_.jointMass*skew1(parent.joint_.jointCom_P - composite.joint_.jointCom_P)*skew1(parent.joint_.jointCom_P - composite.joint_.jointCom_P)-children.joint_.jointMass*skew1(children.joint_.jointPosition_P + children.joint_.jointCom_P - composite.joint_.jointCom_P)*skew1(children.joint_.jointPosition_P + children.joint_.jointCom_P - composite.joint_.jointCom_P);
  
  return composite;
}

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
  // set the variable
  Eigen::Matrix3d identity, link_inertia;
  std::vector<Eigen::VectorXd> sA, sV, sP, sF, Ab;
  std::vector<Eigen::MatrixXd> si, ic, sd, X_BP, X_BP_dot, MA;
  
  identity << 1, 0, 0,
              0, 1, 0,
              0, 0, 1;

  link_inertia << 0.001, 0, 0,
                  0, 0.001, 0,
                  0, 0, 0.001;
  
  // set the body
  seokju::Body link1 (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, 0},
                      Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Zero(3,3), Eigen::Matrix3d::Zero(3,3), {0, 0, 0}, 0);

  seokju::Body link2 (seokju::Body::Joint::Type::REVOLUTE, {1, 0, 0}, {0, 0, 0.3},
                      identity, identity, link_inertia, {0, 0, 0.2}, 0);

  seokju::Body link3 (seokju::Body::Joint::Type::REVOLUTE, {1, 0, 0}, {0, 0, 0.3},
                      identity, identity, link_inertia, {0, 0, 0.2}, 1);

  seokju::Body link4 (seokju::Body::Joint::Type::PRISMATIC, {0, 1, 0}, {0, 0, 0.3},
                      identity, identity, link_inertia, {0, 0, 0.2}, 2);

  // set the parent
  link1.setParent({});
  link2.setParent({&link1});
  link3.setParent({&link2});
  link4.setParent({&link3});

  // compute kinematics down the tree
  link4.computeKinematicsDownTheTree(gc, gv, sA, sV, sP, sF, si, ic, MA, Ab, sd, X_BP, X_BP_dot);

  // compute articulated body up the tree
  link4.computeArticulatedBody(gv, gf);

  // compute generalized accleration
  link4.computeGA(gv, gf);

  return link4.joint_.ga;
}

#endif // ME553_2023_SOLUTIONS_FINALEXAM_20233536_HPP_