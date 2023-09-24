//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2023_SOLUTIONS_EXERCISE5_20233536_HPP_
#define ME553_2023_SOLUTIONS_EXERCISE5_20233536_HPP_

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
      joint_.motion_W = Eigen::MatrixXd::Zero(6,6);
      joint_.motion_P_d = Eigen::MatrixXd::Zero(6,1);
      joint_.motion_W_d = Eigen::MatrixXd::Zero(6,6);
      joint_.X = Eigen::MatrixXd::Zero(6,6);
      joint_.Xd = Eigen::MatrixXd::Zero(6,6);
      
      // Floating base type
      if (joint_.type == Joint::Type::FLOATING) {
        // set the joint postion
        joint_.jointPosition_W = gc.head(3);
        // convert quaternion to rotation matrix
        joint_.rotation_W << 1-2*pow(gc[5],2)-2*pow(gc[6],2), 2*gc[4]*gc[5]-2*gc[3]*gc[6], 2*gc[4]*gc[6]+2*gc[3]*gc[5],
                              2*gc[4]*gc[5]+2*gc[3]*gc[6], 1-2*pow(gc[4],2)-2*pow(gc[6],2), 2*gc[5]*gc[6]-2*gc[3]*gc[4],
                              2*gc[4]*gc[6]-2*gc[3]*gc[5], 2*gc[5]*gc[6]+2*gc[3]*gc[4], 1-2*pow(gc[4],2)-2*pow(gc[5],2);
        // set the joint linear and angular velocity/acceleration
        joint_.jointLinearvel_W = gv.head(3);
        joint_.jointAngularvel_W = gv.segment(3,3);
        joint_.jointAngularacc_W = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearacc_W << 0, 0, 9.81;
        // compute the postion of the center of the mass using rotation
        joint_.jointCom_P = joint_.rotation_W*joint_.jointCom_P;

        // set the spatial accleration, fictitious force, and spatial inertia
        joint_.v << joint_.jointLinearvel_W, joint_.jointAngularvel_W;
        joint_.a << joint_.jointLinearacc_W, joint_.jointAngularacc_W;
        joint_.ff << joint_.jointMass*skew(joint_.jointAngularvel_W)*skew(joint_.jointAngularvel_W)*joint_.jointCom_P, skew(joint_.jointAngularvel_W)*(joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P))*joint_.jointAngularvel_W;
        joint_.jointI << joint_.jointMass, -joint_.jointMass*skew(joint_.jointCom_P),
                        joint_.jointMass*skew(joint_.jointCom_P), joint_.rotation_W*joint_.jointInertia*joint_.rotation_W.transpose()-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P);
        // set the motion subspace for floating base
        joint_.motion_W = Eigen::MatrixXd::Identity(6,6);
        joint_.motion_W_d = Eigen::MatrixXd::Zero(6,6);
        // set the transformation matrix for computing ABA
        joint_.X = Eigen::MatrixXd::Identity(6,6);
        joint_.Xd = Eigen::MatrixXd::Zero(6,6);

        // save the set for computing ABA
        // three legs case without first leg
        if (sA.size()!=0){
          joint_.sA.assign(sA.begin(), sA.end());
          joint_.sV.assign(sV.begin(), sV.end());
          joint_.sF.assign(sF.begin(), sF.end());
          joint_.sP.assign(sP.begin(), sP.end());
          joint_.Si.assign(si.begin(), si.end());
          joint_.Sd.assign(sd.begin(), sd.end());
          joint_.X_BP.assign(X_BP.begin(), X_BP.end());
          joint_.X_BP_dot.assign(X_BP_dot.begin(), X_BP_dot.end());
          joint_.inertia_C.assign(inertia_C.begin(), inertia_C.end());
          joint_.MA.assign(MA.begin(), MA.end());
          joint_.Ab.assign(Ab.begin(), Ab.end());
        }
        // when the compute down the tree for first leg
        else {
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
        // set the velocity and accleration in body frame
        joint_.jointAngularvel_P = gv[joint_.i-1]*joint_.jointAxis_P;
        joint_.jointLinearvel_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointAngularacc_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearacc_P = Eigen::Vector3d::Zero(3,1);
        // set the motion subspace
        joint_.motion_P << 0, 0, 0, joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2];
      }
      else if(joint_.type == Joint::Type::PRISMATIC) {
        // set the velocity and accleration in body frame
        joint_.jointAngularvel_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearvel_P = gv[joint_.i-1]*joint_.jointAxis_P;
        joint_.jointAngularacc_P = Eigen::Vector3d::Zero(3,1);
        joint_.jointLinearacc_P = Eigen::Vector3d::Zero(3,1);
        // in the case of the prismatic, joint position will be changed
        joint_.jointPosition_P = joint_.jointPosition_P + gc[joint_.i]*joint_.jointAxis_P;
        // set the motion subspace
        joint_.motion_P << joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2], 0, 0, 0;
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
          joint_.jointLinearvel_W = par->joint_.jointLinearvel_W + par->joint_.jointAngularvel_W.cross(par->joint_.rotation_W*joint_.jointPosition_P) + par->joint_.rotation_W*gv[joint_.i-1]*joint_.jointAxis_P;
          joint_.jointLinearacc_W = par->joint_.jointLinearacc_W + par->joint_.jointAngularacc_W.cross(par->joint_.rotation_W*joint_.jointPosition_P)
                                  + skew(par->joint_.jointAngularvel_W)*skew(par->joint_.jointAngularvel_W)*(par->joint_.rotation_W*joint_.jointPosition_P)
                                  + 2*par->joint_.jointAngularvel_W.cross(par->joint_.rotation_W*gv[joint_.i-1]*joint_.jointAxis_P);
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

    // compute nonlinearities using backward recursive
    void computebackrecursive() {
      // set the wrench and wrench of the base
      Eigen::VectorXd wrench, wrench_base;
      // set the wrench of whole leg for computing wrenche of the base
      std::vector<Eigen::VectorXd> wrench_leg;
      // initialize the wrench
      joint_.nonlinear.setZero(18);
      wrench_base.setZero(6);

      // for leg
      for (int i=3; i>=0; i--) {
        wrench.setZero(6);
        // each leg case
        for(int j=3*i+3; j>=3*i+1; j--) {
          if (j+1==3*i+4){
            wrench.tail(3) += skew(-joint_.sP[j])*wrench.head(3);
          }
          else {
            wrench.tail(3) += skew(joint_.sP[j+1]-joint_.sP[j])*wrench.head(3);
          }
          wrench += joint_.sA[j] + joint_.sF[j];
          joint_.nonlinear[j+5] = (joint_.Si[j].transpose()*wrench)[0];
        }
        wrench_base += wrench;
        wrench_leg.push_back(wrench);
      }
      // for base
      for(int i=3; i>=0; i--) {
        wrench_base.tail(3) += skew(joint_.sP[3*i+1]-joint_.sP[0])*wrench_leg[3-i].head(3);
      }
      wrench_base += joint_.sA[0] + joint_.sF[0];
      joint_.nonlinear.head(6) = joint_.Si[0].transpose()*wrench_base;
    }

    // set the articulated body variable
    void computeArticulatedBody(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {
      Eigen::MatrixXd M_A, M_A_;
      Eigen::VectorXd b_A, b_A_, W;
      std::vector<Eigen::MatrixXd> M_A_leg;
      std::vector<Eigen::VectorXd> b_A_leg;
      // for leg
      for (int i=3; i>=0; i--) {
        joint_.MA[3*i+3] = joint_.inertia_C[3*i+3];
        joint_.Ab[3*i+3] = joint_.sF[3*i+3];
        // each leg case
        for(int j=3*i+3; j>=3*i+1; j--) {
          M_A_ = joint_.MA[j];
          b_A_ = joint_.Ab[j];
          if(j%3 == 1) {
            W = joint_.sV[0];
          }
          else {
            W = joint_.sV[j-1];
          }
          // compute articulated inertia
          M_A = joint_.X_BP[j]*M_A_*(-joint_.Si[j]*(joint_.Si[j].transpose()*M_A_*joint_.Si[j]).inverse()*(joint_.Si[j].transpose()*M_A_*joint_.X_BP[j].transpose()) + joint_.X_BP[j].transpose());
          // compute articulated fictitious force
          b_A = joint_.X_BP[j]*(M_A_*(joint_.Si[j]*(joint_.Si[j].transpose()*M_A_*joint_.Si[j]).inverse()*(gf.segment(j+5,1)-joint_.Si[j].transpose()*M_A_*(gv[j+5]*joint_.Sd[j]+joint_.X_BP_dot[j].transpose()*W)
                -joint_.Si[j].transpose()*b_A_)+joint_.Sd[j]*gv[j+5]+joint_.X_BP_dot[j].transpose()*W)+b_A_);
          // update the articulated value from children node
          joint_.MA[j-1] = joint_.inertia_C[j-1] + M_A;
          joint_.Ab[j-1] = joint_.sF[j-1] + b_A;
        }
        M_A_leg.push_back(M_A);
        b_A_leg.push_back(b_A);
      }
      // for base
      joint_.MA[0] = joint_.inertia_C[0];
      joint_.Ab[0] = joint_.sF[0];
      // for base from all legs
      for (int i=0; i<=3; i++) {
        joint_.MA[0] += M_A_leg[i];
        joint_.Ab[0] += b_A_leg[i];
      }
      
    }

    // compute generalized acceleration
    void computeGA(const Eigen::VectorXd& gv, const Eigen::VectorXd& gf){
      Eigen::VectorXd W, Wd;
      W = Eigen::VectorXd::Zero(6,1);
      Wd = Eigen::VectorXd::Zero(6,1);
      // initialize the generalized accleration
      joint_.ga = Eigen::VectorXd::Zero(18,1);
      // compute generalized accleration of base
      joint_.ga.head(6) = (joint_.Si[0].transpose()*joint_.MA[0]*joint_.Si[0]).inverse()*(gf.head(6)-joint_.Si[0].transpose()*joint_.MA[0]*(joint_.Sd[0]*gv.head(6)+joint_.X_BP[0]*joint_.sA[0])-joint_.Si[0].transpose()*joint_.Ab[0]);
      joint_.sA[0] = joint_.Si[0]*joint_.ga.head(6)+joint_.X_BP[0].transpose()*joint_.sA[0];
      // for leg
      for (int i=0; i<=3; i++){
        // each leg case
        for(int j=3*i+1; j<=3*i+3; j++){
          if(j%3==1) {
            // first joint of leg (Because the direction of the computing is downward)
            Wd = joint_.sA[0];
            W = joint_.sV[0];
          }
          else {
            // update the value from parent
            W = joint_.sV[j-1];
            Wd = joint_.sA[j-1];
          }
          // compute generalized accleration of each joint
          joint_.ga.segment(j+5,1) = (joint_.Si[j].transpose()*joint_.MA[j]*joint_.Si[j]).inverse()*(gf.segment(j+5,1)-joint_.Si[j].transpose()*joint_.MA[j]*(joint_.Sd[j]*gv[j+5]+joint_.X_BP_dot[j].transpose()*W+joint_.X_BP[j].transpose()*Wd)-joint_.Si[j].transpose()*joint_.Ab[j]);
          joint_.sA[j] = joint_.Si[j]*joint_.ga.segment(j+5,1)+joint_.Sd[j]*gv[j+5]+joint_.X_BP_dot[j].transpose()*W+joint_.X_BP[j].transpose()*Wd;
        }
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
  // set the inertia parameter
  Eigen::Matrix3d identity, trunk_inertia, imu_inertia, FR_hip_inertia, right_thigh_inertia, calf_inertia, foot_inertia, 
                  FL_hip_inertia, left_thigh_inertia, RR_hip_inertia, RL_hip_inertia;

  std::vector<Eigen::VectorXd> sA, sV, sP, sF, Ab;
  std::vector<Eigen::MatrixXd> si, ic, sd, X_BP, X_BP_dot, MA;
  
  identity << 1, 0, 0,
              0, 1, 0,
              0, 0, 1;

  trunk_inertia << 0.033260231, -0.000451628, 0.000487603,
                  -0.000451628,  0.16117211,  4.8356e-05,
                   0.000487603,  4.8356e-05,  0.17460442;

  imu_inertia << 0.0001, 0, 0,
                 0, 0.000001, 0,
                 0, 0, 0.0001;

  FR_hip_inertia << 0.002903894, 7.185e-05, -1.262e-06,
                    7.185e-05,   0.004907517, 1.75e-06,
                    -1.262e-06,  1.75e-06, 0.005586944;

  right_thigh_inertia << 0.005666803, -3.597e-06, 0.000491446,
                    -3.597e-06,   0.005847229, -1.0086e-05,
                    0.000491446,  -1.0086e-05, 0.000369811;

  calf_inertia << 0.006341369, -3e-09, -8.7951e-05,
                     -3e-09, 0.006355157, -1.336e-06,
                     -8.7951e-05, -1.336e-06, 3.9188e-05;

  foot_inertia << 1.6854e-05, 0.0, 0.0,
                     0.0, 1.6854e-05, 0.0,
                     0.0, 0.0, 1.6854e-05;

  FL_hip_inertia << 0.002903894, -7.185e-05, -1.262e-06,
                    -7.185e-05,   0.004907517, -1.75e-06,
                    -1.262e-06,  -1.75e-06, 0.005586944;

  left_thigh_inertia << 0.005666803, 3.597e-06, 0.000491446,
                    3.597e-06,   0.005847229, 1.0086e-05,
                    0.000491446,  1.0086e-05, 0.000369811;

  RR_hip_inertia << 0.002903894, -7.185e-05, 1.262e-06,
                    -7.185e-05,   0.004907517, 1.75e-06,
                    1.262e-06,  1.75e-06, 0.005586944;

  RL_hip_inertia << 0.002903894, 7.185e-05, 1.262e-06,
                    7.185e-05,   0.004907517, -1.75e-06,
                    1.262e-06,  -1.75e-06, 0.005586944;

  // set the body
  seokju::Body trunk (seokju::Body::Joint::Type::FLOATING, {0, 0, 0}, {gc[0], gc[1], gc[2]},
                      identity, 9.041*identity, trunk_inertia, {0.008465, 0.004045, -0.000763}, 0);

  seokju::Body IMU (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, 0},
                    identity, 0.001*identity, imu_inertia, {0, 0, 0}, 1);

  seokju::Body trunk_IMU (computecompositebody(trunk, IMU));

  seokju::Body FR_hip (seokju::Body::Joint::Type::REVOLUTE, {1, 0, 0}, {0.2399, -0.051, 0},
                      identity, 1.993*identity, FR_hip_inertia, {-0.022191, -0.015144, -1.5e-05}, 7);

  seokju::Body FR_thigh (seokju::Body::Joint::Type::REVOLUTE, {0, 1, 0}, {0, -0.083, 0},
                        identity, 0.639*identity, right_thigh_inertia, {-0.005607, 0.003877, -0.048199}, 8);

  seokju::Body FR_calf (seokju::Body::Joint::Type::REVOLUTE, {0, 1, 0}, {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, {0.002781, 6.3e-05, -0.142518}, 9);

  seokju::Body FR_foot (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, {0, 0, 0}, 0);
  
  seokju::Body FR_calf_foot (computecompositebody(FR_calf, FR_foot));

  seokju::Body FL_hip (seokju::Body::Joint::Type::REVOLUTE, {1, 0, 0}, {0.2399, 0.051, 0},
                      identity, 1.993*identity, FL_hip_inertia, {-0.022191, 0.015144, -1.5e-05}, 10);

  seokju::Body FL_thigh (seokju::Body::Joint::Type::PRISMATIC, {0, 1, 0}, {0, 0.083, 0},
                        identity, 0.639*identity, left_thigh_inertia, {-0.005607, -0.003877, -0.048199}, 11);

  seokju::Body FL_calf (seokju::Body::Joint::Type::REVOLUTE, {0, 1, 0}, {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, {0.002781, 6.3e-05, -0.142518}, 12);

  seokju::Body FL_foot (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, {0, 0, 0}, 0);

  seokju::Body FL_calf_foot (computecompositebody(FL_calf, FL_foot));

  seokju::Body RR_hip (seokju::Body::Joint::Type::REVOLUTE, {1, 0, 0}, {-0.2399, -0.051, 0},
                      identity, 1.993*identity, RR_hip_inertia, {0.022191, -0.015144, -1.5e-05}, 13);

  seokju::Body RR_thigh (seokju::Body::Joint::Type::REVOLUTE, {0, 1, 0}, {0, -0.083, 0},
                        identity, 0.639*identity, right_thigh_inertia, {-0.005607, 0.003877, -0.048199}, 14);

  seokju::Body RR_calf (seokju::Body::Joint::Type::REVOLUTE, {0, 1, 0}, {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, {0.002781, 6.3e-05, -0.142518}, 15);

  seokju::Body RR_foot (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, {0, 0, 0}, 0);

  seokju::Body RR_calf_foot (computecompositebody(RR_calf, RR_foot));

  seokju::Body RL_hip (seokju::Body::Joint::Type::PRISMATIC, {1, 0, 0}, {-0.2399, 0.051, 0},
                      identity, 1.993*identity, RL_hip_inertia, {0.022191, 0.015144, -1.5e-05}, 16);

  seokju::Body RL_thigh (seokju::Body::Joint::Type::PRISMATIC, {0, 1, 0}, 
                        {0, 0.083, 0}, identity, 0.639*identity, left_thigh_inertia, {-0.005607, -0.003877, -0.048199}, 17);

  seokju::Body RL_calf (seokju::Body::Joint::Type::PRISMATIC, {0, 1, 0}, {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, {0.002781, 6.3e-05, -0.142518}, 18);

  seokju::Body RL_foot (seokju::Body::Joint::Type::FIXED, {0, 0, 0}, {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, {0, 0, 0}, 0);

  seokju::Body RL_calf_foot (computecompositebody(RL_calf, RL_foot));

  // set the parent node
  trunk_IMU.setParent({});
  FR_hip.setParent({&trunk_IMU});
  FR_thigh.setParent({&FR_hip});
  FR_calf_foot.setParent({&FR_thigh});
  FL_hip.setParent({&trunk_IMU});
  FL_thigh.setParent({&FL_hip});
  FL_calf_foot.setParent({&FL_thigh});
  RR_hip.setParent({&trunk_IMU});
  RR_thigh.setParent({&RR_hip});
  RR_calf_foot.setParent({&RR_thigh});
  RL_hip.setParent({&trunk_IMU});
  RL_thigh.setParent({&RL_hip});
  RL_calf_foot.setParent({&RL_thigh});
  
  // compute kinematics down the tree
  FR_calf_foot.computeKinematicsDownTheTree(gc, gv, sA, sV, sP, sF, si, ic, MA, Ab, sd, X_BP, X_BP_dot);
  FL_calf_foot.computeKinematicsDownTheTree(gc, gv, FR_calf_foot.joint_.sA, FR_calf_foot.joint_.sV, FR_calf_foot.joint_.sP, FR_calf_foot.joint_.sF, FR_calf_foot.joint_.Si, FR_calf_foot.joint_.inertia_C, FR_calf_foot.joint_.MA, FR_calf_foot.joint_.Ab, FR_calf_foot.joint_.Sd, FR_calf_foot.joint_.X_BP, FR_calf_foot.joint_.X_BP_dot);
  RR_calf_foot.computeKinematicsDownTheTree(gc, gv, FL_calf_foot.joint_.sA, FL_calf_foot.joint_.sV, FL_calf_foot.joint_.sP, FL_calf_foot.joint_.sF, FL_calf_foot.joint_.Si, FL_calf_foot.joint_.inertia_C, FL_calf_foot.joint_.MA, FL_calf_foot.joint_.Ab, FL_calf_foot.joint_.Sd, FL_calf_foot.joint_.X_BP, FL_calf_foot.joint_.X_BP_dot);
  RL_calf_foot.computeKinematicsDownTheTree(gc, gv, RR_calf_foot.joint_.sA, RR_calf_foot.joint_.sV, RR_calf_foot.joint_.sP, RR_calf_foot.joint_.sF, RR_calf_foot.joint_.Si, RR_calf_foot.joint_.inertia_C, RR_calf_foot.joint_.MA, RR_calf_foot.joint_.Ab, RR_calf_foot.joint_.Sd, RR_calf_foot.joint_.X_BP, RR_calf_foot.joint_.X_BP_dot);
  
  // compute articulated body up the tree
  RL_calf_foot.computeArticulatedBody(gv, gf);

  // compute generalized acceration down the tree
  RL_calf_foot.computeGA(gv, gf);

  return RL_calf_foot.joint_.ga;
}

#endif // ME553_2023_SOLUTIONS_EXERCISE5_20233536_HPP_