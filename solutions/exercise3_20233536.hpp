//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_
#define ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_

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
      std::vector<Eigen::MatrixXd> inertia_C, sXi_1, sXi, Si;
      Eigen::MatrixXd mass_matrix, jointI, Xi_1, Xi, motion_W;
      Eigen::VectorXd motion_P;

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
    Eigen::Matrix3d skew(const Eigen::Vector3d& c) {
      Eigen::Matrix3d skew_symmetric;
      skew_symmetric << 0, -c[2], c[1],
                        c[2], 0, -c[0],
                        -c[1], c[0], 0;
      return skew_symmetric;
    } 

    // set the all inertia of link for computing CRB
    void setforCRBA(const Eigen::VectorXd& gc) {
      joint_.motion_P = Eigen::MatrixXd::Zero(6,1);
      joint_.motion_W = Eigen::MatrixXd::Zero(6,6);
      joint_.Xi = Eigen::MatrixXd::Zero(6,6);
      joint_.Xi_1 = Eigen::MatrixXd::Zero(6,6);
      joint_.jointI = Eigen::MatrixXd::Zero(6,6);

      // Floating base type
      if (joint_.type == Joint::Type::FLOATING) {
        // set the joint postion
        joint_.jointPosition_P[0] = gc[0];
        joint_.jointPosition_P[1] = gc[1];
        joint_.jointPosition_P[2] = gc[2];

        // convert quaternion to rotation matrix
        joint_.jointRotation_P  << 1-2*pow(gc[5],2)-2*pow(gc[6],2), 2*gc[4]*gc[5]-2*gc[3]*gc[6], 2*gc[4]*gc[6]+2*gc[3]*gc[5],
                                   2*gc[4]*gc[5]+2*gc[3]*gc[6], 1-2*pow(gc[4],2)-2*pow(gc[6],2), 2*gc[5]*gc[6]-2*gc[3]*gc[4],
                                   2*gc[4]*gc[6]-2*gc[3]*gc[5], 2*gc[5]*gc[6]+2*gc[3]*gc[4], 1-2*pow(gc[4],2)-2*pow(gc[5],2);

        // set the motion space for base
        joint_.motion_W = Eigen::MatrixXd::Identity(6,6);
        joint_.Si.push_back(joint_.motion_W);
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
        // set the motion space for revolute joint
        joint_.motion_P << 0, 0, 0, joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2];
        joint_.Si.push_back(joint_.motion_P);
      }

      else if (joint_.type == Joint::Type::PRISMATIC) {
        // set the motion space for prismatic joint
        joint_.motion_P << joint_.jointAxis_P[0], joint_.jointAxis_P[1], joint_.jointAxis_P[2], 0, 0, 0;
        // set the joint position for prismatic joint using generalized coordinate
        joint_.jointPosition_P = gc[joint_.i]*joint_.jointAxis_P + joint_.jointPosition_P;
        joint_.Si.push_back(joint_.motion_P);
      }

      else if (joint_.type == Joint::Type::FIXED) {
        joint_.motion_P = Eigen::VectorXd::Zero(6,1);
        joint_.Si.push_back(joint_.motion_P);
      }

      // set the insertia and transformation matrix
      joint_.jointI << joint_.jointMass, -joint_.jointMass*skew(joint_.jointCom_P),
                      joint_.jointMass*skew(joint_.jointCom_P), joint_.jointInertia-joint_.jointMass*skew(joint_.jointCom_P)*skew(joint_.jointCom_P);
      
      joint_.Xi << joint_.jointRotation_P, Eigen::Matrix3d::Zero(3,3),
                    skew(joint_.jointPosition_P)*joint_.jointRotation_P, joint_.jointRotation_P;
      joint_.Xi_1 << joint_.jointRotation_P.transpose(), -joint_.jointRotation_P.transpose()*skew(joint_.jointPosition_P),
                     Eigen::Matrix3d::Zero(3,3), joint_.jointRotation_P.transpose();

      // save the inertia and transformation matrix set
      joint_.inertia_C.push_back(joint_.jointI);
      joint_.sXi_1.push_back(joint_.Xi_1);
      joint_.sXi.push_back(joint_.Xi);
      
      // save the inertia and transformation matrix each children node for CRBA
      for (auto& child: children_) {
        child->setforCRBA(gc);
        joint_.inertia_C.insert(joint_.inertia_C.end(), child->joint_.inertia_C.begin(), child->joint_.inertia_C.end());
        joint_.sXi_1.insert(joint_.sXi_1.end(), child->joint_.sXi_1.begin(), child->joint_.sXi_1.end());
        joint_.sXi.insert(joint_.sXi.end(), child->joint_.sXi.begin(), child->joint_.sXi.end());
        joint_.Si.insert(joint_.Si.end(), child->joint_.Si.begin(), child->joint_.Si.end());
      }
    }

    // compute using CRBA
    void computeCRBA(std::vector<Eigen::MatrixXd> Ic, std::vector<Eigen::MatrixXd> Xi_1, std::vector<Eigen::MatrixXd> Xi, std::vector<Eigen::MatrixXd> Si) {
      // initialize mass matrix to zero
      joint_.mass_matrix = Eigen::MatrixXd::Zero(18, 18);
      Eigen::MatrixXd F, M, M_, M__;
      int j, k, cnt;
      k = 17;
      cnt = 1;
      for (int i=0; i<=2; i++){
        for (int p=3; p<=5; p++){
          Xi[0](p,i) = 0;
          Xi_1[0](i,p) = 0;
        }
      }

      // compute CRB for four leg
      while (k>=5){
        // compute from children to root in body coordinate
        for (int i=k; i>=k-3; i--) {
          if (i!=k-3) {
            // transform the inertia matrix
            Ic[i-1] += Xi[i]*Ic[i]*Xi_1[i];
          }
          // calculate diagonal component
          F = Ic[i]*Si[i];
          M = Si[i].transpose()*F;
          if (i%4!=1) {
            joint_.mass_matrix(i+cnt,i+cnt) = M(0,0);
            j = i;
            // calculate off-diagonal component
            while (j!=k-3){
              F = Xi[j]*F;
              j = j - 1;
              M_ = F.transpose()*Si[j];
              joint_.mass_matrix(i+cnt,j+cnt) = M_(0,0);
              joint_.mass_matrix(j+cnt,i+cnt) = M_(0,0);
            }
            // compute the component of mass matrix about each whole leg to base
            F = Xi[0]*Xi[k-3]*F;
            M__ = F.transpose();
            for (int p=0; p<=5; p++) {
              joint_.mass_matrix(i+cnt, p) = M__(0,p);
              joint_.mass_matrix(p, i+cnt) = M__(0,p);
            }
          } 
        }
        k = k-4;
        cnt += 1;  
      }   
      // compute CRB for base (four leg to trunk)
      for (int i=3; i>=0; i--) {
        Ic[0] += Xi[2+4*i]*Ic[2+4*i]*Xi_1[2+4*i];
      }

      //compute CRB for base (IMU to trunk)
      Ic[0] += Xi[1]*Ic[1]*Xi_1[1];
      M = Xi[0]*Eigen::MatrixXd::Identity(6,6)*Ic[0]*Eigen::MatrixXd::Identity(6,6)*Xi_1[0];
      for (int i = 0; i<=5; i++) {
        for(int p = 0; p<=5; p++) {
          joint_.mass_matrix(i,p) = M(i,p);
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


/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  // set the inertia parameter
  Eigen::Matrix3d identity, trunk_inertia, imu_inertia, FR_hip_inertia, right_thigh_inertia, calf_inertia, foot_inertia, 
                  FL_hip_inertia, left_thigh_inertia, RR_hip_inertia, RL_hip_inertia;
  
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
  seokju::Body trunk (seokju::Body::Joint::Type::FLOATING,
                      {0, 0, 0}, 
                      {gc[0], gc[1], gc[2]},
                      identity, 9.041*identity, trunk_inertia, 
                      {0.008465, 0.004045, -0.000763}, 0);

  seokju::Body IMU (seokju::Body::Joint::Type::FIXED,
                      {0, 0, 0}, 
                      {0, 0, 0},
                      identity, 0.001*identity, imu_inertia, 
                      {0, 0, 0}, 1);

  seokju::Body FR_hip (seokju::Body::Joint::Type::REVOLUTE,
                      {1, 0, 0}, 
                      {0.2399, -0.051, 0},
                      identity, 1.993*identity, FR_hip_inertia, 
                      {-0.022191, -0.015144, -1.5e-05}, 7);

  seokju::Body FR_thigh (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0}, 
                        {0, -0.083, 0},
                        identity, 0.639*identity, right_thigh_inertia, 
                        {-0.005607, 0.003877, -0.048199}, 8);

  seokju::Body FR_calf (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0}, 
                        {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, 
                        {0.002781, 6.3e-05, -0.142518}, 9);

  seokju::Body FR_foot (seokju::Body::Joint::Type::FIXED,
                        {0, 0, 0}, 
                        {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, 
                        {0, 0, 0}, 0);

  seokju::Body FL_hip (seokju::Body::Joint::Type::REVOLUTE,
                      {1, 0, 0}, 
                      {0.2399, 0.051, 0},
                      identity, 1.993*identity, FL_hip_inertia,
                      {-0.022191, 0.015144, -1.5e-05}, 10);

  seokju::Body FL_thigh (seokju::Body::Joint::Type::PRISMATIC,
                        {0, 1, 0}, 
                        {0, 0.083, 0},
                        identity, 0.639*identity, left_thigh_inertia, 
                        {-0.005607, -0.003877, -0.048199}, 11);

  seokju::Body FL_calf (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0}, 
                        {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, 
                        {0.002781, 6.3e-05, -0.142518}, 12);

  seokju::Body FL_foot (seokju::Body::Joint::Type::FIXED,
                        {0, 0, 0}, 
                        {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, 
                        {0, 0, 0}, 0);

  seokju::Body RR_hip (seokju::Body::Joint::Type::REVOLUTE,
                      {1, 0, 0}, 
                      {-0.2399, -0.051, 0},
                      identity, 1.993*identity, RR_hip_inertia, 
                      {0.022191, -0.015144, -1.5e-05}, 13);

  seokju::Body RR_thigh (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0}, 
                        {0, -0.083, 0},
                        identity, 0.639*identity, right_thigh_inertia, 
                        {-0.005607, 0.003877, -0.048199}, 14);

  seokju::Body RR_calf (seokju::Body::Joint::Type::REVOLUTE,
                        {0, 1, 0}, 
                        {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, 
                        {0.002781, 6.3e-05, -0.142518}, 15);

  seokju::Body RR_foot (seokju::Body::Joint::Type::FIXED,
                        {0, 0, 0}, 
                        {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, 
                        {0, 0, 0}, 0);

  seokju::Body RL_hip (seokju::Body::Joint::Type::PRISMATIC,
                      {1, 0, 0}, 
                      {-0.2399, 0.051, 0},
                      identity, 1.993*identity, RL_hip_inertia, 
                      {0.022191, 0.015144, -1.5e-05}, 16);

  seokju::Body RL_thigh (seokju::Body::Joint::Type::PRISMATIC,
                        {0, 1, 0}, 
                        {0, 0.083, 0},
                        identity, 0.639*identity, left_thigh_inertia, 
                        {-0.005607, -0.003877, -0.048199}, 17);

  seokju::Body RL_calf (seokju::Body::Joint::Type::PRISMATIC,
                        {0, 1, 0}, 
                        {0, 0, -0.25},
                        identity, 0.207*identity, calf_inertia, 
                        {0.002781, 6.3e-05, -0.142518}, 18);

  seokju::Body RL_foot (seokju::Body::Joint::Type::FIXED,
                        {0, 0, 0}, 
                        {0, 0, -0.25},
                        identity, 0.06*identity, foot_inertia, 
                        {0, 0, 0}, 0);

  // set the children node
  trunk.setChildren({&IMU});
  trunk.setChildren({&FR_hip});
  trunk.setChildren({&FL_hip});
  trunk.setChildren({&RR_hip});
  trunk.setChildren({&RL_hip});
  FR_hip.setChildren({&FR_thigh});
  FR_thigh.setChildren({&FR_calf});
  FR_calf.setChildren({&FR_foot});
  FL_hip.setChildren({&FL_thigh});
  FL_thigh.setChildren({&FL_calf});
  FL_calf.setChildren({&FL_foot});
  RR_hip.setChildren({&RR_thigh});
  RR_thigh.setChildren({&RR_calf});
  RR_calf.setChildren({&RR_foot});
  RL_hip.setChildren({&RL_thigh});
  RL_thigh.setChildren({&RL_calf});
  RL_calf.setChildren({&RL_foot});

  // save the set for CRBA
  trunk.setforCRBA(gc);

  // computing mass matrix using CRBA
  trunk.computeCRBA(trunk.joint_.inertia_C, trunk.joint_.sXi_1, trunk.joint_.sXi, trunk.joint_.Si);

  return trunk.joint_.mass_matrix;
}

#endif // ME553_2023_SOLUTIONS_EXERCISE3_20233536_HPP_