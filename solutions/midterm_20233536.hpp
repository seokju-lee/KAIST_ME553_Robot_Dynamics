//
// Created by jemin on 23. 4. 20.
//

#ifndef ME553_2022_SOLUTIONS_MIDTERM_20233536_HPP_
#define ME553_2022_SOLUTIONS_MIDTERM_20233536_HPP_

namespace seokju {
class Body {
  public:
    struct Joint{
      /// variables
      Eigen::Vector3d jointPosition_W;
      Eigen::Vector3d jointAxis_W;
      Eigen::Matrix3d rotation_W;
      Eigen::Matrix2d mass;
      Eigen::MatrixXd p_j, a_j;

      /// robot definition
      Eigen::Vector3d jointAxis_P;
      Eigen::Vector3d jointPosition_P;
      enum class Type {
        FIXED,
        FLOATING,
        REVOLUTE,
        PRISMATC
      } type = Type::REVOLUTE;
      Eigen::Matrix3d jointRotation_P, jointMass, jointInertia;
    };

    // set the body for saving the defined value
    Body(const Joint::Type type,
         const Eigen::Vector3d& jointAxis_P, 
         const Eigen::Vector3d& jointPosition_P,
         const Eigen::Matrix3d& jointRotation_P,
         const Eigen::Matrix3d& jointMass,
         const Eigen::Matrix3d& jointInertia) {
      joint_.type = type;
      joint_.jointAxis_P = jointAxis_P;
      joint_.jointPosition_P = jointPosition_P;
      joint_.jointRotation_P = jointRotation_P;
      joint_.jointMass = jointMass;
      joint_.jointInertia = jointInertia;
    }

    void computeKinematicsDownTheTree(const Eigen::VectorXd& gc) {
      // initialize the jacobian
      joint_.p_j = Eigen::MatrixXd::Zero(3,2);
      joint_.a_j = Eigen::MatrixXd::Zero(3,2);

      // Prismatic joint type
      if (joint_.type == Joint::Type::PRISMATC) {
        // set the joint postion
        joint_.jointPosition_W[0] = gc[0];
        joint_.jointPosition_W[1] = 0;
        joint_.jointPosition_W[2] = 0;

        joint_.rotation_W << 1, 0, 0,
                             0, 1, 0,
                             0, 0, 1;

        for (int i = 0; i <= 2; i++) {
          joint_.p_j(i,0) = joint_.jointAxis_P[i];
        }
      }
      // Revolute joint type
      else if (joint_.type == Joint::Type::REVOLUTE) {
        // Rotation about y-axis
        if (joint_.jointAxis_P[0] == 0 && joint_.jointAxis_P[1] == 1 && joint_.jointAxis_P[2] == 0){
          joint_.jointRotation_P << cos(gc[1]), 0, sin(gc[1]),
                                       0,       1,      0,
                                   -sin(gc[1]), 0, cos(gc[1]);
        }
      }
      // compute the mass matrix using the projected newton-euler method (summation of the first joint)
      joint_.mass = joint_.p_j.transpose()*joint_.jointMass*joint_.p_j + joint_.a_j.transpose()*joint_.jointInertia*joint_.a_j;

      // using recursive method
      for (auto& par: parent_) { 
        par->computeKinematicsDownTheTree(gc);
        // compute the kinematics from the parent value
        joint_.rotation_W = par->joint_.rotation_W*joint_.jointRotation_P;
        joint_.jointPosition_W = joint_.rotation_W*joint_.jointPosition_P + par->joint_.jointPosition_W;

        // set the positional, angular jacobian
        joint_.p_j = par->joint_.p_j;
        joint_.a_j = par->joint_.a_j;

        for (int i = 0; i <= 2; i++) {
          joint_.p_j(i,1) = joint_.jointAxis_P.cross(joint_.jointPosition_W-par->joint_.jointPosition_W)[i];
          joint_.a_j(i,1) = joint_.jointAxis_P[i];
        }
        
        // compute the mass matrix using the projected newton-euler method (summation of the second joint)
        joint_.mass = par->joint_.mass + joint_.p_j.transpose()*joint_.jointMass*joint_.p_j + joint_.a_j.transpose()*joint_.jointInertia*joint_.a_j;
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

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  // set matrix
  Eigen::Matrix3d identity, slider_inertia, rod_inertia;
  Eigen::Matrix2d mass_matrix;
  
  // define the matrix (identity, inertia matrix)
  identity << 1, 0, 0,
              0, 1, 0,
              0, 0, 1;

  slider_inertia << 2, 0, 0,
                    0, 1, 0,
                    0, 0, 2;

  rod_inertia << 1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;
  
  // set the body
  seokju::Body slider (seokju::Body::Joint::Type::PRISMATC,
                      {1, 0, 0}, 
                      {gc[0], 0, 0},
                      identity, 2*identity, slider_inertia);
  seokju::Body rod (seokju::Body::Joint::Type::REVOLUTE,
                      {0, 1, 0}, 
                      {0, 0, 0.5},
                      identity, 5*identity, rod_inertia);
  
  // set the parent body
  slider.setParent({});
  rod.setParent({&slider});
  
  // compute the kinematics to obtain the mass matrix
  rod.computeKinematicsDownTheTree(gc);
  
  return rod.joint_.mass;
}

#endif //ME553_2022_SOLUTIONS_MIDTERM_20233536_HPP_
