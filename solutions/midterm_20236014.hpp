#ifndef ME553_2022_SOLUTIONS_MIDTERM_20236014_HPP_
#define ME553_2022_SOLUTIONS_MIDTERM_20236014_HPP_

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  // masses from URDF
  double m_sliderbar, m_slider, m_rod;
  m_sliderbar = 0.0;
  m_slider = 2.0;
  m_slider = m_slider + m_sliderbar;
  m_rod = 5.0;

  // inertia matrices from URDF
  Eigen::Matrix3d I_sliderbar, I_slider, I_rod;
  double ixx, ixy, ixz, iyy, iyz, izz;
  I_sliderbar = Eigen::Matrix3d::Identity();
  ixx=2.0; ixy=0.0; ixz=0.0; iyy=1.0; iyz=0.0; izz=2.0;
  I_slider << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
  I_rod = Eigen::Matrix3d::Identity();
  // composite inertia matrix according to Huygens-Steiner theorem
  I_slider = I_sliderbar + I_slider;

  // rotation axes from URDF
  Eigen::Vector3d axis_slider, axis_bar;
  axis_slider << 1, 0 ,0;
  axis_bar << 0, 1, 0;

  // angle axis to rotation matrix conversion
  Eigen::Matrix3d R_11p;
  Eigen::AngleAxisd angAx;
  angAx.angle() = gc[1]; angAx.axis() = axis_bar;
  R_11p = angAx.toRotationMatrix();

  Eigen::Vector3d r_rodCOM, r_slider2rodCOM;
  r_rodCOM << 0.0, 0.0, 0.5;

  // position vector from slider to rod center of mass
  r_slider2rodCOM << R_11p*r_rodCOM;

  // jacobian matrices for each body
  Eigen::MatrixXd Jp_slider(3,2), Jp_bar(3,2), Ja_slider(3,2), Ja_bar(3,2);
  Jp_slider << axis_slider, Eigen::Vector3d::Zero();
  Jp_bar << axis_slider, -r_slider2rodCOM.cross(axis_bar);
  Ja_slider << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero();
  Ja_bar << Eigen::Vector3d::Zero(), axis_bar;

  // mass matrix computed from PNE
  Eigen::Matrix2d M;
  M = Jp_slider.transpose()*m_slider*Jp_slider + Ja_slider.transpose()*I_slider*Ja_slider
          + Jp_bar.transpose()*m_rod*Jp_bar + Ja_bar.transpose()*I_rod*Ja_bar;
  //std::cout<<M<<std::endl;
  return M;
}

#endif //ME553_2022_SOLUTIONS_MIDTERM_20236014_HPP_
