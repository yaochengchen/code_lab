 #include <TVector3.h>
  
  const double MagneticInclination = -4.7155;
  const double PI = TMath::Pi();
  const double E_conversion_factor = 2.99792458E4;//1./sqrt(4*PI*Epsilon);
  
  const double antenna_1_theta = 89.327*TMath::DegToRad();// zenith, 0 means upward, 90 means horizontal
  const double antenna_2_theta = 88.891*TMath::DegToRad();
  const double antenna_3_theta = 88.613*TMath::DegToRad();
  const double antenna_4_theta = 90.275*TMath::DegToRad();
  
  const double antenna_1_phi = (90 - 118.282)*TMath::DegToRad();// 0 means geo-east, nagative south, positive north
  const double antenna_2_phi = (90 - 117.096)*TMath::DegToRad();
  const double antenna_3_phi = (90 - 115.572)*TMath::DegToRad();
  const double antenna_4_phi = (90 - 107.281)*TMath::DegToRad();
  
  TVector3 antenna_direction[4][3]; // 4 antennas, 0: east(x) 1: north(y) upward(z)
  //TVector3 antenna_upward[4];
  
  // coreas simulation
  TVector3 coreas_x, coreas_y, coreas_z;
  
  
  double antenna_theta[4] = {0.};
  double antenna_phi[4] = {0.};
  TVector3 E_theta_station[4];
  TVector3 E_phi_station[4];
  
  
