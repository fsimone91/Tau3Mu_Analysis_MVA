m3m[1.62,2.0];
sig_m0_A1_2017[1.777, 1.72, 1.84];
sig_sigma_A1_2017[0.02, 0.0, 0.1];
sig_alpha_A1_2017[ 1, -5., 5.]; 
sig_n_A1_2017[15, 0.0, 30.0]; 
sig_gaus_sigma_A1_2017[0.01,0.0,0.06];
cb_fraction_A1_2017[0.5,0,1];

t3m_sig_CBshape_A1_2017  = CBShape(m3m, sig_m0_A1_2017, sig_sigma_A1_2017, sig_alpha_A1_2017, sig_n_A1_2017);
t3m_sig_GSshape_A1_2017  = Gaussian(m3m,sig_m0_A1_2017,sig_gaus_sigma_A1_2017);


sig_m0_B1_2017[1.777, 1.72, 1.84];
sig_sigma_B1_2017[0.025, 0.001, 0.07];
sig_alpha_B1_2017[ 1, -20., 20.]; 
sig_n_B1_2017[2, 0.0, 5.0]; 
sig_gaus_sigma_B1_2017[0.025,0.001,0.1];
cb_fraction_B1_2017[0.5,0,1];

t3m_sig_CBshape_B1_2017  = CBShape(m3m, sig_m0_B1_2017, sig_sigma_B1_2017, sig_alpha_B1_2017, sig_n_B1_2017);
t3m_sig_GSshape_B1_2017  = Gaussian(m3m,sig_m0_B1_2017,sig_gaus_sigma_B1_2017);


sig_m0_C1_2017[1.777, 1.72, 1.84];
sig_sigma_C1_2017[0.02, 0.0, 0.05];
sig_alpha_C1_2017[ 1, -20., 20.]; 
sig_n_C1_2017[2, 0.0, 5.0]; 
sig_gaus_sigma_C1_2017[0.02,0.0,0.1];
cb_fraction_C1_2017[0.5,0,1];

t3m_sig_CBshape_C1_2017  = CBShape(m3m, sig_m0_C1_2017, sig_sigma_C1_2017, sig_alpha_C1_2017, sig_n_C1_2017);
t3m_sig_GSshape_C1_2017  = Gaussian(m3m,sig_m0_C1_2017,sig_gaus_sigma_C1_2017);


sig_m0_A2_2017[1.777, 1.72, 1.84];
sig_sigma_A2_2017[0.02, 0.0, 0.05];
sig_alpha_A2_2017[ 1, -20., 20.]; 
sig_n_A2_2017[2, 0.0, 5.0]; 
sig_gaus_sigma_A2_2017[0.02,0.0,0.1];
cb_fraction_A2_2017[0.5,0,1];

t3m_sig_CBshape_A2_2017  = CBShape(m3m, sig_m0_A2_2017, sig_sigma_A2_2017, sig_alpha_A2_2017, sig_n_A2_2017);
t3m_sig_GSshape_A2_2017  = Gaussian(m3m,sig_m0_A2_2017,sig_gaus_sigma_A2_2017);


sig_m0_B2_2017[1.777, 1.72, 1.84];
sig_sigma_B2_2017[0.02, 0.0, 0.05];
sig_alpha_B2_2017[ 1, -20., 20.]; 
sig_n_B2_2017[2, 0.0, 5.0]; 
sig_gaus_sigma_B2_2017[0.02,0.0,0.1];
cb_fraction_B2_2017[0.5,0,1];

t3m_sig_CBshape_B2_2017  = CBShape(m3m, sig_m0_B2_2017, sig_sigma_B2_2017, sig_alpha_B2_2017, sig_n_B2_2017);
t3m_sig_GSshape_B2_2017  = Gaussian(m3m,sig_m0_B2_2017,sig_gaus_sigma_B2_2017);


sig_m0_C2_2017[1.777, 1.72, 1.84];
sig_sigma_C2_2017[0.02, 0.0, 0.05];
sig_alpha_C2_2017[ 1, -20., 20.]; 
sig_n_C2_2017[2, 0.0, 5.0]; 
sig_gaus_sigma_C2_2017[0.02,0.0,0.1];
cb_fraction_C2_2017[0.5,0,1];

t3m_sig_CBshape_C2_2017  = CBShape(m3m, sig_m0_C2_2017, sig_sigma_C2_2017, sig_alpha_C2_2017, sig_n_C2_2017);
t3m_sig_GSshape_C2_2017  = Gaussian(m3m,sig_m0_C2_2017,sig_gaus_sigma_C2_2017);


bkg_exp_slope_A1_2017[-1.0,-6.0, 1.0];
bkg_exp_slope_B1_2017[-5.0,-6.0,-0.0];
bkg_exp_slope_C1_2017[-5.0,-6.0,-0.0];
bkg_exp_slope_A2_2017[-5.0,-6.0,-0.0];
bkg_exp_slope_B2_2017[-5.0,-6.0,-0.0];
bkg_exp_slope_C2_2017[-5.0,-6.0,-0.2];


bkg_exp_offset_A1_2017[0.0,-20.0,20.0];
bkg_exp_offset_A2_2017[0.0,-10.0,10.0];
bkg_exp_offset_B1_2017[0.0,-10.0,10.0];
bkg_exp_offset_B2_2017[0.0,-10.0,10.0];
bkg_exp_offset_C1_2017[0.0,-10.0,10.0];
bkg_exp_offset_C2_2017[0.0,-10.0,10.0];


//bkg_exp_shape_A1_2017 = RooExponential(m3m,bkg_exp_slope_A1_2017, bkg_exp_offset_A1_2017);
//bkg_exp_shape_A2_2017 = RooExponential(m3m,bkg_exp_slope_A2_2017, bkg_exp_offset_A2_2017);
//bkg_exp_shape_B1_2017 = RooExponential(m3m,bkg_exp_slope_B1_2017, bkg_exp_offset_B1_2017);
//bkg_exp_shape_B2_2017 = RooExponential(m3m,bkg_exp_slope_B2_2017, bkg_exp_offset_B2_2017);
//bkg_exp_shape_C1_2017 = RooExponential(m3m,bkg_exp_slope_C1_2017, bkg_exp_offset_C1_2017);
//bkg_exp_shape_C2_2017 = RooExponential(m3m,bkg_exp_slope_C2_2017, bkg_exp_offset_C2_2017);


//bkg_exp_slope_A1_2017[-2.0,-6.0,-0.001];
//bkg_exp_slope_A2_2017[-2.0,-6.0,-0.001];
//bkg_exp_slope_B1_2017[-3.0,-6.0,-0.001];
//bkg_exp_slope_B2_2017[-2.0,-6.0,-0.001];
//bkg_exp_slope_C1_2017[-2.0,-6.0,-0.001];
//bkg_exp_slope_C2_2017[-2.0,-6.0,-0.001];
//
//bkg_exp_offset_A1_2017[0.0,0.0,5.0];
//bkg_exp_offset_A2_2017[0.0,0.0,5.0];
//bkg_exp_offset_B1_2017[0.0,0.0,5.0];
//bkg_exp_offset_B2_2017[0.0,0.0,5.0];
//bkg_exp_offset_C1_2017[0.0,0.0,5.0];
//bkg_exp_offset_C2_2017[0.0,0.0,5.0];


sig_m0_A1_2018[1.777, 1.72, 1.84];
sig_sigma_A1_2018[0.02, 0.0, 0.1];
sig_alpha_A1_2018[ 1, -5., 5.]; 
sig_n_A1_2018[15, 0.0, 30.0]; 
sig_gaus_sigma_A1_2018[0.01,0.0,0.06];
cb_fraction_A1_2018[0.5,0,1];

t3m_sig_CBshape_A1_2018  = CBShape(m3m, sig_m0_A1_2018, sig_sigma_A1_2018, sig_alpha_A1_2018, sig_n_A1_2018);
t3m_sig_GSshape_A1_2018  = Gaussian(m3m,sig_m0_A1_2018,sig_gaus_sigma_A1_2018);


sig_m0_B1_2018[1.777, 1.72, 1.84];
sig_sigma_B1_2018[0.02, 0.0, 0.05];
sig_alpha_B1_2018[ 1, -20., 20.]; 
sig_n_B1_2018[2, 0.0, 5.0]; 
sig_gaus_sigma_B1_2018[0.02,0.0,0.1];
cb_fraction_B1_2018[0.5,0,1];

t3m_sig_CBshape_B1_2018  = CBShape(m3m, sig_m0_B1_2018, sig_sigma_B1_2018, sig_alpha_B1_2018, sig_n_B1_2018);
t3m_sig_GSshape_B1_2018  = Gaussian(m3m,sig_m0_B1_2018,sig_gaus_sigma_B1_2018);


sig_m0_C1_2018[1.777, 1.72, 1.84];
sig_sigma_C1_2018[0.02, 0.0, 0.05];
sig_alpha_C1_2018[ 1, -20., 20.]; 
sig_n_C1_2018[2, 0.0, 5.0]; 
sig_gaus_sigma_C1_2018[0.02,0.0,0.1];
cb_fraction_C1_2018[0.5,0,1];

t3m_sig_CBshape_C1_2018  = CBShape(m3m, sig_m0_C1_2018, sig_sigma_C1_2018, sig_alpha_C1_2018, sig_n_C1_2018);
t3m_sig_GSshape_C1_2018  = Gaussian(m3m,sig_m0_C1_2018,sig_gaus_sigma_C1_2018);


sig_m0_A2_2018[1.777, 1.72, 1.84];
sig_sigma_A2_2018[0.02, 0.0, 0.05];
sig_alpha_A2_2018[ 1, -20., 20.]; 
sig_n_A2_2018[2, 0.0, 5.0]; 
sig_gaus_sigma_A2_2018[0.02,0.0,0.1];
cb_fraction_A2_2018[0.5,0,1];

t3m_sig_CBshape_A2_2018  = CBShape(m3m, sig_m0_A2_2018, sig_sigma_A2_2018, sig_alpha_A2_2018, sig_n_A2_2018);
t3m_sig_GSshape_A2_2018  = Gaussian(m3m,sig_m0_A2_2018,sig_gaus_sigma_A2_2018);


sig_m0_B2_2018[1.777, 1.72, 1.84];
sig_sigma_B2_2018[0.02, 0.0, 0.05];
sig_alpha_B2_2018[ 1, -20., 20.]; 
sig_n_B2_2018[2, 0.0, 5.0]; 
sig_gaus_sigma_B2_2018[0.02,0.0,0.1];
cb_fraction_B2_2018[0.5,0,1];

t3m_sig_CBshape_B2_2018  = CBShape(m3m, sig_m0_B2_2018, sig_sigma_B2_2018, sig_alpha_B2_2018, sig_n_B2_2018);
t3m_sig_GSshape_B2_2018  = Gaussian(m3m,sig_m0_B2_2018,sig_gaus_sigma_B2_2018);


sig_m0_C2_2018[1.777, 1.72, 1.84];
sig_sigma_C2_2018[0.02, 0.0, 0.05];
sig_alpha_C2_2018[ 1, -20., 20.]; 
sig_n_C2_2018[2, 0.0, 5.0]; 
sig_gaus_sigma_C2_2018[0.02,0.0,0.1];
cb_fraction_C2_2018[0.5,0,1];

t3m_sig_CBshape_C2_2018  = CBShape(m3m, sig_m0_C2_2018, sig_sigma_C2_2018, sig_alpha_C2_2018, sig_n_C2_2018);
t3m_sig_GSshape_C2_2018  = Gaussian(m3m,sig_m0_C2_2018,sig_gaus_sigma_C2_2018);


bkg_exp_slope_A1_2018[-1.0,-6.0, 1.0];
bkg_exp_slope_B1_2018[-5.0,-6.0,-0.0];
bkg_exp_slope_C1_2018[-5.0,-6.0,-0.0];
bkg_exp_slope_A2_2018[-5.0,-6.0,-0.0];
bkg_exp_slope_B2_2018[-5.0,-6.0,-0.0];
bkg_exp_slope_C2_2018[-5.0,-6.0,-0.2];


bkg_exp_offset_A1_2018[0.0,-20.0,20.0];
bkg_exp_offset_A2_2018[0.0,-10.0,10.0];
bkg_exp_offset_B1_2018[0.0,-10.0,10.0];
bkg_exp_offset_B2_2018[0.0,-10.0,10.0];
bkg_exp_offset_C1_2018[0.0,-10.0,10.0];
bkg_exp_offset_C2_2018[0.0,-10.0,10.0];


//bkg_exp_shape_A1_2018 = RooExponential(m3m,bkg_exp_slope_A1_2018, bkg_exp_offset_A1_2018);
//bkg_exp_shape_A2_2018 = RooExponential(m3m,bkg_exp_slope_A2_2018, bkg_exp_offset_A2_2018);
//bkg_exp_shape_B1_2018 = RooExponential(m3m,bkg_exp_slope_B1_2018, bkg_exp_offset_B1_2018);
//bkg_exp_shape_B2_2018 = RooExponential(m3m,bkg_exp_slope_B2_2018, bkg_exp_offset_B2_2018);
//bkg_exp_shape_C1_2018 = RooExponential(m3m,bkg_exp_slope_C1_2018, bkg_exp_offset_C1_2018);
//bkg_exp_shape_C2_2018 = RooExponential(m3m,bkg_exp_slope_C2_2018, bkg_exp_offset_C2_2018);


//bkg_exp_slope_A1_2018[-2.0,-6.0,-0.001];
//bkg_exp_slope_A2_2018[-2.0,-6.0,-0.001];
//bkg_exp_slope_B1_2018[-3.0,-6.0,-0.001];
//bkg_exp_slope_B2_2018[-2.0,-6.0,-0.001];
//bkg_exp_slope_C1_2018[-2.0,-6.0,-0.001];
//bkg_exp_slope_C2_2018[-2.0,-6.0,-0.001];
//
//bkg_exp_offset_A1_2018[0.0,0.0,5.0];
//bkg_exp_offset_A2_2018[0.0,0.0,5.0];
//bkg_exp_offset_B1_2018[0.0,0.0,5.0];
//bkg_exp_offset_B2_2018[0.0,0.0,5.0];
//bkg_exp_offset_C1_2018[0.0,0.0,5.0];
//bkg_exp_offset_C2_2018[0.0,0.0,5.0];

sqrtS[13000., 13000., 13000.]

//bkg_exp_offset[0.0,-10.0,10.0];
//bkg_exp_shape  = RooExponential(m3m,bkg_exp_slope);



