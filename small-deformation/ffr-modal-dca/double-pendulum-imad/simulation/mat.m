function matrices = mat(q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t)
                     
%% Body 1
sigma1 = 2700;
U1  = 9e-4;
m21 = 1;
l1  = 0.545;
Y1  = 7.3e10;
II1 = 1.6875e-008;

binv = bodyinertia(1, sigma1, U1, m21, l1, Y1, II1, q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
matrices.body1_zeta11 = binv.body_zeta11;
matrices.body1_zeta12 = binv.body_zeta12;
matrices.body1_zeta13 = binv.body_zeta13;
matrices.body1_zeta21 = binv.body_zeta21;
matrices.body1_zeta22 = binv.body_zeta22;
matrices.body1_zeta23 = binv.body_zeta23;


matrices.Lambda1RR = binv.LambdaRR;
matrices.Lambda1RF = binv.LambdaRF;
matrices.Lambda1FR = binv.LambdaFR;
matrices.Lambda1FF = binv.LambdaFF;
matrices.beta1R    = binv.betaR; 
matrices.beta1F    = binv.betaF;
matrices.Gamma1R1  = binv.GammaR1;
matrices.Gamma1F1  = binv.GammaF1;
matrices.Gamma1R2  = binv.GammaR2;
matrices.Gamma1F2  = binv.GammaF2;
matrices.A2hat1    = binv.A2hat1;
matrices.phi2k1    = binv.phi2k1;
matrices.ExA1      = binv.ExA1;
matrices.S1        = binv.S1;
matrices.Xj01      = binv.Xj01;
%% Body 2
sigma2 = 2700;
U2  = 4e-4;
m22 = 3;
l2  = 0.675;
Y2  = 7.3e10;
II2 = 3.3333e-9;


binv = bodyinertia(2, sigma2, U2, m22, l2, Y2, II2, q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
matrices.body2_zeta11 = binv.body_zeta11;
matrices.body2_zeta12 = binv.body_zeta12;
matrices.body2_zeta13 = binv.body_zeta13;
matrices.body2_zeta21 = binv.body_zeta21;
matrices.body2_zeta22 = binv.body_zeta22;
matrices.body2_zeta23 = binv.body_zeta23;

matrices.Lambda2RR = binv.LambdaRR;
matrices.Lambda2RF = binv.LambdaRF;
matrices.Lambda2FR = binv.LambdaFR;
matrices.Lambda2FF = binv.LambdaFF;
matrices.beta2R    = binv.betaR; 
matrices.beta2F    = binv.betaF;
matrices.Gamma2R1  = binv.GammaR1;
matrices.Gamma2F1  = binv.GammaF1;
matrices.Gamma2R2  = binv.GammaR2;
matrices.Gamma2F2  = binv.GammaF2;
matrices.A2hat2    = binv.A2hat2;
matrices.phi2k2    = binv.phi2k2;
matrices.ExA2      = binv.ExA2;
matrices.S2        = binv.S2;
matrices.Xj02      = binv.Xj02;