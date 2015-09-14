function binv = bodyinertia(i, sigma, U, m, l, Y, II, q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t)
                    
% Coeff of shape function
c1 = 1.5;
c2 = 0.5;
%% Body 1
if (i == 1)
   
r1      = y1;
r2      = y2;
rd1     = yd1;
rd2     = yd2;
theta   = q1;
omega3  = qd1;

Lflex  =  LFlex(c1,c2,l,sigma,U,r1,r2,rd1,rd2,omega3,theta);
Lpoint =  LPoint(c1,c2,l,m,r1,r2,rd1,rd2,omega3,theta);
gamma  =  Gamma(c1,c2,l,r1,r2,theta);
beta   =  Beta(c1,c2,l,Y,II,sigma,U,m,r1,r2,theta);

A2hat = [ 0;
          cos(theta) * (-omega3 ^ 2 * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) - 0.2e1 * omega3 * (c1 - c2) * rd2) - sin(theta) * (-omega3 ^ 2 * (c1 - c2) * r2 + 0.2e1 * omega3 * (rd1 - r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) * rd2));
          sin(theta) * (-omega3 ^ 2 * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) - 0.2e1 * omega3 * (c1 - c2) * rd2) + cos(theta) * (-omega3 ^ 2 * (c1 - c2) * r2 + 0.2e1 * omega3 * (rd1 - r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) * rd2));
         ];

r =    [(l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) * cos(theta) - sin(theta) * (c1 - c2) * r2;
        (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) * sin(theta) + cos(theta) * (c1 - c2) * r2;
         0];

ExA =  [0
       -0.1200000000e1 * cos(theta) * rd2 ^ 2/l
       -0.1200000000e1 * sin(theta) * rd2 ^ 2/l;];

binv.A2hat1 = A2hat;


Sk1 = [eye(3), cross(r); zeros(3,3),eye(3)];
Sk1k2 = Sk1(3:5,3:5);

phi2k = [0,   cos(theta), sin(theta);
         0,  -0.1200000000e1 * r2 / l * cos(theta) - 0.10e1 * sin(theta),  -0.1200000000e1 * r2 / l * sin(theta) + 0.10e1 * cos(theta)];

binv.ExA1   = ExA;
binv.phi2k1 = phi2k;
binv.S1     = Sk1k2; 

Xrot        = Xrotz(-theta);
Xj0         = Xrot(3:5,3:5);
binv.Xj01   = Xj0;
%% Body 2
elseif ( i == 2)

r1      = z1;
r2      = z2;
rd1     = zd1;
rd2     = zd2;
theta   = q1+q2;
omega3  = qd1+qd2;
Lflex  =  LFlex(c1,c2,l,sigma,U,r1,r2,rd1,rd2,omega3,theta);
Lpoint =  LPoint(c1,c2,l,m,r1,r2,rd1,rd2,omega3,theta);
gamma  =  Gamma(c1,c2,l,r1,r2,theta);
beta   =  Beta(c1,c2,l,Y,II,sigma,U,m,r1,r2,theta);

A2hat = [  0;
           cos(theta) * (-omega3 ^ 2 * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) - 0.2e1 * omega3 * (c1 - c2) * rd2) - sin(theta) * (-omega3 ^ 2 * (c1 - c2) * r2 + 0.2e1 * omega3 * (rd1 - r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) * rd2));
           sin(theta) * (-omega3 ^ 2 * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) - 0.2e1 * omega3 * (c1 - c2) * rd2) + cos(theta) * (-omega3 ^ 2 * (c1 - c2) * r2 + 0.2e1 * omega3 * (rd1 - r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) * rd2));
         ];

r =    [(l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) * cos(theta) - sin(theta) * (c1 - c2) * r2;
        (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1) * sin(theta) + cos(theta) * (c1 - c2) * r2;
         0];

ExA =  [0
       -0.1200000000e1 * cos(theta) * rd2 ^ 2/l
       -0.1200000000e1 * sin(theta) * rd2 ^ 2/l;];

binv.A2hat2 = A2hat;


Sk1 = [eye(3), cross(r); zeros(3,3),eye(3)];
Sk1k2 = Sk1(3:5,3:5);

phi2k = [0,   cos(theta), sin(theta);
         0,  -0.1200000000e1 * r2 / l * cos(theta) - 0.10e1 * sin(theta),  -0.1200000000e1 * r2 / l * sin(theta) + 0.10e1 * cos(theta)];

binv.ExA2 = ExA;
binv.phi2k2 = phi2k;
binv.S2     = Sk1k2;

Xrot        = Xrotz(-theta);
Xj0         = Xrot(3:5,3:5);
binv.Xj02   = Xj0;
end
%% Inverse Inertias

LambdaRR = Lflex.LRR+Lpoint.lRR;
LambdaRF = Lflex.LRF+Lpoint.lRF;
LambdaFR = Lflex.LFR+Lpoint.lFR;
LambdaFF = Lflex.LFF+Lpoint.lFF;

betaR  =  Lflex.BR  - beta.BR + Lpoint.BR; 
betaF  =  Lflex.BF  - beta.BF + Lpoint.BF;

GammaR1 = gamma.k1gr;
GammaF1 = gamma.k1gf;

GammaR2 = gamma.k2gr;
GammaF2 = gamma.k2gf;

%% To use in the ydd and zdd
binv.LambdaRR = LambdaRR;
binv.LambdaRF = LambdaRF;
binv.LambdaFR = LambdaFR;
binv.LambdaFF = LambdaFF;

binv.betaR  =  betaR; 
binv.betaF  =  betaF;

binv.GammaR1 = GammaR1;
binv.GammaF1 = GammaF1;

binv.GammaR2 = GammaR2;
binv.GammaF2 = GammaF2;

%% Inverse Inertia handle-1
invLambdaFF = inv(LambdaFF);

invmat1 = inv(LambdaRR - LambdaRF*inv(LambdaFF)*LambdaFR);
binv.body_zeta11 = (invmat1*(GammaR1 - LambdaRF*inv(LambdaFF)*GammaF1));
binv.body_zeta12 = (invmat1*(GammaR2 - LambdaRF*inv(LambdaFF)*GammaF2));
binv.body_zeta13 = (invmat1*(-betaR + LambdaRF*inv(LambdaFF)*betaF));
%% eta matrices from Eq 5.49 Rudra
eta1 = (Sk1k2' - phi2k'*invLambdaFF*LambdaFR);   
eta2 = (phi2k'*invLambdaFF*GammaF1);
eta3 = (phi2k'*invLambdaFF*GammaF2); 
eta4 = (A2hat-phi2k'*invLambdaFF*betaF+ExA);

%% Inverse inertia handle-2
binv.body_zeta21 = (eta2+eta1*binv.body_zeta11); %Eq 5.51 Rudra
binv.body_zeta22 = (eta3+eta1*binv.body_zeta12); %Eq 5.51 Rudra
binv.body_zeta23 = (eta4+eta1*binv.body_zeta13); %Eq 5.51 Rudra

