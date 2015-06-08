
clear
clc
ne  =  12;
L   =  1.2;
Y   =  0.7e6;
I   =  1.215e-8;
A   =  0.0018;
rho =  5540;
l   =  L/ne;
g   =  9.81;
II  =  I;
ipe = zeros(8,ne);
ive = zeros(8,ne);
theta1_i =  0;
for i  = 1:ne
    e1 = (i*l-l)*cos(theta1_i);
    e2 = (i*l-l)*sin(theta1_i);
    e3 = cos(theta1_i);
    e4 = sin(theta1_i);
    
    e5 = i*l*cos(theta1_i);
    e6 = i*l*sin(theta1_i);
    e7 = cos(theta1_i);
    e8 = sin(theta1_i);
    
    ipe(:,i) = [e1 e2 e3 e4 e5 e6 e7 e8]';
end
size(ipe)


[M11 M12 M21 M22 eta1 iM22]  = MassMatrix(l,rho,A);
M    =  [M11, M12; M21, M22];
iM   =  inv(M);
dt   =  0.0001;
endT =  1;
Time =  0:dt:endT;
tic
[pe ve]  =  ODEDCAANCFRK4(dt,Time,ipe,ive,ne,iM,l,rho,A,Y,I);
toc

visualize(l,Time,ne,pe,ve,rho,A,g,Y,II,M,L)

