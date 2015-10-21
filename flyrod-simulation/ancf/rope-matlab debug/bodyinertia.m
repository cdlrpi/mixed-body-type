function binv = bodyinertia(i, ne, e, iM, l, rho, A, Y, I)

e1 = e(1);
e2 = e(2);
e3 = e(3);
e4 = e(4);
e5 = e(5);
e6 = e(6);
e7 = e(7);
e8 = e(8);



% gamma  =  Gamma(e3,e4,e7,e8);
beta   =  Beta(e1,e2,e3,e4,e5,e6,e7,e8,l,rho,A,Y,I);

% A1 = iM*[gamma.L11;zeros(4,4)];
% A2 = iM*[zeros(4,4);gamma.L22];
A3 = iM*[beta.L13;beta.L23];

A1 = iM*[eye(4);zeros(4,4)];
A2 = iM*[zeros(4,4);eye(4)];

binv.body_zeta11 = A1(1:4,:);
binv.body_zeta12 = A2(1:4,:);
binv.body_zeta13 = A3(1:4,:);

binv.body_zeta21 = A1(5:8,:);
binv.body_zeta22 = A2(5:8,:);
binv.body_zeta23 = A3(5:8,:);
keyboard
