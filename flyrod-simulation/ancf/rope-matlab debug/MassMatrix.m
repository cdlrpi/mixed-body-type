function  [M11 M12 M21 M22 eta1 iM22] = MassMatrix(l,rho,A)

M1 = [0.13e2 / 0.35e2 * l * rho * A 0 0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 0.9e1 / 0.70e2 * l * rho * A 0 -0.13e2 / 0.420e3 * l ^ 2 * rho * A 0;];
M2 = [0 0.13e2 / 0.35e2 * l * rho * A 0 0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 0.9e1 / 0.70e2 * l * rho * A 0 -0.13e2 / 0.420e3 * l ^ 2 * rho * A;];
M3 = [0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 l ^ 3 * rho * A / 0.105e3 0 0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 -l ^ 3 * rho * A / 0.140e3 0;];
M4 = [0 0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 l ^ 3 * rho * A / 0.105e3 0 0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 -l ^ 3 * rho * A / 0.140e3;];
M5 = [0.9e1 / 0.70e2 * l * rho * A 0 0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 0.13e2 / 0.35e2 * l * rho * A 0 -0.11e2 / 0.210e3 * l ^ 2 * rho * A 0;];
M6 = [0 0.9e1 / 0.70e2 * l * rho * A 0 0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 0.13e2 / 0.35e2 * l * rho * A 0 -0.11e2 / 0.210e3 * l ^ 2 * rho * A;];
M7 = [-0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 -l ^ 3 * rho * A / 0.140e3 0 -0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 l ^ 3 * rho * A / 0.105e3 0;];
M8 = [0 -0.13e2 / 0.420e3 * l ^ 2 * rho * A 0 -l ^ 3 * rho * A / 0.140e3 0 -0.11e2 / 0.210e3 * l ^ 2 * rho * A 0 l ^ 3 * rho * A / 0.105e3;];

M   = [M1;M2;M3;M4;M5;M6;M7;M8];
M11 = M(1:4,1:4);
M12 = M(1:4,5:8);
M21 = M(5:8,1:4);
M22 = M(5:8,5:8);

eta1   = inv(M11 - M12*inv(M22)*M21);
iM22   = inv(M22);
