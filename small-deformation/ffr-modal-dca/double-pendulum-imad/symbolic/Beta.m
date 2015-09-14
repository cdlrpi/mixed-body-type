function beta = Beta(c1,c2,l,Y,II,sigma,U,m,r1,r2,theta)
                    
g = 9.81;

beta213 = 0.3e1 / 0.20e2 * cos(theta) * r2 ^ 2 * c2 ^ 2 * sigma * U * g - 0.3e1 / 0.10e2 * cos(theta) * r2 ^ 2 * c1 * c2 * sigma * U * g - (sin(theta) * c2 / l ^ 3 * r2 - 0.2e1 / 0.3e1 * cos(theta) * r2 ^ 2 * c1 ^ 2 / l ^ 4) * sigma * U * g * l ^ 4 / 0.4e1 - (-sin(theta) * c1 / l ^ 2 * r2 + cos(theta) / l ^ 2 * r1) * sigma * U * g * l ^ 3 / 0.3e1 - cos(theta) * sigma * U * g * l ^ 2 / 0.2e1 - (-sin(theta) * (c1 - c2) * r2 + cos(theta) * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1)) * m * g;
beta221 = 0;
beta222 = -sigma * U * g * l - m * g;
beta231 = -sin(theta) * l * sigma * U * g / 0.3e1 - 0.4e1 / 0.3e1 * cos(theta) ^ 2 * Y * U / l * r1 + sin(theta) * (-0.4e1 / 0.3e1 * sin(theta) * Y * U / l * r1 - m * g);
beta232 = 0.3e1 / 0.10e2 * sin(theta) * r2 * c2 ^ 2 * sigma * U * g - 0.3e1 / 0.5e1 * sin(theta) * r2 * c1 * c2 * sigma * U * g - (-0.4e1 / 0.3e1 * sin(theta) * r2 * c1 ^ 2 / l ^ 4 - cos(theta) * c2 / l ^ 3) * sigma * U * g * l ^ 4 / 0.4e1 - cos(theta) * c1 * l * sigma * U * g / 0.3e1 - (-cos(theta) * r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) - sin(theta) * (c1 - c2)) * sin(theta) * (-0.12e2 * Y * II * c2 ^ 2 / l ^ 3 * r2 + 0.12e2 * Y * II * c1 / l ^ 3 * c2 * r2 - 0.4e1 * Y * II * c1 ^ 2 / l ^ 3 * r2) + (-sin(theta) * r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) + cos(theta) * (c1 - c2)) * (cos(theta) * (-0.12e2 * Y * II * c2 ^ 2 / l ^ 3 * r2 + 0.12e2 * Y * II * c1 / l ^ 3 * c2 * r2 - 0.4e1 * Y * II * c1 ^ 2 / l ^ 3 * r2) - m * g);


B21 =  beta213;
B22 = [beta221;beta222];
B23 = [beta231;beta232];

beta.BR = [B21;B22];
beta.BF =  B23;


