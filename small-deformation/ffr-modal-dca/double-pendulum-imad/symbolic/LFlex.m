function  Lflex = LFlex(c1,c2,l,sigma,U,r1,r2,rd1,rd2,omega3,theta)
                       
                       
L33 = -0.2250000000e0 * sigma * U * c2 ^ 2 * r2 ^ 2 * r1 + 0.1666666667e0 * sigma * U * l * c2 * c1 * r2 ^ 2 - 0.2222222222e0 * sigma * U * r2 ^ 2 * c1 ^ 2 * r1 - 0.2700000000e0 * sigma * U / l * r2 ^ 4 * c1 * c2 ^ 3 + 0.3833333333e0 * sigma * U / l * r2 ^ 4 * c2 ^ 2 * c1 ^ 2 - 0.2500000000e0 * sigma * U / l * c2 * r2 ^ 4 * c1 ^ 3 + 0.7363636364e-1 * sigma * U / l * r2 ^ 4 * c2 ^ 4 + 0.6349206352e-1 * sigma * U / l * r2 ^ 4 * c1 ^ 4 - 0.1142857143e0 * sigma * U * l * r2 ^ 2 * c2 ^ 2 - 0.6666666660e-1 * sigma * U * l * r2 ^ 2 * c1 ^ 2 + 0.3333333333e0 * sigma * U * l ^ 3 + 0.5000000000e0 * sigma * U * l ^ 2 * r1 + 0.2000000000e0 * sigma * U * l * r1 ^ 2 - 0.6666666668e-10 * sigma * U * sin(theta) * r2 ^ 3 * c1 ^ 3 * cos(theta) + 0.4285714288e0 * sigma * U * c1 * c2 * r1 * r2 ^ 2;
L34 = 0.1500000000e0 * sigma * U * sin(theta) * r2 ^ 2 * c2 ^ 2 + 0.1666666667e0 * sigma * U * r2 ^ 2 * sin(theta) * c1 ^ 2 + 0.2500000000e0 * sigma * U * l * r2 * cos(theta) * c2 - 0.3333333333e0 * sigma * U * l * cos(theta) * c1 * r2 - 0.3333333333e0 * sigma * U * l * sin(theta) * r1 - 0.5000000000e0 * sigma * U * l ^ 2 * sin(theta) - 0.3000000000e0 * sigma * U * sin(theta) * r2 ^ 2 * c1 * c2;
L35 = -0.1500000000e0 * sigma * U * cos(theta) * r2 ^ 2 * c2 ^ 2 - 0.1666666667e0 * sigma * U * r2 ^ 2 * cos(theta) * c1 ^ 2 + 0.2500000000e0 * sigma * U * l * r2 * sin(theta) * c2 - 0.3333333333e0 * sigma * U * l * sin(theta) * c1 * r2 + 0.3333333333e0 * sigma * U * l * cos(theta) * r1 + 0.5000000000e0 * sigma * U * l ^ 2 * cos(theta) + 0.3000000000e0 * sigma * U * cos(theta) * r2 ^ 2 * c1 * c2;
L43 = 0.1500000000e0 * sigma * U * sin(theta) * r2 ^ 2 * c2 ^ 2 + 0.1666666667e0 * sigma * U * r2 ^ 2 * sin(theta) * c1 ^ 2 + 0.2500000000e0 * sigma * U * l * r2 * cos(theta) * c2 - 0.3333333333e0 * sigma * U * l * cos(theta) * c1 * r2 - 0.3333333333e0 * sigma * U * l * sin(theta) * r1 - 0.5000000000e0 * sigma * U * l ^ 2 * sin(theta) - 0.3000000000e0 * sigma * U * sin(theta) * r2 ^ 2 * c1 * c2;
L44 = 0.1000000000e1 * sigma * U * l;
L45 = 0;
L53 = -0.1500000000e0 * sigma * U * cos(theta) * r2 ^ 2 * c2 ^ 2 - 0.1666666667e0 * sigma * U * r2 ^ 2 * cos(theta) * c1 ^ 2 + 0.2500000000e0 * sigma * U * l * r2 * sin(theta) * c2 - 0.3333333333e0 * sigma * U * l * sin(theta) * c1 * r2 + 0.3333333333e0 * sigma * U * l * cos(theta) * r1 + 0.5000000000e0 * sigma * U * l ^ 2 * cos(theta) + 0.3000000000e0 * sigma * U * cos(theta) * r2 ^ 2 * c1 * c2;
L54 = 0;
L55 = 0.1000000000e1 * sigma * U * l;
L37 = 0.1666666667e0 * sigma * U * l * r2 * c2 - 0.2000000000e0 * sigma * U * l * c1 * r2;
L38 = 0.3000000000e0 * sigma * U * c2 ^ 2 * r2 ^ 2 * c1 - 0.3095238095e0 * sigma * U * c1 ^ 2 * c2 * r2 ^ 2 - 0.1000000000e0 * sigma * U * r2 ^ 2 * c2 ^ 3 + 0.1111111111e0 * sigma * U * r2 ^ 2 * c1 ^ 3 - 0.1666666667e0 * sigma * U * l * r1 * c2 + 0.2000000000e0 * sigma * U * l * r1 * c1 + 0.2500000000e0 * sigma * U * l ^ 2 * c1 - 0.2000000000e0 * sigma * U * l ^ 2 * c2;
L47 = 0.3333333333e0 * sigma * U * l * cos(theta);
L48 = 0.6000000000e0 * sigma * U * c1 * c2 * cos(theta) * r2 + 0.2500000000e0 * sigma * U * sin(theta) * c2 * l - 0.3333333332e0 * sigma * U * cos(theta) * r2 * c1 ^ 2 - 0.3000000001e0 * sigma * U * c2 ^ 2 * cos(theta) * r2 - 0.3333333333e0 * sigma * U * l * sin(theta) * c1;
L57 = 0.3333333333e0 * sigma * U * l * sin(theta);
L58 = -0.2500000000e0 * sigma * U * cos(theta) * c2 * l - 0.3333333332e0 * sigma * U * sin(theta) * r2 * c1 ^ 2 + 0.6000000000e0 * sigma * U * c1 * c2 * sin(theta) * r2 - 0.3000000001e0 * sigma * U * c2 ^ 2 * sin(theta) * r2 + 0.3333333333e0 * sigma * U * l * cos(theta) * c1;
L73 = 0.1666666667e0 * sigma * U * l * r2 * c2 - 0.2000000000e0 * sigma * U * l * c1 * r2;
L74 = 0.3333333333e0 * sigma * U * l * cos(theta);
L75 = 0.3333333333e0 * sigma * U * l * sin(theta);
L83 = 0.1111111111e0 * sigma * U * r2 ^ 2 * c1 ^ 3 - 0.1666666667e0 * sigma * U * l * r1 * c2 - 0.9999999999e-1 * sigma * U * r2 ^ 2 * c2 ^ 3 + 0.2000000000e0 * sigma * U * l * r1 * c1 - 0.2000000000e0 * sigma * U * l ^ 2 * c2 + 0.2500000000e0 * sigma * U * l ^ 2 * c1 + 0.3000000000e0 * sigma * U * c2 ^ 2 * r2 ^ 2 * c1 - 0.3095238095e0 * sigma * U * c1 ^ 2 * c2 * r2 ^ 2;
L84 = 0.6000000000e0 * sigma * U * c1 * c2 * cos(theta) * r2 + 0.2500000000e0 * sigma * U * sin(theta) * c2 * l - 0.3333333332e0 * sigma * U * cos(theta) * r2 * c1 ^ 2 - 0.3000000001e0 * sigma * U * c2 ^ 2 * cos(theta) * r2 - 0.3333333333e0 * sigma * U * l * sin(theta) * c1;
L85 = -0.2500000000e0 * sigma * U * cos(theta) * c2 * l - 0.3333333332e0 * sigma * U * sin(theta) * r2 * c1 ^ 2 + 0.6000000000e0 * sigma * U * c1 * c2 * sin(theta) * r2 - 0.3000000001e0 * sigma * U * c2 ^ 2 * sin(theta) * r2 + 0.3333333333e0 * sigma * U * l * cos(theta) * c1;
L77 = 0.2000000000e0 * sigma * U * l;
L78 = 0.4285714287e0 * sigma * U * r2 * c1 * c2 - 0.2222222222e0 * sigma * U * r2 * c1 ^ 2 - 0.2250000000e0 * sigma * U * r2 * c2 ^ 2;
L87 = 0.4285714287e0 * sigma * U * r2 * c1 * c2 - 0.2222222222e0 * sigma * U * r2 * c1 ^ 2 - 0.2250000000e0 * sigma * U * r2 * c2 ^ 2;
L88 = 0.2945454545e0 * sigma * U / l * r2 ^ 2 * c2 ^ 4 + 0.2539682539e0 * sigma * U / l * r2 ^ 2 * c1 ^ 4 - 0.3333333334e0 * sigma * U * l * c1 * c2 + 0.2000000000e0 * sigma * U * l * c1 ^ 2 + 0.1428571429e0 * sigma * U * l * c2 ^ 2 - 0.9999999998e0 * sigma * U / l * c2 * r2 ^ 2 * c1 ^ 3 - 0.1080000000e1 * sigma * U / l * r2 ^ 2 * c1 * c2 ^ 3 + 0.1533333333e1 * sigma * U / l * r2 ^ 2 * c2 ^ 2 * c1 ^ 2;

B3 = 0.1000000000e-9 * sigma * U * (-0.1080000000e11 * r2 ^ 3 * c1 * c2 ^ 3 * omega3 * rd2 + 0.1533333333e11 * r2 ^ 3 * c2 ^ 2 * omega3 * c1 ^ 2 * rd2 - 0.1000000000e11 * c2 * r2 ^ 3 * c1 ^ 3 * omega3 * rd2 + 0.6000000000e10 * c2 ^ 2 * r2 * c1 * rd2 ^ 2 * l - 0.2250000000e10 * c2 ^ 2 * r2 ^ 2 * omega3 * rd1 * l - 0.2222222222e10 * r2 ^ 2 * c1 ^ 2 * omega3 * rd1 * l - 0.6190476192e10 * c1 ^ 2 * c2 * r2 * rd2 ^ 2 * l - 0.1333333334e10 * l ^ 2 * c1 ^ 2 * r2 * omega3 * rd2 - 0.2285714286e10 * l ^ 2 * c2 ^ 2 * r2 * omega3 * rd2 + 0.3333333334e10 * l ^ 2 * r2 * c2 * omega3 * c1 * rd2 + 0.4285714288e10 * r2 ^ 2 * c1 * c2 * omega3 * rd1 * l - 0.4500000000e10 * c2 ^ 2 * r2 * r1 * omega3 * rd2 * l - 0.4444444448e10 * r1 * omega3 * r2 * c1 ^ 2 * rd2 * l + 0.8571428576e10 * c1 * c2 * r1 * r2 * omega3 * rd2 * l + 0.5000000000e10 * l ^ 3 * omega3 * rd1 + 0.2945454545e10 * r2 ^ 3 * c2 ^ 4 * omega3 * rd2 + 0.2539682541e10 * r2 ^ 3 * c1 ^ 4 * omega3 * rd2 + 0.4000000000e10 * l ^ 2 * r1 * omega3 * rd1 - 0.2000000000e10 * r2 * c2 ^ 3 * rd2 ^ 2 * l + 0.2222222222e10 * c1 ^ 3 * r2 * rd2 ^ 2 * l) / l;
B4 = 0.1000000000e-9 * sigma * U * (-0.3000000001e10 * c2 ^ 2 * cos(theta) * rd2 ^ 2 - 0.3333333332e10 * cos(theta) * rd2 ^ 2 * c1 ^ 2 + 0.1666666667e10 * cos(theta) * omega3 ^ 2 * r2 ^ 2 * c1 ^ 2 + 0.1500000000e10 * c2 ^ 2 * cos(theta) * omega3 ^ 2 * r2 ^ 2 + 0.6000000000e10 * c1 * c2 * cos(theta) * rd2 ^ 2 - 0.6666666666e10 * l * sin(theta) * omega3 * rd1 - 0.3333333333e10 * l * cos(theta) * omega3 ^ 2 * r1 + 0.6666666668e10 * sin(theta) * omega3 * r2 * c1 ^ 2 * rd2 - 0.2500000000e10 * sin(theta) * omega3 ^ 2 * c2 * r2 * l + 0.5000000000e10 * cos(theta) * omega3 * c2 * rd2 * l + 0.6000000001e10 * c2 ^ 2 * sin(theta) * r2 * omega3 * rd2 - 0.3000000000e10 * c1 * c2 * cos(theta) * omega3 ^ 2 * r2 ^ 2 + 0.3333333333e10 * l * sin(theta) * omega3 ^ 2 * c1 * r2 - 0.6666666666e10 * l * cos(theta) * omega3 * c1 * rd2 - 0.1200000000e11 * c1 * c2 * sin(theta) * r2 * omega3 * rd2 - 0.5000000000e10 * l ^ 2 * cos(theta) * omega3 ^ 2);
B5 = 0.1000000000e-9 * sigma * U * (-0.6666666668e10 * cos(theta) * omega3 * r2 * c1 ^ 2 * rd2 - 0.6666666666e10 * l * sin(theta) * omega3 * c1 * rd2 - 0.3000000000e10 * c1 * c2 * sin(theta) * omega3 ^ 2 * r2 ^ 2 - 0.6000000001e10 * c2 ^ 2 * cos(theta) * r2 * omega3 * rd2 + 0.5000000000e10 * sin(theta) * omega3 * c2 * rd2 * l + 0.2500000000e10 * cos(theta) * omega3 ^ 2 * c2 * r2 * l - 0.3333333333e10 * l * cos(theta) * omega3 ^ 2 * c1 * r2 + 0.1200000000e11 * c1 * c2 * cos(theta) * r2 * omega3 * rd2 - 0.3000000001e10 * c2 ^ 2 * sin(theta) * rd2 ^ 2 - 0.5000000000e10 * l ^ 2 * sin(theta) * omega3 ^ 2 - 0.3333333332e10 * sin(theta) * rd2 ^ 2 * c1 ^ 2 + 0.1500000000e10 * c2 ^ 2 * sin(theta) * omega3 ^ 2 * r2 ^ 2 + 0.6000000000e10 * c1 * c2 * sin(theta) * rd2 ^ 2 + 0.6666666666e10 * l * cos(theta) * omega3 * rd1 + 0.1666666667e10 * sin(theta) * omega3 ^ 2 * r2 ^ 2 * c1 ^ 2 - 0.3333333333e10 * l * sin(theta) * omega3 ^ 2 * r1);
B7 = 0.3333333334e0 * sigma * U * omega3 * c2 * rd2 * l - 0.2142857144e0 * sigma * U * omega3 ^ 2 * r2 ^ 2 * c1 * c2 + 0.4285714287e0 * sigma * U * rd2 ^ 2 * c1 * c2 - 0.2222222222e0 * sigma * U * rd2 ^ 2 * c1 ^ 2 + 0.1111111111e0 * sigma * U * omega3 ^ 2 * r2 ^ 2 * c1 ^ 2 - 0.2250000000e0 * sigma * U * rd2 ^ 2 * c2 ^ 2 + 0.1125000000e0 * sigma * U * omega3 ^ 2 * r2 ^ 2 * c2 ^ 2 - 0.4000000000e0 * sigma * U * l * omega3 * c1 * rd2 - 0.2000000000e0 * sigma * U * l * omega3 ^ 2 * r1 - 0.2500000000e0 * sigma * U * l ^ 2 * omega3 ^ 2;
B8 = -0.1000000000e-19 * sigma * U * (-0.5000000000e20 * c2 * r2 ^ 3 * omega3 ^ 2 * c1 ^ 3 + 0.9999999998e20 * c2 * r2 * rd2 ^ 2 * c1 ^ 3 + 0.1666666667e11 * l * c1 ^ 3 * omega3 * r2 * rd2 - 0.2222222222e20 * l * r2 * c1 ^ 2 * omega3 ^ 2 * r1 + 0.1666666667e20 * c1 * omega3 ^ 2 * c2 * r2 * l ^ 2 + 0.4285714287e20 * r2 * c1 * c2 * l * omega3 ^ 2 * r1 - 0.1428571429e11 * c1 ^ 2 * c2 * l * r2 * omega3 * rd2 - 0.2250000000e20 * c2 ^ 2 * r2 * l * omega3 ^ 2 * r1 + 0.1472727273e20 * r2 ^ 3 * c2 ^ 4 * omega3 ^ 2 - 0.2945454545e20 * r2 * c2 ^ 4 * rd2 ^ 2 - 0.5400000000e20 * r2 ^ 3 * c1 * c2 ^ 3 * omega3 ^ 2 + 0.1080000000e21 * r2 * c1 * c2 ^ 3 * rd2 ^ 2 + 0.1269841270e20 * r2 ^ 3 * c1 ^ 4 * omega3 ^ 2 - 0.2539682539e20 * r2 * c1 ^ 4 * rd2 ^ 2 - 0.1142857143e20 * c2 ^ 2 * l ^ 2 * omega3 ^ 2 * r2 + 0.5000000001e10 * l * cos(theta) * c1 ^ 3 * sin(theta) * omega3 ^ 2 * r2 ^ 2 + 0.3333333334e20 * c2 * l ^ 2 * omega3 * rd1 + 0.7666666666e20 * r2 ^ 3 * c2 ^ 2 * omega3 ^ 2 * c1 ^ 2 - 0.1533333333e21 * r2 * c2 ^ 2 * rd2 ^ 2 * c1 ^ 2 - 0.6666666660e19 * l ^ 2 * r2 * c1 ^ 2 * omega3 ^ 2 - 0.4000000000e20 * c1 * l ^ 2 * omega3 * rd1) / l;



%% Construct Matrices
Lflex.LRR = [L33 L34 L35;
             L43 L44 L45;
             L53 L54 L55];
             
Lflex.LRF = [L37, L38;
             L47, L48;
             L57, L58];

Lflex.LFR = [L73, L74, L75;
             L83, L84, L85];

Lflex.LFF = [L77, L78; L87, L88];   
 
B11 =  B3;
B12 = [B4;B5];
B13 = [B7;B8];

Lflex.BR = [B11;B12];
Lflex.BF =  B13;


