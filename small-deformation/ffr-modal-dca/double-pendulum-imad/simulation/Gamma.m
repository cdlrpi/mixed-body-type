function gamma = Gamma(c1,c2,l,r1,r2,theta)
                      
f1g33 = 1;
f1g34 = 0;
f1g35 = 0;
f1g43 = 0;
f1g44 = 1;
f1g45 = 0;
f1g53 = 0;
f1g54 = 0;
f1g55 = 1;
f1g73 = 0;
f1g74 = 0;
f1g75 = 0;
f1g83 = 0;
f1g84 = 0;
f1g85 = 0;

f2g33 = 1;
f2g34 = -cos(theta) * (c1 - c2) * r2 - sin(theta) * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1);
f2g35 = -sin(theta) * (c1 - c2) * r2 + cos(theta) * (l + r1 - r2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) / 0.2e1);
f2g43 = 0;
f2g44 = 1;
f2g45 = 0;
f2g53 = 0;
f2g54 = 0;
f2g55 = 1;
f2g73 = 0;
f2g74 = cos(theta);
f2g75 = sin(theta);
f2g83 = 2 * c1 / l - 3 * c2 / l;
f2g84 = -cos(theta) * r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) - sin(theta) * (c1 - c2);
f2g85 = -sin(theta) * r2 * (0.9e1 / 0.5e1 * c2 ^ 2 / l - 0.3e1 * c1 / l * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / l) + cos(theta) * (c1 - c2);



f1g1 = [f1g33, f1g34, f1g35];
f1g2 = [f1g43, f1g44, f1g45;
        f1g53, f1g54, f1g55];
f1g3 = [f1g73, f1g74, f1g75;
        f1g83, f1g84, f1g85];
    
f2g1 = [f2g33, f2g34, f2g35];
f2g2 = [f2g43, f2g44, f2g45;
        f2g53, f2g54, f2g55];
f2g3 = [f2g73, f2g74, f2g75;
        f2g83, f2g84, f2g85];


gamma.k1gr = [f1g1;f1g2];
gamma.k1gf = f1g3;

gamma.k2gr = [f2g1;f2g2];
gamma.k2gf = f2g3;


    
