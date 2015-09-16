function  [outacc] = mainfunc(t,dt,p,v)
persistent fdd gdd ydd zdd
c1 = 1.5;
c2 = 0.5;
q1 = p(1);
y1 = p(2);
y2 = p(3);
q2 = p(4);
z1 = p(5);
z2 = p(6);

qd1 = v(1);
yd1 = v(2);
yd2 = v(3);
qd2 = v(4);
zd1 = v(5);
zd2 = v(6);


S  = [1;0;0];
D2 = eye(3);
phid1   = 2 * c1 / 0.545 - 3 * c2 / 0.545;
phid2   = 2 * c1 / 0.675 - 3 * c2 / 0.675;

if (t == 0+dt)

qdd(1) =  0;
qdd(2) =  0;
ydd    =  [0;0];
zdd    =  [0;0];
 
fdd = 0;
gdd = 0;
elseif (t>0 && t<1/6)

a(1)   =  pi * (-1 + 72 * t ^ 3) / 0.4e1;
ad(1)  =  0.54e2 * pi * t ^ 2;
add(1) =  0.108e3 * pi * t;
a(2)   =  -a(1);
ad(2)  =  -ad(1);
add(2) =  -add(1);

matrices = mat(q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
invphi   = mainpass(matrices,add,fdd,gdd,D2);    %Main pass
back     = disassembly(invphi,add,fdd,gdd); %Back Substitution pass

ydd    = inv(-matrices.Lambda1FF)*(matrices.Lambda1FR*back.A1{1,1}-matrices.Gamma1F1*back.f1{1,1}-matrices.Gamma1F2*back.f2{1,1}+matrices.beta1F);
zdd    = inv(-matrices.Lambda2FF)*(matrices.Lambda2FR*back.A1{1,2}-matrices.Gamma2F1*back.f1{1,2}-matrices.Gamma2F2*back.f2{1,2}+matrices.beta2F);

fdd    = phid1*ydd(2);
gdd    = phid2*zdd(2);

qdd(1)     = S'*(back.A1{1,1});
qdd(2)     = S'*(back.A1{1,2}-back.A2{1,1})+fdd;

elseif (t>=1/6 && t<1/3)

a(1)   =  pi * (-18 * t + 108 * t ^ 2 - 144 * t ^ 3) / 0.4e1;
ad(1)  =  pi * (-18 + 216 * t - 432 * t ^ 2) / 0.4e1;
add(1) =  pi * (216 - 864 * t) / 0.4e1;
a(2)   =  -a(1);
ad(2)  =  -ad(1);
add(2) =  -add(1);

matrices = mat(q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
invphi   = mainpass(matrices,add,fdd,gdd,D2);    %Main pass
back     = disassembly(invphi,add,fdd,gdd); %Back Substitution pass

ydd    = inv(-matrices.Lambda1FF)*(matrices.Lambda1FR*back.A1{1,1}-matrices.Gamma1F1*back.f1{1,1}-matrices.Gamma1F2*back.f2{1,1}+matrices.beta1F);
zdd    = inv(-matrices.Lambda2FF)*(matrices.Lambda2FR*back.A1{1,2}-matrices.Gamma2F1*back.f1{1,2}-matrices.Gamma2F2*back.f2{1,2}+matrices.beta2F);

fdd    = phid1*ydd(2);
gdd    = phid2*zdd(2);

qdd(1)     = S'*(back.A1{1,1});
qdd(2)     = S'*(back.A1{1,2}-back.A2{1,1})+fdd;

elseif (t>=1/3 && t<1/2)

a(1)   =  pi * (-8 + 54 * t - 108 * t ^ 2 + 72 * t ^ 3) / 0.4e1;
ad(1)  =  pi * (54 - 216 * t + 216 * t ^ 2) / 0.4e1;
add(1) =  pi * (-216 + 432 * t) / 0.4e1;
a(2)   =  -a(1);
ad(2)  =  -ad(1);
add(2) =  -add(1);

matrices = mat(q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
invphi   = mainpass(matrices,add,fdd,gdd,D2);    %Main pass
back     = disassembly(invphi,add,fdd,gdd); %Back Substitution pass

ydd    = inv(-matrices.Lambda1FF)*(matrices.Lambda1FR*back.A1{1,1}-matrices.Gamma1F1*back.f1{1,1}-matrices.Gamma1F2*back.f2{1,1}+matrices.beta1F);
zdd    = inv(-matrices.Lambda2FF)*(matrices.Lambda2FR*back.A1{1,2}-matrices.Gamma2F1*back.f1{1,2}-matrices.Gamma2F2*back.f2{1,2}+matrices.beta2F);

fdd    = phid1*ydd(2);
gdd    = phid2*zdd(2);

qdd(1)     = S'*(back.A1{1,1});
qdd(2)     = S'*(back.A1{1,2}-back.A2{1,1})+fdd;

elseif (t >= 1/2)

a(1)   =  pi/4;
ad(1)  =  0;
add(1) =  0;
a(2)   =  -a(1);
ad(2)  =  -ad(1);
add(2) =  -add(1);

matrices = mat(q1, qd1, y1, yd1, y2, yd2, q2, qd2, z1, zd1, z2, zd2, t);
invphi   = mainpass(matrices,add,fdd,gdd,D2);    %Main pass
back     = disassembly(invphi,add,fdd,gdd);          %Back Substitution pass

ydd    = inv(-matrices.Lambda1FF)*(matrices.Lambda1FR*back.A1{1,1}-matrices.Gamma1F1*back.f1{1,1}-matrices.Gamma1F2*back.f2{1,1}+matrices.beta1F);
zdd    = inv(-matrices.Lambda2FF)*(matrices.Lambda2FR*back.A1{1,2}-matrices.Gamma2F1*back.f1{1,2}-matrices.Gamma2F2*back.f2{1,2}+matrices.beta2F);

fdd    = phid1*ydd(2);
gdd    = phid2*zdd(2);

qdd(1)     = S'*(back.A1{1,1});
qdd(2)     = S'*(back.A1{1,2}-back.A2{1,1})+fdd;

end

outacc = [qdd(1) ydd' qdd(2) zdd'];