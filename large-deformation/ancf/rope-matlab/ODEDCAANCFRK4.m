function  [pe ve] = ODEDCAANCFRK4(dt,Time,ipe,ive,ne,iM,l,rho,A,Y,I)
pe(:,:,1)  = ipe;
ve(:,:,1)  = ive;
p          = ipe;
v          = ive;


for i = 2:length(Time)
t   = Time(i)

    
k1  = MainFunction(p, ne, iM, l, rho, A, Y, I);
kp1 = v;

k2  = MainFunction(p + dt/2 * (kp1), ne, iM, l, rho, A, Y, I);
kp2 = v + k1*dt/2;

k3  = MainFunction(p + dt/2 * (kp2), ne, iM, l, rho, A, Y, I);
kp3 = v + k2*dt/2;

k4  = MainFunction(p + dt   * (kp3), ne, iM, l, rho, A, Y, I);
kp4 = v + k3*dt;


v = v + dt/6*(k1+2*k2+2*k3+k4);
p = p + dt/6*(kp1+2*kp2+2*kp3+kp4);

pe(:,:,i)  = p;
ve(:,:,i)  = v;
end
keyboard
