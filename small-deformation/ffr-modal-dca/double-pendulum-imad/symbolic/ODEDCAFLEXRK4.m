function  x = ODEDCAFLEXRK4(x,dt,Time,IC)

q1  = IC(1);
qd1 = IC(2);
y1  = IC(3);
yd1 = IC(4);
y2  = IC(5);
yd2 = IC(6);
q2  = IC(7);
qd2 = IC(8);
z1  = IC(9);
zd1 = IC(10);
z2  = IC(11);
zd2 = IC(12);


x(1,:)  = IC;

p = [];
v = [];
i = 1;
for j = 1:size(x,2)/2
    p = [p,x(1,i)];
    v = [v,x(1,i+1)];
    i = i+2;
end

for (i = 2:length(Time))
t   = Time(i)    
    
k1  = mainfunc(t, dt, p, v);
kp1 = v;

k2  = mainfunc(t, dt, p + dt/2 * (kp1), v + k1*dt/2);
kp2 = v + k1*dt/2;

k3  = mainfunc(t, dt, p + dt/2 * (kp2), v + k2*dt/2);
kp3 = v + k2*dt/2;

k4  = mainfunc(t, dt, p + dt   * (kp3), v + dt*k3);
kp4 = v + k3*dt;

v = v + dt/6*(k1+2*k2+2*k3+k4);
p = p + dt/6*(kp1+2*kp2+2*kp3+kp4);

k = 1;
for j    = 1:size(p,2)
x(i,k)   = p(j);
x(i,k+1) = v(j);
k = k+2;
end
end

