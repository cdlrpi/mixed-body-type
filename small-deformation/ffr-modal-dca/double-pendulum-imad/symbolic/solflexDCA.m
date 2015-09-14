clear 
clc
q1    = -pi/4;
qd1   =  0;
y1    =  2.231386235071618e-007;
yd1   =  0;
y2    = -0.004257625936800;
yd2   =  0;
q2    =  0.773676534118413;
qd2   =  0;
z1    =  6.465113539578467e-009;
zd1   =  0;
z2    = -0.013524057996059;
zd2   =  0;
IC   =  [q1 qd1 y1 yd1 y2 yd2 q2 qd2 z1 zd1 z2 zd2];
dt   =  0.0002;
endT =  2;
Time =  0:dt:endT;
x    =  zeros(length(Time),length(IC));                 
tic
x    =  ODEDCAFLEXRK4(x,dt,Time,IC);
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c1    = 1.5;
c2    = 0.5;
la    = 0.545;
lb    = 0.675;
x1    = 0:0.001:0.545;
x2    = 0:0.001:0.675;
phid1 = 2 * c1 / la - 3 * c2 / la;
T     = 0:dt:endT;
count = 1;
rxx2  = [];
ryy2  = [];
fix   = [];

figure(1)
h = figure(1);
jf     = get(h,'JavaFrame');
set(jf,'Maximized',1);

for i = 1:20:length(T)

p1 = [];
p2 = [];
Tt(count) = T(i);

q1  = x(i,1);
qd1 = x(i,2);
u1  = x(i,3);
ud1 = x(i,4);
u2  = x(i,5);
ud2 = x(i,6);

q2  = x(i,7);
qd2 = x(i,8);
v1  = x(i,9);
vd1 = x(i,10);
v2  = x(i,11);
vd2 = x(i,12);

f1  = u1 - u2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / la - 0.3e1 * c1 / la * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / la) / 0.2e1;
f2  = u2;
g1  = v1 - v2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / lb - 0.3e1 * c1 / lb * c2 + 0.4e1 / 0.3e1 * c1 ^ 2 / lb) / 0.2e1;
g2  = v2;

theta = q1;
the1  = theta;

Rot1 = [cos(theta) , -sin(theta), 0;
        sin(theta) ,  cos(theta), 0;
        0, 0, 1];
r1 =   Rot1*[(la+f1); f2; 0];
rxx1(count)  = r1(1);
ryy1(count)  = r1(2);

theta = q1+q2;
the2  = theta;
Rot2 = [cos(theta) , -sin(theta), 0;
        sin(theta) ,  cos(theta), 0;
        0, 0, 1];
r2 =   r1 + Rot2*[(lb+g1); g2; 0];
rxx2(count)  = r2(1);
ryy2(count)  = r2(2);

for j = 1:length(x1)
dx1(j) = x1(j) ^ 2 / la ^ 2 * u1 - u2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / la ^ 6 * x1(j) ^ 5 - 0.3e1 * c1 / la ^ 5 * c2 * x1(j) ^ 4 + 0.4e1 / 0.3e1 * c1 ^ 2 / la ^ 4 * x1(j) ^ 3) / 0.2e1; (c1 * x1(j) ^ 2 / la ^ 2 - c2 * x1(j) ^ 3 / la ^ 3) * u2;
dy1(j) = (1.5*(x1(j)/la)^2 - 0.5*(x1(j)/la)^3)*u2;
p11    = Rot1*[x1(j)+dx1(j) ; dy1(j); 0];
p1     = [p1,p11];
end

for j = 1:length(x2)
dx2(j) = x2(j) ^ 2 / lb ^ 2 * v1 - v2 ^ 2 * (0.9e1 / 0.5e1 * c2 ^ 2 / lb ^ 6 * x2(j) ^ 5 - 0.3e1 * c1 / lb ^ 5 * c2 * x2(j) ^ 4 + 0.4e1 / 0.3e1 * c1 ^ 2 / lb ^ 4 * x2(j) ^ 3) / 0.2e1; (c1 * x2(j) ^ 2 / lb ^ 2 - c2 * x2(j) ^ 3 / lb ^ 3) * v2;
dy2(j) = (1.5*(x2(j)/lb)^2 - 0.5*(x2(j)/lb)^3)*v2;
p22   =  p1(:,length(x1)) + Rot2*[x2(j)+dx2(j) ; dy2(j); 0];
p2     = [p2,p22];
end

subplot(2,2,1)
plot(p1(1,:),p1(2,:),'black','LineWidth',3);
hold on
plot(0,0,'v','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',10);
plot(p2(1,:),p2(2,:),'black','LineWidth',1.5);
plot(rxx1(count),ryy1(count),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.9 .9 .9],'MarkerSize',10);
plot(rxx2(count),ryy2(count),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[.5 .5 .5],'MarkerSize',13);
hold off
text(0.1,0.5,'\downarrow g','FontSize',16)
xlabel('X coordinate [m]','FontSize',12)
ylabel('Y coordinate [m]','FontSize',12)
text(-0,-.4,{'Imad M Khan','Computational Dynamics Lab', 'Rensselaer Polytechnic Institute'},'FontSize',6,'EdgeColor','k','BackgroundColor',[1 1 1])
grid on
axis square
axis([-0.04 1.25 -0.5 0.7]);

subplot(2,2,2)
if (Tt<=0.57)
plot(rxx2(1:count),ryy2(1:count),'black','LineWidth',1.5);
fix = count;
else
plot(rxx2(1:fix),ryy2(1:fix),'black','LineWidth',1.5);
end
title('End point trajectory, t = 0 .. 0.57 secs','FontSize',12);
xlabel('X coordinate [m]','FontSize',12)
ylabel('Y coordinate [m]','FontSize',12)
grid on
axis square
axis([1 1.25 -0.6 0.7]);


subplot(2,2,3)
plot(Tt(1:count),rxx2(1:count),'black','LineWidth',1.5);
xlabel('Time [sec]','FontSize',12)
ylabel('X coordinate [m]','FontSize',12)
grid on
axis([0 2 1 1.25]);

subplot(2,2,4)
plot(Tt(1:count),ryy2(1:count),'black','LineWidth',1.5);
xlabel('Time [sec]','FontSize',12)
ylabel('Y coordinate [m]','FontSize',12)
grid on
axis([0 2 -0.6 0.6]);
pause(0.0001);
count = count+1;
end


