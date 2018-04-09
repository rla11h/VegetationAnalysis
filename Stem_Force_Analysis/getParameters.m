close all; clear all; clc

dt = 1/1000;
ftd = dlmread('02_14_18/test1.txt');
t = ftd(:,1) - ftd(1,1);
fx = ftd(:,2);
fy = ftd(:,3);
fz = -(ftd(:,4) - ftd(1,4));
tx = ftd(:,5);
ty = ftd(:,6);
tz = ftd(:,7); 
ex = ftd(:,8);
ey = ftd(:,9);
ez = -ftd(:,10);
dx = ex - ex(1);
C = smoothdata(fz,'gaussian',1000);


subplot(3,1,1)
%plot(t,fz); hold on;
plot(t,C,'LineWidth',1);
% ylim([-5 10]);
xlim([t(1) t(end)])
xlabel('Time (sec)','interpreter','latex');
ylabel('Force (N)','interpreter','latex');
title('Force Z-Direction','interpreter','latex');
legend('Signal');
grid on;

subplot(3,1,2)
plot(t,dx); hold on;
xlim([t(1) t(end)])
xlabel('Time (sec)','interpreter','latex');
ylabel('Displacement (m)','interpreter','latex');
title('Gripped Displacment X-Direction','interpreter','latex');
grid on;

gca = subplot(3,1,3);
theta = atan2(dx,ez);
plot(t,theta)
xlim([t(1) t(end)]);

[ti, thetaI] = ginput(1);
[tf, thetaF] = ginput(1);
%[tfs, theatFS] = ginput(1);

n = find(t >= tf,1);
F = C(n);
H = ez(n);
D = dx(n);
theta = abs(theta(n));
thetadot = abs((thetaF - thetaI) / (tf-ti));
K = F*H/theta;
B = F*H/thetadot;
J = 0.12333;

K = F * H / (cos(theta)^2 * theta);
B = F*H/thetadot;
J = 0.12333;

t=0:0.01:10;
y0 =[pi/4 0];
[T,Y] = ode45(@(t,theta) stem(t,theta,J,B,K),t,y0);
plot(T,Y(:,1))

