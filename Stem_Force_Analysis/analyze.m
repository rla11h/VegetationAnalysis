% ***************** EE POSITION IS ONLY PUBLISHING AT 11HZ  ********************

close all; clear all; clc

dt = 1/1000;
ftd = dlmread('Datasets/Control_Stem/test3.txt');

%% Extract Data
t = ftd(:,1) - ftd(1,1);
Fx = ftd(:,2) - ftd(1,2);
Fy = ftd(:,3) - ftd(1,3);
Fz = -(ftd(:,4) - ftd(1,4));
Tx = ftd(:,5);
Ty = ftd(:,6);
Tz = ftd(:,7);
Ex = ftd(:,8);
Ey = ftd(:,9);
Ez = -ftd(:,10);

% EE height only changes by a max of about 0.5cm. We are assuming this is
% negligible and hold height of the end-effector constant
h_ = Ez.*0 + Ez(1);
h = Ez(1);

%% Process Data Streams
FzLP = filter(0.01,[1 -0.99],Fz);   % Pass Z forces through a low-pass filter
FxLP = filter(0.01,[1 -0.99],Fx);   % Pass X forces through a low-pass filter
FyLP = filter(0.01,[1 -0.99],Fy);   % Pass Y forces through a low-pass filter

Fzs = smoothdata(Fz,'gaussian',1000); % Run Z force information through a gaussian smoothing algorithm
Fxs = smoothdata(Fx,'gaussian',1000); % Run X force information through a gaussian smoothing algorithm
Fys = smoothdata(Fy,'gaussian',1000); % Run Y force information through a gaussian smoothing algorithm

sum=0;
for i=1:1000              %First second of data stream
    sum = sum + Ex(i);
    Ex_offset = sum / i;  %Offset position data by the average position over the first second of data stream
end
dx = Ex - Ex_offset;
theta = atan2(dx,h_);  % stem angle
theta2 = asin(dx/h_);

%% Estimate Properties
e = 0.002;
lastup = 0; lastdn = 0; moving = 0; a=1; j=1; Fzsum =0; num=0; thsum=0; k = []; bsum = 0; prev_move_state = 0;
stiffnessPts = []; measuringPts = []; measuringDPts = []; measuringBPts = []; movingPts = []; find_k = 0;
for i=2:length(t)  % run through all data points causally
    
    if i==15000
        disp('here');
    end
    if dx(i) > dx(i-1)  % because the sampling rate is so low (FUCK!) we need to save the value of the last step to compare with
        lastup = dx(i-1);
        lastup_i = i;
        last_th = theta(i-1);
    elseif dx(i) < dx(i-1)  % because the sampling rate is so low (FUCK!) we need to save the value of the last step to compare with
        lastdn = dx(i-1);
        lastdn_i = i;
    %    last_th = theta(i-1);
    end
    if (abs(dx(i)-lastup) < e) || (abs(dx(i)-lastdn) < e)% Checks for no movement with a sampling rate of 11Hz and publishing rate of 1kHz
   % if i > 100 && (abs(dx(i)-dx(i-100)) < e)
          moving = 0;
    else
        moving = 1;
      %  e = dx(i)/25;
        if prev_move_state == 0     % Log the time at the start of a forward motion
            t_start = t(i);
            H_start = theta(i);
        end
        movingPts = [movingPts; [t(i), dx(i)]];
    end
    
    % Average the z forces and angles during each hold period
    if dx(i) > 0.05 && moving == 0   % if arm is not moving and not in rest position
        num=num+1;
        Fzsum = Fzsum + FzLP(i);
        Fzavg = Fzsum/num;
        thsum = thsum+theta(i);
        thavg = thsum/num;
        find_k = 1; % Throw the flag to estimate the stiffness
        measuringPts = [measuringPts; [t(i) FzLP(i)]];
    else
        if find_k == 1          %Calculate an approximate stiffness value once forces have been averaged over a holding period
            find_k = 0; % reset flag
            k(a) = (Fzavg*h/thavg)/cos(thavg)^2; % approximate stiffness
            rad(a) = thavg;
            stiffnessPts = [stiffnessPts; [t(i) FzLP(i)]];
            a=a+1;
        end
        num=0;  % Clear sum and counter
        Fzsum=0;
    end
    
     masum = 0; maLength = 100; 
    if length(k) == 2    % If the two hold phases have been characterized begin looking for the start of the thrust       
        if dx(i) > lastup && moving == 1
            k_current = StiffnessProfile(theta(i), k, rad);        % Calculate stiffness based on the profile of the stem
            F_current = 1/h*k_current*theta(i)*cos(theta(i))^2;    % Calculate expected force due to stiffness alone
            F_dam = FzLP(i) - F_current;                           % Residual force due to damping component
            if F_dam < 0                                           % Prevent impossible values
                F_dam = 0;
            end
            for m=-maLength/2:maLength/2   % This is non causal. Cant do real-time
         %   for m=-maLength:0               % This is causal
                masum = masum + theta(i+m);
            end
            MAtheta = masum/maLength;
            if MAtheta < H_start      % Fixes some weirdness caused by adding the MA filter 
                MAtheta = H_start;
            end
           % thdot = (theta(i)-H_start)/(t(i)-t_start);      %  -------------> ** WARNING: ROUGH APPROXIMATION MAY LEAD TO ERROR ** <--------------
            thdot = (MAtheta-H_start)/(t(i)-t_start);        %  -------------> ** WARNING: NONCAUSAL CALCULATION IMPOSSIBLE IN REALTIME ** <--------------
            %thdot = abs(thdot); % maybe(?)
            bcurr = (F_dam*h/thdot)/cos(MAtheta)^2;
            if ~(isnan(bcurr) || bcurr == Inf)      % Log damping when damping has real value (bad values came from t(i)-t_start = 0 on first iteration
                b(j) = bcurr; 
                Hd(j) = thdot;
                bsum = bsum + b(j);
                bavg = bsum/j;
                j=j+1;  
                measuringBPts = [measuringBPts; [t(i) FzLP(i)]];
            end
            measuringDPts = [measuringDPts; [t(i) FzLP(i)]];
        end
       
    end
    
    prev_move_state = moving;
end

%%
figure(1)
subplot(3,1,1)
%plot(t,fz); hold on;
plot(t,FzLP,'r--','LineWidth',1);
hold on;
for i=1:length(k)
    text(stiffnessPts(i,1)-4,stiffnessPts(i,2)+0.1,['K = ' num2str(k(i)) ' N/rad']);
end

text(stiffnessPts(3,1)-10, stiffnessPts(3,2)+0.1, ['B = ' num2str(bavg) ' N*s/rad']);

for i=1:length(measuringDPts)
    hold on;
    plot(measuringDPts(i,1),measuringDPts(i,2),'o');
end
% ylim([-5 10]);
xlim([t(1) t(end)])
xlabel('Time (sec)','interpreter','latex');
ylabel('Force (N)','interpreter','latex');
title('Force Z-Direction','interpreter','latex');
grid on;

subplot(3,1,2)
plot(t,dx); hold on;
for i=1:length(movingPts)
    hold on;
    plot(movingPts(i,1),movingPts(i,2),'o');
end
%plot(fs(:,1),fs(:,2))
xlim([t(1) t(end)])
xlabel('Time (sec)','interpreter','latex');
ylabel('Displacement (m)','interpreter','latex');
title('Gripper Displacment X-Direction','interpreter','latex');
grid on;

subplot(3,1,3)
theta = atan2(dx,Ez);
plot(t,radtodeg(theta))
ylabel('Displacement (deg)','interpreter','latex');
yyaxis right
plot(t,theta)
xlim([t(1) t(end)]);
xlabel('Time (sec)','interpreter','latex');
ylabel('Displacement (rad)');
title('Gripper Displacment Radians','interpreter','latex');
grid on;

figure(2);
plot(measuringBPts(:,1),Hd)
%{
figure(2)
hold on
plot(t,FzLP)
plot(t,FxLP)
plot(t,FyLP)
xlim([t(1) t(end)])
xlabel('Time (sec)','interpreter','latex');
ylabel('Force (N)','interpreter','latex');
title('Force MultiDirectional','interpreter','latex');
legend('Z axis','X axis','Y axis');
grid on;
%}

