%% Code for Part G - LQG Simulation

clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

% Intitial Conditions
x0 = [0, 0, deg2rad(3), 0, deg2rad(-5), 0, 0, 0, 0, 0, 0, 0]; % Intitial Conditions

t_a = linspace(0, 60); %time vector

A = [0 1 0 0 0 0; 
    0 0 -m1*g/M 0 -m2*g/M 0; 
    0 0 0 1 0 0; 
    0 0 (-g/l1)*((m1/M)+1) 0 (-m2*g)/(M*l1) 0; 
    0 0 0 0 0 1; 
    0 0 (-m1*g)/(M*l2) 0 (-g/l2)*((m2/M)+1) 0]; % A matrix

B = [0; 
    1/M;
    0; 
    1/(M*l1); 
    0; 
    1/(M*l2)]; % B matrix

% C matrices (output matrices)
C = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0]; % (x)

D = 0; % D matrix

noise1 = 0.0001; % Variance
G = [noise1 0 0 0 0 0;
    0 noise1 0 0 0 0;
    0 0 noise1 0 0 0;
    0 0 0 noise1 0 0;
    0 0 0 0 noise1 0;
    0 0 0 0 0 noise1]; % process noise

QN = G; %  Process and sensor noise covariance matrices

noise2 = 0.001; % variance
RN = [noise2 0 0;
    0 noise2 0;
    0 0 noise2]; %  Process and sensor noise covariance matrices

[Lun1,P1,E1] = lqe(A, G, C, QN, RN); % L = luenberger observer

% Noise signals
W = sqrt(QN) * randn(6, length(t_a));
V = sqrt(RN) * randn(size(C, 1), length(t_a));

% Linear System - LQR

% Gain matrices
Q = [1500 0 0 0 0 0;
    0 0 0 0 0 0; 
    0 0 25000 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 20000 0; 
    0 0 0 0 0 0];

R = 0.001;

sys = ss(A, B, C, D);
[K,S,P] = lqr(sys,Q,R); % outputs K controller

% Non-Linear System

% function defining system at each time instance in ode45
function nonL = nonLinearModel1(t, t_a, a, K, M, m1, m2, l1, l2, g, L, W, V) 

% a is 12x6 state variable matrix, with first 6x6 part being the system state and the send 6x6 part being the state estimator/observer

x = a(1);
theta1 = a(3);
theta2 = a(5);

% F = U = -K*X
Uk = -K*[a(1);
    a(2);
    a(3);
    a(4);
    a(5);
    a(6);];

% Noise only in x state variable
Bd = [1;
    0;
    0;
    0;
    0;
    0];

% Signals over time
Wt = interp1(t_a, W', t)';  % Interpolated process noise
Vt = interp1(t_a, V', t)';  % Interpolated measurement noise

% Bd * Wt
BdWt = Bd.*Wt;

% Combined state matrices
nonL = zeros(12,1);

% non linear system (X)
b = (Uk - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l2*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2) + BdWt(1);
nonL(2) = b + BdWt(2);
nonL(3) = a(4) + BdWt(2);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3)) + BdWt(4);
nonL(5) = a(6) + BdWt(5);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5)) + BdWt(6);

Y = [a(1) + Vt(1);
    Vt(2);
    Vt(3)];

Xhat = [a(7);
    0;
    0];

obs = L*(Y - Xhat);

% Observer (X^)
nonL(7) = a(2) + obs(1);
nonL(8) = b + obs(2);
nonL(9) = a(4) + obs(3);
nonL(10) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3)) + obs(4);
nonL(11) = a(6) + obs(5);
nonL(12) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5)) + obs(6);
end

% Plots 
[t1,y1] = ode45(@(t,a)nonLinearModel1(t,t_a,a,K,M,m1,m2,l1,l2,g,Lun1,W,V),t_a,x0);

figure(1)
tiled = tiledlayout(3,1);

% Plot x
nexttile
plot(t1,y1(:,1))
yline(0)
title('X')
xlabel('Time (s)') 
ylabel('Meters')
grid on

% Plot theta1
nexttile
plot(t1,180*y1(:,3)/pi)
yline(0)
title('Theta1')
xlabel('Time (s)') 
ylabel('Degrees') 
grid on

% Plot theta1
nexttile
plot(t1,180*y1(:,5)/pi)
yline(0)
title('Theta2')
xlabel('Time (s)') 
ylabel('Degrees')
grid on

title(tiled, 'LQG Controller Results for Response to Initial Conditions and Noise Simulated on Nonlinear System')
