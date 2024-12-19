%% Code for Part F - Linear and Nonlinear Observers

clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

% Intitial Conditions
x0 = [1, 0, deg2rad(3), 0, deg2rad(-5), 0, 0, 0, 0, 0, 0, 0]; 
t_a = linspace(0, 60); %, 0.5); %time vector

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
C1 = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0]; % (x)

C2 = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 1 0]; % (x, theta2)

C3 = [1 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 0 1 0]; % (x, theta1, theta2)

D = 0; % D matrix

noise1 = 0.1; % Covariance
G = [noise1 0 0 0 0 0;
    0 noise1 0 0 0 0;
    0 0 noise1 0 0 0;
    0 0 0 noise1 0 0;
    0 0 0 0 noise1 0;
    0 0 0 0 0 noise1]; % process noise

QN = G; %  Process and sensor noise covariance matrices

noise2 = 0.0001; % Covariance
RN = [noise2 0 0;
    0 noise2 0;
    0 0 noise2]; %  Process and sensor noise covariance matrices

% L = luenberger observer
[Lun1,P1,E1] = lqe(A, G, C1, QN, RN); 
[Lun2,P2,E2] = lqe(A, G, C2, QN, RN);
[Lun3,P3,E3] = lqe(A, G, C3, QN, RN);

% Non-Linear System Observers

% a = [X, X^] = [state, observer]
% Function for first observer state
function nonL = nonLinearModel1(a, F, M, m1, m2, l1, l2, g, L) % observer matrix is just x
x = a(1);
theta1 = a(3);
theta2 = a(5);

C = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0];

nonL = zeros(6,1);

% non linear system - nonL = X'(t)
b = (F - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l2*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2);
nonL(2) = b;
nonL(3) = a(4);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3));
nonL(5) = a(6);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5));

X = [a(1);
    a(2);
    a(3);
    a(4);
    a(5);
    a(6);];

Y = [a(1);
    0;
    0];

Xhat = [a(7);
    0;
    0];

obs = L*(Y - Xhat);

% X^' = AX^ + Bk*Uk + L(Y-CX^)
nonL(7) = a(2) + obs(1);
nonL(8) = b + obs(2);
nonL(9) = a(4) + obs(3);
nonL(10) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3)) + obs(4);
nonL(11) = a(6) + obs(5);
nonL(12) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5)) + obs(6);
end

% a = [X, X^] = [state, observer]
% Function for second observer state
function nonL = nonLinearModel2(a, F, M, m1, m2, l1, l2, g, L) % observer matrix is x, theta 2
x = a(1);
theta1 = a(3);
theta2 = a(5);

C = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 1 0];

nonL = zeros(6,1);

% non linear system - nonL = X'(t)
b = (F - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l2*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2);
nonL(2) = b;
nonL(3) = a(4);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3));
nonL(5) = a(6);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5));

X = [a(1);
    a(2);
    a(3);
    a(4);
    a(5);
    a(6);];

Y = [a(1);
    0;
    a(5)];

Xhat = [a(7);
    0;
    a(11)];

obs = L*(Y - Xhat);

% X^' = AX^ + Bk*Uk + L(Y-CX^)
nonL(7) = a(2) + obs(1);
nonL(8) = b + obs(2);
nonL(9) = a(4) + obs(3);
nonL(10) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3)) + obs(4);
nonL(11) = a(6) + obs(5);
nonL(12) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5)) + obs(6);
end

% a = [X, X^] = [state, observer]
% Function for third observer state
function nonL = nonLinearModel3(a, F, M, m1, m2, l1, l2, g, L) % observer matrix is x, theta 1, theta 2
x = a(1);
theta1 = a(3);
theta2 = a(5);

C = [1 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 0 1 0];

nonL = zeros(6,1);

% non linear system - nonL = X'(t)
b = (F - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l2*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2);
nonL(2) = b;
nonL(3) = a(4);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3));
nonL(5) = a(6);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5));

X = [a(1);
    a(2);
    a(3);
    a(4);
    a(5);
    a(6);];

Y = [a(1);
    a(3);
    a(5)];

Xhat = [a(7);
    a(9);
    a(11)];

obs = L*(Y - Xhat);

% X^' = AX^ + Bk*Uk + L(Y-CX^)
nonL(7) = a(2) + obs(1);
nonL(8) = b + obs(2);
nonL(9) = a(4) + obs(3);
nonL(10) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3)) + obs(4);
nonL(11) = a(6) + obs(5);
nonL(12) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5)) + obs(6);
end

% Plots 
[t1,y1] = ode45(@(t,a)nonLinearModel1(a,0,M,m1,m2,l1,l2,g,Lun1),t_a,x0);

figure(1)
tiled1 = tiledlayout(3,1);

% Plot x
nexttile
plot(t1,y1(:,7))
yline(0)
title('X')
xlabel('Time (s)') 
ylabel('Meters')
grid on

% Plot theta1
nexttile
plot(t1,180*y1(:,9)/pi)
yline(0)
title('Theta1')
xlabel('Time (s)') 
ylabel('Degrees') 
grid on

% Plot theta1
nexttile
plot(t1,180*y1(:,11)/pi)
yline(0)
title('Theta2')
xlabel('Time (s)') 
ylabel('Degrees')
grid on

title(tiled1, '[X(t)] Observer for Response to Initial Conditions for Nonlinear System')

[t2,y2] = ode45(@(t,a)nonLinearModel2(a,0,M,m1,m2,l1,l2,g,Lun2),t_a,x0);

figure(2)
tiled2 = tiledlayout(3,1);

% Plot x
nexttile
plot(t2,y2(:,7))
yline(0)
title('X')
xlabel('Time (s)') 
ylabel('Meters')
grid on

% Plot theta1
nexttile
plot(t2,180*y2(:,9)/pi)
yline(0)
title('Theta1')
xlabel('Time (s)') 
ylabel('Degrees') 
grid on

% Plot theta1
nexttile
plot(t2,180*y2(:,11)/pi)
yline(0)
title('Theta2')
xlabel('Time (s)') 
ylabel('DegreesS')
grid on

title(tiled2, '[X(t), Theta2(t)] Observer for Response to Initial Conditions for Nonlinear System')

[t3,y3] = ode45(@(t,a)nonLinearModel3(a,0,M,m1,m2,l1,l2,g,Lun3),t_a,x0);

figure(3)
tiled3 = tiledlayout(3,1);

% Plot x
nexttile
plot(t3,y3(:,7))
yline(0)
title('X')
xlabel('Time (s)') 
ylabel('Meters')
grid on

% Plot theta1
nexttile
plot(t3,180*y3(:,9)/pi)
yline(0)
title('Theta1')
xlabel('Time (s)') 
ylabel('Degrees') 
grid on

% Plot theta1
nexttile
plot(t3,180*y3(:,11)/pi)
yline(0)
title('Theta2')
xlabel('Time (s)') 
ylabel('Degrees')
grid on

title(tiled3, '[X(t), Theta1(t), Theta2(t)] Observer for Response to Initial Conditions for Nonlinear System')
