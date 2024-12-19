clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

syms x dx theta1 dtheta1 theta2 dtheta2
t_a = linspace(0, 45); %time vector

x0 = [0, 0, deg2rad(3), 0, deg2rad(-5), 0]; % Intitial Conditions
x0Deg = [0, 0, 3, 0, -5, 0]; % Intitial Conditions


% Linear System - LQR
A = [0 1 0 0 0 0; 
    0 0 -m1*g/M 0 -m2*g/M 0; 
    0 0 0 1 0 0; 
    0 0 (-g/l1)*((m1/M)+1) 0 (-m2*g)/(M*l1) 0; 
    0 0 0 0 0 1; 
    0 0 (-m1*g)/(M*l2) 0 (-g/l2)*((m2/M)+1) 0];

B = [0; 
    1/M;
    0; 
    1/(M*l1); 
    0; 
    1/(M*l2)];

C = [1 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 0 1 0];

D=0;

% Gain matrices
Q = [20 0 0 0 0 0;
    0 0 0 0 0 0; 
    0 0 7500 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 5000 0; 
    0 0 0 0 0 0];

R = 0.00001;

% Simulate LQR on linear system
figure(1)
sys = ss(A, B, C, D);
[K,S,P] = lqr(sys,Q,R);
sys1 = ss(A-B*K, B, C, D);
ip = initialplot(sys1, x0);
title('Response to Initial Conditions (Linear System)')
ip.OutputNames = ["X (Meters)"; 'Theta1 (Radians)'; 'Theta2 (Radians)'];
grid on

% Lyapunovs indirect method stability test
A_BK = A+B*K
eigenValues = eig(A-B*K)

% Non-Linear System - LQR

function nonL = nonLinearModel(a, F, M, m1, m2, l1, l2, g)
% x = a(1);
% dx = a(2);
% theta1 = a(3);
% dtheta1 = a(4)
% theta2 = a(5);
% dtheta2 = a(6);
nonL = zeros(6,1);

% non linear system - nonL = X'(t)
b = (F - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l2*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2);
nonL(2) = b;
nonL(3) = a(4);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3));
nonL(5) = a(6);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5));
end

[t,y] = ode45(@(t,a)nonLinearModel(a,-K*a,M,m1,m2,l1,l2,g),t_a,x0);

figure(2);
tiled = tiledlayout(3,1);

% Plot x
nexttile
plot(t,y(:,1))
yline(0)
title('x')
xlabel('Time (s)') 
ylabel('X (Meters)')
grid on

% Plot theta1
nexttile
plot(t,y(:,3))
yline(0)
title('Theta1')
xlabel('Time (s)') 
ylabel('Theta1 (Radians)') 
grid on

% Plot theta1
nexttile
plot(t,y(:,5))
yline(0)
title('Theta2')
xlabel('Time (s)') 
ylabel('Theta2 (Radians)')
grid on

title(tiled,'Response to Initial Conditions (Non-Linear System)')
