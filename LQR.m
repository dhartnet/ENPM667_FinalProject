M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

syms x dx theta1 dtheta1 theta2 dtheta2

% Initial conditions
x0 = [0, 0, 3, 0, -5, 0]; % Intitial Conditions
t_a = linspace(0, 0.5, 30); %time vector

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

D = [0;
0;
0]; 

% Gain matrices
Q = [20 0 0 0 0 0; 
0 0 0 0 0 0; 
0 0 7500 0 0 0; 
0 0 0 0 0 0; 
0 0 0 0 5000 0; 
0 0 0 0 0 0];

R = 0.00001;

% Evaluate system, get K controller and plot system with U = -KX
sys = ss(A, B, C, D);
[K,S,P] = lqr(sys,Q,R);
sys1 = ss(A-B*K, B, C, D);
initial(sys1, x0)

% Eigen values for stability test
eig(A-B*K) % Lyapunov's indirect stability test

% Non-Linear System - LQR

a = [x dx theta1 dtheta1 theta2 dtheta2];

% Funcion defining system
function nonL = nonLinearModel(a, F, M, m1, m2, l1, l2, g)

% x = a(1);
% dx = a(2);
% theta1 = a(3);
% dtheta1 = a(4)
% theta2 = a(5);
% dtheta2 = a(6);

nonL = zeros(6,1);

b = (F - m1*(g*cos(a(3))*sin(a(3)) + l1*sin(a(3))*(a(4)^2)) - m2*(g*cos(a(5))*sin(a(5)) + l1*sin(a(5))*(a(6)^2))) / (M + m1*(sin(a(3))^2) + m2*(sin(a(5))^2));
nonL(1) = a(2);
nonL(2) = b;
nonL(3) = a(4);
nonL(4) = (1/l1)*cos(a(3))*b - (g/l1)*sin(a(3));
nonL(5) = a(6);
nonL(6) = (1/l2)*cos(a(5))*b - (g/l2)*sin(a(5));
end

% Simulate and plot
[~,y] = ode45(@(t,a)nonLinearModel(a,K*a,M,m1,m2,l1,l2,g),t_a,x0);

figure(1);
hold on
plot(t_a,y(:,1))
plot(t_a,y(:,3))
plot(t_a,y(:,5))
ylabel('state variables')
xlabel('time (sec)')
title('Non-Linear system using LQR controller')
legend('x','theta1','theta2')
