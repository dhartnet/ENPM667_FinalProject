clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

syms x dx theta1 dtheta1 theta2 dtheta2

x0 = [0, 0, deg2rad(3), 0, deg2rad(-5), 0]; % Intitial Conditions

t_a = linspace(0, 30); %, 0.5); %time vector

% Linear System LQR to get K matrix

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

Q = [20 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 7500 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 5000 0; 
    0 0 0 0 0 0];

R = 0.00001;

sys = ss(A, B, C, D);
[K,S,P] = lqr(sys,Q,R);

