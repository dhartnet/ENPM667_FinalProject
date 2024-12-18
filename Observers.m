clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

% Initial conditions (in degrees)
x0Deg = [0, 0, 3, 0, -5, 0];

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

noise1 = 0.1; % covariance
G = [noise1 0 0 0 0 0;
    0 noise1 0 0 0 0;
    0 0 noise1 0 0 0;
    0 0 0 noise1 0 0;
    0 0 0 0 noise1 0;
    0 0 0 0 0 noise1]; % process noise

QN = G; %  Process and sensor noise covariance matrices

noise2 = 0.0001; % covariance
RN = [noise2 0 0;
    0 noise2 0;
    0 0 noise2]; %  Process and sensor noise covariance matrices

% Evaluate disturbances and system and get L = luenberger observer
[L1,P1,E1] = lqe(A, G, C1, QN, RN);
[L2,P2,E2] = lqe(A, G, C2, QN, RN);
[L3,P3,E3] = lqe(A, G, C3, QN, RN);

% Define Ac as (A - LC)
Ac1 = A-L1*C1;
Ac2 = A-L2*C2;
Ac3 = A-L3*C3;

% Observer equation = (A-LC)X + Bk*Uk
% Y = CX
% Initial conditions for C1
figure(1)
sys1 = ss(Ac1, B, C1, D); % State estimation
ip1 = initialplot(sys1, x0Deg);
title('[X(t)] Observer for Response to Initial Conditions')
ip1.OutputNames = ['X (Meters)'; 'Theta1 (Degrees)'; 'Theta2 (Degrees)'];
grid on

% Initial conditions for C2
figure(2)
sys2 = ss(Ac2, B, C2, D); % State estimation
ip2 = initialplot(sys2, x0Deg);
title('[X(t), Theta2(t)] Observer for Response to Initial Conditions')
ip2.OutputNames = ['X (Meters)'; 'Theta1 (Degrees)'; 'Theta2 (Degrees)'];
grid on

% Initial conditions for C3
figure(3)
sys3 = ss(Ac3, B, C3, D); % State estimation
ip3 = initialplot(sys3, x0Deg);
title('[X(t), Theta1(t), Theta2(t)] Observer for Response to Initial Conditions')
ip3.OutputNames = ['X (Meters)'; 'Theta1 (Degrees)'; 'Theta2 (Degrees)'];
grid on

% Step for C1
figure(4)
st1 = stepplot(sys1);
output_names = get(st1, 'OutputName');
output_names(1) = 'X (Meters)';
output_names(2) = 'Theta1 (Degrees)';
output_names(3) = 'Theta2 (Degrees)';
set(st1, 'OutputName', output_names)
title('[X(t)] Observer for Response to Step Input')
grid on

% Step for C2
figure(5)
st2 = stepplot(sys2);
output_names2 = get(st2, 'OutputName');
output_names2(1) = 'X (Meters)';
output_names2(2) = 'Theta1 (Degrees)';
output_names2(3) = 'Theta2 (Degrees)';
set(st2, 'OutputName', output_names2)
title('[X(t), Theta2(t)] Observer for Response to Step Input')
grid on

% Step for C3
figure(6)
st3 = stepplot(sys3);
output_names3 = get(st3, 'OutputName');
output_names3(1) = 'X (Meters)';
output_names3(2) = 'Theta1 (Degrees)';
output_names3(3) = 'Theta2 (Degrees)';
set(st3, 'OutputName', output_names3)
title('[X(t), Theta1(t), Theta2(t)] Observer for Response to Step Input')
grid on
