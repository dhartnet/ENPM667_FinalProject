%% Code for Part E - Output Vector Observability

clear
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.81;

% A matrix
A = [0 1 0 0 0 0; 
    0 0 -m1*g/M 0 -m2*g/M 0; 
    0 0 0 1 0 0; 
    0 0 (-g/l1)*((m1/M)+1) 0 (-m2*g)/(M*l1) 0; 
    0 0 0 0 0 1; 
    0 0 (-m1*g)/(M*l2) 0 (-g/l2)*((m2/M)+1) 0];

% C matrices (outputs)
C1 = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 0 0]; % (x)

C2 = [0 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 0 1 0]; % (theta1, theta2)

C3 = [1 0 0 0 0 0; 
    0 0 0 0 0 0; 
    0 0 0 0 1 0]; % (x, theta2)

C4 = [1 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 0 0 0 1 0]; % (x, theta1, theta2)

% Take transpose
A = A.';
C1 = C1.';
C2 = C2.';
C3 = C3.';
C4 = C4.';

% if n = 6, system is observable
n1 = rank([C1 A*C1 (A^2)*C1 (A^3)*C1 (A^4)*C1 (A^5)*C1])
n2 = rank([C2 A*C2 (A^2)*C2 (A^3)*C2 (A^4)*C2 (A^5)*C2])
n3 = rank([C3 A*C3 (A^2)*C3 (A^3)*C3 (A^4)*C3 (A^5)*C3])
n4 = rank([C4 A*C4 (A^2)*C4 (A^3)*C4 (A^4)*C4 (A^5)*C4])
