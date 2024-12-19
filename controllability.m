%% Code for Part C - Controllability Conditions

clear
syms m1;
syms m2;
syms M;
syms l1;
syms l2;
syms g;

% State Space representation matrices
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

% defining controllability matrix columns
e2 = A*B;
e3 = A^(2)*B;
e4 = A^(3)*B;
e5 = A^(4)*B;
e6 = A^(5)*B;

% Controllability matrix
E = [B(1,1) e2(1,1) e3(1, 1) e4(1, 1) e5(1, 1) e6(1, 1); 
    B(2,1) e2(2,1) e3(2, 1) e4(2, 1) e5(2, 1) e6(2, 1); 
    B(3,1) e2(3,1) e3(3, 1) e4(3, 1) e5(3, 1) e6(3, 1); 
    B(4,1) e2(4,1) e3(4, 1) e4(4, 1) e5(4, 1) e6(4, 1); 
    B(5,1) e2(5,1) e3(5, 1) e4(5, 1) e5(5, 1) e6(5, 1); 
    B(6,1) e2(6,1) e3(6, 1) e4(6, 1) e5(6, 1) e6(6, 1);]

% determinant to then see conditions for determinant not equal to 0
determinant = det(E)
