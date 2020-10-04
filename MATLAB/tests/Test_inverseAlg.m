%% Test SDA
%  GG* non singular

clc
clear all 

N = 200;
A = diag([-.5, -.3]);
B = [1,1]';
C = [.1, .2];
D = 1;
x0 = [0, 0]';
[G, d] = get_G(A, B, C, D, x0, N);


%m-input l-output system; in supervector description m*(N+1) - input
%l*(N+1)-output
m = length(B(1,:));
l = length(C(:,1));

u0 = repmat(0.1,m*(N+1),1);
r = 5*rand(l*(N+1),1);

R = eye(m*(N+1));
Q = eye(l*(N+1));

beta = .1;

L = eye(length(G)) - beta*G*G';


do_plot = 1;
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = IA(G,d, beta,r, u0, R, Q,do_plot)
title('Right inverse model algorithm')
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA(G,d, beta,r, u0, R, Q, do_plot)
title('Steepest Descent Algorithm')
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
title('Suppression of evs')

