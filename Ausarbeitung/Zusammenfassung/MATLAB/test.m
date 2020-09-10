%% Test IA, SDA, suspession of evs, beta as evs, resuced system
% Q, R = eye
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


%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));

u0 = repmat(0.1,l*(N+1),1);
r = 5*rand(m*(N+1),1);

R = eye(l*(N+1));
Q = eye(m*(N+1));

beta = .1;

L = eye(length(G)) - beta*G*G';


do_plot = 1;

[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = IA(G,d, beta,r, u0, R, Q,do_plot)
title('Inverse model algorithm')
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA(G,d, beta,r, u0, R, Q, do_plot)
title('Steepest Descent Algorithm')
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
title('Suppression of evs')
[u_inf4, e_inf4, y_inf4, impr4, iteration_number4] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)
title('beta = evs')
%Reduced system :( 
figure
order = 1;
[u_inf5, e_inf5, y_inf5, impr5, iteration_number5] = GA_reducedSystem(A, B, C, D, x0, N, order, 1, beta, r, u0, R, Q,do_plot)
title('Reduced System')                         

%% Test Test IA, SDA, suspession of evs, beta as evs, resuced system
% (Alg. 7.3, 7.4, 7.5) -> list isn't full
% Q, R = eye
%   GG* singular. Mostly: error doesn't converge to 0, because the initial
%   input signal u0 doesn't lie in image(GG*)
N = 10;
A = diag([-.2, -.14]);
B = [1 0 3; 0 2 1]; 
C = [.1 .5; -.4 .1; .3 -.8];
D = diag([1,4,0]);
x0 = [0, 0]';

%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));


[G, d] = get_G(A, B, C, D, x0, N);
u0 = repmat(0.1,l*(N+1),1);
%r = G*G'*5*rand(l*(N+1),1) + G*u0 + d; %to ensure that e0 is in image(GG*)
r = 5*rand(m*(N+1),1);

R = eye(l*(N+1));
Q = eye(m*(N+1));
beta = .1


do_plot = 1;

[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = IA(G,d, beta,r, u0, R, Q,do_plot)
title('Inverse model algorithm')
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA(G,d, beta,r, u0, R, Q, do_plot)
title('Steepest Descent Algorithm')
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
title('Suppression of evs')
[u_inf4, e_inf4, y_inf4, impr4, iteration_number4] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)
title('beta = evs')
order = 1;
[u_inf5, e_inf5, y_inf5, impr5, iteration_number5] = GA_reducedSystem(A, B, C, D, x0, N, order, 1, beta, r, u0, R, Q,do_plot)
title('Reduced System')  


%% Test SDA (Alg. 7.3, 7.4, 7.5)
%GG* singular
% R =kron(eye(N+1),diag([5,3,2]));
% Q = kron(eye(N+1),diag([.2,1,1.8]));


clear all
clc
N = 20;
A = diag([-.2, -.2]);
B = [1 0 3; 0 2 1]; 
C = [.1 .5; -.4 .1; .3 -.8];
D = diag([1,4,0]);
x0 = [0, 0]';


%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));


[G, d] = get_G(A, B, C, D, x0, N);
u0 = repmat(0.1,l*(N+1),1);
r = G*G'*5*rand(m*(N+1),1) + G*u0 + d;%to ensure that e0 is in image(GG*)
%r = 5*rand(l*(N+1),1);
%r = [0 0 -1, zeros(1,9)]';

%R = eye(m*(N+1));
%Q = eye(l*(N+1));

R =kron(eye(N+1),diag([5,3,2]));
Q = kron(eye(N+1),diag([.2,1,1.8]));
beta = .2;


do_plot = 1;

[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = IA(G,d, beta,r, u0, R, Q,do_plot)
title('Inverse model algorithm')
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA(G,d, beta,r, u0, R, Q, do_plot)
title('Steepest Descent Algorithm')
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
title('Suppression of evs')
[u_inf4, e_inf4, y_inf4, impr4, iteration_number4] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)
title('beta = evs')
order = 1;
[u_inf5, e_inf5, y_inf5, impr5, iteration_number5] = GA_reducedSystem(A, B, C, D, x0, N, order, 1, beta, r, u0, R, Q,do_plot)
title('Reduced System')  
%balred diverges 
[u_inf6, e_inf6, y_inf6, imlpr6, iteration_number6] = epsSDA(G,d,N,l,m,beta,r, u0, R, Q,do_plot)

%% Test Model random G
%Alg 7.8, 7.9 
clc
clear all


rand_sys = drss(10, 5, 10);
rand_sys = prescale(rand_sys);
rand_sys = minreal(rand_sys);
[A, B, C, D, ts] = ssdata(rand_sys);


%rank(D, Dk) must be min{m,l}
x0 = zeros(length(A), 1);
N = 5;



%l-input m-output system; in supervector description l*(N+1) - input
%m*(N+1)-output
l = length(B(1,:));
m = length(C(:,1));


u0 = repmat(0.1,l*(N+1),1);
%r = G*G'*5*rand(12,1) + G*u0 + d;%to ensure that e0 is in image(GG*)
r = 5*rand(m*(N+1),1);

R = eye(l*(N+1));

Q = eye(m*(N+1));
[G, d] = get_G(A, B, C, D, x0, N);


beta = .2

do_plot = 1;

%Possibly needs a lot of iterations -> reduce number of states/N for shorter
%calculating times, eg number of states = 10 and N = 5
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = IA(G,d, beta,r, u0, R, Q,do_plot)
title('Inverse model algorithm')
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2] = SDA(G,d, beta,r, u0, R, Q, do_plot)
title('Steepest Descent Algorithm')
[u_inf3, e_inf3, y_inf3, impr3, iteration_number3] = SDA_suppression_of_evs(G,d, r, u0, R, Q,do_plot)
title('Suppression of evs')
[u_inf4, e_inf4, y_inf4, impr4, iteration_number4] = SDA_beta_as_evs(G,d, r, u0, R, Q,do_plot)
title('beta = evs')
order = 1;%Order for reduced system with balred 
[u_inf5, e_inf5, y_inf5, impr5, iteration_number5] = GA_reducedSystem(A, B, C, D, x0, N, order, 1, beta, r, u0, R, Q,do_plot)
title('Reduced System')  

%%
% Model reduction algoriths are in scripts
% MrPZS % Model reduction Pole Zero Simplification (->minreal)
% MrModeSelection % Model reduction Mode selection
% Balred is already implemented in GA_reducedSystem (doesnt work well :( ) 




