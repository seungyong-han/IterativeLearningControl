%% Satellite Example from Robust Lecture

clc 
clear all


s = zpk('s');
Ts = 1e-3;

G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));

Gd = c2d(G,Ts);
[A,B,C,D] = ssdata(Gd);


systemnames='Gd';

time = .1; %[s]
per = 3; 
%% define time horizon and time values
N = time/Ts;
time = 0:Ts:time;

%% Solve with LQR
Ex3_LQR;

r = r_vec(2:end)
u = u_sv(2:end)
[An,Bn,Cn,Dn, Nnew] = get_non0D_system(A - B*F, B, C-D*F, D, N);
[An,Bn,Cn,Dn, Nnew] = get_non0D_system(A, B, C, D, N);

beta = .0001;[G, d] = get_G(An, Bn, Cn, Dn, x0, Nnew);

[u_inf, e_inf, y_inf, impr,iteration_number, error_history] = RIA(G,d, beta,r, u, 0);%SDA(G,d, beta,r, u, 1, 1, 0);


%% Define start values for the algorithms
R = 1;
Q = 1;
Nmax = 100;


%%
do_plot = 1;
%Ex3_IA;
cont = 10; 

Ex3_IA; 




