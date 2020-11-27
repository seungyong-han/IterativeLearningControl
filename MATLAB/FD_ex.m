%% Set paths & clear all
addpath('help_fcns');
addpath('Algorithms');
addpath('scripts');
%%

clc
clear all
A = [2 1; 4 3];
B = [1; 2];
C = [0 1];
D = 2;

N = 430;
per = 5;
b = .7;  
sample_time = 1e-1;
% Find controller using separation principle
Ex1_LQR;
A = A - B*F;
C = C - D*F; 
%%
%----------------------------------------------
sys = ss(A,B,C,D,sample_time); 
systemnames='sys';
inputvar='[r;u]';
outputvar='[sys-r; u; r-sys]';
input_to_sys='[u]';
P=sysic;

% Slightly adjusted weight from RC lecture transformed to discrete-time 
s = tf('s');
M = 10;
wb = 0.75;
a = 1e-3;

Wp = c2d((s/sqrt(M) + wb)^2/(s + wb*sqrt(a))^2, sample_time);

% H infty design for weighted plant
[Kc, ~, ga] = hinfsyn(blkdiag(Wp, 0.001, 1) * P, 1, 1);

cloop = lft(P, Kc);
S = -cloop(1, 1); %sensitivity here
T = 1 - S; 
K = Kc*S; 
%% Define the supermatrices 
[A,B,C,D] = ssdata(sys); 
[G, d] = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(T); 
T_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(S); 
S_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(K); 
K_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 

%% FB with matrices 
clc 
u = u_sv;
r = r_vec; 
b = .7; 

e = r - G*u - d; 
cont = 1; 
error_history = norm(e); 
iteration_number = 0; 

while cont
   
    iteration_number = iteration_number + 1; 
    u_new = u + b*K_sv*T_sv'*e;
    e_new = r -  G*u_new - d; 
    
    
    if norm(e - e_new)<1e-6
        cont = 0; 
    end
    
    if mod(iteration_number, 1000) == 0
         disp(['curr_error_difference: ', num2str(norm(e_new - e))]); 
    end
    error_history = [error_history; norm(e)]; 

    e = e_new;
    u = u_new;
end

 


