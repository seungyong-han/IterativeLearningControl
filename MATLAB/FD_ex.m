%% Set paths & clear all
addpath('help_fcns');
addpath('Algorithms');
addpath('scripts');

clc 
clear all

%% Define system and controller
% Define system (from Lecture)
s = zpk('s');
sample_time = 1e-3;
sys = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
N = 10;

%Define plant 
sys = c2d(sys,sample_time);
systemnames='sys';
inputvar='[r;u]';
outputvar='[sys-r; r-sys]';
input_to_sys='[u]';
P=sysic;

% Define Controller (from Lecture)
K = ss(7.9212*(s+0.1818)*(s^2-0.2244*s+0.8981)/((s^2+3.899*s+4.745)*(s^2+1.039*s+3.395)));
Kd = c2d(K,sample_time);
Kc = Kd;

% Define closed loop and sensitivity/compl. sensitivily 
cloop = lft(P, Kc);

S = -cloop;
T = 1 - S; 

%% Get supervector form

Nmin = N; % The minimal possible time horizon over all systems
[A, B, C, D] = ssdata(K);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
K_sv = get_G(A, B, C, D, zeros(length(A)), N1); 
if N1<Nmin
    Nmin = N1;
end
    
[A,B,C,D] = ssdata(T);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
T_sv = get_G(A,B,C,D,zeros(length(A)), N1);
if N1<Nmin
    Nmin = N1;
end
[A,B,C,D] = ssdata(sys);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
G = get_G(A,B,C,D,zeros(length(A)), N1);
if N1<Nmin
    Nmin = N1;
end

%% Feedback design 
cont = 1; %while-loop criterion 

% Compensate the matrices, sth all of them have the same dimension 
T_sv = T_sv(length(T_sv) - Nmin + 1:end, length(T_sv) - Nmin + 1:end);
K_sv = K_sv(length(K_sv) - Nmin + 1:end, length(K_sv) - Nmin + 1:end);
G    = G   (length(G)    - Nmin + 1:end, length(G)    - Nmin + 1:end);
T_sv_star = T_sv';

r = 3*ones(length(T_sv),1);%tracking value
e = (eye(length(T_sv)) - T_sv)*r;%initial error
u = zeros(length(T_sv),1);%initial input = 0

%Calculate beta
beta = .2;
beta = beta/norm(T_sv)/norm(T_sv_star);

%Preparation for iteration
M = beta*K_sv*T_sv_star; %Matrix for input calculation
Mt = eye(length(T_sv)) - beta*T_sv*T_sv_star; %Matrix for error calculation
iteration_number = 0;
error_history = [];

%Iteration 
while cont
    %Calculaate new input and errorr
    iteration_number = iteration_number + 1;
    u_new = u + M*e;
    e_new = Mt*e; 
    
    %save the curent error norm
    error_history = [error_history; norm(e_new)];

    %termination criteria 
    if(norm(e - e_new)<1e-6)
        cont = 0;
    end
    
    %print norm of current error difference
    if mod(iteration_number, 100) == 0
        disp(['curr_error_diff=', num2str(norm(e_new - e))]);
    end
    
    %reassign the variables for further iterations 
    e = e_new;
    u_bef = u; %u before for further calculations? 
    u = u_new; 
end

%% Display the results 
disp(['e_inf = ', num2str(error_history(end))])

plot(0:iteration_number-1, error_history);
