%% Set paths & clear all
addpath('help_fcns');
addpath('Algorithms');
addpath('scripts');

clc 
clear all

%% Define system and controller
% Define system (from Lecture)
s = zpk('s');
sample_time = 0.5*1e-1;
sys = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
N = 431;

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


%% LQR
[A, B, C, D] = ssdata(sys);
[A,B,C,D, N] = get_non0D_system(A,B,C,D, N);
N
RP = 0.001; % Weights...
QP = diag([10, 1, 10, 10]);
RQ = 0.001;
QQ = 10*eye(size(A, 1));
P = idare(A, B, QP, RP);   
Q = idare(A', C', QQ, RQ);

% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));

%Calculate the matrices for separation principle based controller
F = (B' * P * B + RP) \ (B' * P * A);
J = ((C * Q * C' + RQ) \ (C * Q * A'))';
M = 5;
rng default;
r_vec = kron(rand(M, 1), ones(N/M, 1));

eig(A - B * F)
eig(A - J*C)
x0 = zeros(size(A, 1), 1);
%x0 = randn(size(x0));

% tracking with LQR
u_sv = [];%input supervector
y_sv = [];%output supervector

% Calculate matrix for reference tracking
TMP = (eye(length(A)) - (A - B*F));
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M = M\eye(length(M));


x_hat = x0;
x = x0;
for i = 1:N
    r = r_vec(i);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end

% Plot
close all
hold on
plot(0:N-1, y_sv);
plot(0:N-1, r_vec);
legend('output $y(t)$', 'reference signal', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off

%% Inverse Algorithm
u0 = u_sv;
r = r_vec;
R = 1;
Q = 1;
do_plot = 1;
b = .1; 
[G, d] = get_G(A-B*F, B, C - D*F, D, x0, N-1);
close all
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] = IA(G,D, d, b,r, u0, do_plot, 1);
xmax = iteration_number1;
if(xmax>0)
    xlim([0, xmax])
end
legend('$e_k$', 'interpreter', 'latex');

% Plot results LQR vs IA 
close all
hold on
plot(0:N-1, y_sv);
plot(0:N-1, r, 'LineWidth', 1.8);
plot(0:N-1, y_inf1{end}, '--');
plot(0:N-1, y_inf1{1}, '--');
plot(0:N-1, y_inf1{2}, '--');
legend('LQR', 'Reference', 'IA');
xlim([0, N-1])
hold off


%% Get supervector form
% Ob D = 0 ist sollte hier doch immer egal sein?
% Ich glaube auch er definiert T als (I + GK)^{-1}GK
% wobei G und K die supervektor matrizen sind..
% Das sollte aber wohl auch im direkten zusammenhang zum den supervekoren
% matrizen von T stehen. (hoffentlich)

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
    if(norm(e - e_new)<1e-5)
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

%%


