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
r_vec = [r_vec; r_vec(end)]; 

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
for i = 1:N+1
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
plot(0:N, y_sv);
plot(0:N, r_vec);
legend('output $y(t)$', 'reference signal', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off

%% Define new system 
A = A - B*F; 
C = C -D*F; 
%%
sys = ss(A,B,C,D,sample_time); 
systemnames='sys';
inputvar='[r;u]';
outputvar='[sys-r; r-sys]';
input_to_sys='[u]';
P=sysic;

Kc = hinfsyn(sys, 1, 1);
cloop = lft(P, Kc);
S = -cloop; %sensitivity
T = 1 - S; 
K = Kc*S; 
%% Define the supermatrices 
A_sys = A;  %save the matrices for the case they will be needed later 
B_sys = B;
C_sys = C; 
D_sys = D;


[A,B,C,D] = ssdata(sys); 
[G, d] = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(T); 
T_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(S); 
S_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 

[A,B,C,D] = ssdata(K); 
K_sv = get_G(A,B,C,D,zeros(length(A), 1), N); 




%% FB Design using reverse time simulation
clc 

[A,B,C,D] = ssdata(T); 
n = length(A);
b = .1;

e = S_sv*r; %error
p = zeros(n, N+1); %for inverse time simulation
v = zeros(N+1, 1); %output of inverse time simulation 
v(N+1) = D'*e(N+1);%p(:,N+1) = 0;

error_history = norm(e);
cont = 1;
iteration_number = 0;

while cont 
    
    iteration_number = iteration_number + 1; 
    %First without weights: R, Q = I -- should work in some way...
    
    %Compute the output of T'*e_k
    for t = N:-1:1
        p(:, t) = A'*p(:, t+1) + C'*e(t+1);%e_k
        v(t) = B'*p(:, t) + D'*e(t);%e_k
    end
    
    u_new = u + b*K_sv*v; 
    e_new = r - G*u_new - d; 
    error_history = [error_history, norm(e_new)]; 
    if norm(e_new - e)<10^-6
        cont = 0;
    end
    
    if mod(iteration_number, 1000) == 0
         disp(['curr_error_difference: ', num2str(norm(e_new - e))]); 
    end

    e = e_new;
    u = u_new;

end






%% Display the results 
disp(['e_inf = ', num2str(error_history(end))])

plot(0:iteration_number-1, error_history);

%%


