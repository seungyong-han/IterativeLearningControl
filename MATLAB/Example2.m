%% Example 2: 
clc
clear all
%% System Setup

%COMPleib('REA1') example: 
A = [1.3800   -0.2077    6.7150   -5.6760;
   -0.5814   -4.2900         0    0.6750;
    1.0670    4.2730   -6.6540    5.8930;
    0.0480    4.2730    1.3430   -2.1040];
 
 B =[0             0;
    5.6790         0;
    1.1360   -3.1460;
    1.1360         0];
 
 C = [1     0     1    -1
     0     1     0     0
     0     0     1    -1];
  
ts = 1e-3;
l = length(B(1,:));
m = length(C(:,1));
D = zeros(m,l);
x0 = zeros(4,1);%[1, .1, .2, -.3]';
sys = ss(A, B, C, D);
dsys = c2d(sys, ts);
[A,B,C,D] = ssdata(dsys);

N = 50;
b = .1; 


%% Find controller using separation principle
Ex2_LQR;
%% Define start values for the algorithms

rel_deg = get_relDeg(A,B,C,D);
[A,B,C,D, N] = get_non0D_system(A,B,C,D, N)
[G, d] = get_G(A-B*F, B, C - D*F, D, x0, N-1);

u0 = u_sv(rel_deg*l+1:end);
r = r_vec(rel_deg*m+1:end);
y_sv = y_sv(rel_deg*m+1:end);
R = 1;
Q = 1;
do_plot = 1;
%% LIA
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] = IA(G,d, b,r, u0, do_plot,1);
%% Plot results

%reshape solutions (3-dim output)
y_sv = reshape(y_sv, [length(y_sv)/m, m]);
r = reshape(r, [length(r)/m, m]);
y = reshape(y_inf1{end}, [length(y_inf1{end})/m, m]);
%%
close all
hold on
plot(0:N-1, y_sv(:,3));
plot(0:N-1, y(:,3));
plot(0:N-1, r(:,3), 'LineWidth', 1.8);
legend('LQR', 'IA', 'Reference');
hold off

%% reference as constant signal
%Load e0 for feasible r
e0 = load('Workspace/LIA_feas.mat');
e0 = e0.v0;

r = e0 + G*u0 + d;
[u_inf1, e_inf1, y_inf1, impr1, iteration_number1, error_history1] = IA(G,d, b,r, u0, do_plot,1);
r = reshape(r, [length(r)/m, m]);
y = y_inf1{end};
y = reshape(y, [length(y)/m, m]);
%%

close all
hold on
plot(0:N-1, y(:,1));
plot(0:N-1, r(:,1), 'LineWidth', 1.8);
legend('IA', 'Reference');
hold off