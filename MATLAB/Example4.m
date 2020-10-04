clc
clear all
ex = 'REA4';
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib(ex);



A = [2e-7 14; .1 4];
B = [100 0; 2 1];
C = [1 10];
D = [1e-4 1e4];

l = length(B(1,:));
m =length(C(:, 1));
%D = zeros(m, l);



N = 40;
b = .1;  

% Find controller using separation principle
Ex4_LQR;

% Define start values for the algorithms
u0 = u_sv;
r = r_vec;
R = 1;
Q = 1;
do_plot = 1;
[A,B,C,D, N] = get_non0D_system(A,B,C,D, N)

[G, d] = get_G(A-B*F, B, C - D*F, D, x0, N-1);
%save_plot(fig, 'Ex1_LQR');




rank(G)
disp(['cond number = ',num2str(cond(G))]);
% IA vs LQR example
Ex4_IA;
%save_plot(fig, 'IA_N40');
%% SDA example
Ex1_SDA;
%save_plot(fig, 'SDA_N40');
%%
% Plot SDA output
close all
fig = figure; 
hold on
plot(0:N-1, y_inf2{end});
plot(0:N-1, r_vec);
legend('output $y(t)$ with SDA', 'reference signal', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off
%save_plot(fig, 'SDA_N40_output');

%% SDA example vs SE example
Ex1_SDAvsES;
save_plot(fig, 'Ex1_SDAvsES');
%% Uncertaiy example for SDA
Ex1_uncDemo_SDA; 
save_plot(fig, 'Ex1_uncSDA');
%% Uncertainty proof for SDA 
Ex1_uncProof_SDA;
%% Uncertainty proof for IA beta = .05 vs beta = .4
N = 20;
Ex1_uncProof_IA;