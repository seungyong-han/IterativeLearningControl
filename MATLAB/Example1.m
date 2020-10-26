clc
clear all
A = [2 1; 4 3];
B = [1; 2];
C = [0 1];
D = 2;

N = 50;
b = .1;  

% Find controller using separation principle
Ex1_LQR;

% Define start values for the algorithms
u0 = u_sv;
r = r_vec;

do_plot = 1;
[A,B,C,D, N] = get_non0D_system(A,B,C,D, N)

[G, d] = get_redG(A-B*F, B, C - D*F, D, x0, N);
%save_plot(fig, 'Ex1_LQR');
disp(['cond number = ',num2str(cond(G))]);
%%
% IA vs LQR example
Ex1_IA;
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

Q = 1;
R = 1;
b = .1; 
Ex1_SDAvsES;
%% Uncertaiy example for SDA
Ex1_uncDemo_SDA; 
%% Uncertainty proof for SDA 
Ex1_uncProof_SDA;
%% Uncertainty proof for IA beta = .05 vs beta = .4
N = 20;
Ex1_uncProof_IA;
%%
