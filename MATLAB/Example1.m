clc
clear all
A = [2 1; 4 3];
B = [1; 2];
C = [0 1];
D = 2;

N = 20;
b = .1;  

% Find controller using separation principle
Ex1_LQR;

% Define start values for the algorithms
u0 = u_sv;
r = r_vec;
R = 1;
Q = 1;
do_plot = 1;
[A,B,C,D, N] = get_non0D_system(A,B,C,D, N)

[G, d] = get_G(A-B*F, B, C - D*F, D, x0, N-1);
%% IA vs LQR example
Ex1_IA;
%% SDA vs SE example
Ex1_SDA;
%% Uncertaiy example for SDA
Ex1_uncDemo_SDA; 
%% Uncertainty proof for SDA 
Ex1_uncProof_SDA;
%% Uncertainty proof for IA beta = .05 vs beta = .4
N = 20;
Ex1_uncProof_IA;