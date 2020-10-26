%% Satellite Example from Robust Lecture

clc 
clear all

s = zpk('s');
Ts = 1e-3;

G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
Gd = c2d(G,Ts);
[A,B,C,D] = ssdata(Gd);

systemnames='Gd';

time = 30; %[s]
per = 6; 
%% define time horizon and time values
N = time/Ts + 1; 
time = 0:Ts:time;

%% Solve with LQR
Ex3_LQR;
%% Define start values for the algorithms
u = u_sv;
R = 1;
Q = 1;
Nmax = 3e4;

%%
do_plot = 1;
%Ex3_IA;
Ex3_IA; 




