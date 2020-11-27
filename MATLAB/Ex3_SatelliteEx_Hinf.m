%Satellite Example
clc 
clear all

%Simulate solution iwth Hinf 
s = zpk('s');
sample_time = 1e-3;

G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
Gd = c2d(G,sample_time);
systemnames='Gd';

% inputvar='[d;n;r;u]';

% outputvar='[Gd+d-r;u;r-n-d-Gd]';
inputvar='[r;u]';
outputvar='[Gd - r;r-Gd]';
input_to_Gd='[u]';
P=sysic;

K = ss(7.9212*(s+0.1818)*(s^2-0.2244*s+0.8981)/((s^2+3.899*s+4.745)*(s^2+1.039*s+3.395)));
Kd = c2d(K,sample_time);
%Kd = hinfsyn(P, 1, 1);
P_star_K=ss(lft(P,Kd));
isstable(P_star_K)


addpath('models&livescripts')
ampl = 10;

time = '1'; 
simulation = sim('myModel','SimulationMode','normal', 'StopTime', time);
data_satellite = simulation.get('data_satellite');
u = simulation.get('u');
time = data_satellite.Time(:,1);
tracking_val = data_satellite.Data(:,1);
ouput_val = data_satellite.Data(:,2);


u = u.Data(:,1);
x = simulation.xout{1}.Values.Data;

plot(time,ouput_val);

%%
cont = 50; 
u_sv = u; 
r_vec = tracking_val; 
[A, B, C, D,] = get_stable_ss(Gd, P, Kd);

Nmax = length(u) - 1; 
x0 = x(1:length(A))'; 
Ex3_IA; 



