addpath('Matlab_COMPlib_r1_1\COMPlib_r1_1');
[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib('AC1');
% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));
Ts = 1e-3;

sys = ss(A,B,C,zeros(m,l)); 
sys_d = c2d(sys, Ts);

N = 3e4;
b = .1;  
l
m
C*B

%%
Ex4_LQR;
%%
Ex4_IA; 