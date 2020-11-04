addpath('Matlab_COMPlib_r1_1\COMPlib_r1_1');
%%
clc 
clear all
%  s = zpk('s');
%  sys = (s + 1)/(s^3 + 3*s^2 + 3*s + 2);
%  [A,B,C,D] = ssdata(sys);

[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib('AC1');%REA1:1e4 iterations

if(length(A)>15)
    keyboard;
end
% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));
Ts = 1e-3;

sys = ss(A,B,C,zeros(m,l)); 
sys_d = c2d(sys, Ts);
[A,B,C,D] = ssdata(sys_d); 
N = 3e4;
Nmax = N; 
b = .1;  
l
m
C*B
f = @(p)norm((C*(A^p - A^N)*(eye(length(A)) - A)^-1*B));

Ex4_LQR;

abs(eig(A - B*F))
max(abs(A - B*F))

A = A - B*F;
C = C - D*F;
A1 = (eye(length(A)) - A)^-1; 

tol = 10^-6/2/norm(C)/norm(B)/norm(A1)

for i = 1:5
    p = 10^i; 
    if(tol > norm(A^p))
        i
        break
    end
end
f = @(p)(2*norm(C)*norm(B)*norm((eye(length(A)) - A)^-1)*norm(A^p));%

%%
Ex4_IA; 
