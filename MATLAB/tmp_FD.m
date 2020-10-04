clc 
clear all
%% Predefined
N = 10;

s = zpk('s');
sample_time = 1e-3;

sys = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));

[A, B, C, D, ts] = ssdata(sys);
l = length(B(1, :));
m = length(C(:, 1));


sys = c2d(sys,sample_time);
% A = [2 1; 4 3];
% B = [1; 2];
% C = [0 1];
% D = 2;
% sys = ss(A,B,C,D, sample_time);


systemnames='sys';

inputvar='[r;u]';
outputvar='[sys-r; r-sys]';
input_to_sys='[u]';
P=sysic;

Kc = hinfsyn(P, m, l);

K = ss(7.9212*(s+0.1818)*(s^2-0.2244*s+0.8981)/((s^2+3.899*s+4.745)*(s^2+1.039*s+3.395)));
Kd = c2d(K,sample_time);
Kc = Kd;
% 
% sys = zpk(sys);
% Kc = zpk(Kc);
cloop = lft(P, Kc);

sens = loopsens(sys,Kc);

T = sens.Ti;
S = 1 - T; 

S = -cloop;
T = eye(length(S)) - S; 
 
K =    Kc/(eye(length(sys)) + sys*Kc); 
T_star = conj(T); 

%% Get supervector form

Nmin = N;
[A, B, C, D] = ssdata(K);
u0 = zeros(length(l));
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
% [A,B,C,D] = ssdata(T_star);
% T_sv_star = get_G(A,B,C,D,zeros(length(A)), N);
T_sv_star = T_sv'; 
% 
[A,B,C,D] = ssdata(Kc);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
Kc_sv = get_G(A,B,C,D,zeros(length(A)), N1);
if N1<Nmin
    Nmin = N1;
end
% 
% [A,B,C,D] = ssdata(S);
% [A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
% S_sv = get_G(A,B,C,D,zeros(length(A)), N1);
% 
[A,B,C,D] = ssdata(sys);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
G = get_G(A,B,C,D,zeros(length(A)), N1);
if N1<Nmin
    Nmin = N1;
end

[A,B,C,D] = ssdata(cloop);
[A,B,C,D,N1] = get_non0D_system(A,B,C,D,N);
cl_sv = get_G(A,B,C,D,zeros(length(A)), N1);
if N1<Nmin
    Nmin = N1;
end

%%
%FD - Algorithm



cont = 1;
r = 3;
beta = .2;
beta = beta/norm(T_sv)/norm(T_sv_star);
T_sv = T_sv(length(T_sv) - Nmin + 1:end,length(T_sv) - Nmin + 1:end);
K_sv = K_sv(length(K_sv) - Nmin + 1:end,length(K_sv) - Nmin +1:end);

%Kc_sv = Kc_sv(length(Kc_sv) - length(T_sv) + 1:end, length(Kc_sv) - length(T_sv) + 1:end);
%K_sv = K_sv(length(K_sv) - length(T_sv) + 1:end, length(K_sv) - length(T_sv) + 1:end);
%S_sv = S_sv(length(S_sv) - length(T_sv) + 1:end, length(S_sv) - length(T_sv) + 1:end);
%T_sv = T_sv(length(T_sv) - length(Kc_sv) + 1:end, length(T_sv) - length(Kc_sv) + 1:end);
T_sv_star = T_sv';

%r = get_pulse(length(T_sv), floor(length(T_sv)/4));
r = 3*ones(length(T_sv),1);
e = (eye(length(T_sv)) - T_sv)*r;

u = u0*zeros(length(T_sv),1);

%M = beta*Kc_sv*S_sv*T_sv_star;
M = beta*K_sv*T_sv_star; 
Mt = eye(length(T_sv)) - beta*T_sv*T_sv_star;
iteration_number = 0;
error_history = [];
while cont
    iteration_number = iteration_number + 1;
    u_new = u + M*e;
    e_new = Mt*e; 
    error_history = [error_history; norm(e_new)];

    if(norm(e - e_new)<1e-6)
        cont = 0;
    end
    
    if mod(iteration_number, 100) == 0
        disp(['curr_error_diff=', num2str(norm(e_new - e))]);
    end
    e = e_new;
    u_bef = u; 
    u = u_new; 

end

disp(['e_inf = ', num2str(error_history(end))])

    plot(0:iteration_number-1, error_history);
    G = G(2:end, 2:end);
%% IA&SDA
% 
% [A,B,C,D] = ssdata(lft(P, Kc));
% [G,d] = get_G(A,B,C,D,zeros(length(A),1),N);
% beta = .1;
% r= 1*ones(N+1, 1);
% 
% [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = IA(G,d,beta,r, u0, 1, 1);
% [u_inf, e_inf, y_inf, impr,iteration_number, error_history] = SDA(G,d,beta,r, u0, 1, 1, 1);
% 
% 
% disp(['e_inf = ', num2str(error_history(end))]); 



