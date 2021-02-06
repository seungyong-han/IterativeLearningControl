% %% Satellite Example from Robust Control Lecture
% 
% clc 
% clear all
% 
% 
% s = zpk('s');
% Ts = 1e-2;
% 
% G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
% %[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib('AC1');
% %l = length  (B(1,:));
% %m = length(C(:,1));
% %G = ss(A, B, C, zeros(m, l)); 
% Gd = c2d(G,Ts);
% [A,B,C,D] = ssdata(Gd);
% 
% n = length(A); 
% 
% systemnames='Gd';
% 
% time = 15; %[s]
% per = 3; 
% %% define time horizon and time values
% N = time/Ts;
% time = 0:Ts:time;

% %% Solve with LQR
% Ex3_LQR;
% 
% count = 50;
% [A,B,C,D,count_new] = get_non0D_system(A - B*F, B, C-D*F, D, count); 
% Nnew = N - (count - count_new); 
% r = r_vec(count - count_new+1:end, :)';
% y_sv = y_sv(:,count - count_new+1:end); 
% %u_sv = zeros(N+1,1); %to achieve minimum energy solution 
% u_sv = u_sv(:,count-count_new+1:end)'; 
% 
% x_sv = x_sv(:, count - count_new + 1:end); 
% 
% count = count_new; 

%% solve with controller from the lecture

clc 
clear all

s = zpk('s');
sample_time = 1e-3;

G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
Gd = c2d(G,sample_time);
systemnames='Gd';

% inputvar='[d;n;r;u]';
% outputvar='[Gd+d-r;u;r-n-d-Gd]';
inputvar='[r;u]';
outputvar='[Gd+d-r;r-Gd]';                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ;r-Gd]';
input_to_Gd='[u]';
P=sysic;

K = ss(7.9212*(s+0.1818)*(s^2-0.2244*s+0.8981)/((s^2+3.899*s+4.745)*(s^2+1.039*s+3.395)));
Kd = c2d(K,sample_time);
%Kd = hinfsyn(P, 1, 1);
P_star_K=ss(lft(P,Kd));
isstable(P_star_K)

addpath('models&livescripts')
ampl = 2;
simulation = sim('myModel','SimulationMode','normal', 'StopTime', '5');
data_satellite = simulation.get('data_satellite');
u = simulation.get('u');
time = data_satellite.Time(:,1);
tracking_val = data_satellite.Data(:,1);
output_val = data_satellite.Data(:,2);


u = u.Data(:,1);
x = simulation.xout{1}.Values.Data;

%plot(time,output_val);

[A,B,C,D] = get_stable_ss(Gd, P, Kd); %ssdata(P_star_K); 

N = length(u); 
[A,B,C,D,Nnew] = get_non0D_system(A, B, C, D, N); 
l = length(B(1, :)); 
m = length(C(:, 1)); 
diff = N - Nnew; 

%%
y_sv = output_val;
u_sv = u; 
%%

x0 = zeros(length(A),1); 
r_vec = tracking_val; 
u = zeros(length(u), 1); %begin with zero input => only the input from initial controller 
r = zeros(length(y_sv), 1);%want the error to be zero  
e = r - (y_sv - tracking_val); %zero - act error


Gst_line = [D']; 
G_line = [D]; 
b = 1e-3;
Amult = eye(length(A)); 
d = C*x0; 
for i = 1:Nnew
    Gst_line = [Gst_line, (C*Amult*B)'];
    G_line = [(C*Amult*B), G_line]; 
    Amult = Amult*A;
    d = [d; C*Amult*x0]; 
end
beta = b*1/norm(Gst_line, 1)/norm(G_line, 1);
Gst_line = beta*Gst_line; 


cont = true;
iter_nb = 0; 
N = Nnew;
while cont

iter_nb = iter_nb + 1; 
for i = 1:N+1
    u(i) = u(i) + Gst_line(:, 1:N-i+1)*e(1:N-i+1);
end

% for i = 1:N+1
%     y(i, 1) = G_line(:, 1:l*i)*u(1:i); 
% end
% y = y + d;

y = lsim(Gd, u + u_sv, time); 


e_new = -y + r_vec; 
  




if norm(e - e_new)<1e-6
    cont = false;
end

if mod(iter_nb, 100) == 0
    disp('curr error diff:')
    norm(e - e_new)
    disp('curr error:')
    norm(e_new)
end
e = e_new; 

end
%%
u_alg = u + u_sv; 

%%
Astr{1} = eye(length(A)); 
AB = B; 
for i = 2:count+1
    Astr{i} = A*Astr{i-1}; 
    AB = [AB, Astr{i}*B]; 
end
A_N = Astr{end};

%%
N = 0;
%Nmax = length(u_sv);
Nmax = 1500; 
p = 0; 


R = eye(l*(count+1));
Q = eye(m*(count+1));

glb_iteration_number = 0;
glb_error_history = {}; 
glb_u_history = {}; 
glb_y_history = {}; 

max_iteration_number = 0; 

while N < Nmax
    
    glb_iteration_number = glb_iteration_number + 1; 
    N = N + count;
    
    %Define x0
    if glb_iteration_number == 1
        x0 = x_sv(:, N+1);    
    else
        x0 = A_N*x0 + AB*flip(u_inf(:,end));
        p = p+1; 
    end
    
    %Define u0
    u0 = u_sv(N - count + 1:N+1, :);
    u0 = reshape(u0, [(count+1)*l, 1]);
    
    %Supervector notation
    [G, d] = get_G(A,B, C, D, x0, count);
    %tracking value and options 
    r = r_vec(N - count + 1:N+1, :);
    r = reshape(r, [(count+1)*m, 1]);
    do_plot = 0;
    do_print = 0;
       
    %run algorithm 
%     if N<200
%         beta = 1e-8;
%     end
%     if N>=200 && N<500
%         beta = 1e-5;
%     end
    beta = 1e-4; 
    Satellite_ex_SDAOffline; 
    
    %track global values
    glb_error_history{glb_iteration_number} = error_history; 
    glb_u_history{glb_iteration_number}     = u_inf; 
    glb_y_history{glb_iteration_number}     = y_inf; 
    
    
    if max_iteration_number<iteration_number
        max_iteration_number = iteration_number; 
    end
    
    %print where we are 
    if mod(N, 100)==0
    fprintf('Iteration %d of %d\n', N+1, Nmax);
    end

end 
%%
for i = 1:glb_iteration_number - 1
    tmp = glb_y_history{i};
    y = [y; tmp(:,end)]; 
    tmp = glb_u_history{i};
    u = [u; tmp(:,end)]; 
end
