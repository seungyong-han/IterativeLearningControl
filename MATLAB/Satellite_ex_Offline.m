%% Satellite Example from Robust Control Lecture

clc 
clear all


s = zpk('s');
Ts = 1e-2;

G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
%[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib('AC1');
%l = length  (B(1,:));
%m = length(C(:,1));
%G = ss(A, B, C, zeros(m, l)); 
Gd = c2d(G,Ts);
[A,B,C,D] = ssdata(Gd);

n = length(A); 

systemnames='Gd';

time = 15; %[s]
per = 3; 
%% define time horizon and time values
N = time/Ts;
time = 0:Ts:time;

%% Solve with LQR
Ex3_LQR;

count = 50;
[A,B,C,D,count_new] = get_non0D_system(A - B*F, B, C-D*F, D, count); 

r = r_vec(count - count_new+1:end, :);

%u_sv = zeros(N+1,1); %to achieve minimum energy solution 
u_sv = u_sv(:,count-count_new+1:end)'; 

x_sv = x_sv(:, count - count_new + 1:end); 

count = count_new; 

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
