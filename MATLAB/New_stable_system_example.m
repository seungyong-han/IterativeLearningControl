clear all 
s = zpk('s');
sample_time = 1e-2;
t = 0:sample_time:10; 
N = length(t) - 1; 
G = ss(0.036*(s+25.28)/(s^2*(s^2+0.0396*s+1)));
Gd = c2d(G,sample_time);
[A,B,C,D] = ssdata(Gd); 


[A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny]=COMPleib('AC1');

l = length(B(1,:)); 
m = length(C(:,1)); 
n = length(A); 
D = zeros(m,l);

sys = ss(A, B, C, D); 
sys = c2d(sys, sample_time); 
[A,B,C,D] = ssdata(sys); 

Qw = eye(length(A)); 
Rw = eye(l); 

P = dare(A, B, Qw , Rw);
Q = dare(A', C', Qw , Rw);


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m); 
J = J*C*Q*A';
J = J'; 

TMP = (eye(length(A)) - A + B*F);
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M =pinv(M);
disp('Done')
%% Calculate result with LQR 
unitstep = t<=4;
unitstep1 = t>=8;
unitstep3 = t>=2; 
r_vec = [unitstep + unitstep1; unitstep3 - unitstep1;unitstep];%size [m, N+1];

%r_vec = ones(m,N+1); 
x = zeros(length(A),1); 
x_hat = x; 
u_sv = [];
y_sv = []; 
for i = 1:N+1
    r = r_vec(:, i);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end
%r_vec = reshape(r_vec, [m*(N+1), 1]); 

disp('Done')
%%
A = A - B*F;
C = C - D*F; 

[A,B,C,D, Nnew] = get_non0D_system(A, B, C,D, N); 
sys = ss(A,B,C,D); 
diff = N - Nnew; 
disp('Done')
%%
tracking_val = r_vec; 



%% Old method: full G matrix 
Nmax = N; 
N = 0;
iteration_number = 0;
p = 0; 
count = 20; 

u_inf = [];
e_inf = [];
y_inf = [];
impr = [];
act_trValue = [];
tic
colNb = 1;
plots = {};
plots{1} = [];
plots{2} = [];
plots{3} = [];
plots{4} = [];
plots{5} = [];

l = length  (B(1,:));
m = length(C(:,1));
R = eye(l*(count+1));
Q = eye(m*(count+1));
large_u = 0; 

Amult = A;
AB = B; 
for i = 1:count
    AB = [AB Amult*B];
    Amult = A*Amult; 
end
A_N = Amult; 

while N < Nmax
    
   
    N = N + count;
    
    if N == count
        x0 = x(:, 1);    
    else
        x0 = A_N*x0 + AB*flip(u_inf(p*count*(m) + 1 + p*m:count*(p+1)*m + (p+1)*m));
        p = p+1; 
    end
   
    u0 = zeros((count+1)*m, 1); %to achieve minimum energy solution 
    [G, d] = get_G(A,B, C, D, x0, count);

    beta = 1e-4;
    r = tracking_val(:, N - count + 1:N+1);
    r = reshape(r, [length(r)*m,1]); 
    do_plot = 0;
    do_print = 0;
    
    [u_inf1, e_inf1, y_inf1, impr1, iteration_number1] = SDA(G,d, beta,r, u0,R,1, do_plot);%RIA(G, d, beta, r, u0, 0); %
   
    u_inf = [u_inf; u_inf1{length(u_inf1)}];
    y_inf = [y_inf; y_inf1{length(y_inf1)}];
    e_inf = [e_inf; e_inf1{length(e_inf1)}];
    act_trValue = [act_trValue; r];
    impr = [impr; impr1];
    iteration_number = iteration_number + iteration_number1;
    if mod(N, 100)==0
    fprintf('Iteration %d of %d\n', N+1, Nmax);
    end
    q = length(u_inf1);

end 
ell_time = toc;
%% Sparse method (No full G matrix, ony G_vec)
disp("STOP HERE")
keyboard;
x0 = zeros(length(A),1); 

u = 100*zeros((N+1)*l, 1); 
r = r_vec;

Gst_line = [D']; 
G_line = [D]; 
b = .1;
Amult = eye(length(A)); 
d = C*x0; 
Mr = M*r(1:m,:); 

for i = 1:Nnew
    Gst_line = [Gst_line, (C*Amult*B)'];
    G_line = [(C*Amult*B), G_line]; 
    Amult = Amult*A;
    d = [d; C*Amult*x0]; 
    Mr = [Mr; M*r(i*m+1:(i+1)*m)]; 
end
beta = b/norm(Gst_line, inf)/norm(G_line, inf);
Gst_line = beta*Gst_line; 

e = r; 

u = u(diff*l+1:end); 
r = r(diff*m+1:end); 
e = e(diff*m+1:end); 

cont = true;
iter_nb = 0; 
N = Nnew;


u_sv = u; 

cell_nb = 1; 
while cont

iter_nb = iter_nb + 1; 
%u_sv = u_sv + beta*K0*e;
%y_sv = G*u_sv; 

for i = 0:N
     u(i*l+1:(i+1)*l) = u(i*l+1:(i+1)*l) + Gst_line(:, 1:(N-i+1)*m)*e(i*m+1:end);%SDA
     y(i*m+1:(i+1)*m,1) = G_line(:,(N-i)*m+1:end)*(u(1:(i+1)*l));% + Mr(1:(i+1)*m)); 
end

if any(isnan(y))
    keyboard; 
end


e_new = r - y; 

if norm(e_new)>norm(e)
        keyboard; 
end


if norm(e - e_new)<1e-6
    cont = false;
end

% if mod(iter_nb, 100) == 0
%     disp('curr error diff:')
%     norm(e - e_new)
%     disp('curr error:')
%     norm(e_new)
%     if norm(e_new)>norm(e)
%         keyboard; 
%     end
% end
try
if mod(iter_nb, 1000)==0
    load save_outputs; 
    save_outputs{cell_nb}{1} = iter_nb; 
    save_outputs{cell_nb}{2} = y; 
    save('save_outputs.mat', 'save_outputs');
    clear save_outputs; 
end
catch
    disp('Couldnt save the output values'); 
end
e = e_new; 
end



