%% Uncertain Example - Theory 
clc
clear all
%certain system
A = [2 1; 1 3];
B = [1; 2];
C = [0 1];
D = 2;

P = dare(A, B, 1, 1);
Q = dare(A', C', 1,1);

l = length(B(1,:));
m = length(C(:,1));


F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 


N = 20;

if(N<5)
    r_vec = ones(N, 1);
else
    r_vec = get_pulse(N, floor(N/4));
end


x0 = [0 ; 0 ];


u_sv = [];
y_sv = [];

TMP = (eye(length(A)) - A + B*F);
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M = M\eye(length(M));

x_hat = x0;
x = x0;
for i = 1:N
    r = r_vec(i);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end


A = A - B*F;
C = C - D*F; 
N = N-1;




[G,d] = get_G(A,B,C,D,x0, N);
K0 = G^-1;


do_plot = 1;
b = .05;
u0 = u_sv;
r = r_vec;
R = 1;
Q = 1;
do_plot = 0;



[u_inf1, e_inf1, y_inf1, impr1,iteration_number1, error_history1] = RIA_unc(G, K0,d, b,r, u0, R, Q, do_plot)
A = [2 1; 0 3];
A = A - B*F;
[G,d] = get_G(A, B, C, D, x0, N);
[u_inf2, e_inf2, y_inf2, impr,iteration_number2, error_history2] = RIA_unc(G, K0,d, b,r, u0, R, Q, do_plot)

hold on
plot(0:N, y_inf1{end},'LineWidth', 2);
plot(0:N, y_sv, 'LineWidth', 1);
plot(0:N, r_vec,'LineWidth', 1);
legend('Output $y(t)$ without ILC', 'Output $y(t)$ with RIA', 'reference signal', 'interpreter', 'latex');
xlabel('time t');
ylabel('output and reference');
xlim([0, N-1]);
hold off


%% Yalmip
a = ureal('a', 1, 'PlusMinus', [-1, 1]);
A = [2 1; a 3];
A = A - B*F;
C = C - D*F;
ts = 1e-3;

usys = ss(A, B, C, D, ts);
[A, B1, B2, C1, D11, D12, C2, D21, D22, Delta] = get_expand_system(usys);

get_G_large(usys, x0, N);
[G1, G2, G3, G4, d1, d2, Delta] = get_G_large(usys, x0, N);

[m1,n1] = size(G1);
[m,n]   = size(G);

G = G((m-m1+1):end, (n-n1+1):end);

b = .05
K = 2*b*G^-1; 
%K = 1/norm(G)^2*G';

AL = eye(length(G1)) + G1;
BL = G2;
CL = G3*K;
DL = G4;

X = sdpvar(length(AL));
%X = [2 0; 0 .5];

%X = eye(2);
X = kron([-1,0; 0 1], X);

P = sdpvar(length(DL));
Q = sdpvar(length(DL));
P = kron([1 0; 0 -1], P) + kron([0 1; 1 0], Q);


M1 = [eye(length(AL)), zeros(length(AL), length(BL(1,:))); AL, BL]
M2 = [CL, DL; zeros(length(DL),length(CL(1,:))), eye(length(DL))];
cnstr = [M1'*X*M1  + M2'*P*M2<=0]
cnstr = cnstr + [Q' + Q == 0]


%cnstr = cnstr + [-X + A'*X*A <= 0, A;*X*B]

sol = solvesdp(cnstr) 


