%%


% Calculate stabilizing solutions for discrete algebraic Riccati equation

%Rw = diag([1, .1, 2]);
%Qw = diag([.1, .2,1, 1, 1]);
Qw = eye(length(A)); 
Rw = eye(l); 

P = dare(A, B, 1 , 1);
Q = dare(A', C', 1 , 1);


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 
%%



x0 = zeros(length(A), 1); 

% tracking with LQR
u_sv = [];%input supervector
y_sv = [];%output supervector

% Calculate matrix for reference tracking
TMP = (eye(length(A)) - A + B*F);
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M =pinv(M);

x_hat = x0;
x = x0;

%r_vec = reshape(r_vec, [m, N+1]);
r_vec = .1*ones(m*(N+1), 1); 
r_vec = [ones(1, N+1); zeros(1, N+1); zeros(1, N+1)]

for i = 1:N+1
    r = r_vec(:, i);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end

%%
%Plot
close all
fig = figure; 
hold on
y_out = reshape(y_sv, [m, N+1]); 
u_out = reshape(u_sv, [l, N+1]);
r_out = reshape(r_vec,[m, N+1]);

plot(0:N, r_vec(1, :));
plot(0:N, y_out(1, :));
legend('Reference', 'LQR', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off

%%

