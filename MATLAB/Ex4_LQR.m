%%


% Calculate stabilizing solutions for discrete algebraic Riccati equation

Rw = diag([1, .1, 2]);
Qw = diag([.1, .2,1, 1, 1]);
P = dare(A, B, Qw , Rw);
Q = dare(A', C', Qw , Rw);


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 
%%

r_vec = ones(m*(N+1), 1); 


x0 = zeros(length(A), 1); 

% tracking with LQR
u_sv = [];%input supervector
y_sv = [];%output supervector

% Calculate matrix for reference tracking
TMP = (eye(length(A)) - A + B*F);
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M = M\eye(length(M));

x_hat = x0;
x = x0;

r_vec = reshape(r_vec, [m, N+1]);
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
plot(0:N, r_vec(3, :));
plot(0:N, y_out(3, :));
legend('Reference', 'LQR', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off

%%

