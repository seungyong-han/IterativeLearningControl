
% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));

% Calculate stabilizing solutions for discrete algebraic Riccati equation
P = dare(A, B, eye(m) , eye(l));
Q = dare(A', C', eye(m) , eye(l));


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 


% tracking signal as pulse 
if(N<5)
    r_vec = ones(N+1, 1);
else
    r_vec = get_pulse(N+1, 5);
end


x0 = [0 ; 0 ];

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
for i = 1:N+1
    r = r_vec(i);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end

%Plot
close all
fig = figure; 
hold on
plot(0:N, y_sv);
plot(0:N, r_vec);
legend('LQR', 'reference signal', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off
