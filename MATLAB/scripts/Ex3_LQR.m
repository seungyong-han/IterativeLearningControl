
% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));

% Calculate stabilizing solutions for discrete algebraic Riccati equation
P = dare(A, B, 1,1);
Q = dare(A', C', 1,1);


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 


r_vec = 2*ones(length(time),1);%log(time/300+1)'; 
%r_vec = [r_vec zeros(length(time),1)]; 
%r_vec = [r_vec zeros(length(time),1)]; 

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
x_sv = x; 
for i = 1:N+1
    r = r_vec(i, :);
    u = -F*x_hat + M*r'; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv u];
    y_sv = [y_sv y];
    
    x_sv = [x_sv x]; 
end


x_sv = x_sv(:, 1:N+1); 


%%
% Plot

close all
hold on
plot(time, r_vec(:,1));
plot(time, y_sv(1,:));
%grid off
%plot(time, zeros(length(time),1), '--', 'color', 'black', 'LineWidth', 2)
legend('reference', 'output $y(t)$', 'interpreter', 'latex');
xlabel('time $t$ [s]', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, time(end)]);
%ylim([-.1, 1.2])
hold off
