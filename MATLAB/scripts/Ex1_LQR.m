
% l-input m-output system
l = length(B(1,:));
m = length(C(:,1));

% Calculate stabilizing solutions for discrete algebraic Riccati equation

Rw = 2*eye(l); 
Qw = 3*eye(m); 
P = dare(A, B, Qw , Rw);
Q = dare(A', C', Qw , Rw);


%Calculate the matrices for separation principle based controller
F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 


if(N == 45 || N == 30 || N == 81)
    r_vec = .3*ones(floor(N/3),1); 
    r_vec = [r_vec; .5*ones(floor(N/3),1)]; 
    r_vec = [r_vec; ones(floor(N/3),1)]; 
    r_vec = [r_vec; 1];
end
if(N == 100)
    r_vec = .3*ones(floor(N/3),1); 
    r_vec = [r_vec; .5*ones(floor(N/3),1)]; 
    r_vec = [r_vec; ones(floor(N/3),1)]; 
    r_vec = [r_vec; 1; 1];
end

if(N == 50)
    r_vec = .3*ones(17,1); 
    r_vec = [r_vec; .5*ones(17,1)]; 
    r_vec = [r_vec; ones(17,1)]; 
end

if(N == 40)
    r_vec = .5*ones(N/2,1); 
    r_vec = [r_vec; ones(N/2,1)]; 
    %r_vec = [r_vec; ones(13,1)]; 
    r_vec = [r_vec; 1];
end

if(N == 20)
    r_vec = .5*ones(10,1); 
    r_vec = [r_vec; ones(10,1)]; 
    r_vec = [r_vec; 1];
end

if(N == 430)
    r_vec = .5*ones(100, 1); 
    r_vec = [r_vec; ones(120, 1)]; 
    r_vec = [r_vec; 1.2*ones(110, 1)];
    r_vec = [r_vec; .7*ones(100,1)]; 
    r_vec = [r_vec; r_vec(end)]; 
end

if(N == 1000)
    r_vec = .5*ones(150, 1); 
    r_vec = [r_vec; ones(200, 1)]; 
    r_vec = [r_vec; 1.2*ones(350, 1)];
    r_vec = [r_vec; .7*ones(150,1)]; 
    r_vec = [r_vec; -.2*ones(150,1)]; 
    r_vec = [r_vec; r_vec(end)]; 
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
plot(0:N, r_vec);
plot(0:N, y_sv);
legend('Reference', 'LQR', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off
