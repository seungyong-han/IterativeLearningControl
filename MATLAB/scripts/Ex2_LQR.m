P = dare(A, B, 1, 1);
Q = dare(A', C', 1,1);

l = length(B(1,:));
m = length(C(:,1));


F = (eye(l) + B'*P*B)\eye(l);
F = F*B'*P*A;

J = (eye(m) + C*Q*C')\eye(m);
J = J*C*Q*A';
J = J'; 

ts = 1e-3;

if(N<5)
    r_vec = ones(3*(N+1), 1);
else
    r_vec1 = get_pulse(N, floor(N/5));
    r_vec2 = get_pulse(N, floor(N/15));
    r_vec3 = get_pulse(N, floor(N/20));
    
    r_vec = [r_vec1, r_vec2, r_vec3];
    r_vec= reshape(r_vec, [length(r_vec)*3, 1]);
end
r_vec = .2*r_vec;




u_sv = [];
y_sv = [];

TMP = (eye(length(A)) - A + B*F);
TMP = TMP\eye(length(A));
M = (C - D*F)*TMP*B + D;
M = M\eye(length(M));

x_hat = x0;
x = x0;
for i = 1:N
    r = r_vec(i*3 - 2:i*3);
    u = -F*x_hat + M*r; 
    y = C*x + D*u;
    x_hat = (A - J*C)*x_hat + (B - J*D)*u + J*y;
    
    x = A*x + B*u;   
    
    u_sv = [u_sv; u];
    y_sv = [y_sv; y];
end

% hold on
% plot(0:N-1, y_sv);
% plot(0:N-1, r_vec);
% legend('output $y(t)$', 'reference signal', 'interpreter', 'latex');
% xlabel('time $t$', 'interpreter', 'latex');
% ylabel('output and reference');
% xlim([0, N-1]);
% hold off
