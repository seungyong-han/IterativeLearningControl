%% Compare SDA and SE
R = 1;
Q = 1;


hold on
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] =
SDA(G,d, b,r, u0, R, Q,do_plot);%For Q = 1

Q = 1e-2*eye(m*(N+1));

Q(1,1) = 1*Q(1,1);
Q(2,2) = 100*Q(2,2);
Q(3,3) = 100*Q(3,3);
Q(4,4) = 100*Q(4,4);

%R = eye(l*(N+1)); 

R = diag([ones(4, 1); .1*ones(47,1)]);


b = .1; 
[u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] = SDA(G,d, b,r_vec, u0, R, Q,0);%For Q = 1/2^(2t)I



hold off
xmax =iteration_number2;
xlim([0, xmax])
% legend('$e_k$', 'interpreter', 'latex');
% set(gca, 'YScale', 'log')
% close all 
% fig = figure; 
% p1 = plotresults(error_history2, iteration_number2, 'SDA', 1, -50, .03, 'Red');
% hold on
% plot(0:xmax-1, zeros(xmax,1), '--', 'color', 'black');
% hold off
% legend(p1);

%disp(['original error = ', num2str(one_error)]);
%disp(['final error = ', num2str(error_history2(end))]);
%disp(['diff: ',num2str(-error_history2(end)+one_error)]);


%Plot
close all
plot(0:N, r_vec);
hold on
plot(0:N, y_sv);
plot(0:N, y_inf2{end}, '--');
legend('Reference','LQR', 'SDA');
xlim([0, N-1])
hold off