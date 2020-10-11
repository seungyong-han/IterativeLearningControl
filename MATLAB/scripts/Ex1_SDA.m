%% Compare SDA and SE
R = 1;
Q = 1;


hold on
% [u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] =
% SDA(G,d, b,r, u0, R, Q,do_plot);%For Q = 1

Q = eye(length(G));
tend = 5;
for t = 1:tend
    Q(t,t) = 10; 
end

[u_inf2, e_inf2, y_inf2, impr2, iteration_number2, error_history2] = SDA(G,d, b,r, u0, R, Q,0);%For Q = 1/2^(2t)I



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

disp(['original error = ', num2str(one_error)]);
disp(['final error = ', num2str(error_history2(end))]);
disp(['diff: ',num2str(-error_history2(end)+one_error)]);


%Plot
close all
fig = figure; 
hold on
plot(0:N-1, y_inf2{end});
plot(0:N-1, r_vec);
legend('SDA', 'reference signal', 'interpreter', 'latex');
xlabel('time $t$', 'interpreter', 'latex');
ylabel('output and reference');
xlim([0, N-1]);
hold off