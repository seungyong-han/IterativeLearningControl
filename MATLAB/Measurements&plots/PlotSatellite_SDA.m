%% Plot satellite example for SDA 
close all
clear all
load('data\Satellite_1e3.mat'); 
fig = figure;
hold on

len = length(time);


%All
% plot(time(1:len),tracking_val(1:len), 'DisplayName', 'tracking value');
% plot(time(1:len),ouput_val(1:len),'DisplayName', 'k = 0');
% 
% for t = 1:5
%     plot(time(1:len), plots{t}(1:len), 'DisplayName', ['k = ',num2str(t*5)]);
% end

%First, tracking value, last
plot(time(1:len),act_trValue(1:len), 'DisplayName', 'tracking value');
plot(time(1:len),ouput_val(1:len),'DisplayName', 'k = 0');
plot(time(1:len), plots{1}(1:len), 'DisplayName', ['SDA k = ',num2str(5)]);
plot(time(1:len), plots{4}(1:len), 'DisplayName', ['SDA k = ',num2str(20)]);



legend();
fprintf('TOTAL TIME: %f', ell_time);
fprintf('TOTAL ITERATION NUMBER: %d', iteration_number);
fprintf('AVERGAGE IMPROVEMENT: %f', sum(impr)/length(impr));
ylim([-.8 1.7])
xlabel('time')
ylabel('output value')
title('Satellite Example')

cd ..;
saveas(fig, 'Figures\Satellite_SDA.jpg');
cd Measurements&plots;

hold off