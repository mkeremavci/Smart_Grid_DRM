clear
clc

%% Initilization
consumption = 7.5:0.001:10.5;
N = 100;
L = mean(consumption) * N;

%% Pricing 1
Ah = N/2;

%% Pricing 2
omega = [-2, 0, 2]; omega = (omega+4) / 4;
Bh = N/2 * omega;

%% Plots
% Pricing 0
price = zeros(3,length(consumption));
price(1,:) = 1.2 * (L + omega(1))/N;
price(2,:) = 1.2 * (L + omega(2))/N;
price(3,:) = 1.2 * (L + omega(3))/N;
figure;
hold on
plot(consumption, price(1,:),'LineWidth',2)
plot(consumption, price(2,:),'LineWidth',2)
plot(consumption, price(3,:),'LineWidth',2)
xline(9)
yline(price(2,1501))
grid on
title('Pricing Model - 0')
ylabel('Unit Price per Hour (')
xlabel('Individual Consumption')
lll = legend('$\omega_h$ = -2','$\omega_h$ = 0','$\omega_h$ = 2');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

% Pricing 1
price = zeros(3,length(consumption));
price(1,:) = 0.8 .* (Ah.*consumption + L + omega(1))./N;
price(2,:) = 0.8 .* (Ah.*consumption + L + omega(2))./N;
price(3,:) = 0.8 .* (Ah.*consumption + L + omega(3))./N;
figure;
hold on
plot(consumption, price(1,:),'LineWidth',2)
plot(consumption, price(2,:),'LineWidth',2)
plot(consumption, price(3,:),'LineWidth',2)
xline(9)
yline(price(2,1501))
grid on
title('Pricing Model - 1')
ylabel('Unit Price per Hour')
xlabel('Individual Consumption')
lll = legend('$\omega_h$ = -2','$\omega_h$ = 0','$\omega_h$ = 2');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off

% Pricing 0
price = zeros(3,length(consumption));
price(1,:) = 0.8 .* (Bh(1).*consumption + L + omega(1))./N;
price(2,:) = 0.8 .* (Bh(2).*consumption + L + omega(2))./N;
price(3,:) = 0.8 .* (Bh(3).*consumption + L + omega(3))./N;
figure;
hold on
plot(consumption, price(1,:),'LineWidth',2)
plot(consumption, price(2,:),'LineWidth',2)
plot(consumption, price(3,:),'LineWidth',2)
xline(9)
yline(price(2,1501))
grid on
title('Pricing Model - 2')
ylabel('Unit Price per Hour')
xlabel('Individual Consumption')
lll = legend('$\omega_h$ = -2','$\omega_h$ = 0','$\omega_h$ = 2');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',14);
hold off



