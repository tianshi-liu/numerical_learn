clear all;
clc;
%%%% set domain and mesh %%%
xmin = -10;
xmax = 10;
CFL = 0.5;
dx = 0.05;
dt = dx * CFL;
x = (xmin + dx : dx : xmax)';
Nx = size(x, 1);
%%%% set initial value %%%%%%
recmin = -2;
recmax = 2;
u0 = (1 + cos((x - (recmin + recmax) / 2) / (recmax - recmin) * 2 * pi)) / 2;
u0((x < recmin) | (x > recmax)) = 0;
%%%% time marching %%%%%
u = u0;
t = 0;
tmax = 15 * dt;
while (t < tmax)
    u = u - CFL * ([u(2 : end, 1); u(1, 1)] - u);
    t = t + dt;
end
f = figure(1);
plot(x, u0, 'r', 'LineWidth', 2);
hold on
plot(x, u, 'b', 'LineWidth', 2);
set(gca, 'FontSize', 15);
ylim([-0.2 1.2]);
xlim([xmin xmax]);
print(f, 'advection_1d_downwind.png', '-dpng');
%ylim([0, 1]);