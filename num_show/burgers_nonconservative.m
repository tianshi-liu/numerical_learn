clear all;
clc;
%%%% set domain and mesh %%%
xmin = -0.1;
xmax = 1.1;
CFL = 0.5;
dx = 0.01;
dt = dx * CFL;
x = (xmin + dx : dx : xmax)';
Nx = size(x, 1);
%%%% set initial value %%%%%%
ul = 1.2;
ur = 0.4;
u0 = zeros(Nx, 1);
u0(x < 0) = ul;
u0(x >= 0) = ur;
%%%% time marching %%%%%
u = u0;
t = 0;
tmax = 1;
while (t < tmax)
    u = u - CFL * u .* (u - [u(1, 1); u(1 : end - 1,1)]);
    t = t + dt;
end
%%%% accurate solution %%%
u_true = zeros(Nx, 1);
shock_position = 0.8; % given by analytic solution
u_true(x < shock_position) = ul;
u_true(x >= shock_position) = ur;
f = figure(1);
plot(x, u0, 'r', 'LineWidth', 2);
hold on
plot(x, u_true, 'g', 'LineWidth', 2);
hold on
plot(x, u, 'bo', 'LineWidth', 2);
set(gca, 'FontSize', 15);
xlim([xmin xmax]);
ylim([0.3 1.3]);
print(f, 'burgers_nonconservative.png', '-dpng');