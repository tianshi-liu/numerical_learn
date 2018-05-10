clear all;
clc;
%%%% set domain and mesh %%%
xmin = -0.1;
xmax = 1.1;
% xmin = -10;
% xmax = 10;
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
%%%% smooth initial value
% recmin = -8;
% recmax = -4;
% u0 = (1 + cos((x - (recmin + recmax) / 2) / (recmax - recmin) * 2 * pi)) / 2;
% u0((x < recmin) | (x > recmax)) = 0;
%%%% time marching %%%%%
u = u0;
t = 0;
tmax = 1;
% tmax = 25;
% tstep = 0;
% drawevery = 0.5;
% isfirst = true;
while (t < tmax)
    u = u - CFL * (u .^ 2 - [u(1, 1); u(1 : end - 1,1)] .^ 2) / 2;
    t = t + dt;
%     if (mod(tstep, floor(drawevery / dt)) == 0)
%         f=figure(1);
%         set(gcf,'position',[200,100,450,200]); %left bottom width height
%         set(gcf,'Color',[1,1,1]);
%         plot(x, u0, 'r', 'LineWidth', 2);
%         hold on
%         plot(x, u, 'b', 'LineWidth', 2);
%         set(gca, 'FontSize', 15);
%         ylim([-0.2 1.2]);
%         xlim([xmin xmax]);
%         %set(gca,'PlotBoxAspectRatio',[2,1,1]);
%         %set(gca,'DataAspectRatio',[2,1,1]);  
%         generate_gif(f, 'burgers_show.gif', isfirst);
%         isfirst = false;
%     end
%     tstep = tstep + 1;
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
print(f, 'burgers_conservative.png', '-dpng');

function generate_gif(f, filename, isfirst)
    drawnow;
    frame=getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if isfirst
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    close(f);
end