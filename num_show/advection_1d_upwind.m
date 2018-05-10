clear all;
clc;
%%%% set domain and mesh %%%
xmin = -10;
xmax = 10;
CFL = 0.01 * 4;
dx = 0.1 / 4;
dt = dx * CFL;
x = (xmin + dx : dx : xmax)';
Nx = size(x, 1);
%%%% set initial value %%%%%%
%%%%%% cosine
recmin = -2;
recmax = 2;
u0 = (1 + cos((x - (recmin + recmax) / 2) / (recmax - recmin) * 2 * pi)) / 2;
u0((x < recmin) | (x > recmax)) = 0;
%%%% Ricker's wavelet
% center = 0;
% sigma = 1;
% u0 = 2 / (sqrt(3 * sigma) * sqrt(sqrt(pi))) * (1 - (x / sigma) .^ 2) .* exp(- (x / sigma) .^ 2 ./ 2);
%%%% time marching %%%%%
u = u0;
t = 0;
tmax = 20;
tstep = 0;
drawevery = 0.5;
isfirst = true;
while (t <= tmax + dt)
    u = u - CFL * (u - [u(end, 1); u(1 : end - 1,1)]);
    t = t + dt;
    if (mod(tstep, floor(drawevery / dt)) == 0)
        f=figure(1);
        set(gcf,'Color',[1,1,1]);
        plot(x, u0, 'r', 'LineWidth', 2);
        hold on
        plot(x, u, 'b', 'LineWidth', 2);
        set(gca, 'FontSize', 15);
        ylim([-0.2 1.2]);
        xlim([xmin xmax]);
        %set(gca,'DataAspectRatio',[98,100,1]);  
        generate_gif(f, 'advection_1d_upwind.gif', isfirst);
        isfirst = false;
    end
    tstep = tstep + 1;
end
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
% figure(1);
% plot(x, u0, 'r');
% hold on
% plot(x, u, 'b');
%ylim([0, 1]);