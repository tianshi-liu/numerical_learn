clear all
clc;
xmin = -10;
xmax = 10;
dx = 0.1;
x = (xmin:dx:xmax)';
center = 0;
width = 2;
u0 = 1 - abs(x - center) / width;
u0(u0 < 0) = 0;
tstep = 0;
maxtstep = 80;
drawevery = 5;
isfirst = true;
while tstep < maxtstep
    u = zeros(size(u0));
    u(1:end - tstep) = u(1:end - tstep) + u0(tstep + 1:end) / 2;
    u(tstep + 1:end) = u(tstep + 1:end) + u0(1:end - tstep) / 2;
    if (mod(tstep, floor(drawevery)) == 0)
        f=figure(1);
        set(gcf,'position',[200,100,450,200]); %left bottom width height
        set(gcf,'Color',[1,1,1]);
        plot(x, u0, 'r', 'LineWidth', 2);
        hold on
        plot(x, u, 'b', 'LineWidth', 2);
        set(gca, 'FontSize', 15);
        ylim([-0.2 1.2]);
        xlim([xmin xmax]);
        %set(gca,'PlotBoxAspectRatio',[2,1,1]);
        %set(gca,'DataAspectRatio',[2,1,1]);  
        generate_gif(f, 'tent_show.gif', isfirst);
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