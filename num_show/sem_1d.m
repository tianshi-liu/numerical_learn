clear all
clc;
%%%% mesh discretization %%%%%
xmin = -10;
xmax = 10;
dx = 0.1;
x = xmin : dx : xmax;
Nx = size(x, 2);
vmax = 1;
CFL = 0.01;
dt = CFL * dx / vmax;
%%%% Gauss points and Lagrange derivatives%%%%%
order = 4;
[xi, w, ~] = lglnodes(order);
xi = xi(end: -1: 1, 1);
w = w(end: -1: 1, 1);
wd = zeros(order + 1, order + 1);
for i = 1:1:order + 1
    if (i == 1)
        xi_drop = xi(2:end, 1);
    elseif (i == order + 1)
        xi_drop = xi(1:end - 1, 1);
    else
        xi_drop = [xi(1:i - 1, 1); xi(i + 1:end, 1)];
    end
    p = poly(xi_drop);
    pd = polyder(p);
    wd(:, i) = polyval(pd, xi) / polyval(p, xi(i, 1));
end
%%%% define nodes %%%%%%
x_elements = bsxfun(@plus, (xi + 1) / 2 * (x(1, 2:end) - x(1, 1:end - 1)), x(1, 1:end - 1)); %% order+1 * elements
x_nodes = [reshape(x_elements(1:end - 1, :), [order * (Nx - 1), 1]); x(1, end)];
N_nodes = size(x_nodes, 1);
%%%% set media parameter
rho = ones(size(x_elements));
mu = ones(size(x_elements));
%%%% compute element mass and stiffness matrix
Me = bsxfun(@times, rho, w);
Ke = bsxfun(@times, reshape(wd', [order + 1, 1, order + 1]), reshape(wd', [1, order + 1, order + 1]));
Ke = bsxfun(@times, reshape(Ke, [order + 1, order + 1, order + 1, 1]), reshape(mu, [1, 1, order + 1, Nx - 1]));
Ke = bsxfun(@times, Ke, reshape(w, [1, 1, order + 1, 1]));
Ke = squeeze(sum(Ke, 3));
Me = bsxfun(@times, Me, x(1, 2:end) - x(1, 1:end - 1)) / 2;
Ke = bsxfun(@rdivide, Ke, reshape(x(1, 2:end) - x(1, 1:end - 1), [1, 1, Nx - 1])) * 2;
Me(1, 2:end) = Me(1, 2:end) + Me(end, 1:end - 1);
Me = [reshape(Me(1:order, :), [N_nodes - 1, 1]); Me(end, end)];
%%%% set node-to-element mapping
node2element = bsxfun(@plus, (1:1:order)', (0:Nx - 2) * order);
node2element = [node2element; [node2element(1, 2:end), order * (Nx - 1) + 1]];
%%%% set boundary condition
boundary_type = 'Dirichlet'; % Dirichlet, Neumann, Periodic
if strcmp(boundary_type, 'Dirichlet')
    Me = Me(2:end - 1, 1);
elseif strcmp(boundary_type, 'Periodic')
    Me(1, 1) = Me(1, 1) + Me(end, 1);
    Me = Me(1:end - 1, 1);
end
%%%% set initial value
%%%%%%%%cosine
recmin = -2;
recmax = 2;
displacement0 = (1 + cos((x_nodes - (recmin + recmax) / 2) / (recmax - recmin) * 2 * pi)) / 2;
displacement0((x_nodes < recmin) | (x_nodes > recmax)) = 0;
%%%%%%%%%%%%%%
%%% Ricker's wavelet
% center = 0;
% sigma = 1;
% displacement0 = 2 / (sqrt(3 * sigma) * sqrt(sqrt(pi))) * (1 - (x_nodes / sigma) .^ 2) .* exp(- (x_nodes / sigma) .^ 2 ./ 2);
%%%%%%%%%%%%%%%
displacement = displacement0;
velocity = zeros(N_nodes, 1);
%%%%time marching%%%%%%%%%%
tmax = 40;
t = 0;
tstep = 0;
drawevery = 0.5;
isfirst = true;
while (t <= tmax + dt)
    displacement_pred = displacement + dt * velocity;
    acceleration = get_acceleration(displacement, order, Nx, Ke, Me, node2element, N_nodes, boundary_type);
    acceleration_pred = get_acceleration(displacement_pred, order, Nx, Ke, Me, node2element, N_nodes, boundary_type);
    velocity_pred = velocity + dt * acceleration;
    velocity = velocity + dt * (acceleration + acceleration_pred) / 2;
    displacement = displacement + dt * (velocity + velocity_pred) / 2;
    %%%%create animation%%%%%
%     if (mod(tstep, floor(drawevery / dt)) == 0)
%         f=figure(1);
%         set(gcf,'Color',[1,1,1]);
%         plot(x_nodes, displacement0, 'r', 'LineWidth', 2);
%         hold on
%         plot(x_nodes, displacement, 'b', 'LineWidth', 2);
%         set(gca, 'FontSize', 15);
%         ylim([-0.2 1.2]);
%         xlim([xmin xmax]);
%         annotation('textbox',[0.66 0.8 0.3 0.06],...
%             'String',{['t =' num2str(t)]},...
%             'HorizontalAlignment','center',...
%             'LineStyle','None',...
%             'FontSize',15);
%         %set(gca,'DataAspectRatio',[98,100,1]);  
%         generate_gif(f, 'sem_1d.gif', isfirst);
%         isfirst = false;
%     end
    t = t + dt;
    tstep = tstep + 1;    
end
figure(1);
plot(x_nodes, displacement, 'b', 'LineWidth', 2);
hold on
plot(x_nodes, displacement0, 'r', 'LineWidth', 2);

function acceleration = get_acceleration(displacement, order, Nx, Ke, Me, node2element, N_nodes, boundary_type)
    acceleration = zeros (N_nodes, 1);
    f_internal = bsxfun(@times, Ke, reshape(displacement(node2element), [order + 1, 1, Nx - 1]));
    f_internal = squeeze(sum(f_internal, 1));
    f_internal(1, 2:end) = f_internal(1, 2:end) + f_internal(end, 1:end - 1);
    f_internal = [reshape(f_internal(1:order, :), [N_nodes - 1, 1]); f_internal(end, end)];
    if strcmp(boundary_type, 'Dirichlet')
        acceleration(2:end - 1, 1) = - f_internal(2:end - 1, 1) ./ Me;        
    elseif strcmp(boundary_type, 'Periodic')
        f_internal(1, 1) = f_internal(1, 1) + f_internal(end, 1);
        acceleration(1:end - 1, 1) = - f_internal(1:end - 1, 1) ./ Me;
    else
        acceleration = - f_internal ./ Me;
    end
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
    