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
    wd(:, i) = polyval(pd, xi) / polyval(p, xi(i, 1)); %samples*functions
end
%%%% define nodes %%%%%%
x_elements = bsxfun(@plus, (xi + 1) / 2 * (x(1, 2:end) - x(1, 1:end - 1)), x(1, 1:end - 1)); %% order+1 * elements
%x_nodes = [reshape(x_elements(1:end - 1, :), [order * (Nx - 1), 1]); x(1, end)];
%N_nodes = size(x_nodes, 1);
%%%% set media parameter
rho = ones(size(x_elements));
mu = ones(size(x_elements));
%%%% compute element mass and stiffness matrix
Me = bsxfun(@times, rho, w);
Ke = bsxfun(@times, wd, w);
Ke = bsxfun(@times, reshape(Ke, [order + 1, order + 1, 1]), reshape(mu, [order + 1, 1, Nx - 1]));
Me = bsxfun(@times, Me, x(1, 2:end) - x(1, 1:end - 1)) / 2;
% Ke = bsxfun(@times, reshape(wd', [order + 1, 1, order + 1]), reshape(wd', [1, order + 1, order + 1]));
% Ke = bsxfun(@times, reshape(Ke, [order + 1, order + 1, order + 1, 1]), reshape(mu, [1, 1, order + 1, Nx - 1]));
% Ke = bsxfun(@times, Ke, reshape(w, [1, 1, order + 1, 1]));
% Ke = squeeze(sum(Ke, 3));
% Me = bsxfun(@times, Me, x(1, 2:end) - x(1, 1:end - 1)) / 2;
% Ke = bsxfun(@rdivide, Ke, reshape(x(1, 2:end) - x(1, 1:end - 1), [1, 1, Nx - 1])) * 2;
% Me(1, 2:end) = Me(1, 2:end) + Me(end, 1:end - 1);
% Me = [reshape(Me(1:order, :), [N_nodes - 1, 1]); Me(end, end)];
%%%% set node-to-element mapping
% node2element = bsxfun(@plus, (1:1:order)', (0:Nx - 2) * order);
% node2element = [node2element; [node2element(1, 2:end), order * (Nx - 1) + 1]];
%%%% set boundary condition
boundary_type = 'Outward'; % Dirichlet, Neumann, Outward
% if strcmp(boundary_type, 'Dirichlet')
%     Me = Me(2:end - 1, 1);
% elseif strcmp(boundary_type, 'Periodic')
%     Me(1, 1) = Me(1, 1) + Me(end, 1);
%     Me = Me(1:end - 1, 1);
% end
% suppose Neumann boundary
%%%% set initial value
%%%%%%%%zero
displacement0 = zeros(size(x_elements));
%%%%%%%%cosine
% recmin = -2;
% recmax = 2;
% displacement0 = (1 + cos((x_elements - (recmin + recmax) / 2) / (recmax - recmin) * 2 * pi)) / 2;
% displacement0((x_elements < recmin) | (x_elements > recmax)) = 0;
%%%%%%%%%%%%%%
%%% Ricker's wavelet
% center = 0;
% sigma = 1;
% displacement0 = 2 / (sqrt(3 * sigma) * sqrt(sqrt(pi))) * (1 - (x_nodes / sigma) .^ 2) .* exp(- (x_nodes / sigma) .^ 2 ./ 2);
%%%%%%%%%%%%%%%
displacement = displacement0;
velocity = zeros(size(x_elements));
strain = squeeze(sum(bsxfun(@times, Ke, reshape(displacement, [1, order + 1, Nx - 1])), 2))./ Me;
strain = strain .* rho;
%%%%%%%%%%%%%define source type%%%%%%%%%
source_type = 'dislocation';
dislocation_position = floor(Nx / 2);  %% dislocation between element position & position+1
sourcetime_func = @(t) pulse(t, 5);
%%%%time marching%%%%%%%%%%
tmax = 15;
t = 0;
tstep = 0;
drawevery = 0.2;
isfirst = true;
while (t <= tmax + dt)
    %displacement_half = displacement + dt * velocity;
    %acceleration = get_acceleration(displacement, order, Nx, Ke, Me, node2element, N_nodes, boundary_type);
    %acceleration_half = get_acceleration(displacement_half, order, Nx, Ke, Me, node2element, N_nodes, boundary_type);
    if strcmp(source_type, 'dislocation')
        mu_source_left = mu(order + 1, dislocation_position);
        mu_source_right = mu(1, dislocation_position + 1);
        wave_speed_source_left = sqrt(mu_source_left / rho(order + 1, dislocation_position));
        wave_speed_source_right = sqrt(mu_source_right / rho(1, dislocation_position + 1));
        dislocation_flux_1 = - mu_source_left * mu_source_right / ...
            (mu_source_left * wave_speed_source_right + mu_source_right * wave_speed_source_left);
        dislocation_flux_2_left = - wave_speed_source_left * mu_source_right / ...
            (mu_source_left * wave_speed_source_right + mu_source_right * wave_speed_source_left);
        dislocation_flux_2_right = wave_speed_source_right * mu_source_left / ...
            (mu_source_left * wave_speed_source_right + mu_source_right * wave_speed_source_left);
    end    
    [F1, F2] = get_flux(velocity, strain, order, Nx, rho, mu, boundary_type);
    edge1 = [F1(1, 1:end - 1); F1(1, 2:end)];
    edge2 = [F2(1, 1:end - 1); F2(1, 2:end)];
    if strcmp(source_type, 'dislocation')
        vd = sourcetime_func(t);
        edge1(1, dislocation_position + 1) = edge1(1, dislocation_position + 1) + dislocation_flux_1 * vd;            
        edge1(2, dislocation_position) = edge1(2, dislocation_position) + dislocation_flux_1 * vd;
        edge2(1, dislocation_position + 1) = edge2(1, dislocation_position + 1) + dislocation_flux_2_right * vd;
        edge2(2, dislocation_position) = edge2(2, dislocation_position) + dislocation_flux_2_left * vd;        
    end
    [acceleration, strain_rate] = get_rate(velocity, strain, order, Nx, Ke, Me, edge1, edge2, mu);
    velocity_pred = velocity + dt * acceleration;
    strain_pred = strain + dt * strain_rate;
    [F1, F2] = get_flux(velocity_pred, strain_pred, order, Nx, rho, mu, boundary_type);
    edge1 = [F1(1, 1:end - 1); F1(1, 2:end)];
    edge2 = [F2(1, 1:end - 1); F2(1, 2:end)];
    if strcmp(source_type, 'dislocation')
        vd = sourcetime_func(t);
        edge1(1, dislocation_position + 1) = edge1(1, dislocation_position + 1) + dislocation_flux_1 * vd;            
        edge1(2, dislocation_position) = edge1(2, dislocation_position) + dislocation_flux_1 * vd;
        edge2(1, dislocation_position + 1) = edge2(1, dislocation_position + 1) + dislocation_flux_2_right * vd;
        edge2(2, dislocation_position) = edge2(2, dislocation_position) + dislocation_flux_2_left * vd;        
    end
    [acceleration_pred, strain_rate_pred] = get_rate(velocity_pred, strain_pred, order, Nx, Ke, Me, edge1, edge2, mu);
    displacement = displacement + dt * (velocity + velocity_pred) / 2;
    velocity = velocity + dt * (acceleration + acceleration_pred) / 2;
    strain = strain + dt * (strain_rate + strain_rate_pred) / 2;
%%%%create animation %%%%%%%%%%   
    if (mod(tstep, floor(drawevery / dt)) == 0)
        f=figure(1);
        set(gcf,'Color',[1,1,1]);
        plot(reshape(x_elements, [], 1), reshape(displacement0, [], 1), 'r', 'LineWidth', 2);
        hold on
        plot(reshape(x_elements, [], 1), reshape(displacement, [], 1), 'b', 'LineWidth', 2);
        set(gca, 'FontSize', 15);
        ylim([-1.2 1.2]);
        xlim([xmin xmax]);
        annotation('textbox',[0.66 0.8 0.3 0.06],...
            'String',{['t =' num2str(t)]},...
            'HorizontalAlignment','center',...
            'LineStyle','None',...
            'FontSize',15);
        %set(gca,'DataAspectRatio',[98,100,1]);  
        generate_gif(f, 'dg_1d_dislocation.gif', isfirst);
        isfirst = false;
    end
    t = t + dt;
    tstep = tstep + 1;    
end
% f=figure(1);
% plot(reshape(x_elements, [], 1), reshape(displacement0, [], 1), 'r', 'LineWidth', 2);
% hold on
% plot(reshape(x_elements, [], 1), reshape(displacement, [], 1), 'b', 'LineWidth', 2);

function [acceleration, strain_rate] = get_rate(velocity, strain, order, Nx, Ke, Me, edge1, edge2, mu)
acceleration = -squeeze(sum(bsxfun(@times, Ke, reshape(strain, [order + 1, 1, Nx - 1])), 1));
acceleration(1, :) = acceleration(1, :) - edge1(1, :);
acceleration(end, :) = acceleration(end, :) + edge1(2, :);
acceleration = acceleration ./ Me;
strain_rate = squeeze(sum(bsxfun(@times, Ke, reshape(velocity, [1, order + 1, Nx - 1])), 2));
strain_rate(1, :) = strain_rate(1, :) - mu(1, :) .* (edge2(1, :) - velocity(1, :));
strain_rate(end, :) = strain_rate(end, :) + mu(end, :) .* (edge2(2, :) - velocity(end, :));
strain_rate = strain_rate ./ Me;
end

function [F1, F2] = get_flux(velocity, strain, order, Nx, rho, mu, boundary_type)
wave_speed = sqrt(mu ./ rho);
if strcmp(boundary_type, 'Neumann')
    velocity_neg = [velocity(1, 1), velocity(order + 1, 1:Nx - 1)];
    velocity_pos = [velocity(1, 1:Nx - 1), velocity(order + 1, Nx - 1)];
    strain_neg = [-strain(1, 1), strain(order + 1, 1:Nx - 1)];
    strain_pos = [strain(1, 1:Nx - 1), -strain(order + 1, Nx - 1)];
elseif strcmp(boundary_type, 'Dirichlet')
    velocity_neg = [-velocity(1, 1), velocity(order + 1, 1:Nx - 1)];
    velocity_pos = [velocity(1, 1:Nx - 1), -velocity(order + 1, Nx - 1)]; 
    strain_neg = [strain(1, 1), strain(order + 1, 1:Nx - 1)];
    strain_pos = [strain(1, 1:Nx - 1), strain(order + 1, Nx - 1)];
elseif strcmp(boundary_type, 'Outward')
    velocity_neg = [0, velocity(order + 1, 1:Nx - 1)];
    velocity_pos = [velocity(1, 1:Nx - 1), 0]; 
    strain_neg = [0, strain(order + 1, 1:Nx - 1)];
    strain_pos = [strain(1, 1:Nx - 1), 0];   
end
wave_speed_neg = [wave_speed(1, 1), wave_speed(order + 1, 1:Nx - 1)];
wave_speed_pos = [wave_speed(1, 1:Nx - 1), wave_speed(order + 1, Nx - 1)];
mu_neg = [mu(1, 1), mu(order + 1, 1:Nx - 1)];
mu_pos = [mu(1, 1:Nx - 1), mu(order + 1, Nx - 1)];
F1 = mu_neg .* mu_pos .* ((velocity_pos - velocity_neg) + (wave_speed_neg .* strain_neg + wave_speed_pos .* strain_pos))...
    ./ (mu_neg .* wave_speed_pos + mu_pos .* wave_speed_neg);
F2 = ((mu_pos .* wave_speed_neg .* velocity_pos + mu_neg .* wave_speed_pos .* velocity_neg) - ...
    wave_speed_pos .* wave_speed_neg .* (mu_neg .* strain_neg - mu_pos .* strain_pos))...
    ./ (mu_neg .* wave_speed_pos + mu_pos .* wave_speed_neg);
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

function dislocation = pulse(t, tmax)
%%%%cosine
%     if t < tmax
%         dislocation = (1 + cos((t - tmax / 2) / tmax * 2 * pi)) / 2;
%     else
%         dislocation = 0;
%     end
%%%%%sine
    if t < tmax
        dislocation = sin((t - tmax) / tmax * 2 * pi);
    else
        dislocation = 0;
    end
end