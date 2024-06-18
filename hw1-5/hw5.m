Nx = 101;
Ny = 101;
L = 1;
x = linspace(0,L,Nx);
y = linspace(0,L,Ny);
xx = linspace(0,L,21);
yy = linspace(0,L,21);
dx = x(2)-x(1);
dy = y(2)-y(1);
h = dx;
nu = 0.001;
dt = 0.002;

grd=meshgrid(x,y);
ggrd=meshgrid(xx,yy);

u = zeros(size(grd));
v = zeros(size(grd));
uu = zeros(size(ggrd));
vv = zeros(size(ggrd));

v(end,2:end-1)=1;
%for i=2:Ny-1
%    v(end,i)=y(i)*y(i)*(1-y(i))*(1-y(i))*16;
%end

omega = zeros(size(grd));
omega = vel_to_omega(omega,u,v,dx,dy);

norm(omega,inf)
norm(u,inf)

[psi,~] = get_stream(omega,dx,dy);

flag_end = false;
out_iter = 0;

figure;
set(gcf, 'Position', [100, 100, 1500, 400]);

while flag_end == false
    omega = get_bound_omega(psi,u,v,omega,dx,dy);
    omega = forward_omega(psi,omega,nu,dx,dy,dt);
    [psi,flag_end] = get_stream(omega,dx,dy,1e-7,psi);
    [u,v] = get_u(psi,dx,dy,u,v);
    out_iter = out_iter + 1;

    if mod(out_iter, 100) == 0
        current_time = out_iter * dt;

        %subplot(1,2,1);
        subplot('Position', [0.03, 0.1, 0.3, 0.8]);
        %caxis([0,1]);
        %colormap('jet');
        contour(x, y, psi, 100);
        title(['Stream Function Visualization (Iteration: ', num2str(out_iter), ', Time: ', num2str(current_time), ' s)']);
        xlabel('x');
        ylabel('y');
        colorbar;  

        %subplot(1,2,2);
        for i = 1:21
            for j = 1:21
                uu(i,j)=u((i-1)*5+1,(j-1)*5+1);
                vv(i,j)=v((i-1)*5+1,(j-1)*5+1);
            end
        end
                

        subplot('Position', [0.37, 0.1, 0.3, 0.8]);
        %quiver(x,y,v,u,'LineWidth', 0.5);
        quiver(xx,yy,vv,uu,'LineWidth', 0.5);
        title('Velocity Field Visualization');
        xlabel('x');
        ylabel('y');

        subplot('Position', [0.70, 0.1, 0.3, 0.8]);
        contour(x, y, omega, 100);
        title('Vorticity Visualization');
        xlabel('x');
        ylabel('y');
        colorbar;  
        
        pause(0.05); 
    end
end



function [new_u,new_v] = get_u(psi,dx,dy,u,v)
    [Nx, Ny] = size(psi);
    new_u = u;
    new_v = v;
    new_u(2:end-1,2:end-1)=(psi(2:end-1,3:end)-psi(2:end-1,1:end-2))/2/dy;
    new_v(2:end-1,2:end-1)=-(psi(3:end,2:end-1)-psi(1:end-2,2:end-1))/2/dx;
end

function omega_new = vel_to_omega(omega,u,v,dx,dy)
    omega_new = omega;
    omega_new(:,2:end-1) = - (u(:,3:end)-u(:,1:end-2))/2/dy;
end

function [psi, end_signal] = get_stream(omega_int, dx, dy, tol, psi_init)
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5 || isempty(psi_init)
        [Nx, Ny] = size(omega_int);
        psi_init = zeros(Nx,Ny);
    end
    
    change = 1.0;
    iter = 0;
    psi = psi_init;
    
    while change > tol
        psi_new = psi;
        psi_new(2:end-1, 2:end-1) = (omega_int(2:end-1, 2:end-1) ...
            + (psi(3:end, 2:end-1) + psi(1:end-2, 2:end-1)) / dy / dy ...
            + (psi(2:end-1, 3:end) + psi(2:end-1, 1:end-2)) / dx / dx) ...
            / (2 / dx / dx + 2 / dy / dy);
        
        change = max(max(abs(psi_new - psi)));
        psi = psi_new;
        iter = iter + 1;
    end
    
    end_signal = (iter == 1);
end

function omega_new = get_bound_omega(psi,u,v,omega,dx,dy)
    [Nx,Ny] = size(omega);
    omega_new=zeros(Nx,Ny);
    omega_new = omega;
    omega_new(1,:) = 2/dx/dx*(psi(1,:)-psi(2,:))-2/dx*v(1,:);
    omega_new(end, :) = 2 / dy / dy * (psi(end, :) - psi(end-1, :)) +2/dx*v(end,:);
    omega_new(:, 1) = 2 / dx / dx * (psi(:, 1) - psi(:, 2)) + 2 / dx *u(:,1);
    omega_new(:, end) = 2 / dx / dx * (psi(:, end) - psi(:, end-1)) - 2 / dx * u(:, end);
end

function omega_new = forward_omega(psi,omega,nu,dx,dy,dt)
    uw_x = -(psi(2:end-1,3:end) - psi(2:end-1,1:end-2)) / 4 / dx / dy .* (omega(3:end,2:end-1) - omega(1:end-2,2:end-1));
    vw_y = +(psi(3:end,2:end-1) - psi(1:end-2,2:end-1)) / 4 / dx / dy .* (omega(2:end-1,3:end) - omega(2:end-1,1:end-2));
    diff_x = nu * (omega(3:end, 2:end-1) + omega(1:end-2, 2:end-1) - 2 * omega(2:end-1, 2:end-1)) / dy / dy;
    diff_y = nu * (omega(2:end-1, 3:end) + omega(2:end-1, 1:end-2) - 2 * omega(2:end-1, 2:end-1)) / dx / dx;
    omega_new = omega;
    omega_new(2:end-1,2:end-1) = omega(2:end-1, 2:end-1) + dt * (uw_x + vw_y + diff_y + diff_x);
end