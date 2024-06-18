figure;
ax0 = axes;
hold(ax0,'on');
grid(ax0,'on');
axis(ax0,[0 3 -1.5 1.5]);

figure;
ax1 = axes;
hold(ax1,'on');
grid(ax1,'on');
axis(ax1,[0 0.04 0 0.09]);
xlabel(ax1,'dx');
ylabel(ax1,'error');
title(ax1,'error and spatial step size');

figure;
ax2 = axes;
hold(ax2,'on');
grid(ax2,'on');
%axis(ax2,[0 0.11 0 0.11]);
xlabel(ax2,'dt');
ylabel(ax2,'error');
title(ax2,'error and time step size');


u_acc=laxwendroff(ax0,3,1,3000,80000,1);

dx_values = [];
dx2_values = [];
xerrors = [];
x2errors = [];
dt_values = [];
terrors = [];

for nx = 100:2:500
    u_xn = laxwendroff(ax0,3,1,nx,3000,0);
    err = u_xn-u_acc;
    dx_n = 3 / nx;

    dx_values(end+1) = dx_n;
    xerrors(end+1) = abs(err);
end

plot(ax1,dx_values,xerrors);

for nx = 100:2:500
    u_xn = newlaxwendroff(ax0,3,1,nx,0);
    err = u_xn-u_acc;
    dx_n = 3 / nx;

    dx2_values(end+1) = dx_n;
    x2errors(end+1) = abs(err);
end

plot(ax1,dx2_values,x2errors);

for nt = 100:2:1000
    u_xn = laxwendroff(ax0,3,1,100,nt,2);
    err = u_xn-u_acc;
    dt_n = 1 / nt;

    dt_values(end+1) = dt_n;
    terrors(end+1) = abs(err);
end
plot(ax2,dt_values,terrors);

function u_mid = laxwendroff(ax,L,T,Nx,Nt,ifacc)

    dx = L / Nx;     
    dt = T / Nt;     
    a = 1;
    nu = 0.04;
    c = a * dt / dx;
    d = nu * dt / dx / dx;

    x = linspace(0, L, Nx+1);
    x = x(1:end-1);
    u = zeros(Nx, Nt+1);

    u(:, 1) = sin(2 * pi * x);

    for n = 1:Nt
        un = u(:, n);
        unp1 = u(:, n+1);
        for i = 2:Nx-1
            unp1(i) = (d - c/2) * un(i+1) + (1 - 2*d) * un(i) + (d + c/2) * un(i-1);
        end

        unp1(1) = unp1(Nx-1);
        unp1(Nx) = unp1(2);
        u(:, n+1) = unp1;
    end
    
    if ifacc == 1
        plot(ax,x,u(:,end));
        xlabel(ax,'x');
        ylabel(ax,'u(x,t)');
        title(ax,sprintf('Time = 1,dx=%.7f,dt=%.7f',dx,dt));
    end
  
    if ifacc == 2
        u_mid = u(floor(Nx/4)+1,end);
    else
        u_mid = u(floor(Nx/2)+1,end);
    end
end

function u_mid = newlaxwendroff(ax,L,T,Nx,ifacc)

    dx = L / Nx;     
    %dt = T / Nt;
    dt = (sqrt(36*0.04*0.04+4*dx*dx)-6*0.04)/2;
    a = 1;
    Nt = ceil(1/dt);
    nu = 0.04;
    c = a * dt / dx;
    d = nu * dt / dx / dx;
  
    x = linspace(0, L, Nx+1);
    x = x(1:end-1);
    u = zeros(Nx, Nt+1);

    u(:, 1) = sin(2 * pi * x);

    for n = 1:Nt
        un = u(:, n);
        unp1 = u(:, n+1);
        for i = 2:Nx-1
            unp1(i) = (d - c/2) * un(i+1) + (1 - 2*d) * un(i) + (d + c/2) * un(i-1);
        end

        unp1(1) = unp1(Nx-1);
        unp1(Nx) = unp1(2);
        u(:, n+1) = unp1;
    end
    
    if ifacc == 1
        plot(ax,x,u(:,end));
        xlabel(ax,'x');
        ylabel(ax,'u(x,t)');
        title(ax,sprintf('Time = 1,dx=%.7f,dt=%.7f',dx,dt));
    end
  
    %u_mid = u(floor(Nx/2)+1,end);
    u_mid = u(1,end);
    %u_mid = u(end,end);
end