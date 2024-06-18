clear; %close all; clc;
global gamma

%% Parameters
CFL     = 0.75;	% CFL number
tFinal	= 0.80;	% Final time
nEx      = 80;  % Number of cells in x
nEy      = 50;  % Number of cells in y
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
plot_fig= 1;

% Discretize spatial domain
a=0; b=4; dx=(b-a)/nEx; nx=nEx; x=linspace(a,b,nx);
c=0; d=2.4; dy=(d-c)/nEy; ny=nEy; y=linspace(c,d,ny);
[Y,X] = meshgrid(y,x);
[X0,Y0,psi_x,psi_y,eta_x,eta_y,J] = Mesh0();
%disp(Y(end,:));
%disp(eta_x(end,:));
%disp(J(end,:));

[rho0,u0,v0,p0,tFinal,CFL] = Euler_IC2d(X,Y);
E0 = p0./((gamma-1)*rho0)+0.5*(u0.^2+v0.^2);  % Total Energy density
a0 = sqrt(gamma*p0./rho0);            % Speed of sound
q0=reshape([rho0 rho0.*u0 rho0.*v0 rho0.*E0],nx,ny,4);        % vec. of conserved properties


% Discretize time domain
lambda0=max(max(abs(u0)+abs(v0)+a0)); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

%% Solver Loop
figure('Position', [100, 100, 1400, 400]);
%fg1=figure;
%fg2=figure;
%fg3=figure;
% h=createButton;
while t<tFinal

    % RK Initial step
    %qo = q;
    qreal = q;
    for i=1:nx
        for j=1:ny
            for k=1:4
                q(i,j,k)=q(i,j,k)/J(i,j);
            end
        end
    end
    qo = q;
    
    q = GetBound(q,J);
    % 1st stage
    dFx=WENO5LF2d(lambda,q,dx,1,eta_x,eta_y,J);
    dFy=WENO5LF2d(lambda,q,dy,2,eta_x,eta_y,J);
    dF = dFx+dFy;
    q = qo-dt*dF;
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);

    q = GetBound(q,J);
    % 2nd Stage
    dFx=WENO5LF2d(lambda,q,dx,1,eta_x,eta_y,J);
    dFy=WENO5LF2d(lambda,q,dy,2,eta_x,eta_y,J);
    dF = dFx+dFy;
    q = 0.75*qo+0.25*(q-dt*dF);
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);

    q = GetBound(q,J);
    % 3rd stage
    dFx=WENO5LF2d(lambda,q,dx,1,eta_x,eta_y,J);
    dFy=WENO5LF2d(lambda,q,dy,2,eta_x,eta_y,J);
    dF = dFx+dFy;
    q = (qo+2*(q-dt*dF))/3;
    q(1,:,:)=q0(1,:,:); q(end,:,:)=q0(end,:,:); q(:,1,:)=q0(:,1,:); q(:,end,:)=q0(:,end,:);
    q(1,:,:)=q(2,:,:); q(end,:,:)=q(end-1,:,:); q(:,1,:)=q(:,2,:); q(:,end,:)=q(:,end-1,:);

    q = GetBound(q,J);

    for i=1:nx
        for j=1:ny
            for k=1:4
                q(i,j,k)=q(i,j,k)*J(i,j);
            end
        end
    end

    % compute primary properties
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho; p=(gamma-1)*rho.*(E-0.5*(u.^2+v.^2));
    a=sqrt(gamma*p./rho); if min(p)<0; error('negative pressure found!'); end

    % Update dt and time
    lambda=max(max(abs(u)+abs(v)+a)); dt=CFL*dx/lambda; if t+dt>tFinal; dt=tFinal-t; end

    % Update time and iteration counter
	t=t+dt; it=it+1;

    %figure;
    
    %streamline(X0, Y0, u, v, X0, Y0);
    %title(['Streamlines, t=' num2str(t)]);
    %xlabel('X');
    %ylabel('Y');
    %cla;quiver(X0, Y0, u, v);
    %title(['Velocity Field, t=' num2str(t)]);
    %xlabel('X');
    %ylabel('Y');
    %cla;plot(x,u); grid on; hold on; scatter(x,u);
    %figure(fg1);
    %cla;fdisplay(X0,Y0,rho);title(['rho,t=' num2str(t)]);
    %colorbar;
    %figure(fg2);
    %cla;fdisplay(X0,Y0,p);title(['p,t=' num2str(t)]);
    %colorbar;
    %figure(fg3);
    %cla;fdisplay(X0,Y0,(sqrt(u.^2+v.^2))./a);title(['Ma,t=' num2str(t)]);
    %colorbar;  % 添加颜色条显示量级


    subplot(1,3,1); % Create subplot 1
    fdisplay(X0,Y0,rho); title(['Density \rho, t=' num2str(t)]);
    colorbar;
    
    subplot(1,3,2); % Create subplot 2
    fdisplay(X0,Y0,p); title(['Pressure p, t=' num2str(t)]);
    colorbar;
    
    subplot(1,3,3); % Create subplot 3
    fdisplay(X0,Y0,(sqrt(u.^2+v.^2))./a); title(['Mach number Ma, t=' num2str(t)]);
    colorbar;

    pause(0.001)

end
