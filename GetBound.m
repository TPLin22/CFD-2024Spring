function [q_b] = GetBound(q,jacob)

%
% q_b = q;

%r0(:,:) = 1.185;
%u0(:,:) = 0;
%v0(:,:) = 99719;
%p0(:,:) = 686.47;
gamma = 1.4;

rho=q(:,:,1).*jacob(:,:); 
u=q(:,:,2)./rho.*jacob(:,:); 
v=q(:,:,3)./rho.*jacob(:,:); 
E=q(:,:,4)./rho.*jacob(:,:);
p=(gamma-1)*rho.*(E-0.5*(u.^2 + v.^2));

rho_b = rho;
u_b = u;
v_b = v;
E_b = E;
p_b = p;

% left
u_b(1,:) = 686.47;
rho_b(1,:) = 1.185;
v_b(1,:) = 0;
p_b(1,:) = 99719;

%right
u_b(end,:) = u_b(end-1,:);
rho_b(end,:) = rho_b(end-1,:);
v_b(end,:) = v_b(end-1,:);
p_b(end,:) = p_b(end-1,:);

%bottom
rho_b(:,1) = rho_b(:,2);
p_b(:,1) = p_b(:,2);
u_b(1:20,1) = 686.47;
u_b(21:end,1) = 663.08;
v_b(1:20,1) = 0;
v_b(21:end,1) = 177.67;

%top
rho_b(:,end) = rho_b(:,end-1);
p_b(:,end) = p_b(:,end-1);
u_b(:,end) = u_b(:,end-1);
v_b(:,end) = v_b(:,end-1);

q_b=reshape([rho_b./jacob rho_b.*u_b./jacob rho_b.*v_b./jacob rho_b.*E_b./jacob],80,50,4);
end