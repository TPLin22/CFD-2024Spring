function res = WENO5LF2d(a,w,dx,dir,e_x,e_y,jacob)

if dir==1
    turn=[1 0 0];
    G = F(w,1,e_x,e_y,jacob);
else
    turn=[0 1 0];
    G = F(w,2,e_x,e_y,jacob);
end


v=0.5*(G+a*w); u=circshift(0.5*(G-a*w),-1*turn);

%% Right Flux
vmm = circshift(v,2*turn);
vm  = circshift(v,turn);
vp  = circshift(v,-1*turn);
vpp = circshift(v,-2*turn);

% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2;
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux
umm = circshift(u,2*turn);
um  = circshift(u,turn);
up  = circshift(u,-1*turn);
upp = circshift(u,-2*turn);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2;
B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

%% Compute finite volume residual term, df/dx.
res = (hp-circshift(hp,turn)+hn-circshift(hn,turn))/dx;
end

% Compute flux vector
function flux = F(q,i,e_x,e_y,jacob)
    global gamma

    % primary properties
    rho=q(:,:,1).*jacob(:,:); u=q(:,:,2)./rho.*jacob(:,:); v=q(:,:,3)./rho.*jacob(:,:); E=q(:,:,4)./rho.*jacob(:,:);
    p=(gamma-1)*rho.*(E-0.5*(u.^2 + v.^2));

    % flux vector of conserved properties
    if(i==1)
        %flux=reshape([rho.*u rho.*u.^2+p rho.*u.*v u.*(rho.*E+p)],size(q));
        flux=reshape([rho.*u./jacob (rho.*u.^2+p)./jacob rho.*u.*v./jacob u.*(rho.*E+p)./jacob],size(q));
    else
        %flux=reshape([rho.*v rho.*u.*v rho.*v.^2+p v.*(rho.*E+p)],size(q));
        flux=reshape([(e_x.*rho.*u+e_y.*rho.*v)./jacob (e_x.*(rho.*u.^2+p)+e_y.*rho.*u.*v)./jacob (e_y.*(rho.*v.^2+p)+e_x.*rho.*u.*v)./jacob (e_y.*v.*(rho.*E+p)+e_x.*u.*(rho.*E+p))./jacob],size(q));
    end
end
