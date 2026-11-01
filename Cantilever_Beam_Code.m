%--------------------------------------------------------------------------
%---- Natural frequency and mode shape of a beam with ROTATIONAL SPRING ---
%---- support at the root (x = 0) and a FREE tip (x = L)                ---
%---- Rayleigh-Ritz (assumed-modes) solution, following the Hamilton's   --
%---- principle derivation in the "Rotational_Spring_Support" notes      --
%---- Modified below to replace the rigid clamp at x = 0 with an         --
%---- elastic rotational restraint (spring, k_theta); the tip at x = L   --
%---- remains free. Changes are flagged in the comments.                 --
%--------------------------------------------------------------------------
clc; clear;
xl = 0.0;  xu = 1.0; ng = 48; nf = 101;
nw = 3; % Number of functions for displacement field
xi = zeros(1,ng);    xi = gauleg(ng);     xn = linspace(0,1,nf);
%--------------------------------------------------------------------------

%---------------------------------Beam geometry----------------------------
L = 0.5;  % Length of the beam (m)
b = 0.03; % width (m)
h = 0.01; % thickness (m)
A = b*h; % cross-sectional area (m^2)
I = (1/12)*b*h^3; % Area moment of inertia (m^4)
%--------------------------------------------------------------------------

%----------------------------Material properties---------------------------
E = 210e9; % Young's modulus (N/m^2)
rho = 7800; % Material density (kg/m^3)
%--------------------------------------------------------------------------

%--------------------[NEW] Rotational spring support at x = 0--------------
k_theta = 5000; % Rotational spring stiffness at the root (N.m/rad)
% Non-dimensional spring parameter K* = k_theta*L/(E*I), matching the
% (k_theta/L^2) term that appears in the non-dimensional energy eqn. (13):
%   K* -> 0        : spring offers no rotational resistance -> pinned-free beam
%   K* -> infinity : spring is rigid                        -> clamped-free cantilever
Kstar = k_theta*L/(E*I);
%--------------------------------------------------------------------------

%--------------------[CHANGED] Set of coordinate (trial) functions---------
% The ONLY essential (geometric) boundary condition here is v(0) = 0
% (translation is still prevented at the root); the ROTATION at x = 0 is
% now elastically restrained by the spring, i.e. it is a NATURAL boundary
% condition that is satisfied weakly through the extra spring-stiffness
% term added to K below (see eqns. 13-19), NOT built into phi_j.
%
% phi_j therefore only needs phi_j(0) = 0, e.g. (eqn. 15):
for j = 1:nw
    phi(j,:) = sin(j*pi*xi);
end
phid  = rd(phi);
phidd = rd(phid);

% The Gauss-Legendre nodes returned by gauleg() are strictly interior to
% (0,1), so phid/phidd cannot be sampled exactly at xi = 0. The slope
% needed for the spring term is instead taken analytically from
% phi_j = sin(j*pi*xi):
%   d(phi_j)/d(xi) |{xi=0} = j*pi*cos(j*pi*xi) |{xi=0} = j*pi
dphi0 = pi*(1:nw);
%--------------------------------------------------------------------------

%--------------------------------Computation-------------------------------
% [CHANGED] Stiffness matrix = beam-bending term + rotational-spring term
% (eqn. 19):  K(k,j) = (EI/L^3)*int(phi_j'' phi_k'') + (k_theta/L^2)*phi_j'(0)*phi_k'(0)
for k = 1:nw
    for j = 1:nw
        for l = 1:ng
            dIsm(l) = phidd(k,l)*phidd(j,l);
        end
        Kbeam  = (E*I/L^3)*mqg(dIsm,0,1);
        Kspr   = (k_theta/L^2)*dphi0(k)*dphi0(j);
        K(k,j) = Kbeam + Kspr;
    end
end

% Mass matrix is unaffected by the spring (eqn. 19)
for k = 1:nw
    for j = 1:nw
        for l = 1:ng
            dImm(l) = phi(k,l)*phi(j,l);
        end
        M(k,j) = (rho*A*L)*mqg(dImm,0,1);
    end
end

[eigvctrm,omegasqm] = eig(K,M);
omegasq = diag(omegasqm);
omega = sqrt(omegasq);

omega1 = omega(1)
omega2 = omega(2)
omega3 = omega(3)

eigvctr1 = eigvctrm(:,1);
eigvctr2 = eigvctrm(:,2);
eigvctr3 = eigvctrm(:,3);

mdsp1 = displacement(nw,eigvctr1,phi);
mdsp2 = displacement(nw,eigvctr2,phi);
mdsp3 = displacement(nw,eigvctr3,phi);

mdsp1nf = spline(xi,mdsp1,xn);
mdsp2nf = spline(xi,mdsp2,xn);
mdsp3nf = spline(xi,mdsp3,xn);
%--------------------------------------------------------------------------

%----------------------[CHANGED] Exact (transcendental) solution-----------
% For a uniform Euler-Bernoulli beam that is translation-fixed but
% elastically restrained in rotation (spring k_theta) at x = 0, and free
% at x = L, the exact eigenvalue equation (derived from the same boundary
% conditions used above) is:
%
%   lambda*(sinh(lambda)*cos(lambda) - sin(lambda)*cosh(lambda)) ...
%                              + Kstar*(1 + cos(lambda)*cosh(lambda)) = 0
%
% with   omega = lambda^2*sqrt(E*I/(rho*A*L^4))
%
% Limiting cases (useful checks):
%   Kstar -> 0        :  tan(lambda) = tanh(lambda)            (pinned-free beam)
%   Kstar -> infinity :  cos(lambda)*cosh(lambda) = -1          (clamped-free cantilever,
%                                                    lambda = 1.875, 4.694, 7.855, ...)
freqeqn = @(lam) lam.*(sinh(lam).*cos(lam) - sin(lam).*cosh(lam)) ...
                 + Kstar*(1 + cos(lam).*cosh(lam));

lam_scan  = 0.05:0.001:15;
f_scan    = freqeqn(lam_scan);
lam_roots = [];
for m = 1:length(lam_scan)-1
    if f_scan(m)*f_scan(m+1) < 0
        r = fzero(freqeqn,[lam_scan(m) lam_scan(m+1)]);
        if isempty(lam_roots) || min(abs(lam_roots-r)) > 1e-3
            lam_roots(end+1) = r; %#ok<AGROW>
        end
        if length(lam_roots) >= 3
            break
        end
    end
end

omega1an = lam_roots(1)^2*sqrt(E*I/(rho*A*L^4))
omega2an = lam_roots(2)^2*sqrt(E*I/(rho*A*L^4))
omega3an = lam_roots(3)^2*sqrt(E*I/(rho*A*L^4))
%--------------------------------------------------------------------------

%----------------------------Mode shape plots------------------------------
figure(1)
subplot(2,2,1)
plot(xn,mdsp1nf/max(mdsp1nf));
title('Mode 1'); xlabel('\xi = x/L'); ylabel('Normalized v');
subplot(2,2,2)
plot(xn,mdsp2nf/(max(mdsp2nf)-min(mdsp2nf)));
title('Mode 2'); xlabel('\xi = x/L'); ylabel('Normalized v');
subplot(2,2,3)
plot(xn,mdsp3nf/(max(mdsp3nf)-min(mdsp3nf)));
title('Mode 3'); xlabel('\xi = x/L'); ylabel('Normalized v');
%--------------------------------------------------------------------------
%------------------------------------End-----------------------------------
%--------------------------------------------------------------------------
