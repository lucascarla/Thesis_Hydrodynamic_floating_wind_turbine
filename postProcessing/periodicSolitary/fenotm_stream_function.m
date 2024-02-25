function [eta,B,Ubar,k,Q,R, surface_elevation] = fenton(H,d,LorTandU,N,Hsteps,Npos,plt)

% This is an implementation of the Fourier approximation method for finding
% stream function waves described in section 3.1.2 of the paper by J.D.
% Fenton, "Numerical Methods for Nonlinear Waves" (1999), Numerical Methods
% for Nonlinear Waves, Advances in Coastal and Ocean Engineering, Vol. 5,
% pp241-324. Equation numbers in the code refer to this paper. At the time
% of writing the paper could be downloaded here:
%
% http://www.johndfenton.com/Papers/Fenton99Liu-Numerical-methods-for-nonlinear-waves.pdf
%
% [eta,B,Ubar,k,Q,R] = fenton(H,d,LorTandU,N,Hsteps)
%
% INPUT
% H: Wave height.
% d: Mean water depth.
% LorTandU: 1) Wave length or
% 2) A vector where
% a) LorTandU(1) is the wave period.
% b) LorTandU(2) is the current.
% c) LorTandU(3) is the type of current with
% LorTandU(3) = 1 corresponding to LorTandU(2) = u1 and
% LorTandU(3) = 2 corresponding to LorTandU(2) = u2.
% See eq. (3.13) in Fentons paper for details.
% N: Number of Fourier modes and bins along a half wave length.
% Hsteps: Number of steps in which to gradually increase the wave height.
% Npos: Number of points along half wave length in which to return
% surface elevation.
% plt: 0 by default, set to 1 to plot progress.
%
%
% OUTPUT
% eta: Double array containing Npos positions along x-axis and the
% corredponding surface elevation above seabed.
% B: Fourier coefficients in stream function expansion (3.5)
% Ubar: Mean horizontal velocity as seen from co-moving frame (3.5)
% k: Wave number = [wave length]/(2*pi)
% Q: Constant in kinematic boundary condition (3.7)
% R: Constant in dynamic boundary condition (3.8)
%
% J. Roenby, DHI, 28 June 2012
%
% Notes to self:
% Compared results with those of NGJ's routine. The eta's and k differ by
% around 1e-12. The B's, Q and R do not match due to the fact that NGJ's
% code does not work with nondimensionalized variables.
%
% Note that for very long waves almost all of the "action" takes place near
% the crest. Therefore N must be chosen very large to resolve the wave.
% Tried to change code to work with varying x step size, but jacobian
% became ill-conditioned.
% Also tried to implement a version where the number of points along the
% half wave length is larger that one plus the number of Fourier
% coefficients. The system is then overdetermined, and matlab automatically
% uses its QR factorisation algorithm when the "\" operator is called in
% the Fmin subfunction. This, however, resulted in a much worse convergence
% with F(x0) only reaching ~1e-5, with the original version reaching
% F(x0)~1e-12.
%
% A challenging case with H 98% of maximum allowed from page 3 in
% http://johndfenton.com/Steady-waves/Instructions.pdf
% is [eta,B,Ubar,k,Q,R] = fenton(0.786,1,50,140,20,500,1)

if nargin < 7
    plt = 1; %plot progress
end
if nargin < 5
    Hsteps = 1;
end
if nargin < 4
    N = 100; %N=Npos-1, so the number of wavegauges per half wavelength minus 1
end
if nargin < 3
    LorTandU = [14, 0,0]; %wavelength
end
if nargin < 2
    d = 20;
end
if nargin < 1
    H = 20; %Wave height
end
g = 9.81;

if numel(LorTandU) == 1 %wave length speficied by user
    lambda = LorTandU;
    k = 2*pi/lambda;
    %Nondimensionalized wave length
    LorTandU = lambda/d;
elseif numel(LorTandU) == 3 %wave period and current specified by user
    T = LorTandU(1);
    u = LorTandU(2);
    sigma = 2*pi/T;
    %Initial guess of wave number from linear theory (3.22)
    if sigma.^2*d/g >= 1
        kh0 = sigma.^2*d/g;
    else
        kh0 = sigma*sqrt(d/g);
    end
    k = fzero(@(kh) sigma.^2*d/g - kh.*tanh(kh),kh0)/d;
    %Nondimensionalized waver period and current
    LorTandU(1) = T*sqrt(g/d);
    LorTandU(2) = u/sqrt(g*d);
end

H0 = H/Hsteps;
kd = k*d;
keta = kd + 1/2*k*H0*cos([0:N]*pi/N);
u = sqrt(tanh(kd));
B = 1/2*k*H0/sqrt(tanh(kd))*[1 zeros(1,N-1)];
q = kd*sqrt(tanh(kd));
r = kd + 1/2*u^2;

if plt
    close all
    figure(1); clf;
    xx = [0:N]/N*pi/k;
    plot(xx,keta/k,'.-r')
% axis equal
    hold on
end

x = [keta, B, u, kd, q, r];
h = H/d;
% for n = 1:Hsteps
%     [x,F0] = Fmin(@conditions,x,n*h/Hsteps,LorTandU);
%     k = x(2*N+3)/d;
%     xx = [0:N]/N*pi/k;
%     eta = x(1:N+1)/k;
%     %surface_elevation=eta;
%     if plt
%         figure(2)
%         plot(xx,eta,'.r')
%         drawnow
%     end
% end

k = x(2*N+3)/d;
B = x(N+2:2*N+1);
eta = x(1:N+1)/k;
Ubar = x(2*N+2)/sqrt(k/g);
Q = x(2*N+4)/sqrt(k^3/g);
R = x(2*N+5)/(k/g);

%Cosine transform of surface elevation
%(3.26) with j running from 0 to N
[m,j] = meshgrid((0:N),(0:N));
etapp = [1/2*eta(1); eta(2:end-1)'; 1/2*eta(end)];
E = cos(j.*m*pi/N)*etapp;

%(3.27) with with division by N as required in the inverse cosine transform
E(1) = 1/2*E(1); E(end) = 1/2*E(end);
if nargin < 6
    Npos = N+1;
end
pos = [0:Npos-1]/(Npos-1)*pi/k;
[j,X] = meshgrid((0:N),pos);
ETA = 2/N*cos(j*k.*X)*E;
eta = [pos; ETA(:)'];

if plt
    figure(2)
    plot(pos,ETA,'.k')
    xlabel('x')
    ylabel('\eta')
end
function F = conditions(x,h,LorTandU)

%Returning a vector of length 2N+5 that will be the zero vector if the
%input corresponds to a wave solution.

N = (numel(x)-5)/2;

keta = x(1:N+1); %Height function eta nondimensionalized by wavenumber at the points kx = [0:N]*pi/N.
B = x(N+2:2*N+1); %The N Fourier coefficients in the stream function expansion.
u = x(2*N+2); %Nondimensionalized mean horizontal velocity in co-moving frame, Ubar*sqrt(k/g).
kd = x(2*N+3); %Mean depth normalized by wavenumber, k*d.
q = x(2*N+4); %Nondimensionalized volume flux in co-moving frame, q = Q*sqrt(k^3/g).
r = x(2*N+5); %Nondimensionalized constant in Bernoulli equation, R*k/g.

F = zeros(1,2*N+5);

%Creating meshes for efficient evaluation of boundary conditions
[ke,b] = meshgrid(keta,B);
[kx,j] = meshgrid([0:N]*pi/N,[1:N]);

%Kinematic boundary condition (3.7)
F(1:N+1) = - u*keta + sum(b.*sinh(j.*ke)./cosh(j*kd).*cos(j.*kx),1) + q;

%Dynamic boundary condition (3.8)
U = - u + sum(j.*b.*cosh(j.*ke)./cosh(j*kd).*cos(j.*kx),1);
V = sum(j.*b.*sinh(j.*ke)./cosh(j*kd).*sin(j.*kx),1);
F(N+2:2*N+2) = 1/2*U.^2 + 1/2*V.^2 + keta - r;

clear ke b kx j

%Mean depth equation (3.9)
F(2*N+3) = 1/N*(1/2*(keta(1) + keta(end)) + sum(keta(2:end-1))) - kd;
    
%Wave height equation (3.10)
F(2*N+4) = keta(1) - keta(end) - kd*h;

if numel(LorTandU) == 1
    l = LorTandU;
    %Relative wavelength equation (3.11)
    F(2*N+5) = kd - 2*pi/l;
elseif numel(LorTandU) == 3
    t = LorTandU(1);
    if LorTandU(3) == 1
        u1 = LorTandU(2);
        %Equation for mean current in laboratory frame (3.15)
        F(2*N+5) = sqrt(kd)*u + kd*u1 - 2*pi/t;
    elseif LorTandU(3) == 2
        u2 = LorTandU(2);
        %Equation for mass transport velocity (3.16)
        F(2*N+5) = q/sqrt(kd) + kd*u2 - 2*pi/t;
    end
end

function dF = jac(F,x,varargin)

%Calculates the jacobian

h = 1e-10; %HARDCODED STEP SIZE
dF = zeros(length(x));
for j = 1:length(x)
    dx = zeros(size(x));
    dx(j) = h;
    dF(:,j) = (F(x+dx,varargin{:})-F(x,varargin{:}))/h;
end

function [x0,F0] = Fmin(F,x0,varargin)

%Newton method for finding x0 = (x1,...,xN) such that F(x0) = (0,...,0).

tol = 1e-5; %HARDCODED TOLERANCE
err = 2*tol;
while err > tol
    F0 = F(x0,varargin{:}).';
dx = - jac(F,x0,varargin{:})\F0;
x0 = x0(:).' + dx(:).';
    err = max(abs(F0))
end
if nargout > 1
    F0 = F(x0,varargin{:});
end