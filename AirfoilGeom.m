% Airfoil curve fitting 
% from Kulfan 'Universal Parametric Geometry Representation' 2008

clc; clear; close all;

af = readmatrix("naca2412.txt");
N = length(af);
Np = 20;

%% Kulfan
psi = flip(af(1:ceil(N/2),1));
ZETAu = flip(af(1:ceil(N/2),2));
ZETAl = af(ceil(N/2):N,2);


N1 = 0.5;
N2 = 1.0;

C = psi.^(N1) .* (1-psi).^(N2);

Delta_zeta_u = af(1,2);
Delta_zeta_l = af(N,2);

% Bernstein Poly
n = 5;
i = linspace(0,n,n+1);

K = factorial(n)./(factorial(i).*factorial(n-i));
S = K.*psi.^i .*(1-psi).^(n-i);
BigS = C.*S;

Bu = ZETAu - psi*Delta_zeta_u;
Bl = ZETAl - psi*Delta_zeta_l;

Au = lsqr(BigS,Bu, 1e-8);
Al = lsqr(BigS,Bl, 1e-8);
zetau = BigS*Au + psi*Delta_zeta_u;
zetal = BigS*Al + psi*Delta_zeta_l;

% use weights for larger number of x spacing
psi = 0.5*(1- cos(linspace(0,pi,Np)))';
C = psi.^(N1) .* (1-psi).^(N2);

Sf = K.*psi.^i .*(1-psi).^(n-i);
BigSf = C.*Sf;


zetau = BigSf*Au + psi*Delta_zeta_u;
zetal = BigSf*Al + psi*Delta_zeta_l;

coords = [[flip(psi); psi(2:Np)] [flip(zetal); zetau(2:Np)]];
coords2 = AFGeom(af,Np,n);

% Plot
figure
plot(flip(af(1:ceil(N/2),1)), ZETAu,'o',flip(af(1:ceil(N/2),1)),ZETAl,'o')
hold on
plot(psi,zetau)
plot(psi,zetal)
plot(coords(:,1),coords(:,2),'o')
xlabel('x/c'); ylabel('y/c');
axis equal