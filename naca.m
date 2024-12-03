% NACA Plotting Scratch to test naca4 function
clear; clc; close all;

m = .00;
p = 0.0;
t = .12;

n = 51;
x = linspace(0,1,n);
if p ==0
    y =0.*x;   
else
    y = m/p^2 * (2*p.*x - x.^2);
end
yc = (m/(1-p)^2)*((1-2*p) + 2*p.*x - x.^2);
px = p*(n-1)+1;

yc(1:px) = y(1:px);

yt = (t/0.2)*(0.2969.*x.^0.5 - 0.126.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1015.*x.^4);

theta = atan2(diff(yc), diff(x));
theta = [theta theta(end)];
xu = x - yt.*sin(theta);
xl = x + yt.*sin(theta);
yu = yc + yt.*cos(theta);
yl = yc - yt.*cos(theta);

% figure
% plot(x,yu,x,yl)
% 
% X = [flip(x) x]';
% Y = [flip(yu) yl]';
% 
% figure
% plot(X,Y)

[X,Y] = naca4('2412',51);

figure
plot(X,Y)
axis('equal')





