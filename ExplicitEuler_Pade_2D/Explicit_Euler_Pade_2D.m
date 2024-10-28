% =========================================================================
% A Matlab R2021b program for conduction heat transfer
% 
% Created by:
% Theodoret Putra Agatho
% University of Atma Jaya Yogyakarta
% Department of Informatics
% 08/10/2024
% =========================================================================
% -------------------------------------------------------------------------
% Main Program
% -------------------------------------------------------------------------
clear all; clc; close all;

% Hyperparameters
alpha = 0.1;

Global;

% Initialize the domain's matrix
% dy = 0.03; dx = 0.03;
Nx = 32; xmin = 0; xmax = 4; % xmax = Nx*dy;
Ny = Nx; ymin = xmin; ymax = xmax; % ymax = Ny*dx;
dx = (xmax - xmin)/(Nx-1);
dy = (ymax - ymin)/(Nx-1);

x = linspace(xmin,xmax,Nx); y = linspace(ymin,ymax,Ny);
[X,Y] = meshgrid(x,y);


% Initialize the time and steps
tmin = 0; tmax = 1.0; time_step = 100;
dt = (tmax - tmin)/time_step;


% Initialize the initial condition
U = 100.*ones(Ny,Nx); % consider "U" as current y and "y" as next y

% Initialize the boundary condition
left = find(X == xmin); right = find(X == xmax);
top = find(Y == ymin); bottom = find(Y == ymax);
U = dirichlet(U,left,right); U = neumann(U,top,bottom);


% Time loop
loop_start = tic; % loop stopwatch
for i = 1:time_step
%     if (mod(i,time_step/10) == 0 || i == 1)
%         disp("time step = " + i);
%         figure(1); surf(Y,X,U); shading interp; view(0,90);
%     end
    [hyy,hxx] = derv2(U);

    U = U + dt.*alpha.*(hyy + hxx);

    U = dirichlet(U,left,right); 
    U = neumann(U,top,bottom);
end
loop_end = toc(loop_start);


% Analytical (/exact) Solution
L = ymax;
Uetemp = zeros(Ny,Nx);

mmax = 201;
for m = 1:mmax
    Uetemp = Uetemp + exp(-(m*pi/L).^2.*alpha.*tmax).*((1-(-1).^m)./(m.*pi)).*sin(m*pi*y/L);
end
Ue = 300 + 2*(100-300).*Uetemp;


% Visualisation
figure(1); surf(Y,X,U); shading interp;
figure(2); surf(Y,X,Ue); shading interp;


% Error
error = mean(sqrt((U-Ue).^2),'all')

% -------------------------------------------------------------------------
% Functions
% -------------------------------------------------------------------------
% Boundary Condition Function =============================================
function u = dirichlet(u,left,right)
    u(left) = 300; u(right) = 300;
end

function u = neumann(u,top,bottom)
    u(top) = u(top+1); u(bottom) = u(bottom-1);
end