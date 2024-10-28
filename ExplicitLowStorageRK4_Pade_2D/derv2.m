function [qxx,qyy]   = derv2(q)
% Compact 4th order finite differnce for second order derivative
Global;
alpha = 0.1;
ak = 1.2;

qxx = zeros(Nx,Ny); qyy = zeros(Nx,Ny); 

vert = 1:Ny; 
horz = 1:Nx;

a = zeros(Nx,1); b = zeros(Nx,1); c = zeros(Nx,1); d = zeros(Nx,1);
     a(horz) =  0.1;   b(horz) = 1.0;  c(horz) =  0.1; 
     a(1)    =  0.0;   b(1)    = 1.0;  c(1)    = 11.0;
     a(Nx)   = 11.0;   b(Nx)   = 1.0;  c(Nx)   =  0.0;


 for j = 1:Ny
     d(1)      = (13*q(1,j)-27*q(2,j)+15*q(3,j)-q(4,j))/dx/dx;
     d(2:Nx-1) = ak*(q(3:Nx,j)-2*q(2:Nx-1,j)+q(1:Nx-2,j))/dx/dx;
     d(Nx)     = (13*q(Nx,j)-27*q(Nx-1,j)+15*q(Nx-2,j)-q(Nx-3,j))/dx/dx;

     d = tridg(a,b,c,d,Nx);
      qxx(horz,j) = d(horz);
 end
 
a = zeros(Ny,1); b = zeros(Ny,1); c = zeros(Ny,1); d = zeros(Ny,1);
     a(vert) =  0.1;   b(vert) = 1.0;  c(vert) =  0.1; 
     a(1)    =  0.0;   b(1)    = 1.0;  c(1)    = 11.0;
     a(Ny)   = 11.0;   b(Ny)   = 1.0;  c(Ny)   =  0.0;

 for i = 1:Nx
     d(1)      = (13*q(i,1)-27*q(i,2)+15*q(i,3)-q(i,4))/dy/dy;
     d(2:Ny-1) = ak*(q(i,3:Ny)-2*q(i,2:Ny-1)+q(i,1:Ny-2))/dy/dy;
     d(Ny)     = (13*q(i,Ny)-27*q(i,Ny-1)+15*q(i,Ny-2)-q(i,Ny-3))/dy/dy;
     d = tridg(a,b,c,d,Ny);
     qyy(i,vert) = d(vert);
 end

 return




