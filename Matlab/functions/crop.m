function IMout = crop(IM,Ny,Nx)
[Nx0,Ny0] = size(IM);
x0 = round(Nx0/2-Nx/2);
y0 = round(Ny0/2-Ny/2);

IMout = IM(y0:y0+Ny-1,x0:x0+Nx-1);












