%% Import the line-by-line Green's functions, numerically calculated by
%% Mathematica, to build the 2D Green's functions.
%% Then, the 3D heat source density is calculated, and then convoluted with
%% the Green's function to compute thte 3D temperature profile.
%% Computes the data for the Figure 6g-i, i.e., the temperature maps
%% as a function of h and lambda.

%Authors: Guillaume Baffou
%affiliation: CNRS, Institut Fresnel 
%Date: Jan 2021

% Used to make Figure 6ghi of ref https://doi.org/10.1063/5.0057185

% GitHub: https://github.com/baffou/ThermoNanoLayer
% Please cite https://doi.org/10.1063/5.0057185 if you use this code in a publication. Thx.


clear
addpath(genpath('functions'))

%% Parameters

h2 = [1 5:5:100]*1e-9; % metal layer thickness [m]
Rmax  = 10e-6 ;   % max radial coordinate for the temperature calculation [m]
Nr    = 200 ;       % number of calculated points along r
lambda= (400:25:1200)*1e-9;     % wavelength [m]
NA    = 1.0;        % Numerical aperture of the illumination
mat1  = 'BK7';      % Bottom medium, from which light is coming
mat2  = 'Au';       % middle, absorbing layer
mat3  = 'Water';    % top medium
Plaser=1e-3;    % laser power [W]

%% calculations
c0 = 299792458;
eps0 = 8.85418782e-12;

% materials' properties
n1 = indexRead(lambda,mat1);
n2 = indexRead(lambda,mat2);
n3 = indexRead(lambda,mat3);

% illumination properties
Rlaser    = 0.61*lambda/NA ;    % laser radius
I0 = Plaser./(2*pi*Rlaser.^2);    % laser irradiance at the center of the Gaussian shape
E0 = sqrt(2*I0./(n1*c0*eps0));    % electric field amplitude at the center of the beam, before entering the material 

skinDepth = lambda./(4*pi*imag(n2));

Nl = length(lambda);
Nh = length(h2);

t12 = 2*n1  ./(n1+n2);
r12 = (n1-n2)./(n1+n2);
r23 = (n2-n3)./(n2+n3);


%% Import the 1D Green's functions G{zs}(z,r) calculated with Mathematica
%% and build the 2D Green's functions

rList = dlmread('rList.txt');

Nr0 = length(rList);
Nz0 = 10; %Number of Green's function along z
GreenProfiles0 = cell(Nz0,1);

for izs = 1:Nz0
    GreenProfiles0{izs} = zeros(Nz0,Nr0);
    for iz = 1:Nz0
        GreenProfiles0{izs}(iz,:) = dlmread(['GreensFunctions10_GlassGoldWater/G' num2str(izs) 't' num2str(iz) '.txt'])';
    end
end


%% loop over wavelength lambda and layer thickness h
Tmax = zeros(Nl,Nh);
P0 = zeros(Nl,Nh);
HeatmapXZ = cell(Nl,Nh);
TmapXZ = cell(Nl,Nh);
HeatmapXY = cell(Nl,Nh);
Qtotal = zeros(Nl,Nh);
Ptotal = zeros(Nl,Nh);
figure
subplot(1,2,1)
imagesc(Tmax);
colorbar
subplot(1,2,2)
imagesc(P0);
colorbar
drawnow
for ih = 1:Nh
    % Mathematica non-uniform mesh grid, conversion in meter unit
    rListm = rList*h2(ih);
    zListm = ((1:Nz0)-0.5)/Nz0*h2(ih);
    dz = 1/Nz0*h2(ih);
    % Definition of a uniform mesh grid
    Rm2 = (0:2*Nr-1)/(2*Nr-1)*2*Rmax; % factor of 2 for the further convolution
    [Rgrid,Zgrid] = meshgrid(Rm2,zListm);

    %% defining the coordinates vectors

    Rm = (0:Nr-1)/(Nr-1)*Rmax;    % radial coordinates [0:Rmax]
    dr = 1/(Nr-1)*Rmax;
    Rms = [flip(Rm(:,2:end),2),Rm];    % radial coordinates [-Rmax:Rmax]
    [Xgrid3D,Ygrid3D,Zgrid3D] = meshgrid(Rms,Rms,zListm); % 3D matrix of xyz coordinates
    Rgrid3D2 = Xgrid3D.^2+Ygrid3D.^2;
    Xgrid2D = (-(2*Nr-1):2*Nr-1)/(2*Nr-1)*2*Rmax; % factor of 2 for the subsequent convolution
    [X2D,Y2D] = meshgrid(Xgrid2D,Xgrid2D);
    RR = sqrt(X2D.^2+Y2D.^2);


    %% Calculation of the x-y Green's function G(zs,0,0;z,x,y) for convolution

    GreenXY = cell(Nz0,Nz0);

    for iz = 1:Nz0
        for izs = 1:Nz0
            Gfunc = GreenProfiles0{izs}(iz,:);
            GreenXY{izs}{iz} = interp1(rListm/h2(ih),Gfunc/h2(ih),RR/h2(ih));
        end
    end

    Green = cell(Nz0,1);
    GreenBig = cell(Nz0,1);

    for izs = 1:Nz0
        Green{izs}  =  interp2( rListm , zListm , GreenProfiles0{izs} , Rgrid , Zgrid) ;
        GreenBig{izs} = [flip(Green{izs}(:,2:end),2),Green{izs}];
    end

%    parfor il = 1:Nl
    for il = 1:Nl
        k0 = 2*pi./lambda(il);
        omega = k0*c0;

        %% defining the heat source density
        % 3F fields of Ez and q
        Ez = E0(il)*t12(il)*(exp(1i*k0*n2(il)*Zgrid3D)+r23(il)*exp(1i*k0*n2(il)*(2*h2(ih)-Zgrid3D)))./(1+r12(il)*r23(il)*exp(2*1i*k0*n2(il)*h2(ih)));
        dV = dr*dr*dz;
        Qtotal(il,ih) = omega*eps0*imag(n2(il)^2)*sum(abs(Ez(:)).^2)/2*dV; % total power integrated/absorbed within the layer
        Ptotal(il,ih) = I0(il)*(2*Rmax)^2;
        
        % Gaussian heating:
        Heat = exp(-Rgrid3D2./(2*Rlaser(il).^2)).*(sqrt(Rgrid3D2)<=2.5*Rlaser(il)).*omega*eps0.*imag(n2(il)^2).*abs(Ez).^2/2;
        P0(il,ih) = sum(Heat(:))*dV;
        
        Heat = P0(il,ih)*Heat/(sum(Heat(:)*dV));
        
        disp([num2str(ih) '/' num2str(Nh)])
        

        HeatmapXZ{il,ih} = permute(Heat(:,round(end/2),:),[1,3,2]);
        HeatmapXY{il,ih} = Heat(:,:,1);

        
        %% Convolution, plane by plane

        Tmap = Heat*0;
        for iz = 1:Nz0
            TmapBig = 0;
            for izs = 1:Nz0
                TmapBig = TmapBig+conv2(Heat(:,:,izs),GreenXY{izs}{iz});
            end
            Tmap(:,:,iz) = crop(TmapBig,2*Nr-1,2*Nr-1)*dV;
        end
        TmapXZ{il,ih} = permute(Tmap(:,round(end/2),:),[1,3,2]);
        Tmax(il,ih) = max(Tmap(:));

        drawnow
    end %Delta
end %lambda


%%

figure
subplot(1,2,1)
imagesc(1e9*h2,1e9*lambda,100*P0/Plaser);
title('Absorption (%)')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')


subplot(1,2,2)
imagesc(1e9*h2,1e9*lambda,Tmax);
title('temperature')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')


%%
hindexmin = 1;
P0big = imresize(P0(:,hindexmin:end),8);
TmaxBig = imresize(Tmax(:,hindexmin:end),8);

figure
subplot(1,2,1)
him = imagesc(1e9*h2(hindexmin:end),1e9*lambda,100*P0big/Plaser);
title('Absorption (%)')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')


subplot(1,2,2)
him2 = imagesc(1e9*h2(hindexmin:end),1e9*lambda,TmaxBig);
title('temperature')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')
drawnow


%% Save data

% dlmwrite('P0.txt',100*P0/Plaser)
% dlmwrite('Tmax.txt',Tmax)
% dlmwrite('P0biground.txt',round(100*P0big/Plaser))
% dlmwrite('TmaxBiground.txt',round(TmaxBig/2)*2)
% dlmwrite('P0big.txt',100*P0big/Plaser)
% dlmwrite('TmaxBig.txt',TmaxBig)
% dlmwrite('logP0big.txt',log10(100*P0big/Plaser))
% dlmwrite('logTmaxBig.txt',log10(TmaxBig))


