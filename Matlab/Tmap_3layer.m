%% Computes the 3D temperature profile (Tmap) and the 3D heat source density (Heat)
% within a finite-thickness absorbing layer, embedded in two semi-infinite media
% when the layer is illuminated/heated by a laser focused on the layer by a
% microscope objective lens of given NA.
% For this purpose, the program imports the line-by-line Green's functions,
% numerically calculated using Mathematica, to build the 2D Green's functions.
% Then, the 3D heat source density is calculated, and convoluted with
% the Green's function to compute the 3D temperature profile.
% Also compute the absorbance A of the layer via absorbed power calculations.

%Benoit Rogez, Guillaume Baffou
%CNRS, Institut Fresnel 
%Jan 2021

% Used to make Figure 4de of ref https://doi.org/10.1063/5.0057185

% GitHub: https://github.com/baffou/ThermoNanoLayer
% Please cite https://doi.org/10.1063/5.0057185 if you use this code in a publication. Thx.



clear
addpath(genpath('functions'))

%% Parameters

Delta  = 20e-9;  % metal layer thickness [m]
Rmax   = 10e-6;  % max radial coordinate for the temperature calculation [m]
Nr     = 200;    % number of calculated points along r
lambda = 800e-9;% wavelength [m]
NA     = 0.7;    % Numerical aperture of the illumination
Plaser = 1e-3;   % laser power [W]

% Materials. Must correspond to the thermal conductivities used in the Mathematica notebook to compute the Green's functions

mat1  = 'BK7';      % Bottom medium, from which light is coming
mat2  = 'Au';       % middle, absorbing layer
mat3  = 'Water';    % top medium

%% calculations
c0  = 299792458;
eps0= 8.85418782e-12;

% materials' properties
n1  = indexRead(lambda,mat1);
n2  = indexRead(lambda,mat2);
n3  = indexRead(lambda,mat3);

% illumination properties
Rlaser = 0.61*lambda/NA;        % laser radius
I0 = Plaser./(2*pi*Rlaser.^2);  % laser irradiance at the center of the Gaussian shape
E0 = sqrt(2*I0./(n1*c0*eps0));  % electric field amplitude at the center of the Gaussian shape

% skinDepth=lambda./(4*pi*imag(n2));

Nl = length(lambda);
Nh = length(Delta);

t12 =  2*n1  ./(n1+n2);
r12 = (n1-n2)./(n1+n2);
r23 = (n2-n3)./(n2+n3);

%% import the 1D Green's functions G{zs}(z,r) calculated with Mathematica
%% and build the 2D Green's functions

rList=dlmread('functions/GreensFunctions20_GlassGoldWater/rList.txt');

Nr0 = length(rList);
Nz0 = 20; % Number of Green's function along z
GreenProfiles0 = cell(Nz0,1);

%for izs=1:Nz0
%    GreenProfiles0{izs}=zeros(Nz0,Nr0);
%    for iz=1:Nz0
%        GreenProfiles0{izs}(iz,:)=dlmread(['functions/GreensFunctions20_GlassGoldWater/G' num2str(izs) 't' num2str(iz) '.txt'])';
%    end
%end

for izs = 1:Nz0
    GreenProfiles0{izs} = zeros(Nz0,Nr0);
    for iz = 1:Nz0
        filename = ['functions/GreensFunctions20_GlassGoldWater/G' num2str(izs) 't' num2str(iz) '.txt'];
        fileID = fopen(filename);
        C = textscan(fileID,'%s');
        fclose(fileID);
        Nc = numel(C{1});
        for ic = 1:Nc
            if contains(C{1}{ic},'*') % values that contain '*^-6' in the exported Mathematica file that Matlab cannot read
                pos = find(C{1}{ic}=='*');
                num = str2double(C{1}{ic}(1:pos-1));
                pos2 = find(C{1}{ic}=='-');
                pow = str2double(C{1}{ic}(pos2+1:end));
                GreenProfiles0{izs}(iz,ic) = num*10^-pow;
            else
                GreenProfiles0{izs}(iz,ic) = str2double(C{1}{ic});
            end
        end  
    end
end



%% loop over lambda and Delta
Tmax     = zeros(Nl,Nh);
HeatmapXY= cell(Nl,Nh);
Qtotal   = zeros(Nl,Nh);
Ptotal   = zeros(Nl,Nh);
Abs      = zeros(Nl,Nh);

for ih = 1:Nh
    % Mathematica non-uniform mesh grid, dimensionment of length scales
    rListm = rList*Delta(ih);
    Zm     = ((1:Nz0)-0.5)/Nz0*Delta(ih);
    dz     = 1/Nz0*Delta(ih);
    % Definition of a uniform mesh grid
    Rm2 = (0:2*Nr-1)/(2*Nr-1)*2*Rmax; % factor of 2 for the further convolution
    [Rgrid,Zgrid] = meshgrid(Rm2,Zm);

    %% defining the coordinates vectors

    Rm = (0:Nr-1)/(Nr-1)*Rmax;    % radial coordinates [0:Rmax]
    dr = 1/(Nr-1)*Rmax;
    Rms= [flip(Rm(:,2:end),2),Rm];    % radial coordinates [-Rmax:Rmax]
    [Xgrid3D,Ygrid3D,Zgrid3D] = meshgrid(Rms,Rms,Zm); % 3D matrix of xyz coordinates
    Rgrid3D2 = Xgrid3D.^2+Ygrid3D.^2;
    Xgrid2D  = (-(2*Nr-1):2*Nr-1)/(2*Nr-1)*2*Rmax; % factor of 2 for the subsequent convolution
    [X2D,Y2D]= meshgrid(Xgrid2D,Xgrid2D);
    RR = sqrt(X2D.^2+Y2D.^2);

    %% Calculation of the x-y Green's function G(zs,0,0;z,x,y) for convolution

    GreenXY = cell(Nz0,Nz0);

    for iz = 1:Nz0
        for izs = 1:Nz0
            Gfunc = GreenProfiles0{izs}(iz,:);
            GreenXY{izs}{iz} = interp1(rListm/Delta(ih),Gfunc/Delta(ih),RR/Delta(ih));
        end
    end

    for il = 1:Nl
        k0 = 2*pi./lambda(il);
        omega = k0*c0;

        %% defining the heat source density
        % 3F fields of Ez and q
        Ez = E0(il)*t12(il)*(exp(1i*k0*n2(il)*Zgrid3D)+r23(il)*exp(1i*k0*n2(il)*(2*Delta(ih)-Zgrid3D)))./(1+r12(il)*r23(il)*exp(2*1i*k0*n2(il)*Zgrid3D));
        dV = dr*dr*dz;
        Qtotal(il,ih) = omega*eps0*imag(n2(il)^2)*sum(abs(Ez(:)).^2)/2*dV; % total power integrated/absorbed within the layer
        Ptotal(il,ih) = I0(il)*(2*Rmax)^2; %total incident power
        Abs(il,ih) = Qtotal(il,ih)/Ptotal(il,ih); % Absorbance of the layer, can be compared with the result of Fig2.m
        % Gaussian heating:
        Heat0 = exp(-Rgrid3D2./(2*Rlaser(il).^2)).*(sqrt(Rgrid3D2)<=2.5*Rlaser(il)).*omega*eps0.*imag(n2(il)^2).*abs(Ez).^2/2;
        P0 = sum(Heat0(:))*dV;
        
        % Uniform heating in x-y-z
        %Heat=Rgrid3D2<=Rlaser(il)^2;
        
        Heat = P0*Heat0/(sum(Heat0(:)*dV));
        
        disp([num2str(ih) '/' num2str(Nh)])
        
        figure,plot(Zm,permute(Heat(round(end/2),round(end/2),:),[3,1,2]))%,ylim([0,1])
        title 'heat source density z-profile'
        xlabel 'z [nm]'
        ylabel 'heat source (arb units)'
        pause(0.5)

        HeatmapXY{il,ih} = Heat(:,:,1);

        
        %% Convolution, plane by plane
        count = 0;    
        Tmap = Heat*0;
        for iz = 1:Nz0
            TmapBig = 0;
            for izs = 1:Nz0
                TmapBig = TmapBig+conv2(Heat(:,:,izs),GreenXY{izs}{iz});
            end
            Tmap(:,:,iz) = crop(TmapBig,2*Nr-1,2*Nr-1)*dV;
        end
        Tmax(il,ih) = max(Tmap(:));

    end %Delta
end %lambda
        
%% plot heat and temperature maps in XY

close all
figure,
subplot(1,2,1)
imagesc(1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),HeatmapXY{Nl,1}) % in mW/µm^3
title('heat source density at z=0')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'colormap',parula(1024))
hc1 = colorbar;
ylabel(hc1,'W/m^3')
xlabel('µm')
ylabel('µm')
subplot(1,2,2)
imagesc(1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),Tmap(:,:,round(end/2))) % in K
title('temperature distribution at z=0')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'colormap',hot(1024))
hc2 = colorbar;
ylabel(hc2,'K')
xlabel('µm')
ylabel('µm')

%% Plot heat and temperature maps in XZ
hfig = figure;
Pos = get(0, 'Screensize');
hfig.Position(1) = Pos(1);
hfig.Position(3) = Pos(3);

subplot(2,1,1)
[~,c] = contourf(1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),1e9*Zm,permute(Heat(:,round(end/2),:),[3,1,2]),5);
c.LineWidth = 1;
hc2 = colorbar;
ylabel(hc2,'W/m^3')
set(gca,'colormap',parula(1024))
caxis([0 max(Heat(:))])
title('Heat source density XZ map')
xlabel('µm')
ylabel('nm')

subplot(2,1,2)
[~,c] = contourf(1e6*(Rm2(1:end-1)-mean(Rm2(1:end-1))),1e9*Zm,permute(Tmap(:,round(end/2),:),[3,1,2]),20);
c.LineWidth = 1;
hc1 = colorbar;
ylabel(hc1,'°C')
set(gca,'colormap',hot(1024))
caxis([0 max(Tmap(:))])
title('Temperature increase XZ map')
xlabel('µm')
ylabel('nm')

%%
figure
plot(1e6*Rm2(1:end-1),Tmap(:,round(end/2),round(end/2)))
xlabel('µm')
ylabel('K')
title('Temperature increase profile at z=h/2')
set(gca,'YLim',[0 ceil(max(Tmap(:,round(end/2),round(end/2))))])
