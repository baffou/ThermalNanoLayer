%% Computes the Transmission, Reflection and Absorption of a multilayer
%% Multilayer system (up to 5 layers) as a function of the metal layer thickness and incidence angle (Figure 2)
%% Metal layer is layer #3.

% Benoit Rogez, Guillaume Baffou
% CNRS, Institut Fresnel 
% Jan 2021

% Used to make Figure 2 of ref https://doi.org/10.1063/5.0057185

% GitHub: https://github.com/baffou/ThermoNanoLayer
% Please cite https://doi.org/10.1063/5.0057185 if you use this code in a publication. Thx.

clear
addpath(genpath('functions'))

%% Illumination

lambda0 = 800;      % wavelength [nm]
Pol     = 'TM';     % polarization, 'TM' or 'TE'
theta   = 0:0.5:90; % range of incidence angles [deg]

%% Definition of the layer structure

h3=0:1:200;         % range of metal layer thicknesses [nm]

%matList = {'BK7' 'Ti' 'Au' 'H2O' 'H2O'};
%matList = {'BK7' 'Cr' 'Au' 'Air' 'Air'};
%matList = {'Air' 'Air' 'Au' 'Air' 'Air'};
%matList = {'BK7' 'Cr' 'Au' 'Water' 'Water'};
matList = {'BK7' 'Cr' 'Au' 'Water' 'Water'};    %multilayer composition (see the list of available materials in the 'materials' folder)

hList = {'Inf'   0   h3     0  'Inf'};    %layers' thicknesses. 1st and 5th are necessarily 'Inf'.


%% Calculations
Nh = length(hList{3}); Nt=length(theta);

Reflection   = zeros(Nh,Nt);
Transmission = zeros(Nh,Nt);
Absorption   = zeros(Nh,Nt);

n = zeros(5,1);
for imat = 1:5
    n(imat) = indexRead(lambda0,matList{imat});
end

for ih = 1:Nh              
    for it = 1:Nt
        hList0   = hList;
        hList0{3}= hList{3}(ih);
        [T,R,A]  = Fresnel_TRA(lambda0,Pol,theta(it),n,hList0);
        Reflection(ih,it)  = R;
        Transmission(ih,it)= T;
        Absorption(ih,it)  = A;
    end
end

%% Plotting of the results
hfig = figure('Position',[0 0 1024 400]);
Pos = get(0, 'Screensize');
hfig.Position(1) = Pos(1);
hfig.Position(3) = Pos(3);

GT = subplot(1,3,1);
pcolor(theta,hList{3},Transmission)
shading interp
colormap(GT,jet(1024))
caxis([0 1])
title(GT,'Transmission')
xlabel(GT,'Angle of incidence [deg]')
ylabel(GT,'Thickness of medium 3 [nm]')
colorbar(GT)

GR = subplot(1,3,2);
pcolor(theta,hList{3},Reflection)
shading interp
colormap(GR,jet(1024))
caxis([0 1])
title(GR,'Reflection')
xlabel(GR,'Angle of incidence [deg]')
ylabel(GR,'Thickness of medium 3 [nm]')
colorbar(GR)

GA = subplot(1,3,3);
pcolor(theta,hList{3},Absorption)
shading interp
colormap(GA,jet(1024))
caxis([0 1])
title (GA,'Absorption')
xlabel(GA,'Angle of incidence [deg]')
ylabel(GA,'Thickness of medium 3 [nm]')
colorbar(GA)
