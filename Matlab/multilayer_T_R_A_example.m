%% Computes the Reflection, Transmission and Absorption of a multilayer
%% system (up to 5 layers) as a function of the metal layer thickness and
%% the wavelength (Figure 2)

%Authors: Benoit Rogez, Guillaume Baffou
%affiliation: CNRS, Institut Fresnel 
%Date: Jan 2021

% Used to make Figure 6def of ref https://doi.org/10.1063/5.0057185

% GitHub: https://github.com/baffou/ThermoNanoLayer
% Please cite https://doi.org/10.1063/5.0057185 if you use this code in a publication. Thx.



clear
addpath(genpath('functions'))

%% Illumination
lambda = 400:25:1200;   % wavelength [nm]
Pol = 'TM';             % Polarization, 'TM' or 'TE'
theta=0;                % incidence angle [deg]

%% Definition of the layer structure
% h in [nm]
% varying h must be the thrid layer

%Figure='d';
%Figure='e';
Figure='f';

switch Figure           % number of the figure: 'd', 'e' or 'f'
    case 'd'
        h3=5:5:100;     % middle layer thickness [nm]
        nList={'BK7' 'BK7' 'Au' 'H2O' 'H2O'};   % layers' composition, see the list in the 'material' folder
        hList={'inf'    0   h3     0  'Inf'};   % layers' thicknesses
    case 'e'
        h3=0:5:95;
        nList={'BK7' 'Cr' 'Au' 'H2O' 'H2O'};
        hList={'inf'    5   h3     0  'Inf'};
    case 'f'
        h3=1:20;
        nList={'BK7' 'BK7' 'Cr' 'H2O' 'H2O'};
        hList={'inf'    0   h3     0  'Inf'};
end


%% Calculations
Nh=length(hList{3});
Nl=length(lambda);
Reflection=zeros(Nl,Nh);
Transmission=zeros(Nl,Nh);
Absorption=zeros(Nl,Nh);
n=zeros(5,1);


% Compute the different parameters T,R,A = f(h3,lambda)
for ih=1:Nh              
    for il=1:Nl              
        for imat=1:5
            n(imat) = indexRead(lambda(il),nList{imat});
        end
        th = deg2rad(theta);
        hList0=hList;
        hList0{3}=hList{3}(ih);
        [T,R,A]=Fresnel_TRA(lambda(il),Pol,th,n,hList0);
        Transmission(il,ih)=T;
        Reflection(il,ih)=R;
        Absorption(il,ih)=A;
    end
end

%% Save results

%dlmwrite(['Fig' Figure '.txt'],Absorption)

%% Plot results

Tbig=imresize(Transmission,8)';
Rbig=imresize(Reflection,8)';
Abig=imresize(Absorption,8)';

hfig=figure;
subplot(3,1,1)
imagesc(hList{3},lambda,100*Tbig);
set(gca, 'YDir','normal')
title('Transmission (%)')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')
subplot(3,1,2)
imagesc(hList{3},lambda,100*Rbig);
set(gca, 'YDir','normal')
title('Reflection (%)')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')
subplot(3,1,3)
imagesc(hList{3},lambda,100*Abig);
set(gca, 'YDir','normal')
title('Absorption (%)')
colormap('parula(64)')
colorbar
xlabel('height (nm)')
ylabel('lambda (nm)')




