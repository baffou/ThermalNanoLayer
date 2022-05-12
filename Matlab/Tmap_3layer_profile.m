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

% Used to make Figure 4bc of ref https://doi.org/10.1063/5.0057185

% GitHub: https://github.com/baffou/ThermoNanoLayer
% Please cite https://doi.org/10.1063/5.0057185 if you use this code in a publication. Thx.



clear
addpath(genpath('functions'))

%% Parameters

Delta  = 10e-9;  % metal layer thickness [m]
Rmax   = 10e-6;  % max radial coordinate for the temperature calculation [m]
Nr     = 200;    % number of calculated points along r
lambda = 800e-9;% wavelength [m]
NA     = 0.7;    % Numerical aperture of the illumination
Plaser = 1e-3;   % laser power [W]

% Materials. Must correspond to the thermal conductivities used in the Mathematica notebook to compute the Green's functions

mat1  = 'BK7';      % Bottom medium, from which light is coming
mat2  = 'Au';       % middle, absorbing layer
mat3  = 'BK7';    % top medium

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

rList=dlmread('functions/GreenFunctionFig4bc/rList.txt');
GreenProfile=rList*0;

filename=['functions/GreenFunctionFig4bc/G10t10.txt'];
fileID = fopen(filename);
C = textscan(fileID,'%s');
fclose(fileID);
Nc=numel(C{1});
for ic=1:Nc
    if contains(C{1}{ic},'*') % values that contain '*^-6' that Matlab cannot read
        pos=find(C{1}{ic}=='*');
        num=str2double(C{1}{ic}(1:pos-1));
        pos2=find(C{1}{ic}=='-');
        pow=str2double(C{1}{ic}(pos2+1:end));
        GreenProfile(ic)=num*10^-pow;
    else
        GreenProfile(ic)=str2double(C{1}{ic});
    end
end
    
    %%
    Tglass=1./(4*pi*1*rList);
    Tgold=1./(4*pi*310*rList);
    
    %%
    figure
    hold on
    plot(rList,GreenProfile)
    plot(rList,Tglass)
    plot(rList,Tgold)
    set(gca,'YLim',[0 .01])
    
    %%
    %folderOut='Fig4bcresults';
    %mkdir(folderOut)
    %dlmwrite([folderOut '/Tgold.txt'],Tgold,'\n')
    %dlmwrite([folderOut '/Tglass.txt'],Tglass,'\n')
    %dlmwrite([folderOut '/Tlayer.txt'],GreenProfile,'\n')
    
    
