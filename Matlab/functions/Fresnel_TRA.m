%% Compute the Fresnel coefficient in Reflection, Transmission and Absorption in Intensity for a 5-layer system

%Benoit Rogez, Guillaume Baffou
%CNRS, Institut Fresnel 
%March 2019
%%
function [T,R,A] = Fresnel_TRA(lambda0,Pol,theta0,n,hList)
%% parameters
%Inputs
 %lambda0: wavelength in vacuum [nm]
 %Pol:   polarisation of the beam: TE or TM
 %theta0:angle of incidence of the beam [degree]
 %n:     5-vector of the complex refractive indices of the 5 media.
 %hList: 5-cell of the thicknesses of the 5 media [nm] (media 1 and 5 are assumed semi-infinite)
%Outputs
 %T:     transmission of the multilayer
 %R:     reflection of the multilayer
 %A:     absorption of the multilayer 

%The beam is incident on medium 1

%%
theta = deg2rad(theta0);

Nh = numel(n);

h = zeros(5,1);
for il = 2:Nh-1  % loop over the intermediate layers
    h(il) = hList{il};
end  

%Computation of the different wavevectors
k0 = 2*pi/lambda0;
kx = k0*n(1)*sin(theta);               %In-plane component. Same for all media    
kz2 = sqrt(k0*k0*n(2)*n(2)-kx*kx);     %Out-of-plane component
kz3 = sqrt(k0*k0*n(3)*n(3)-kx*kx);
kz4 = sqrt(k0*k0*n(4)*n(4)-kx*kx);

%% Complex Fresnel coefficients for each interface

if strcmp(Pol,'TE') % For TE wave

    %Interface 1-2
    Th(2)= asin(n(1)*sin(theta)/n(2));                                          %Refraction angle
    r12 = (n(1)*cos(theta)-n(2)*cos(Th(2)))/(n(1)*cos(theta)+n(2)*cos(Th(2)));  %Reflection
    t12 = (2*n(1)*cos(theta))/(n(1)*cos(theta)+n(2)*cos(Th(2)));                %Transmission

    %Interface 2-3
    Th(3)= asin(n(1)*sin(theta)/n(3));
    r23 = (n(2)*cos(Th(2))-n(3)*cos(Th(3)))/(n(2)*cos(Th(2))+n(3)*cos(Th(3)));
    t23 = (2*n(2)*cos(Th(2)))/(n(2)*cos(Th(2))+n(3)*cos(Th(3)));

    %Interface 3-4
    Th(4)=asin(n(1)*sin(theta)/n(4));
    r34 = (n(3)*cos(Th(3))-n(4)*cos(Th(4)))/(n(3)*cos(Th(3))+n(4)*cos(Th(4)));
    t34 = (2*n(3)*cos(Th(3)))/(n(3)*cos(Th(3))+n(4)*cos(Th(4)));

    %Interface 4-5
    Th(5)= asin(n(1)*sin(theta)/n(5));
    r45 = (n(4)*cos(Th(4))-n(5)*cos(Th(5)))/(n(4)*cos(Th(4))+n(5)*cos(Th(5)));
    t45 = (2*n(4)*cos(Th(4)))/(n(4)*cos(Th(4))+n(5)*cos(Th(5)));

elseif strcmp(Pol,'TM') % For a TM wave

    %Interface 1-2
    Th(2)= asin(n(1)*sin(theta)/n(2));
    r12 = (n(2)*cos(theta)-n(1)*cos(Th(2)))/(n(2)*cos(theta)+n(1)*cos(Th(2)));
    t12 = (2*n(1)*cos(theta))/(n(2)*cos(theta)+n(1)*cos(Th(2)));

    %Interface 2-3
    Th(3)= asin(n(1)*sin(theta)/n(3));
    r23 = (n(3)*cos(Th(2))-n(2)*cos(Th(3)))/(n(3)*cos(Th(2))+n(2)*cos(Th(3)));
    t23 = (2*n(2)*cos(Th(2)))/(n(3)*cos(Th(2))+n(2)*cos(Th(3)));

    %Interface 3-4
    Th(4)= asin(n(1)*sin(theta)/n(4));
    r34 = (n(4)*cos(Th(3))-n(3)*cos(Th(4)))/(n(4)*cos(Th(3))+n(3)*cos(Th(4)));
    t34 = (2*n(3)*cos(Th(3)))/(n(4)*cos(Th(3))+n(3)*cos(Th(4)));

    %Interface 4-5
    Th(5)= asin(n(1)*sin(theta)/n(5));
    r45 = (n(5)*cos(Th(4))-n(4)*cos(Th(5)))/(n(5)*cos(Th(4))+n(4)*cos(Th(5)));     
    t45 = (2*n(4)*cos(Th(4)))/(n(5)*cos(Th(4))+n(4)*cos(Th(5)));

else
    error('Polarization must be ''TE'' or ''TM''')
end

%% Complex Fresnel coefficients for the whole multilayer
% System 3-5
r35 = (r34+r45*exp(1i*2*h(4)*kz4))/(1+r34*r45*exp(1i*2*h(4)*kz4));
t35 = (t34*t45*exp(1i  *h(4)*kz4))/(1+r34*r45*exp(1i*2*h(4)*kz4));

% System 2-5
r25 = (r23+r35*exp(1i*2*h(3)*kz3))/(1+r23*r35*exp(1i*2*h(3)*kz3));
t25 = (t23*t35*exp(1i  *h(3)*kz3))/(1+r23*r35*exp(1i*2*h(3)*kz3));

% System 1-5
r15 = (r12+r25*exp(1i*2*h(2)*kz2))/(1+r12*r25*exp(1i*2*h(2)*kz2));
t15 = (t12*t25*exp(1i  *h(2)*kz2))/(1+r12*r25*exp(1i*2*h(2)*kz2));

%% Fresnel coefficients for the whole multilayer

%Transmission
T = abs(t15)*abs(t15)*n(5)*real(cos(Th(5)))/(n(1)*cos(theta));
%Reflection
R = abs(r15)*abs(r15);
%Absorption
A = 1-T-R;

end
