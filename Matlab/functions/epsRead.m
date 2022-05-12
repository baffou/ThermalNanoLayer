 % function that returns the complex value of the dielectric constant of a
% material as a function of the wavelength lambda.
% Lambda can be written in [nm] or [m].
%
% fileName is the name of the file that contains a set of permittivies
% values (real and imag values) along with the associated energies.

%Authors: Guillaume Baffou
%affiliation: CNRS, Institut Fresnel 
%Date: Jan 2021


%%
function EPS=epsRead(lambda,fileName)

Nl = length(lambda);
EPS = lambda*0;

for il = 1:Nl
    lambda0 = lambda(il);

    if lambda0<1
        lambda0 = lambda0*1e9;  %restablish lambda in nm
    end
    nMetal = dlmread([fileName '.txt']);

    factor = 1239.8419;  % 1eV -> 1239.8419 (nm)
    eV0 = factor/lambda0;
    if eV0>max(nMetal(:,1)) || eV0<min(nMetal(:,1))
        error(['specified wavelength out of range:' num2str(ceil(factor/max(nMetal(:,1)))) ' - ' num2str(floor(factor/min(nMetal(:,1)))) ' nm'])
    end

    eps1 = nMetal(:,2).^2-nMetal(:,3).^2;
    eps2 = 2.0*nMetal(:,2).*nMetal(:,3);


    eps10 = interp1(nMetal(:,1),eps1,eV0,'spline');
    eps20 = interp1(nMetal(:,1),eps2,eV0,'spline');
    EPS(il) = eps10+1i*eps20;

end
      
      
      
      
