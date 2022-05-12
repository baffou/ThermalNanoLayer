% function that returns the complex value of the dielectric constant of a
% material as a function of the wavelength lambda0.
%
% fileName is the name of the file that contains a set of permittivies
% values (real and imag values) along with the associated energies.

%Authors: Guillaume Baffou
%affiliation: CNRS, Institut Fresnel 
%Date: Jan 2021



%%
function index = indexRead(lambda,fileName)
%possible values of fileName : 'Au'
Nl = length(lambda);
index = lambda*0;

for il = 1:Nl
    
    lambda0 = lambda(il);
    
    if lambda0<1
        lambda0 = lambda0*1e9;  %restablish lambda in nm
    end

    nMetal = dlmread(['materials/' fileName '.txt']);

    factor = 1239.8419;  % 1eV -> 1239.8419 (nm)
    eV0 = factor/lambda0;
    if eV0>max(nMetal(:,1)) || eV0<min(nMetal(:,1))
        error(['specified wavelength out of range:' num2str(ceil(factor/max(nMetal(:,1)))) ' - ' num2str(floor(factor/min(nMetal(:,1)))) ' nm'])
    end

    
    
    n1 = interp1(nMetal(:,1),nMetal(:,2),eV0,'spline');
    n2 = interp1(nMetal(:,1),nMetal(:,3),eV0,'spline');

    index(il) = complex(n1,n2);

end 
      
