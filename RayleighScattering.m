%% SCATTERING BY A SPHERICAL METAL NANOPARTICLE USING RAYLEIGH SCATTERING
%% inputs
%n_bf Index of refraction of the medium
%m_np complex refractive index of the nanoparticle
%lambda wavelength in nm
%D diemeter of the particle in nm
%fv volumetric fraction in %
function [Qext,Qsca,Qabs,m,n_nf]=RayleighScattering(lambda,D,n_bf,m_np,fv)

alfa = pi*(D*10^-9)./(lambda.*10^-9);                               %Size parameter for the nanoparticles
m = m_np./n_bf;                                                     %normalized complex refractive index of the particles
m1 = ((m.^2-1)./(m.^2+2));                                          %calculation aid
m2 = ((m.^4+27*m.^2+38)./(2*m.^2+3));                               %calculation aid
Qabs = 4*alfa.*imag(m1.*(1+((alfa.^2)./15).*m1.*m2));               %Absorption efficiency of the nanoparticle
Qsca = 8/3.*(alfa.^4).*real((m1).^2);                               %Scattering efficiency of the nanoparticle
Qext = Qabs+Qsca;                                                   %Extinction efficiency of the nanoparticle

end