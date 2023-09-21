%% PHOTONS FLUX FROM SOLAR INCIDENT SPECTRUM
%% inputs

function [Ephoton,total_flux_graf,pv_h_abs]=PhotonFlux(heV,c,wl_solar,eqe,rad_out,h)

Ephoton = ((heV*c)./(wl_solar*1e-9));       %create a vector of photon energy in eletron-volts
Absorption_spectrum = subplus(eqe).*rad_out;
pv_h_abs = trapz(wl_solar,Absorption_spectrum);           %Energy absorbed by solar cell

photon_flux_density = (Absorption_spectrum.* wl_solar)/(h*c*1e+9);%incident foton flux as a function of wavelength (ph/m².s.nm)
for m = 1:length(Ephoton)%loop for calculate the total photon flux (ph/m².s)
    total_flux_graf(m) = trapz(wl_solar(1:m),photon_flux_density(1:m),1);%create a vector for each bandgap
end

end