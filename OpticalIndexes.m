%% THE OPTICAL INDEXES OF NANOFLUIDS PV/T SYSTEM
%% inputs

function [g_transm,g1_h_abs,rad_in,rad_nf,nf_h_abs,g2_h_abs,total_abs,rad_out,mean_t_nf,r_g1_g1,r_g_nf,r_g_pv,reflect1,nf_transm,reflect0]= OpticalIndexes(radiation,wl_solar,L_glass,e_glass,SB,e_pv,ng,n_nf,ket,t_nf,kg)

reflect0 = ((ng-1)./(ng+1)).^2;                                     %Air-Glass reflectance
Ka_glass = 4*pi*subplus(kg)./(wl_solar.*1e-9);                      %Glass absorption coefficient
g_transm = exp(((Ka_glass).*-1).*(L_glass));                        %Glass Transmittance
reflect1 = ((subplus(ng)-n_nf)./(subplus(ng)+n_nf)).^2;             %Nanofluid-Glass reflectance
nf_transm = exp((subplus(ket).*-1).*(t_nf));                        %Nanofluido Transmittance
g1_h_abs = trapz(wl_solar,(1-g_transm).*radiation);                 %Heat absorbed by first glass
rad_in = radiation.*g_transm;                                       %Radiation after the first glass cover
rad_nf = rad_in.*nf_transm;                                         %Radiation after the nanofluid
nf_h_abs = trapz(wl_solar,(1-nf_transm).*rad_in);                   %Heat absorbed by nanofluid
g2_h_abs = trapz(wl_solar,(1-g_transm).*rad_nf);                    %Heat absorbed by second glass
total_abs = trapz(wl_solar,(1-nf_transm.*g_transm.^2).*radiation);  %Total heat absorbed
rad_out = rad_nf.*g_transm;                                                     %Radiation after the second glass cover
mean_t_nf = trapz(wl_solar,(nf_transm).*radiation)/trapz(wl_solar,radiation);   %Transmissivity of nanofluid
r_g1_g1 = SB/((1/e_glass)+(1/e_glass)-1);                                       %thermal radioativity resistance for 1º and 2º glass
r_g_nf = (e_glass*SB*(1-mean_t_nf)*(1+(1-e_glass)*mean_t_nf))/((2*e_glass-e_glass^2)*mean_t_nf^2);
r_g_pv = (SB*e_pv*e_glass)/(e_glass+e_pv-e_pv*e_glass);
end