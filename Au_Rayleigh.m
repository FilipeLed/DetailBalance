clc 
clear all
close all
system('taskkill /F /IM EXCEL.EXE');
directoryname ='C:\Users\Filip\Documents\MATLAB';
tic %start to count the time Elapsed time

%% 1º STEP
%GLOBAL CONSTANTS AND PROPERTIES DECLARATION

[h,heV,c,kb,q,SB,Tamb,Tsky,Cp_w,Dvis_w,p_w,k_w,e_glass,h_g_air,n1,n2,n3,Alfa,beta,K,Egap0,e_pv]=PropertiesDeclaration();

%% 2° PVT COLLECTOR CHARACTERISTICS
Cgeo = 1;                       %geometric concentration factor
L_glass = 3*10^-3;              %glass thickness [m]
W_nf = 1;                       %solar collector width [m]
L_nf = 1;                       %solar collector length [m]

%% NANOFLUID SETUP DECLARATION
m_flow = 0.0038;                                    %nanofluid mass flow rate Kg/s
t_nf = 5*10^-3;                                    %nanofluid thickness [m]
fv = 0.0001;                                        %volume fraction vector[%]
D = 10;                                             %nanoparticle size vector[nm]

%% GOLD PROPERTIES DECLARATION
Cp_p = 129;                     %Gold heat capacity J/kg.K at room temperature
p_p = 19320;                    %Gold density kg/m³
k_p = 315;                      %Gold Thermal conductivity W/m.K
       
%% 3º STEP
%DATA EXTRACTION FROM TABLES
[data,wl_water,k_water,n_water,wl_solar,global_values,wl_glass_n,n_glass,wl_glass_k,k_glass,wl_silver,n_silver,k_silver,wl_gold,n_gold,k_gold,wl_copper,n_copper,k_copper,wl_EQE,EQE_si]=DataExtraction();
radiation = global_values*Cgeo;                 %set solar spectrum concetration
solar_constant = trapz(wl_solar, radiation);    %The "Solar constant" for AM1.5 is the sun's total irradiance in W/m²

%% 4º STEP
%DATA INTERPOLATION

kw = fit(wl_water,k_water,'smoothingspline');    %water extinction coefficient interpolation   
nw = fit(wl_water,n_water,'smoothingspline');    %water refractive index interpolation      
ng = fit(wl_glass_n,n_glass,'linearinterp');     %glass refractive index interpolation 
kg = fit(wl_glass_k,k_glass,'linearinterp');     %glass extinction coefficent interpolation 
rad = fit(wl_solar,global_values,'linearinterp');%solar irradiance interpolation
eqe = fit(wl_EQE,EQE_si,'linearinterp');         %Silico solar cell external quantum efficiency interpolation
np = fit(wl_gold,n_gold,'linearinterp');         %metal refractive index interpolation
kp = fit(wl_gold,k_gold,'linearinterp');         %metal extinction coefficient interpolation

%% 5° STEP
%THERMO-PHYSICAL PROPERTIES OF NANOFLUIDS

[A_nf,p_nf,k_nf,Cp_nf,h_g_nf,u_nf]=ThermoPhysical(fv,p_p,p_w,Dvis_w,k_w,k_p,Tamb,D,Cp_p,Cp_w,m_flow,t_nf,W_nf,L_nf);


%% 6º STEP
%NANOFLUID OPTICAL CHARACTERIZATION

[Qext,Qsca,Qabs,m]=RayleighScattering(wl_solar,D,subplus(nw(wl_solar)),complex(np(wl_solar),kp(wl_solar)));

%% 7º STEP
%NANOFLUID FINAL PROPERTIES
n_nf = subplus(nw(wl_solar)) + ((1.5*fv.*subplus(nw(wl_solar)))).*((real(m).^2-1)./(real(m).^2+2));   %nanofluid refractive index
w_abs =(subplus(kw(wl_solar)).*pi*4)./(wl_solar.*10^-9);     %Water absorption coefficient
ke = 3*fv*Qext./(2*D*10^-9);                         %Extinction coefficient of the nanoparticle
ket = ke+w_abs;                                      %Total extinction coefficient (of the nanofluid)

%% 8º STEP
%THE OPTICAL INDEXES OF PV/T SYSTEM
[g_transm,g1_h_abs,rad_in,rad_nf,nf_h_abs,g2_h_abs,total_abs,rad_out,mean_t_nf,r_g1_g1,r_g_nf,r_g_pv,reflect1,nf_transm,reflect0]= OpticalIndexes(radiation,wl_solar,L_glass,e_glass,SB,e_pv,ng(wl_solar),n_nf,ket,t_nf,kg(wl_solar));

%% 9º STEP
%CONVERT DATA FROM SOLAR SPECTRUM SOURCE
[Ephoton,total_flux_graf,pv_h_abs]=PhotonFlux(heV,c,wl_solar,eqe(wl_solar),rad_out,h);
g=fit(Ephoton,total_flux_graf','linearinterp');%Interpolate to get a continuous function (ph/m².s)

%% 10º STEP
%HEAT BALANCE
Tin = Tamb;
F = @(x)[(-x(1)+q*g((Egap0 - (Alfa*x(3)^2)/(x(3) + beta)))-(K*x(3)^3/n2*(exp(-q*(Egap0 - (Alfa*x(3)^2)/(x(3) + beta))/(n3*kb*x(3)))))*(exp((q*(x(2)+x(1)*4.6*1e-08*exp(0.0207*x(3))))/(n1*kb*x(3)))-1)-(x(2)+x(1)*4.6*1e-08*exp(0.0207*x(3)))/153.92*1e-04*exp(799.93/x(3)))
         (q*g(Egap0 - (Alfa*x(3)^2)/(x(3) + beta))+(K*x(3)^3/n2*(exp(-q*(Egap0 - (Alfa*x(3)^2)/(x(3) + beta))/(n3*kb*x(3)))))-((K*x(3)^3/n2*(exp(-q*(Egap0 - (Alfa*x(3)^2)/(x(3) + beta))/(n3*kb*x(3)))))*(1+(q*x(2))/(n1*kb*x(3)))*exp((q*(x(2)+x(1)*(4.6*1e-08*exp(0.0207*x(3)))))/(n1*kb*x(3))))-(2*x(2)+x(1)*(4.6*1e-08*exp(0.0207*x(3))))/(153.92*1e-04*exp(799.93/x(3))))
         (g1_h_abs + r_g_nf*(x(5)^4 - x(4)^4)+ mean_t_nf*r_g1_g1*(x(6)^4 - x(4)^4) + h_g_nf*(x(5) - x(4)) - SB*e_glass*(x(4)^4 - Tsky^4) - h_g_air*(x(4) - Tamb))
         (nf_h_abs - h_g_nf*(2*x(5)-x(4)-x(6)) - r_g_nf*(2*x(5)^4 - x(4)^4- x(6)^4)-(m_flow*Cp_nf/(A_nf))*(x(5) - Tin))
         (g2_h_abs + h_g_nf*(x(5) - x(6)) + r_g_nf*(x(5)^4 - x(6)^4) - mean_t_nf*r_g1_g1*(x(6)^4 - x(4)^4)+ r_g_pv*(x(3)^4 - x(6)^4)- h_g_air*(x(6) - Tamb)-SB*e_glass*(x(6)^4 - Tsky^4))
         (pv_h_abs - x(1)*x(2) - r_g_pv*(x(3)^4 - x(6)^4)- SB*e_glass*(x(3)^4 - Tsky^4)- 2*h_g_air*(x(3) - Tamb))];
         %(total_abs + r_g_pv*(x(3)^4 - x(6)^4) - h_g_air*(x(4) +x(6) - 2*Tamb)- SB*e_glass*(x(4)^4 +x(6)^4 - 2*Tsky^4)-(m_flow*Cp_nf/(A_nf))*(x(7) - Tin))];
  x = fsolve(F,[0;0;Tamb;Tamb;Tamb;Tamb]);
  clc
  P_el = x(1)*x(2); %Electrical power delivered
  P_th = ((m_flow*Cp_nf/(A_nf))*(x(5)-Tin)); %Thermal power delivered
  ef_el = x(1)*x(2)/trapz(wl_solar,radiation);%Electrical efficiency 
  ef_th = ((m_flow*Cp_nf/(A_nf))*(x(5)-Tin))/trapz(wl_solar,radiation);%Thermal efficiency
  Tout = x(5);  %Outlet Temperature
  Jmpp = x(1);   %current of maximum power point
  Vmpp = x(2);   %voltage of maximum power point
  Tpv = x(3);    %PV cell temperature
  Rad_out = trapz(wl_solar,rad_out);
  
%% 11° PLOTS
figure (1)
hold on
yyaxis left
plot(wl_solar,radiation,'-g','LineWidth',1.25)
plot(wl_solar,rad_nf,':b','LineWidth',1.5)
xlabel('Wavelength (nm)')
ylabel('Spectral irradiance  (W/m²/nm)')
xlim([280 1500])
ylim([0 2])
yyaxis right
plot(wl_solar,ket,'-.r','LineWidth',1.25)
ylim([0 14000])
title('NANOPARTICLE USING RAYLEIGH SCATTERING');
xlabel('Wavelength (nm)')
ylabel('Attenuation coefficient (m-¹)')
legend('AM 1.5 Global (ASTM 1273)','transmitted spectrum','Extinction coefficient','Location','northeast','Orientation','vertical')
legend boxoff

figure (2)
plot(wl_solar,Qext*pi*(D/2)^2,'-.r','LineWidth',1.25)
xlim([280 1500])
title('Extinction Cross Sectional using Rayleigh Scattering');
xlabel('Wavelength (nm)')
ylabel('Extinction Cross Sectional Area  (nm²)')
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';

%com.mathworks.mlservices.MLCommandHistoryServices.removeAll;