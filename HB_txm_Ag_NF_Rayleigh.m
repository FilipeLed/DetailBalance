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
fv = 0.0001;                                        %volume fraction vector[%]
D = 10;                                             %nanoparticle size vector[nm]

%% SILVER PROPERTIES DECLARATION
Cp_p = 232;                    %silver heat capacity J/kg.K at room temperature
p_p = 10490;                   %silver density kg/m³
k_p = 429;                     %silver Thermal conductivity W/m.K

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
np = fit(wl_silver,n_silver,'linearinterp');     %metal refractive index interpolation
kp = fit(wl_silver,k_silver,'linearinterp');     %metal extinction coefficient interpolation 

%% 5º STEP
%NANOFLUID OPTICAL CHARACTERIZATION

[Qext,Qsca,Qabs,m]=RayleighScattering(wl_solar,D,subplus(nw(wl_solar)),complex(np(wl_solar),kp(wl_solar)));


%% 6º STEP
%NANOFLUID FINAL PROPERTIES
n_nf = subplus(nw(wl_solar)) + ((1.5*fv.*subplus(nw(wl_solar)))).*((real(m).^2-1)./(real(m).^2+2));   %nanofluid refractive index
w_abs =(subplus(kw(wl_solar)).*pi*4)./(wl_solar.*10^-9);     %Water absorption coefficient
ke = 3*fv*Qext./(2*D*10^-9);                         %Extinction coefficient of the nanoparticle
ket = ke+w_abs;                                      %Total extinction coefficient (of the nanofluid)

%% OUTER LOOP (channel thickness dependence)

t_nf = (0.1:0.1:20)*1e-3;                           %nanofluid thickness vector [m]
m_flow = 0.00005:0.00005:0.02;                      %nanofluid mass flow rate Kg/s vector


for z = 1:length(t_nf)
    
%% INNER LOOP (mass flow rate dependence)    
    for n = 1:length(m_flow)
        
        %% 7° STEP
%THERMO-PHYSICAL PROPERTIES OF NANOFLUIDS

[A_nf,p_nf,k_nf,Cp_nf,h_g_nf,u_nf]=ThermoPhysical(fv,p_p,p_w,Dvis_w,k_w,k_p,Tamb,D,Cp_p,Cp_w,m_flow(n),t_nf(z),W_nf,L_nf);

%% 8º STEP
%THE OPTICAL INDEXES OF PV/T SYSTEM
[g_transm,g1_h_abs,rad_in,rad_nf,nf_h_abs,g2_h_abs,total_abs,rad_out,mean_t_nf,r_g1_g1,r_g_nf,r_g_pv,reflect1,nf_transm,reflect0]= OpticalIndexes(radiation,wl_solar,L_glass,e_glass,SB,e_pv,ng(wl_solar),n_nf,ket,t_nf(z),kg(wl_solar));

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
         (nf_h_abs - h_g_nf*(2*x(5)-x(4)-x(6)) - r_g_nf*(2*x(5)^4 - x(4)^4- x(6)^4)-(m_flow(n)*Cp_nf/(A_nf))*(x(5) - Tin))
         (g2_h_abs + h_g_nf*(x(5) - x(6)) + r_g_nf*(x(5)^4 - x(6)^4) - mean_t_nf*r_g1_g1*(x(6)^4 - x(4)^4)+ r_g_pv*(x(3)^4 - x(6)^4)- h_g_air*(x(6) - Tamb)-SB*e_glass*(x(6)^4 - Tsky^4))
         (pv_h_abs - x(1)*x(2) - r_g_pv*(x(3)^4 - x(6)^4)- SB*e_glass*(x(3)^4 - Tsky^4)- 2*h_g_air*(x(3) - Tamb))];
         %(total_abs + r_g_pv*(x(3)^4 - x(6)^4) - h_g_air*(x(4) +x(6) - 2*Tamb)- SB*e_glass*(x(4)^4 +x(6)^4 - 2*Tsky^4)-(m_flow*Cp_nf/(A_nf))*(x(7) - Tin))];
  x = fsolve(F,[0;0;Tamb;Tamb;Tamb;Tamb]);
  clc
  P_el(z,n) = x(1)*x(2);
  P_th(z,n) = ((m_flow(n)*Cp_nf/(A_nf))*(x(5)-Tin));
  ef_el(z,n) = x(1)*x(2)/trapz(wl_solar,radiation);
  ef_th(z,n) = ((m_flow(n)*Cp_nf/(A_nf))*(x(5)-Tin))/trapz(wl_solar,radiation);
  Tout(z,n) = x(5);
  Jmpp(z,n) = x(1);   %current of maximum power point
  Vmpp(z,n) = x(2);   %voltage of maximum power point
  Tpv(z,n) = x(3);
  
    end
end

%
%% saving data in excel

exer_tm_th = ef_el'+(1 - 298./Tout').*ef_th';
[Y,I]= max(exer_tm_th);
flow_max = m_flow(I);

True = 1;
False = 0;
xlCenter = -4108;
Row0 = {'Silver mxt P_th',NaN,NaN,NaN,NaN,NaN,NaN};
Row1 = {'m/t'};
Row2 = {'Silver mxt P_el',NaN,NaN,NaN,NaN};
Row3 = {'Silver mxt Tout',NaN,NaN,NaN,NaN};
Row4 = {'Silver mxt ef_el',NaN,NaN,NaN,NaN};
Row5 = {'Silver mxt ef_th',NaN,NaN,NaN,NaN};
%directoryname = uigetdir('', 'Escolha uma pasta para salvar o arquivo');
file_name = fullfile(directoryname,'Artigo_Rayleigh.xlsx');
xlswrite(file_name,Row0,'Ag mxt Pth','A1');
xlswrite(file_name,Row1,'Ag mxt Pth','A2');
xlswrite(file_name,m_flow','Ag mxt Pth','A3');
xlswrite(file_name,t_nf,'Ag mxt Pth','B2');
xlswrite(file_name,P_th','Ag mxt Pth','B3');

xlswrite(file_name,Row2,'Ag mxt Pel','A1');
xlswrite(file_name,Row1,'Ag mxt Pel','A2');
xlswrite(file_name,m_flow','Ag mxt Pel','A3');
xlswrite(file_name,t_nf,'Ag mxt Pel','B2');
xlswrite(file_name,P_el','Ag mxt Pel','B3');

xlswrite(file_name,Row3,'Ag mxt Tout','A1');
xlswrite(file_name,Row1,'Ag mxt Tout','A2');
xlswrite(file_name,m_flow','Ag mxt Tout','A3');
xlswrite(file_name,t_nf,'Ag mxt Tout','B2');
xlswrite(file_name,Tout','Ag mxt Tout','B3');

xlswrite(file_name,Row4,'Ag mxt efel','A1');
xlswrite(file_name,Row1,'Ag mxt efel','A2');
xlswrite(file_name,m_flow','Ag mxt efel','A3');
xlswrite(file_name,t_nf,'Ag mxt efel','B2');
xlswrite(file_name,ef_el','Ag mxt efel','B3');

xlswrite(file_name,Row5,'Ag mxt efth','A1');
xlswrite(file_name,Row1,'Ag mxt efth','A2');
xlswrite(file_name,m_flow','Ag mxt efth','A3');
xlswrite(file_name,t_nf,'Ag mxt efth','B2');
xlswrite(file_name,ef_th','Ag mxt efth','B3');

Excel = actxserver('Excel.Application');
Workbooks = Excel.Workbooks;
Excel.Visible = 0; % if you want to see the excel file real time enter = 1;
Workbook = Excel.Workbooks.Open(file_name);

Range = Excel.Range('A1:X1');
Range.Select;
Range.MergeCells = True;
Range.VerticalAlignment = xlCenter;
Range.HorizontalAlignment = xlCenter;
Range.ColumnWidth = 9;

Workbook.Save;
Excel.Quit;
system('taskkill /F /IM EXCEL.EXE');
%
True = 1;
False = 0;
xlCenter = -4108;
Row6 = {'Ag max m_flow',NaN,NaN,NaN,NaN,NaN,NaN};
Row7 = {'t'};
Row8 = {'mass flows'};
file_name = fullfile(directoryname,'OpticalData.xlsx');
xlswrite(file_name,Row6,'Rayleigh','A1');
xlswrite(file_name,Row7,'Rayleigh','A2');
xlswrite(file_name,Row8,'Rayleigh','B2');
xlswrite(file_name,t_nf','Rayleigh','A3');
xlswrite(file_name,flow_max','Rayleigh','B3');

Excel = actxserver('Excel.Application');
Workbooks = Excel.Workbooks;
Excel.Visible = 0; % if you want to see the excel file real time enter = 1;
Workbook = Excel.Workbooks.Open(file_name);

Range = Excel.Range('A1:B1');
Range.Select;
Range.MergeCells = True;
Range.VerticalAlignment = xlCenter;
Range.HorizontalAlignment = xlCenter;
Range.ColumnWidth = 9;

Workbook.Save;
Excel.Quit;

%}