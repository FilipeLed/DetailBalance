%% CONSTANTS AND PROPERTIES DECLARATION
%% inputs
function [h,heV,c,kb,q,SB,Tamb,Tsky,Cp_w,Dvis_w,p_w,k_w,e_glass,h_g_air,n1,n2,n3,Alfa,beta,K,Egap0,e_pv]=PropertiesDeclaration()

%GLOBAL CONSTANTS DECLARATION
h = 6.62607004e-34;         %Plack's constant in J.s
heV = 4.135667662e-15;      %Plack's constant in eV.s
eV = 1.6021766208e-19;      %Units are J
c = 299792458.0;            %speed of light m/s
kb = 1.38064852e-23;        %Boltzmann constant in J/K
q = 1.6021766208e-19;       %Elementary charge in Coulomb
SB = 5.670367e-8;           %Stefan–Boltzmann constant W/(m^2.K^4)

%SYSTEM PROPERTIES DECLARATION

%AIR AND AMBIENT
Tamb = 298.15;                      %Temperature of enviroment in K
Vair = 1;                           %ambiente air velocity in m/s
Lc = 1;                             %characteristic length [m]
Cp_air = 1.00701722725e+3;          %air heat capacity J/kg.K at room temperature
Dvis_air = 18.47269269875e-6;       %air dynamic viscosity Pa.s
p_air = 1.18418609175;              %air density kg/m³
Kvis_air = 1.5599484597434e-5;      %air Kinematic viscosity m²/s
Pr_air = 0.7141254659;              %air Prandtl-Number
k_air = 0.0257;                     %air Thermal conductivity W/m.K
Tsky = ((Tamb^4)*(1-0.261*exp(-7.77e-4*(Tamb-273)^2)))^0.25;    %Temperature of the sky Daguenet's formula
Re_air = (Vair*Lc)/Kvis_air;        %air Reynold number
Nu_air = 0.664*(Re_air^0.5)*(Pr_air^(1/3));%air Nusselt number by Incropera 6th page 410

%WATER
Cp_w = 4.1784e+3;               %water heat capacity J/kg.K at room temperature
Dvis_w = 0.8007e-3;             %water dynamic viscosity Pa.s
p_w = 995.7;                    %water density kg/m³
Kvis_w = 0.801e-6;              %water Kinematic viscosity m²/s
k_w = 0.591;                    %water Thermal conductivity W/m.K

%GLASS
k_g = 0.9;                      %glass Heat conductance W/m.K
e_glass = 0.9;                  %glass emissivity
h_g_air = (Nu_air*k_air)/Lc;    %glass-air/pv-air convective heat transfer coefficient W/m².K

%PV cell
n1 = 1;
n2 =1;
n3 = 1;
Alfa = 7.021*1e-04;
beta = 1108;
K= 500;
Egap0 = 1.1557;                 %bandgap at 0K
e_pv =0.9;                      %pv cell emissivity 
end
