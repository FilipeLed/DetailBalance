%% THE THERMO-PHYSICAL PROPERTIES OF NANOFLUIDS
%% inputs
%fv volume fraction vector[%]
%p_p bulk metal density kg/m³
%p_w water density kg/m³
%Dvis_w water dynamic viscosity Pa.s
%k_w water Thermal conductivity W/m.K
%k_p bulk metal thermal conductivity W/m.K
%Tamb Temperature of enviroment in K
%D nanoparticle diameter in nm
%Cp_p bulk metal heat capacity J/kg.K at room temperature
%Cp_w water heat capacity J/kg.K at room temperature
%m_flow nanofluid mass flow rate Kg/s
%t_nf nanofluid channel thickness in m
%W_nf nanofluid width in m
%L_nf nanofluid channel length in m

function [A_nf,p_nf,k_nf,Cp_nf,h_g_nf,u_nf]=ThermoPhysical(fv,p_p,p_w,Dvis_w,k_w,k_p,Tamb,D,Cp_p,Cp_w,m_flow,t_nf,W_nf,L_nf)

A_nf = W_nf*L_nf;                                   %nanofluid area [m²]
D_h = (4*t_nf*W_nf)/(2*t_nf+2*W_nf);                %hydraulic diameter [m] Dh = 4*Area/Perimetro
p_nf = fv*p_p + (1-fv)*p_w;                         %nanofluid density [Kg/m³] by Pak and Cho
Dvis_nf = (1+2.5*fv)*Dvis_w;                        %nanofluid dynamic viscosity Pa.s by Einstein
k_nf = ((k_p+2*k_w-2*fv*(k_w-k_p))/(k_p+2*k_w+fv*(k_w-k_p)))*k_w; %nanofluid Thermal conductivity W/m.K by Maxwell
%k_nf = k_w*(1+(0.135*(k_w/k_p)^0.273)*(fv^0.467)*((Tamb/20)^0.547)*(100/D)^0.234);%nanofluid Thermal conductivity W/m.K by Patel
Cp_nf = (fv*(p_p*Cp_p)+(1-fv)*(p_w*Cp_w))/p_nf;     %nanofluid heat capacity J/kg.K by Xuan and Roetzel
Pr_nf = (Dvis_nf*Cp_nf)/k_nf;                       %nanofluid Prandtl-Number
u_nf = m_flow/(t_nf*W_nf*p_nf);                     %nanofluid velocity m/s
Re_nf = (p_nf*D_h*u_nf)/Dvis_nf;                    %nanofluid Reynolds number
Nu_nf = 0.023*(Re_nf^0.8)*(Pr_nf^0.4);              %nanofluid Nusselt number by  Dittos-Boelter
h_g_nf = (Nu_nf*k_nf)/D_h;                          %glass-nanofluid convective heat transfer coefficient W/m².K 

end
