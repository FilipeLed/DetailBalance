%% DATA EXTRACTION FROM TABLES
%% inputs
function [data,wl_water,k_water,n_water,wl_solar,global_values,wl_glass_n,n_glass,wl_glass_k,k_glass,wl_silver,n_silver,k_silver,wl_gold,n_gold,k_gold,wl_copper,n_copper,k_copper,wl_aluminun,n_aluminun,k_aluminun,wl_EQE,EQE_si]=DataExtraction()
data = xlsread('OpticalData.xlsx','Data','A3:AF2004');
wl_water = data(1:119,1).*1000;     %nm
k_water = data(1:119,2);            %water extinction coefficient
n_water = data(1:119,3);            %water refractive index 
wl_solar = data(1:2002,7);          %nm
global_values = data(1:2002,10);    %global radiation [W/m²/um]
wl_glass_n = data(1:91,12).*1000;   %um
n_glass = data(1:91,13);            %glass refractive index
wl_glass_k = data(1:106,14).*1000;  %nm
k_glass = data(1:106,15);           %glass extinction coefficent
wl_silver = data(1:519,4)*1000;     %nm
n_silver = data(1:519,5);           %silver refractive index
k_silver = data(1:519,6);           %silver extinction
wl_gold = data(1:69,22)*1000;       %nm
n_gold = data(1:69,23);             %gold refractive index
k_gold = data(1:69,24);             %gold extinction
wl_copper = data(1:74,25)*1000;     %nm
n_copper = data(1:74,26);           %copper refractive index
k_copper = data(1:74,27);           %copper extinction
wl_aluminun = data(1:339,30)*1000;  %nm
n_aluminun = data(1:339,31);        %aluminum refractive index
k_aluminun = data(1:339,32);        %aluminum extinction
wl_EQE = data(1:668,28);            %nm
EQE_si = data(1:668,29);            %Silico solar cell external quantum efficiency
end