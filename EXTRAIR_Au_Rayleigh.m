clc 
system('taskkill /F /IM EXCEL.EXE');
close all
clear all
%DATA EXTRACTION FROM TABLES
data01 = xlsread('Artigo_Rayleigh.xlsx','Au mxt Tout','A3:A2030');%mass flow rate range
data02 = xlsread('Artigo_Rayleigh.xlsx','Au mxt Tout','B2:AAA2');%thickness range
data03 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt Pth','A3:A2030');%volumetric fraction range
data04 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt Pth','B2:AAA2');%thickness range
data = xlsread('Artigo_Rayleigh.xlsx','FV','A3:AAA2030');
data1 = xlsread('Artigo_Rayleigh.xlsx','Diameter','A3:AAA2030');
data2 = xlsread('Artigo_Rayleigh.xlsx','thickness','A3:AAA2030');
data3 = xlsread('Artigo_Rayleigh.xlsx','mass flow','A3:AAA2030');
data4 = xlsread('Artigo_Rayleigh.xlsx','Au mxt Pth','B3:AAA2030');
data5 = xlsread('Artigo_Rayleigh.xlsx','Au mxt Pel','B3:AAA2030');
data6 = xlsread('Artigo_Rayleigh.xlsx','Au mxt Tout','B3:AAA2030');
data7 = xlsread('Artigo_Rayleigh.xlsx','Au mxt efel','B3:AAA2030');
data8 = xlsread('Artigo_Rayleigh.xlsx','Au mxt efth','B3:AAA2030');
data9 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt Pth','B3:AAA2030');
data10 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt Pel','B3:AAA2030');
data11 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt Tout','B3:AAA2030');
data12 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt efel','B3:AAA2030');
data13 = xlsread('Artigo_Rayleigh.xlsx','Au fvxt efth','B3:AAA2030');

D = data1(:,7);
Au_P_el_d = data1(:,8);
Au_P_th_d = data1(:,9);
Au_Tnf_d = data1(:,10);
Au_efel_d = data1(:,11);
Au_efth_d = data1(:,12);

m_flow = data3(:,7);
Au_P_el_m = data3(:,8);
Au_P_th_m = data3(:,9);
Au_Tnf_m = data3(:,10);
Au_efel_m = data3(:,11);
Au_efth_m = data3(:,12);

t_nf = data2(:,7);
Au_P_el_t = data2(:,8);
Au_P_th_t = data2(:,9);
Au_Tnf_t = data2(:,10);
Au_efel_t = data2(:,11);
Au_efth_t = data2(:,12);

fv = data(:,7);
Au_P_el_fv = data(:,8);
Au_P_th_fv = data(:,9);
Au_Tnf_fv = data(:,10);
Au_efel_fv = data(:,11);
Au_efth_fv = data(:,12);

exer_corr = 1/(1-(4*298)/(3*5800)+(298^4)/(3*5800^4)); %exergy correction factor

Au_m_carnot = ((1 - 298./Au_Tnf_m).*Au_efth_m).*exer_corr;

Au_fv_carnot = ((1 - 298./Au_Tnf_fv).*Au_efth_fv).*exer_corr;

Au_d_carnot = ((1 - 298./Au_Tnf_d).*Au_efth_d).*exer_corr;

Au_t_carnot = ((1 - 298./Au_Tnf_t).*Au_efth_t).*exer_corr;


%EXERGY Au
Au_mxt_exer_Th = ((1 - 298./data6).*data8).*exer_corr;
Au_mxt_exer_Tot = data7 + Au_mxt_exer_Th;
Au_fvxt_exer_Th = ((1 - 298./data11).*data13).*exer_corr;
Au_fvxt_exer_Tot = data12 + Au_fvxt_exer_Th;

[A,B] = max(Au_mxt_exer_Th);
flow_max_th = data01(B);
f0 = fit(data02',flow_max_th,'fourier4');
f00 = fit(data02',A','fourier8');

[Y,I]= max(Au_mxt_exer_Tot);
flow_max = data01(I);
f1 = fit(data02',flow_max,'fourier4');
f01 = fit(data02',Y','fourier8');

ener_tfv_tot = data12+data13;
[Z,W]= max(ener_tfv_tot);
fv_max = data03(W);
f2 = fit(data04',fv_max,'fourier8');
f02 = fit(data04',Z','fourier8');

[C,V]= max(data11' - 273.15);
Thic_max = (data04(V))';
f3 = fit(data03(2:length(data03)),Thic_max(2:length(data03)),'exp2');
f03 = fit(data03(2:length(data03)),(C(2:length(data03)))','smoothingspline');

ener_tfv_tot = data12+data13;
[R,U]= min(ener_tfv_tot');
th_min = data04(U);
f4 = fit(data03,th_min','fourier3');
f04 = fit(data03,R','fourier8');

figure(1)
hold on
surf(data02,data01,data7)
shading interp
colormap jet
zlabel('Electrical energetic efficiency,\eta_{el}')
ylabel('Mass flow rate (kg/s)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(2)
hold on
surf(data02,data01,data8)
shading interp
colormap jet
zlabel('Thermal energetic efficiency,\eta_{el}')
ylabel('Mass flow rate (kg/s)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(3)
hold on
surf(data02,data01,data6 - 273.15)
shading interp
colormap jet
zlabel('Outlet temperature (°C)')
ylabel('Mass flow rate (kg/s)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(4)
hold on
surf(data02,data01,Au_mxt_exer_Th)
plot3(data02,f0(data02),f00(data02),'k','LineWidth',1)
shading interp
colormap jet
zlabel('Thermal exergetic efficiency,\eta_{th}')
ylabel('Mass flow rate (kg/s)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(5)
hold on
surf(data02,data01,Au_mxt_exer_Tot)
plot3(data02,f1(data02),f01(data02),'k','LineWidth',1)
shading interp
colormap jet
zlabel('Total exergetic efficiency,\eta_{tot}')
ylabel('Mass flow rate (kg/s)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(6)
plot(data02,f1(data02),'LineWidth',1.5)
xlabel('Nanofluid thickness (m)')
ylabel('Mass flow rate (kg/s)')
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
grid minor, ax.YMinorGrid = 'off';,ax.XMinorGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(7)
hold on
surf(data04,data03,data12)
shading interp
colormap jet
ylim ([0 1/10^3])
zlabel('Electrical energetic efficiency,\eta_{el}')
ylabel('Volumetric fraction(%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(8)
hold on
surf(data04,data03,data13)
shading interp
colormap jet
ylim ([0 1/10^3])
zlabel('Thermal energetic efficiency,\eta_{el}')
ylabel('Volumetric fraction(%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(9)
hold on
surf(data04,data03,Au_fvxt_exer_Th)
shading interp
colormap jet
ylim ([0 1/10^3])
zlabel('Thermal exergy efficiency,\eta_{el}')
ylabel('Volumetric fraction(%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(10)
hold on
surf(data04,data03,Au_fvxt_exer_Tot)
shading interp
colormap jet
ylim ([0 1/10^3])
zlabel('Total exergetic efficiency,\eta_{tot}')
ylabel('Volumetric fraction(%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(11)
hold on
surf(data04,data03,data11 - 273.15)
plot3(f3((data03(2:length(data03)))'),data03(2:length(data03)),f03((data03(2:length(data03)))'),'k','LineWidth',1)
shading interp
colormap jet
zlabel('Outlet temperature (°C)')
ylim ([0 1/10^3])
ylabel('Volumetric fraction(%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(12)
hold on
surf(data04,data03,data12+data13)
plot3(data04,f2(data04),f02(data04),'k','LineWidth',1)
%plot3(f4(data03),data03,f04(data03),'k','LineWidth',1)
shading interp
colormap jet
ylim ([0 1/10^3])
%xlim([0.0005 0.02])
zlabel('Total energetic efficiency,\eta_{tot}')
ylabel('Volumetric fraction (%)', 'Rotation',20)
xlabel('Nanofluid thickness (m)', 'Rotation',-10)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(13)
hold on
plot(data04,f2(data04),'LineWidth',2)
plot(f4(data03),data03,'LineWidth',2)
ylabel('Volumetric fraction (%)')
xlabel('Nanofluid thickness (m)')
xlim([0.0005 0.02])
ylim([1/1e5 1/1e3]);
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
grid minor, ax.YMinorGrid = 'off';,ax.XMinorGrid = 'on';
legend('Maximum efficiency','Minimum efficiency','Location','northeast','Orientation','vertical')
title('Gold Nanofluid','LineWidth',2);
legend boxoff

figure(14)
hold on
plot(data04,f2(data04),'LineWidth',2)
ylabel('Volumetric fraction (%)')
xlabel('Nanofluid thickness (m)')
ylim ([0 1/10^3])
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
grid minor, ax.YMinorGrid = 'off';,ax.XMinorGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(15)
plot(data03,f4(data03),'LineWidth',2)
xlabel('Volumetric fraction (%)')
ylabel('Nanofluid thickness (m)')
ylim([0.0005 0.02])
xlim([1.07/1e5 1/1e3]);
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
grid minor, ax.YMinorGrid = 'off';,ax.XMinorGrid = 'on';
title('Gold Nanofluid','LineWidth',2);

figure(16)
hold on
MarkerIndices1 = 1:5:length(D);
plot(D,Au_efth_d+Au_efel_d,'g','LineWidth',1.5)
plot(D, Au_efth_d, ':r','LineWidth',1.5)
plot(D,Au_efel_d,'--b','LineWidth',1.5)
ylim ([0 max(Au_efth_d+Au_efel_d)*1.33])
xlim([0 max(D)])
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
ylabel('Energetic efficiency,\eta_{}','FontSize',14)
xlabel('Nanoparticle Diameter (nm)','FontSize',14)
legend('Total efficiency','Thermal efficiency','Electrical efficiency','Location','northeast','Orientation','vertical')
legend boxoff
title('Gold Nanofluid','LineWidth',2);

figure(17)
hold on
plot(D,Au_efel_d+Au_d_carnot,'g','LineWidth',1.5)
plot(D,Au_d_carnot,':r','LineWidth',1.5)
plot(D,Au_efel_d,'--b','LineWidth',1.5)
ylim ([0 max(Au_efel_d+Au_d_carnot)*1.33])
xlim([0 max(D)])
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
ylabel('Exergetic efficiency,\eta_{ex}','FontSize',14)
xlabel('Nanoparticle Diameter (nm)','FontSize',14)
legend('Total efficiency','Thermal efficiency','Electrical efficiency','Location','northeast','Orientation','vertical')
legend boxoff
title('Gold Nanofluid','LineWidth',2);

figure(18)
hold on
plot(m_flow,Au_efth_m+Au_efel_m,'g','LineWidth',1.5)
plot(m_flow,Au_efth_m,':r','LineWidth',1.5)
plot(m_flow,Au_efel_m,'--b','LineWidth',1.5)
ylim ([0 max(Au_efth_m+Au_efel_m)*1.33])
xlim([0 max(m_flow)])
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
ylabel('Energetic efficiency,\eta_{}','FontSize',14)
xlabel('Mass flow rate (kg/s)','FontSize',14)
legend('Total efficiency','Thermal efficiency','Electrical efficiency','Location','northeast','Orientation','vertical')
legend boxoff
title('Gold Nanofluid','LineWidth',2);

figure(19)
hold on
plot(m_flow,Au_efel_m+Au_m_carnot,'b','LineWidth',2)
[F G]=max(Au_efel_m+Au_m_carnot);
xlim([0 max(m_flow)])
xline(m_flow(G),'--','LineWidth',1.5);
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
ylabel('Exergetic efficiency,\eta_{ex}','FontSize',14)
xlabel('Mass flow rate (kg/s)','FontSize',14)
title('Gold Nanofluid','LineWidth',2);

figure(20)
hold on
plot(m_flow,Au_Tnf_m - 273.15,'r','LineWidth',2)
grid on, ax = gca;, ax.XGrid = 'on';,ax.YGrid = 'on';
ylabel('Outlet temperature (°C)','FontSize',14)
xlabel('Mass flow rate (kg/s)','FontSize',14)
title('Gold Nanofluid','LineWidth',2);

