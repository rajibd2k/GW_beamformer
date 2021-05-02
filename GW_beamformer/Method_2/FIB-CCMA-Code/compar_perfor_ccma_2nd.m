clear;
close all;
clc;
save_data = 0;
savefigure = 1;
% in this simulation, we design the first-order CDMA with a UCA or UCCA

% Some common parameers
% set some fixed parameters
c               = 340; 
j               = sqrt(-1);
theta_s         = 0/180*pi;
f_vect          = [100:100:500 600:200:8000];
f_vect2         = [100:10:8000];
theta           = (-180:2:180)'/180*pi;

%% 1) design the second-order FIB 

r_1             = 0.03;
r_2             = 0.015;
b_vect          = [0.1035 0.242 0.309 0.242 0.1035]'; % supercardioid

% 1-1) design the first-order FIB with single ring 5 + 0
P               = 1;
RP_vect         = [r_1];
NP_vect         = [2];
MP_vect         = [5];
phi0_vect       = [0]*pi/180;
run main_cdma_ucca_nosync.m    
G_wn_matrix(:,1)       = G_wn;
G_dn_matrix(:,1)       = G_dn;
BP_matrix1             = BP_matrix;

% 1-2) design the first-order FIB with 7 + 1 
P               = 2;
RP_vect         = [r_1; 0];
NP_vect         = [2; 0];
MP_vect         = [7; 1];
phi0_vect       = [0; 0]*pi/180;
run main_cdma_ucca_nosync.m    

G_wn_matrix(:,2)       = G_wn;
G_dn_matrix(:,2)       = G_dn;
BP_matrix2             = BP_matrix;

% 1-3) design the first-order FIB with 7+5
P               = 2;
RP_vect         = [r_1; r_2];
NP_vect         = [2; 2];
MP_vect         = [7; 5];
phi0_vect       = [0; 0]*pi/180;
run main_cdma_ucca_nosync.m
G_wn_matrix(:,3)       = G_wn;
G_dn_matrix(:,3)       = G_dn;
BP_matrix3             = BP_matrix;

% 1-4) design the second-order FIB with 10+5+1 
% now we can show there still has a null caused by the second order
% Bessel functions.
P               = 3;
RP_vect         = [r_1; r_2; 0];
NP_vect         = [2; 2; 0];
MP_vect         = [10; 5; 1];
phi0_vect       = [0; 0; 0]*pi/180;
run main_cdma_ucca_nosync.m

G_wn_matrix(:,4)       = G_wn;
G_dn_matrix(:,4)       = G_dn;
BP_matrix4             = BP_matrix;

%% plot figures

FontSize        = 12;
c_min           = -30;
right_value     = 0.0; 
left_value      = 0.0;
up_value        = 0.03; 
down_value      = 0.01;

f_high          = max(f_vect2);
f_high          = f_high/1;


%% BP
plot_combine_4BP(theta,f_vect,BP_matrix1,BP_matrix2,BP_matrix3,BP_matrix4);

FontSize        = 10;
%%  DF/WNG
fig=figure('position',[0 0 600 700]); linewd = 1.5;linewd2 = 2.5;
H=subplot(211);
PPP=get(H,'pos');%get the position
PPP = PPP + [left_value down_value right_value up_value];
set(H,'pos',PPP)%set new boundary
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;
plot(f_vect2, 10*log10(G_dn_matrix(:,1)),'-k','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_dn_matrix(:,2)),'--m','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_dn_matrix(:,3)),'-.b','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_dn_matrix(:,4)),'-r','LineWidth',linewd2);grid on;

xtickVect  = [0:1000:10000];
xtickLabelVect = {'0','1','2','3','4','5','6','7','8','9','10'};
set(gca, 'XTick', xtickVect, 'XTickLabel', xtickLabelVect);
ytickVect = [-10:2:10];
ytickLabelVect = {'-10','-8','-6','-4','-2','0','2','4','6','8','10'};
set(gca, 'YTick', ytickVect, 'YTickLabel', ytickLabelVect);
My_legend = {'CMA','CCMA-I','CCMA-II','CCMA-III'};
legend(My_legend,'Location','SouthWest','Interpreter','latex')
ylabel('DI (dB)','fontsize',FontSize,'interpreter', 'latex')
xlabel('$f$ (kHz)','fontsize',FontSize,'interpreter', 'latex')

text(3820,-9,'(a)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 f_high -10 10]);box on;grid minor; 

% plot directive factor
H=subplot(212);
PPP=get(H,'pos');%get the position
PPP = PPP + [left_value down_value right_value up_value];
set(H,'pos',PPP)%set new boundary
set(gcf,'DefaultTextInterpreter','latex')
set(fig, 'Color', [1, 1, 1]);hold on;

plot(f_vect2, 10*log10(G_wn_matrix(:,1)),'-k','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_wn_matrix(:,2)),'--m','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_wn_matrix(:,3)),'-.b','LineWidth',linewd);grid on;
plot(f_vect2, 10*log10(G_wn_matrix(:,4)),'-r','LineWidth',linewd2);grid on;
xtickVect  = [0:1000:10000];
xtickLabelVect = {'0','1','2','3','4','5','6','7','8','9','10'};
set(gca, 'XTick', xtickVect, 'XTickLabel', xtickLabelVect);
legend(My_legend,'Location','SouthWest','Interpreter','latex')
ylabel('WNG (dB)','fontsize',FontSize,'interpreter', 'latex')
xlabel('$f$ (kHz)','fontsize',FontSize,'interpreter', 'latex')

text(3820,-47,'(b)','fontsize',FontSize,'interpreter', 'latex')
set(gca,'fontsize',FontSize,'fontname','Times new roman');
axis([0 f_high -50 20]);box on;grid minor; 


