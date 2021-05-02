% Modified Kaiser Window beamformer for Uniform Concentric Circular Array (UCCA)
% Using Gradient Descent Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;

M_1 = 4 ; P = 2 ; central_sensor = 'n' ; %% y/n
Delta_r = 4*10^-2 ; 
c = 340 ; Ts = 1/20000 ; FS = 1/Ts ; 

% % choose frequency
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

theta_d = 90 ; % Evevation - DOA of the SOI between [0,90]
phi_d = 0 ; % Azimuth - DOA of the SOI between [0,90]

[ d, lambda_min, r_1, delta_1, M_all ] = d_UCCA( M_1, central_sensor, P, Delta_r, theta_d, phi_d, f, c, Ts ) ;

% % DS
% %-------------------------------------------------------------------------------------------
% phi_BW = 60 ;
% power_level_diff = 3 ;
% 
% [ d ] = d_UCCA( M_1, central_sensor, P, Delta_r, theta_d, phi_d, f, c, Ts ) ;
% [ h ] = DS_UCCA( d ) ;
% [theta_range, phi_range, B, B_dB, b_phi, b_phi_bb, D, D_bb] = BP_BW_DF_UCCA(h, theta_d, phi_d, power_level_diff, M_1, central_sensor, P, Delta_r, f, c, Ts) ;


% Method 3
%-------------------------------------------------------------------------------------------
b_2N = ones(3,1) / 3 ;
N = (length(b_2N)-1)/2 ; 
phi_BW = 60 ;
power_level_diff = 6 ;
        
[ h ] = Method_3( phi_BW, M_1, central_sensor, P, Delta_r, theta_d, phi_d, f, c, Ts, N ) ;
[theta_range, phi_range, B, B_dB, b_phi, b_phi_bb, D, D_bb] = BP_BW_DF_UCCA(h, theta_d, phi_d, power_level_diff, M_1, central_sensor, P, Delta_r, f, c, Ts) ;
[W, W_bb] = WNG_UCCA(h, d) ;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; close all ;

% % choose frequency
% Ts = 1/8000 ; FS = 1/Ts ; % Hz
% f = [0 : FS/256 : FS/2]' ; % Hz
% f = f/FS ;

Delta_b_phi = b_phi - phi_BW ;


[~,idx_freq] = min(abs(f*FS - 1000)) ;
freq = f(idx_freq)*FS ;
figure();
imagesc( phi_range, theta_range, B(:,:,idx_freq) ) ; colormap('gray') ; hold on ;
plot3(phi_d, theta_d, 1, 'rx') ;
xticks([-180:90:180]) ; yticks([0:45:180]) ;
xlabel('$\phi$  (degrees)') ; ylabel('$\theta$ (degrees)') ;  title(['$| \mathcal{B} (f,\theta,\phi) |, ~ fF_s =~$', num2str(freq), ' Hz']) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');


[~,idx_freq] = min(abs(f*FS - 1000)) ;
freq = f(idx_freq)*FS ;
idx_theta = find( theta_range == theta_d ) ;
B_phi_f = reshape( B(idx_theta,:,idx_freq), length(phi_range), 1) ;
figure();
polarplot( phi_range*pi/180, B_phi_f  ) ; 
title(['$| \mathcal{B} (f,\theta,\phi) |, ~ fF_s =~$', num2str(freq), ' Hz']) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');



figure();
subplot(3,1,1) ; 
plot(f*FS, abs(Delta_b_phi) ) ; title(['$| \Delta b_{\phi}(f) |$']) ;
xlabel('$f F_s$ (Hz)') ; ylabel('degrees') ; 
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(3,1,2) ; 
plot(f*FS, 10*log10( D ) ) ; title(['$\mathcal{D}(f)$']) ;
xlabel('$f F_s$ (Hz)') ; ylabel('dB') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(3,1,3) ; 
plot(f*FS, 10*log10( W ) ) ; title(['$\mathcal{W}(f)$']) ;
xlabel('$f F_s$ (Hz)') ; ylabel('dB') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');


