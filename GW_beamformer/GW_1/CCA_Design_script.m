% Optimized Positions of Sensors in Concentric Circular Array (UCCA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all ; clc ; close all ;

c = 340 ; Ts = 1/16000 ; FS = 1/Ts ; 
f_max = FS / 2 ;
lambda_min = c / f_max ;

% radius of each ring (including the central sensor)
r_p = [0 : 5 : 20]' * 10^(-2) ; % fixed unoptimized rings
P = length(r_p) ;

% Number of sensors per ring as per Nyquist criteria
M_max = pi ./ asin( lambda_min / 4 ./ r_p ) ;
M_max = round( M_max ) ;

for p = 1 : P
    
    if r_p(p) == 0
        M_max(p) = 1 ;
        continue ;
        
    elseif rem( M_max(p), 2) > 0
        M_max(p) = M_max(p) + 1 ;
    
    end
        
end

M_tot_max = sum( M_max ) ;

sensor_fraction = 1 ; % 51/M_tot_max ; %%
M_tot = round( sensor_fraction * M_tot_max ) ; % number of sensors available

if and( sum( r_p == 0 ) , rem( M_tot, 2 ) == 0 )
    M_tot = M_tot + 1 ;
end

% % choose frequency
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

theta_d = 0 ; % Evevation - DOA of the SOI between [0,90]
phi_d = 0 ; % Azimuth - DOA of the SOI between [0,90]

% Sensors angles per ring
%------------------------
P = length(r_p) ;
phi_p_m = cell(1, P) ;
for p = 1 : P
    if r_p(p) ==  0
        phi_p_m{p} = 0 ;
    end
end
init_phi_p_m = phi_p_m ;

F_low = 2000 ; % Hz
F_high = 2000 ; % Hz

[ phi_p_m, D, D_bb, D_bb_low, D_bb_high, D_bb_ratio] = Optimized_Sensor_Positions_CCA( r_p, init_phi_p_m, M_tot, M_max, theta_d, phi_d, f, c, Ts, F_low, F_high ) ;

% checking for rings with no sensors
check = zeros(P,1) ;
for p = 1:P
    if isempty( phi_p_m{p} )
        check(p) = 1 ;
    end
end

phi_p_m( find(check) ) = [] ;
r_p( find(check) ) = [] ;
%save('CCA_design', 'r_p', 'phi_p_m') ;


% Plot the design
%------------------------------------------------
figure();
for p = 1:length(r_p)
    if isempty( phi_p_m{p} )
        continue ;
    end
    polarplot( phi_p_m{p}, r_p(p), 'bo' ) ; hold on ;
end
title(['sensors positions']) ; %axis('tight') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');



% Plot DF vs sensors
%------------------------------------------------
num_sensors = 1 + ([1:length(D_bb)]' - 1)* 2 ;
figure() ;
plot( num_sensors, 10*log10(D_bb) ) ; hold on ;
plot( num_sensors, 10*log10(D_bb_low) , '--') ; hold on ;
plot( num_sensors, 10*log10(D_bb_high) , ':') ; hold on ;
xlim([ min(num_sensors), max(num_sensors)]) ; 
ylim([ 0, 20]) ; 
xticks( round( linspace(min(xlim), max(xlim), 4) ) ) ;
yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
legend({'$\mathcal{D}$','$\mathcal{D}_{\mathrm{lf}}$','$\mathcal{D}_{\mathrm{hf}}$'},'Interpreter','Latex') ; 
xlabel('$M$') ; ylabel('dB') ;  
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');









% num_sensors = 1 + ([1:length(D_bb)]' - 1)* 2 ;
% figure() ;
% subplot(1,3,1) ; plot( num_sensors, 10*log10(D_bb) ) ; title(['DF (broadband)']) ; axis('tight') ;
% xlabel('$M$') ; ylabel('$\mathcal{D}$~dB') ;  
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
% 
% subplot(1,3,2) ; plot( num_sensors, 10*log10(D_bb_low) ) ; title(['DF (bb - lf)']) ; axis('tight') ;
% xlabel('$M$') ; ylabel('$\mathcal{D}$~dB') ;  
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
% 
% subplot(1,3,3) ; plot( num_sensors, 10*log10(D_bb_high) ) ; title(['DF (bb - hf)']) ; axis('tight') ;
% xlabel('$M$') ; ylabel('$\mathcal{D}$~dB') ;   
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');





% subplot(1,5,4) ; 
% plot( num_sensors, D_bb_ratio ) ; hold on ;
% title(['DF (bb - lf/hf)']) ; axis('tight') ;
% xlabel('distribution') ; ylabel('$\mathcal{D}$') ;  
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
% 
% subplot(1,5,5) ; 
% plot( num_sensors, abs([nan ; diff(D_bb_ratio)]) ) ; hold on ;
% title(['Derivative DF (bb - lf/hf)']) ; axis('tight') ;
% xlabel('distribution') ; ylabel('$\mathcal{D}$') ;  
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
