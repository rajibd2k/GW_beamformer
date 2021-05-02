% Modified Kaiser Window beamformer for Uniform Concentric Circular Array (UCCA)
% Using Gradient Descent Algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;

M_1 = 4 ; P = 5 ; central_sensor = 'y' ; %% y/n
Delta_r = 4*10^-2 ; 
c = 340 ; Ts = 1/8000 ; FS = 1/Ts ; 

theta_d = 90 ; % Elevation - DOA of the SOI between [0,90]
phi_d = 0 ; % Azimuth - DOA of the SOI between [0,90]

theta_range = [0:180]' ;
phi_range = [-180:180]' ; 

% Radii and number of sensors'
%----------------------------------------------------------------------------------
% % choose analog frequency
f_vect = [0 : FS/256 : FS/2] ; % Hz
f_vect2 = f_vect ; % for plotting

f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;
r_1 = round( 100* lambda_min / 4 / sin( pi / M_1 ) ) / 100 ; 

if strcmp(central_sensor, 'y') % a central sensor exists
    M_all = [1 ; zeros(P,1) ] ;
    d = cell(1,P+1) ;
    d{1} = ones(1, length(f_vect)) ;
else
    M_all = zeros(P,1) ;
    d = cell(1,P) ;
end

RP_vect = zeros(size(M_all)) ; 
for p = 1 : P

    if p == 1
        M_p = M_1 ;
        r_p = r_1 ;
    else
        r_p = r_1 + Delta_r*(p-1) ;
        M_p =  round( pi / asin( lambda_min / 4 / r_p ) ) ;
        if rem( M_p , 2 ) > 0
            M_p = M_p + 1 ;
        end
    end
    
    if strcmp(central_sensor, 'y') % a central sensor exists
        M_all(p+1) = M_p ;
        RP_vect(p+1) = r_p ;
    else
        M_all(p) = M_p ;
        RP_vect(p) = r_p ;
    end

end
RP_vect = RP_vect(end:-1:1) ;


% Method 3
%-------------------------------------------------------------------------------------------
N = 32 ;
phi_BW = 60 ;
power_level_diff = 6 ;

% Spatial filter
% symmetric hamming window based filter
n = [-N:N]' ;
b_2N = fir1(length(n)-1, (phi_BW/2)/180,'low');
b_2N = b_2N([1:length(n)]) ;
b_2N = b_2N(:) ;

b_vect = b_2N ;

theta_d = theta_d * pi/180 ;
theta_s = phi_d * pi / 180 ;

NP_vect = N*ones(P,1) ;
if strcmp(central_sensor, 'y')
    NP_vect = [NP_vect; 0] ;
end

MP_vect = M_all(end:-1:1) ;

phi0_vect = (0*pi/180)* ones(size(RP_vect)) ;

theta = phi_range*pi/180 ;

P = P + strcmp(central_sensor, 'y') ;

Filter_Para.c = c;
Filter_Para.P = P ;
Filter_Para.theta_s = theta_s ;
Filter_Para.b_vect = b_vect ;
Filter_Para.RP_vect = RP_vect ;
Filter_Para.NP_vect = NP_vect ;
Filter_Para.MP_vect = MP_vect ;
Filter_Para.phi0_vect = phi0_vect ;
Filter_Para.theta_d = theta_d ;

run main_cdma_ucca_nosync.m  

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc ; close all ;

idx_freq = 1+2^6 ;
freq = f_vect(idx_freq) ;
B_phi_f = BP_matrix(:,idx_freq) ;
figure();
polarplot( phi_range*pi/180, B_phi_f  ) ; 
title(['$| \mathcal{B} (f,\theta_d,\phi) |, ~ fF_s =~$', num2str(freq), ' Hz']) ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');



figure();
subplot(1,2,1) ; plot(f_vect, 10*log10( G_dn ) ) ; title(['$\mathcal{D}(f)$']) ;
xlabel('$f F_s$ (Hz)') ; ylabel('dB') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

subplot(1,2,2) ; 
plot(f_vect, 10*log10( G_wn ) ) ; title(['$\mathcal{W}(f)$']) ;
xlabel('$f F_s$ (Hz)') ; ylabel('dB') ;
b=gca;
set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
a=findobj(gcf); % get the handles associated with the current figure
alllines=findall(a,'Type','line');
alltext=findall(a,'Type','text');
set(alllines,'Linewidth',2, 'MarkerSize', 10);
set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');