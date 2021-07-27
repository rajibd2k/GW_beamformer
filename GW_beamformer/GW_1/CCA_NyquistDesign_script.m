% CCA (Nyquist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all ; clc ; close all ;

c = 340 ; Ts = 1/16000 ; FS = 1/Ts ; 
f_max = FS / 2 ;
lambda_min = c / f_max ;

% radius of each ring (including the central sensor)
r_p = [0 : 5 : 20]' * 10^(-2) ;
%r_p = [0 : 10 : 20]' * 10^(-2) ; % 2-ring design
P = length(r_p) ;

% Number of sensors per ring as per radii ratio
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

sensor_fraction = 1 ;
M_tot = round( sensor_fraction * M_tot_max ) ; % number of sensors available

if and( sum( r_p == 0 ) , rem( M_tot, 2 ) == 0 )
    M_tot = M_tot + 1 ;
end

% % choose frequency
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

% Sensors angles per ring
%------------------------
P = length(r_p) ;
phi_p_m = cell(1, P) ;
for p = 1 : P
    if r_p(p) ==  0
        phi_p_m{p} = 0 ;
    else
        K_p = M_max(p)/2 ;
        m = [-K_p + 1 : K_p]' ;
        phi_p_m{p} = m * 2 * pi / M_max(p) ;
    end
end

save('CCA_design_Nyquist', 'r_p', 'phi_p_m') ;

