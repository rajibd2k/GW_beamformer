%%% Sensor Positions for Concentric Circular Array (CCA)
%%% Optimize Directivity of DS-beamformer for Broadside direction (theta_d = 0)

%%% r_p : radii of all the rings (including central sensor) 
%%% phi_p_m : final sensor positions (angles) for all rings 
%%% F_low, F_high : frequency range for calculating D_bb_low, D_bb_high

%%% D : directivity factor of DS beamformer, for all frequencies, for all iterations
%%% D_bb : broadband directivity factor of DS beamformer, for all iterations
%%% D_bb_low : broadband directivity factor (low-frequencies only) of DS beamformer, for all iterations
%%% D_bb_high : broadband directivity factor (high-frequencies only) of DS beamformer, for all iterations
%%% D_bb_ratio : broadband directivity factor ratio (low-frequency / all-frequency) of DS beamformer, for all iterations

%%% theta_d : DOA-Elevation of SOI (in degrees) : 0 (Broadside - Topview)
%%% phi_d : DOA-Azimuth of SOI (in degrees) : inconsequential
%%% f : normalized frequencies [0, 0.5]
%%% c : speed of sound propagation (340 m/s)
%%% Ts : sampling period (= 1/Fs, Fs = 8000 Hz)

%%% init_phi_p_m : initial sensor positions (angles)  
%%% M_tot : total number of sensors to be distributed
%%% M_max : maximum number of sensors allowed per ring

function [ phi_p_m, D, D_bb, D_bb_low, D_bb_high, D_bb_ratio] = Optimized_Sensor_Positions_CCA( r_p, init_phi_p_m, M_tot, M_max, theta_d, phi_d, f, c, Ts, F_low, F_high )

P = length(r_p) ; % number of rings, including the central sensor

init_M = 0 ;
for p = 1 : P
    init_M = init_M + length(init_phi_p_m{p}) ;
end

rem_M = M_tot - init_M ; % remaining number of sensors to be distributed among rings
num_dist = rem_M / 2 ; % sensors are distributed in pairs


% % figure to visualize
% figure(1);
% for p = 1:P
% 
%     if isempty( init_phi_p_m{p} )
%         continue ;
%     end
%     polarplot( init_phi_p_m{p}, r_p(p), 'o' ) ; hold on ;
% 
% end
% 
% title(['sensors positions']) ; %axis('tight') ;
% b=gca;
% set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
% a=findobj(gcf); % get the handles associated with the current figure
% alllines=findall(a,'Type','line');
% alltext=findall(a,'Type','text');
% set(alllines,'Linewidth',2, 'MarkerSize', 10);
% set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
% 
% pause() ;


% Searching for senor positions in each ring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = zeros(num_dist+1,length(f)) ;
D_bb = zeros(num_dist+1, 1) ;
D_bb_low = zeros(num_dist+1, 1) ;
D_bb_high = zeros(num_dist+1, 1) ;
D_bb_ratio = zeros(num_dist+1, 1) ;

% initialize 
% --------------------------------------------------------------
phi_p_m = init_phi_p_m ;
[ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
[ h ] = DS_CCA( d ) ;
% Metrics
[ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m ) ;
[ D(1,:), D_bb(1), D_bb_low(1), D_bb_high(1), D_bb_ratio(1) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;

for idx_dist = 1 : num_dist
    
    optimization_var = zeros(P,1) ;
    M_dist = zeros(P,1) ;
    for p = 1:P
        
        ring_radius = r_p(p) ;
        sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
        M_p = length(sensors_angles) ;
        M_dist(p) = M_p ;
        
        if ring_radius == 0
            optimization_var(p) = inf ;
            continue ;
        end

        tmp_M_p = M_p + 2 ;

        tmp_K_p = tmp_M_p / 2 ;
        tmp_m = [-tmp_K_p+1 : tmp_K_p]' ; % sensors indices
        tmp_phi_p = 2*pi / tmp_M_p ;
        tmp_sensors_angles = tmp_m * tmp_phi_p ;
        
        tmp_phi_p_m = phi_p_m ;
        tmp_phi_p_m{p} = tmp_sensors_angles ;
        
        [ tmp_d ] = d_CCA( r_p, tmp_phi_p_m, theta_d, phi_d, f, c, Ts ) ;
        [ tmp_h ] = DS_CCA( tmp_d ) ;
        % Metrics
        [ tmp_Gamma_distance ] = GammaDistance_CCA( r_p, tmp_phi_p_m ) ;
        [ tmp_D, tmp_D_bb, tmp_D_bb_low, tmp_D_bb_high, tmp_D_bb_ratio ] = DFanalytical_CCA( tmp_d, tmp_h, tmp_Gamma_distance, f, c, Ts, F_low, F_high ) ;

         %optimization_var(p) = -log10(tmp_D_bb_low) ;
         optimization_var(p) = -log10(tmp_D_bb) ;
        % optimization_var(p) = -log10(tmp_D_bb_high) ;
        % optimization_var(p) = -0.7*log10(tmp_D_bb_low) -0.3*log10(tmp_D_bb_high) ;

    end
    
    % check in which rings sensors can be added
    check = (M_dist - M_max) < 0 ;
    check = check + 0 ;
    
    optimization_var = optimization_var .* check ;
    [~, idx_ring] = min( optimization_var ) ;
    M_p = M_dist(idx_ring) + 2 ;
    
    K_p = M_p / 2 ;
    m = [-K_p+1 : K_p]' ; % sensors indices
    phi_p = 2*pi / M_p ;
    sensors_angles = m * phi_p ;
    phi_p_m{idx_ring} = sensors_angles ;
    
    [ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
    [ h ] = DS_CCA( d ) ;
    % Metrics
    [ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m ) ;
    [ D(1+idx_dist,:), D_bb(1+idx_dist), D_bb_low(1+idx_dist), D_bb_high(1+idx_dist), D_bb_ratio(1+idx_dist) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;

    
%     % figure to visualize
%     clf ;
% 
%     figure(1);
%     for p = 1:P
% 
%         if isempty( phi_p_m{p} )
%             continue ;
%         end
%         polarplot( phi_p_m{p}, r_p(p), 'o' ) ; hold on ;
%         
%     end
%     
%     title(['sensors positions']) ; %axis('tight') ;
%     b=gca;
%     set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
%     a=findobj(gcf); % get the handles associated with the current figure
%     alllines=findall(a,'Type','line');
%     alltext=findall(a,'Type','text');
%     set(alllines,'Linewidth',2, 'MarkerSize', 10);
%     set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
%     
%     pause();

    
end
        
end


