%%% Modified Kaiser Window for Concentric Circular Array (UCCA)
%%% M_active : Percentage of sensors, active in each ring 
%%% r_p : radii of all the rings (including central sensor) 
%%% phi_p_m : sensor positions (angles) for all rings 
%%% beta_kaiser_p :  kaiser window width-parameter for each ring, for each frequency
%%% ring_weight_p : weight (emphasis) given to each ring (including central sensor), for each frequency

function [ h ] = Modified_Bidir_Kaiser_CCA( M_active, beta_kaiser_p, ring_weight_p, r_p, phi_p_m, theta_d, phi_d, f, c, Ts )

% 3-D unit vector representing the DOA of the SOI 
cart_d = [sin(theta_d*pi/180)*cos(phi_d*pi/180), sin(theta_d*pi/180)*sin(phi_d*pi/180), cos(theta_d*pi/180)] ; 
cart_opp_d = - cart_d ;

[ d_phase ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ; % theta in degrees

P = length(r_p) ; % number of rings (including the central sensor)
h = cell(1, P) ;   

for p = 1 : P
    
    ring_radius = r_p(p) ;
    sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
    theta_m = pi/2 ; % ring is in horizontal plane
    
    M_p = length(sensors_angles) ;
    
    if ring_radius == 0
        h{p} = ring_weight_p{p} .* d_phase{p} ;
        continue ;
    end
    
    if ring_radius == 0
        m = 0 ;
    else
        K_p = M_p/2 ;
        m = [-K_p+1:K_p]' ;
    end
    
    % Active and dead sensors
    %---------------------------------------------        
    M = round( ( M_active / 100 ) * M_p ) ;
    M = 2 + round( (M-2)/4 )*4 ;
    if M < 2
        M = 2 ;
    elseif M > M_p
        M = M_p ;
    end

    angle_sep = abs( wrapToPi(pi/2 + (theta_d*pi/180) - sensors_angles) ); % theta converted to radians
    [~, dead_sensors] = sort( angle_sep ) ;

    M_dead = M_p - M ;

    sensors_status = ones(M_p,1) ;
    if M_dead > 0
        dead_sensors = dead_sensors( [(1:M_dead/2 ), (M_p-M_dead/2+1):M_p] ) ;
        dead_sensors = sort(dead_sensors) ;
        sensors_status(dead_sensors) = 0 ;
    end      

    % 3-D unit vectors representing positions of sensors in the ring
    % distance between DOA unit vector and sensor unit vector
    %-------------------------------------------------------------------------------
    unit_vec_dist = zeros(M_p, 1) ;
    for idx_m = 1 : length(m)

        cart_p_m = [sin(theta_m)*cos( sensors_angles(idx_m) ), sin(theta_m)*sin( sensors_angles(idx_m) ), cos(theta_m)] ;

        if abs( angdiff( sensors_angles(idx_m), phi_d*pi/180  ) ) <= pi/2

            tmp_dist = dist( [cart_d; cart_p_m]' ) ; 
            tmp_dist = tmp_dist(end,1) ; 

        else

            tmp_dist = dist( [cart_opp_d; cart_p_m]' ) ; 
            tmp_dist = tmp_dist(end,1) ; 

        end

        unit_vec_dist(idx_m) = tmp_dist ;

    end

    if theta_d > 0
        best_aligned_sensor = [sin(theta_m)*cos(phi_d*pi/180), sin(theta_m)*sin(phi_d*pi/180), cos(theta_m)] ; 
        min_poss_dist = [ cart_d; best_aligned_sensor ] ; % always
        min_poss_dist = dist( min_poss_dist' ) ;
        min_poss_dist = min_poss_dist(2,1) ;
        
        orth_phi_d = phi_d + 90 ; 
        worst_aligned_sensor = [sin(theta_m)*cos(orth_phi_d*pi/180), sin(theta_m)*sin(orth_phi_d*pi/180), cos(theta_m)] ; 
        max_poss_dist = [ cart_d; worst_aligned_sensor ] ; % always
        max_poss_dist = dist( max_poss_dist' ) ;
        max_poss_dist = max_poss_dist(2,1) ;
        
        unit_vec_dist = unit_vec_dist - min_poss_dist ;
        unit_vec_dist = unit_vec_dist / (max_poss_dist -  min_poss_dist) ;

%         unit_vec_dist = unit_vec_dist - min_poss_dist ;
%         unit_vec_dist = unit_vec_dist / max(unit_vec_dist) ;
        m_angle_correlation = unit_vec_dist ;
    else
        m_angle_correlation = zeros(size(unit_vec_dist)) ;
    end

    % kaiser window
    %-------------------------------------------------------------------------------
    Kaiser_window_p = zeros(M_p, length(f)) ;
    beta_kaiser = beta_kaiser_p{p} ;
    for idx_f = 1 : length(f)
        Kaiser_window_p( :, idx_f) = besseli( 0 , beta_kaiser(idx_f) * sqrt( 1- (m_angle_correlation).^2 ) ) / besseli( 0 , beta_kaiser(idx_f) ) ;
    end


    % kaiser window beamformer
    %-------------------------------------------------------------------------------
    h_p = Kaiser_window_p .* ( ones(M_p, 1) * ring_weight_p{p} ) ;
    h_p = h_p .* ( sensors_status * ones(1,length(f)) ) ;
    h_p = h_p .* d_phase{p} ;
    h{p} = h_p ;

end


% normalizing 
h_all = [] ; M_all = zeros( size( length(h)+1, 1) ) ;
for p = 1 : length(h)
    h_all = [h_all ; h{p}] ;
    M_all(p+1) = size(h{p},1) ;

    if p < length(h)
        continue ;
    end

    h_all = h_all ./ (   ones( size(h_all,1), 1 ) * sum( abs(h_all) )   ) ;
    M_all = cumsum( M_all ) ;
end

for p = 1 : length(h)
    h{p} = h_all( [M_all(p) + 1 : M_all(p+1)]', :) ;
end 

        
end


