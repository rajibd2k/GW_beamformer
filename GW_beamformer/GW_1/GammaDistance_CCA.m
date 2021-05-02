%%% Find Gamma_distance : Matrix of Sensor-Distances of Concentric Circular Array (UCCA)
%%% r_p : radii of all rings
%%% phi_p_m : angles of sensors in all rings
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Fs : sampling period of the signals (8 kHz)

function [ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m )

    P = length(r_p) ;
    M_all = zeros(P,1) ;
    for p = 1 : P
        M_all(p) = size(phi_p_m{p}, 1) ;
    end
    M_tot = sum(M_all) ;
    
    positions = zeros(2, M_tot) ;
    
    % For each ring
    for p = 1 : P
        
        ring_radius = r_p(p) ;
        sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
        
        % Calculate Gamma
        %----------------------
        positions_x = ring_radius * cos( sensors_angles ) ; % cartesian co-ordinates / speed of propagation
        positions_y = ring_radius * sin( sensors_angles ) ;
        positions_p = [positions_x positions_y]' ;     
        
        % re-organize data
        if p == 1
            beg_sensor = 1 ;
        else
            beg_sensor = sum(M_all(1:p-1)) + 1 ;
        end
        end_sensor = sum(M_all(1:p)) ;
        
        positions( : ,  beg_sensor : end_sensor ) = positions_p ;

    end
        
    Gamma_distance = dist(positions) ;
    
end


