%%% Calculate PowerPattern (PP) for UCCA
%%% h : beamformers at different frequencies, for all concentric rings
%%% r_p : radii of all rings
%%% phi_p_m : angles of sensors in all rings
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Ts = 1 / Fs : sampling period of the signals (8 kHz)

function [theta_range, phi_range, B] = BP_CCA(h, r_p, phi_p_m, f, c, Ts)

P = length(h) ; % number of rings, including the central sensor
M_all = zeros(P,1) ;
for p = 1 : P
        M_all(p) = size(h{p}, 1) ;
end
M_tot = sum(M_all) ;

% re-organize beamformers
h_all = zeros(M_tot, length(f)) ;
for p = 1 : length(h) % = P or P+1

    if p == 1
        beg_sensor = 1 ;
    else
        beg_sensor = sum(M_all(1:p-1)) + 1 ;
    end
    end_sensor = sum(M_all(1:p)) ;

    h_all( beg_sensor : end_sensor , : ) = h{p} ;

end


theta_range = [-180:180]' ; % Elevation angle of incidence of SOI
phi_range = [-180:180]' ; % Azimuth angle of incidence of SOI
B = zeros( length(theta_range), length(phi_range), length(f) ) ;  
for idx_theta_range = 1 : length(theta_range)

    theta = theta_range( idx_theta_range ) ;

    for idx_phi_range = 1 : length(phi_range)

        phi = phi_range( idx_phi_range ) ;

        % evaluate steering vector for DOA (Elevation, Azimuth)
        [ d ] = d_CCA( r_p, phi_p_m, theta, phi, f, c, Ts ) ;

        % re-organize steering vector
        d_all = zeros(M_tot, length(f)) ;
        for p = 1 : length(d) % = P or P+1

            if p == 1
                beg_sensor = 1 ;
            else
                beg_sensor = sum(M_all(1:p-1)) + 1 ;
            end
            end_sensor = sum(M_all(1:p)) ;

            d_all( beg_sensor : end_sensor , : ) = d{p} ;

        end

        B( idx_theta_range, idx_phi_range, : ) = nansum( conj(h_all).* d_all , 1) ;

    end

end


end
