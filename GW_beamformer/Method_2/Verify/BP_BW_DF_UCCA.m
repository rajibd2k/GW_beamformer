%%% Calculate PowerPattern (PP) for UCCA
%%% h : beamformers at different frequencies, for all concentric rings
%%% theta_d, phi_d : Elevation and Azimuth of SOI for which beamformers are obtained
%%% power_level_diff : 6 dB
%%% M1 : number of sensors in first ring (odd)
%%% central_sensor : y/n
%%% P : number of concentric rings in the UCCA
%%% Delta_r : increase in radius from one ring to the next
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Ts = 1 / Fs : sampling period of the signals (8 kHz)

function [theta_range, phi_range, B, B_dB, b_phi, b_phi_bb, D, D_bb] = BP_BW_DF_UCCA(h, theta_d, phi_d, power_level_diff, M_1, central_sensor, P, Delta_r, f, c, Ts)

if strcmp(central_sensor, 'y') % a central sensor exists
    M_all = [1 ; zeros(P,1) ] ;
else
    M_all = zeros(P,1) ;
end

for p = 1 : P
    if strcmp(central_sensor, 'y') % a central sensor exists
        M_all(p+1) = size(h{p+1}, 1) ;
    else
        M_all(p) = size(h{p}, 1) ;
    end
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

theta_range = [0:180]' ; % Elevation angle of incidence of SOI
phi_range = [-180:180]' ; % Azimuth angle of incidence of SOI
B = zeros( length(theta_range), length(phi_range), length(f) ) ;  
for idx_theta_range = 1 : length(theta_range)

    theta = theta_range( idx_theta_range ) ;

    for idx_phi_range = 1 : length(phi_range)

        phi = phi_range( idx_phi_range ) ;

        % evaluate steering vector for DOA (Elevation, Azimuth)
        [ d, ~, ~, ~, M_all ] = d_UCCA( M_1, central_sensor, P, Delta_r, theta, phi, f, c, Ts ) ;
        M_tot = sum(M_all) ;

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

        B( idx_theta_range, idx_phi_range, : ) = abs( nansum( conj(h_all).* d_all ) ) ;

    end

end

B_dB = zeros(size(B)) ;

% Directivity and Beamwidth - Azimuth
%---------------------------------------------------------------------------------------
% beampattern at DOA
[~, idx_theta_DOA] = min( abs( theta_range - theta_d ) ) ;
[~, idx_phi_DOA] = min( abs( phi_range - phi_d ) ) ;

B_d = B(idx_theta_DOA, idx_phi_DOA, :) ;
B_d = reshape( B_d, 1, length(f) ) ;

Delta_theta = mean( abs(diff(theta_range)) ) * pi / 180 ;
Delta_phi = mean( abs(diff(phi_range)) ) * pi / 180 ;

% Directivity
D = zeros( 1, length(f) ) ;
D_den = zeros( 1, length(f) ) ;
% Beamwidth-Azimuth
b_phi = zeros(1, length(f)) ;

for idx_f = 1 : length(f)

    B_d_f = abs( B_d(idx_f) )^2 ;
    B_f = B(:, :, idx_f) ;
    B_f = reshape(B_f, length(theta_range), length(phi_range) ) ;
    B_f = abs( B_f ).^2 ;

    % Directivity
    D_den_f = B_f .* ( sin(theta_range*pi/180) * ones(1, length(phi_range)) ) ;
    D_den_f = nanmean(D_den_f, 2) * length(phi_range) * Delta_phi ;
    D_den_f = nanmean(D_den_f) * length(theta_range) * Delta_theta ;
    D_den_f = D_den_f / (4*pi) ;
    
    D_den(idx_f) = D_den_f ;
    D(idx_f) = B_d_f / D_den_f ;
    
    
    % Beamwidth-Azimuth
    if max( max(B_f) ) > 0
        B_f = B_f / max( max( B_f ) ) ;
    end
    B_dB_f = 10*log10( B_f ) ;
    B_dB(:,:,idx_f) = B_dB_f ;
    
    B_d_dB_f = B_dB_f(idx_theta_DOA, idx_phi_DOA) ;
    B_dB_f_cutoff = B_d_dB_f - power_level_diff ;
    
    tmp = find( B_dB_f(idx_theta_DOA, :) <= B_dB_f_cutoff ) ;
    if isempty(tmp)
        b_phi(idx_f) = 360 ;
        continue ;
    end
    tmp = phi_range(tmp) ;
    
    tmp_minus = tmp( find(tmp < phi_d ) ) ;
    if isempty(tmp_minus)
        tmp_minus = -180 ;
    else
        tmp_minus = max( tmp_minus ) ;
    end
    
    tmp_plus = tmp( find(tmp > phi_d ) ) ;
    if isempty(tmp_plus)
        tmp_plus = 180 ;
    else
        tmp_plus = min( tmp_plus ) ;
    end
    
    b_phi(idx_f) = tmp_plus - tmp_minus ;

end

D_bb = nanmean( abs(B_d).^2  ) / nanmean(D_den) ;

b_phi_bb = nanmean( b_phi ) ; 

end
