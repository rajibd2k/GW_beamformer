%%% Constant Beamwidth Beamforming for Uniform Concentric Circular Array (UCCA)
%%% "Theoretical and Experimental Studies on Broadband Constant Beamwidth Beamforming for Circular Arrays"
%%% Optimize Azimuth only

%%% h : final beamformer, for all frequencies, for all the concentric rings

%%% theta_d : DOA-Elevation of SOI (in degrees)
%%% phi_d : DOA-Azimuth of SOI (in degrees)
%%% f : normalized frequencies [0, 0.5]
%%% c : speed of sound propagation (340 m/s)
%%% Ts : sampling period (= 1/Fs, Fs = 8000 Hz)

%%% phi_BW : desired beamwidth-Azimuth (in degrees)
%%% power_level_diff : power level (in dB) for calculating beamwidth (= 6)

function [ h ] = Method_1( phi_BW, power_level_diff, bessel_epsilon, r_p, phi_p_m, theta_d, phi_d, f, c, Ts )

theta_range = [0:180]' ; % degrees
phi_range = [-180:180]' ; % degrees

[ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;

% DS beamformer
[ h_DS ] = DS_CCA( d ) ;

% reference frequency of each ring - one which has the desired beamwidth
[ f_ref_index ] = BP_BW_fref_ring_CCA(h_DS, theta_d, phi_d, phi_BW, power_level_diff, r_p, phi_p_m, f, c, Ts) ;

% Calculate T_matrix at each frequency of each ring
[ T ] = T_Matrix( r_p, phi_p_m, theta_d, f, c, Ts, bessel_epsilon ) ;

h = cell(1, length(d)) ;
h_all = [] ;
for p = 1 : length(d)
    
    if r_p(p) == 0
        h{p} = h_DS{p} ;
        continue ;
    end
    
    % DS beamformer (at reference frequency) of each ring
    h_p_ref = h_DS{p}(:, f_ref_index(p) ) ;
    
    h_p = zeros( size(d{p}) ) ;
    tmp1 = T{p}(:,:, f_ref_index(p)) ;
    for idx_f = 1 : length(f)
        
        tmp2 = T{p}(:,:, idx_f) ;       
        T_p_f = tmp1 * pinv( tmp2' * tmp2 ) * tmp2' ;
        
        h_p_f = T_p_f' * h_p_ref ;
        
        h_p(:, idx_f) = h_p_f ;
    end
    
    h{p} = h_p ;
    h_all = [h_all ; h{p}] ;
    
end

h_all = sum( abs(h_all), 1 ) ;

% normalizing 
for p = 1 : length(h)
    h{p} = h{p} ./ (   ones( size(h{p},1), 1 ) * h_all ) ;
end

end


