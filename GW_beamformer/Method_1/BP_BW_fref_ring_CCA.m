%%% Calculate PowerPattern (PP) for UCCA, for all concentric rings
%%% h : beamformers at different frequencies, for all concentric rings
%%% theta_d, phi_d : Elevation and Azimuth of SOI for which beamformers are obtained
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Ts = 1 / Fs : sampling period of the signals (8 kHz)
%%% phi_BW : desired beamwidth-Azimuth (in degrees)
%%% power_level_diff : power level (in dB) for calculating beamwidth (= 6)

function [ f_ref_index ] = BP_BW_fref_ring_CCA(h, theta_d, phi_d, phi_BW, power_level_diff, r_p, phi_p_m, f, c, Ts)

f_ref_index = zeros( length(h), 1 ) ;

for p = 1 : length(h) % = P or P+1
    
    [ ~, b_p_phi ] = BW_CCA(h(p), r_p(p), phi_p_m(p), theta_d, phi_d, f, c, Ts, power_level_diff) ;
    
    % Frequency index which satisfied desired beamwidth
    [~, f_des_index_p] = min(abs(b_p_phi - phi_BW) ) ;
    f_ref_index(p) = f_des_index_p ;
   
end

end
