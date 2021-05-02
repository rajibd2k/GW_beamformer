%%% Find steering vectors for each ring of Uniform Concentric Circular Array (UCCA)
%%% r_p : radii of all rings
%%% phi_p_m : angles of sensors in all rings
%%% theta_d : Elevation Angle of incidence of the SOI in degrees [0:180]
%%% phi_d : Azimuth Angle of incidence of the SOI in degrees [-180:180]
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Fs : sampling period of the signals (8 kHz)
%%% d : Steering vector of the SOI at each frequency, for all concentric rings

function [ d, lambda_min ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts )

FS = 1 / Ts ;
f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;

% r_1 = round( 100* lambda_min / 4 / sin( pi / M_1 ) ) / 100 ; 
% delta_1 = 2 * r_1 * sin( pi / M_1 ) ; 

P = length(r_p) ;
d = cell(1, P) ;
for p = 1 : P
    
    ring_radius = r_p(p) ;
    sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
    
    hat_r_p = ring_radius / lambda_min ;

    delay = - hat_r_p * alpha * sin( theta_d*pi/180 ) * cos( phi_d*pi/180 - sensors_angles ) ;

    % steering vector of the UCCA
    [f_mat , delay_mat] = meshgrid(f, delay) ;
    d_p = exp(-1i*2*pi*f_mat.*delay_mat) ;

    d{p} = d_p ;

end 

end


