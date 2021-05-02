%%% Find steering vectors for each ring of Uniform Concentric Circular Array (UCCA)
%%% M1 : number of sensors in first ring
%%% central_sensor : y/n
%%% P : number of concentric rings in the UCCA
%%% Delta_r : increase in radius from one ring to the next
%%% theta_d : Elevation Angle of incidence of the SOI in degrees [0:180]
%%% phi_d : Azimuth Angle of incidence of the SOI in degrees [-180:180]
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Fs : sampling period of the signals (8 kHz)
%%% d : Steering vector of the SOI at each frequency, for all concentric rings

function [ d, lambda_min, r_1, delta_1, M_all ] = d_UCCA( M_1, central_sensor, P, Delta_r, theta_d, phi_d, f, c, Ts )

FS = 1 / Ts ;
f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;

% r_1 = round( 100* lambda_min / 4 / sin( pi / M_1 ) ) / 100 ; 
% delta_1 = 2 * r_1 * sin( pi / M_1 ) ; 

if strcmp(central_sensor, 'y') % a central sensor exists
    M_all = [1 ; zeros(P,1) ] ;
    d = cell(1,P+1) ;
    d{1} = ones(1, length(f)) ;
else
    M_all = zeros(P,1) ;
    d = cell(1,P) ;
end

for p = 1 : P

%     if p == 1
%         M_p = M_1 ;
%         phi_p = 2*pi / M_p ;
%         r_p = r_1 ;
%     else
%         r_p = r_1 + Delta_r*(p-1) ;
%         M_p =  round( pi / asin( lambda_min / 4 / r_p ) ) ;
%         if rem( M_p , 2 ) > 0
%             M_p = M_p + 1 ;
%         end
%         phi_p = 2*pi / M_p ;
%     end

    if p == 1
        M_p = M_1 ; %4
        phi_p = 2*pi / M_p ;
        r_p = 2.2 * 10^(-2) ;
        r_1 = 2.2 * 10^(-2) ;
        delta_1 = 2 * r_1 * sin( pi / M_1 ) ; 
    else
        r_p = 3.0 * 10^(-2) ;
        M_p =  M_1 ; %4
        phi_p = 2*pi / M_p ;
    end
    
    if strcmp(central_sensor, 'y') % a central sensor exists
        M_all(p+1) = M_p ;
    else
        M_all(p) = M_p ;
    end

    %lambda_min = c ; alpha = 1 ; %%
    hat_r_p = r_p / lambda_min ;

    m = [0 : M_p-1]' ;
    phi_p_m = phi_p * m ; % angles for all sensors for p-th ring

    delay = - hat_r_p * alpha * sin( theta_d*pi/180 ) * cos( phi_d*pi/180 - phi_p_m ) ;

    % steering vector of the UCCA
    [f_mat , delay_mat] = meshgrid(f, delay) ;
    d_p = exp(-1i*2*pi*f_mat.*delay_mat) ;

    if strcmp(central_sensor, 'y') % a central sensor exists
        d{p+1} = d_p ;
    else
        d{p} = d_p ;
    end

end 

end


