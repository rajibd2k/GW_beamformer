%%% Calculate the T-Matrix, for each ring, for each frequency, for Uniform Concentric Circular Array (UCCA)

function [ T ] = T_Matrix( r_p, phi_p_m, theta_d, f, c, Ts, bessel_epsilon )

FS = 1 / Ts ;
f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;

% r_1 = round( 100* lambda_min / 4 / sin( pi / M_1 ) ) / 100 ; 
% delta_1 = 2 * r_1 * sin( pi / M_1 ) ; 

P = length(r_p) ;
d = cell(1, P) ;
T = cell(1, P) ;

for p = 1 : P

    ring_radius = r_p(p) ;
    sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
    M_p = length(sensors_angles) ;
    
    if ring_radius == 0
        continue ;
    end
    
    hat_r_p = ring_radius / lambda_min ;

    % Calculate T matrix for each frequency
    eta = 2*pi*f* hat_r_p * alpha * sin( theta_d*pi/180 ) ;
    bessel_order = [-30 : 30]' ;
    bessel_val = zeros(length(bessel_order), length(f)) ;
    for idx_bessel_order = 1 : length(bessel_order)
        bessel_val( idx_bessel_order , : ) = besselj( bessel_order(idx_bessel_order), eta ) ;                
    end
    
    N_p = find( bessel_val(:,end) > bessel_epsilon ) ;
    N_p = max(N_p) - bessel_order(end) ;
    n = [-N_p : N_p]' ;
    
    beg_index = find( bessel_order == -N_p ) ;
    end_index = find( bessel_order == N_p ) ;
    bessel_val = bessel_val( [beg_index : end_index]' , :) ;
    
    [N_p_mat, sensors_angles_mat] = meshgrid(n, sensors_angles) ;
    
    T_p = zeros(M_p, 2*N_p+1, length(f)) ;
    for idx_f = 1 : length(f)
       
        tmp1 = bessel_val( : , idx_f) ;
        tmp1 = ones(length(sensors_angles),1) * tmp1' ;
        tmp1 = ( (1i).^N_p_mat ).* tmp1 ;
        
        tmp2 = -1i* N_p_mat.* sensors_angles_mat ;
        tmp2 = exp( tmp2 ) ;
        
        T_f = tmp1.*tmp2 ;
        T_p(:,:,idx_f) = T_f ;
    end
    
    T{p} = T_p ;

end 
        
end


