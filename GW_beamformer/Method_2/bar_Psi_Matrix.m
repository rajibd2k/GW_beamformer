%%% Calculate the Psi-Matrix, for each ring, for each frequency, for Uniform Concentric Circular Array (UCCA)

function [ bar_Psi, M_all ] = bar_Psi_Matrix( r_p, phi_p_m, theta_d, f, c, Ts, N)

FS = 1 / Ts ;
f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;

num_rings = length(r_p) ;

M_all = zeros(1, num_rings) ;

bar_psi = cell(1, num_rings) ;

n = [-N : N]' ; % bessel order

for p = 1 : num_rings
    
    ring_radius = r_p(p) ;
    sensors_angles = phi_p_m{p} ; % angles for all sensors for p-th ring
    M_p = length( sensors_angles ) ;
    M_all(p) = M_p ;
    
    if ring_radius == 0

        m = 0 ;
        sensors_angles = 0 ; % angles for all sensors for p-th ring

        % Calculate T matrix for each frequency
        eta = zeros(size(f)) ;
        bessel_val = zeros(length(n), length(f)) ;

        tmp = find(n==0) ;
        bessel_val(tmp, : ) = 1 ;

        psi_p = n * sensors_angles' ;
        psi_p = exp( -1i * psi_p ) ;

        bar_psi_p = zeros( length(n), length(m), length(f) ) ;
        for idx_f = 1 : length(f)
            bar_psi_p(:, :, idx_f) = ( bessel_val(:,idx_f) * ones(length(m)) ) .* psi_p ;
        end

        bar_psi{p} = bar_psi_p ;
        
       continue ; 
       
    end
    

    hat_r_p = ring_radius / lambda_min ;
    K_p = M_p / 2 ;
    m = [(-K_p+1) : K_p]' ;

    % Calculate T matrix for each frequency
    eta = 2*pi*f* hat_r_p * alpha * sin( theta_d*pi/180 ) ;

    bessel_val = zeros(length(n), length(f)) ;
    for idx_bessel_order = 1 : length(n)
        bessel_val( idx_bessel_order , : ) = besselj( n(idx_bessel_order), eta ) ;                
    end
    
    psi_p = n * sensors_angles' ;
    psi_p = exp( -1i * psi_p ) ;
    
    bar_psi_p = zeros( length(n), length(m), length(f) ) ;
    for idx_f = 1 : length(f)
        bar_psi_p(:, :, idx_f) = ( bessel_val(:,idx_f) * ones(1,length(m)) ) .* psi_p ;
    end

    bar_psi{p} = bar_psi_p ;

end 

M_tot = sum(M_all) ;

bar_Psi = zeros(length(n), M_tot, length(f)) ;
for p = 1 : num_rings
    
    M = [ (sum(M_all(1:(p-1))) + 1) : sum(M_all(1:p))] ;
    bar_Psi(:,M,:) = bar_psi{p} ;
    
end

bar_Psi = conj( bar_Psi ) ;
        
end


