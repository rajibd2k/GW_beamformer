%%% Calculate the Psi-Matrix, for each ring, for each frequency, for Uniform Concentric Circular Array (UCCA)

function [ bar_Psi, M_all ] = bar_Psi_Matrix( M_1, central_sensor, P, Delta_r, theta_d, f, c, Ts, N)

FS = 1 / Ts ;
f_max = FS / 2 ;
alpha = FS / f_max ; % = 2
lambda_min = c / f_max ;

% r_1 = round( 100* lambda_min / 4 / sin( pi / M_1 ) ) / 100 ; 
% delta_1 = 2 * r_1 * sin( pi / M_1 ) ; 

num_rings = P + strcmp(central_sensor,'y') ;

M_all = zeros(1, num_rings) ;

bar_psi = cell(1, num_rings) ;

n = [-N : N]' ; % bessel order

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
    
    %lambda_min = c ; alpha = 1 ; %%
    hat_r_p = r_p / lambda_min ;
        
    m = [0 : M_p-1]' ;
    phi_p_m = phi_p * m ; % angles for all sensors for p-th ring

    % Calculate psi matrix for each frequency
    eta = 2*pi*f* hat_r_p * alpha * sin( theta_d*pi/180 ) ;

    bessel_val = zeros(length(n), length(f)) ;
    for idx_bessel_order = 1 : length(n)
        bessel_val( idx_bessel_order , : ) = besselj( n(idx_bessel_order), eta ) ;                
    end
    
    psi_p = n * phi_p_m' ;
    psi_p = exp( -1i * psi_p ) ;
    
    bar_psi_p = zeros( length(n), length(m), length(f) ) ;
    for idx_f = 1 : length(f)
        bar_psi_p(:, :, idx_f) = ( bessel_val(:,idx_f) * ones(1,length(m)) ) .* psi_p ;
    end
    
    if strcmp(central_sensor, 'y') % a central sensor exists
        M_all(p+1) = M_p ;
        bar_psi{p+1} = bar_psi_p ;
    else
        M_all(p) = M_p ;
        bar_psi{p} = bar_psi_p ;
    end

end 

if strcmp(central_sensor,'y')
    
    M_p = 1 ;
    M_all(1) = M_p ;
        
    m = [0 : M_p-1]' ;
    phi_p_m = 0 ; % angles for all sensors for p-th ring

    % Calculate psi matrix for each frequency
    eta = zeros(size(f)) ;
    bessel_val = zeros(length(n), length(f)) ;
    
    tmp = find(n==0) ;
    bessel_val(tmp, : ) = 1 ;

    psi_p = n * phi_p_m' ;
    psi_p = exp( -1i * psi_p ) ;
    
    bar_psi_p = zeros( length(n), length(m), length(f) ) ;
    for idx_f = 1 : length(f)
        bar_psi_p(:, :, idx_f) = ( bessel_val(:,idx_f) * ones(length(m)) ) .* psi_p ;
    end
    
    bar_psi{1} = bar_psi_p ;
 
end

M_tot = sum(M_all) ;

bar_Psi = zeros(length(n), M_tot, length(f)) ;
for p = 1 : num_rings
    
    M = [ (sum(M_all(1:(p-1))) + 1) : sum(M_all(1:p))] ;
    bar_Psi(:,M,:) = bar_psi{p} ;
    
end

bar_Psi = conj( bar_Psi ) ;
        
end


