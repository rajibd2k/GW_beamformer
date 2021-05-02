%%% Modified Kaiser Window for Concentric Circular Array (UCCA)
%%% Optimize Azimuth-Directivity

%%% h : final beamformer, for all frequencies, for all the concentric rings
%%% beta_prime : [0,1] : final kaiser window width parameter (log100 scale), for all frequencies, for all rings 
%%% w : [0,1] : final weight (emphasis) given to each ring (including central sensor), for all frequencies
%%% theta_range : [0:180]' (degrees)
%%% phi_range : [-180:180]' (degrees)

%%% Delta_b_theta : beamwidth-Elevation deviation (degree) of final beamformer, for all frequencies, for all iterations
%%% Delta_b_phi : beamwidth-Azimuth deviation (degree) of final beamformer, for all frequencies, for all iterations
%%% D : directivity of final beamformer, for all frequencies, for all iterations
%%% W : WNG of final beamformer, for all frequencies, for all iterations

%%% M_active : Percentage of sensors, active in each ring (= 100)
%%% r_p : radii of all the rings (including central sensor) 
%%% phi_p_m : sensor positions (angles) for all rings 
%%% theta_d : DOA-Elevation of SOI (in degrees)
%%% phi_d : DOA-Azimuth of SOI (in degrees)
%%% f : normalized frequencies [0, 0.5]
%%% c : speed of sound propagation (340 m/s)
%%% Ts : sampling period (= 1/Fs, Fs = 8000 Hz)

%%% GRADIENT DESCENT OPTIMIZATION
%%% init_beta_prime : [0,1] : initialization of kaiser window width parameter (log100 scale), for all frequencies, for all rings 
%%% init_w : [0,1] : initialization of weight (emphasis) given to each ring (including central sensor), for all frequencies
%%% num_iterations : number of iterations 
%%% gradient_stepsize : stepsize (= 0.2) for calculating the gradient 
%%% mu_beta_prime : stepsize (= 0.05) for modifying kaiser window width parameter (log100 scale), for all frequencies 
%%% mu_w : stepsize (= 0.05) for modifying weight given to each ring, for all frequencies 
%%% theta_BW : desired beamwidth-Elevation (in degrees)
%%% phi_BW : desired beamwidth-Azimuth (in degrees)
%%% power_level_diff : power level (in dB) for calculating beamwidth (= 6)

function [ h, beta_prime, w, theta_range, phi_range, Delta_b_theta, Delta_b_phi, D, W] = GradDes_Kaiser_CCA( M_active, init_beta_prime, init_w, num_iterations, gradient_stepsize, mu_beta_prime, mu_w, theta_BW, phi_BW, power_level_diff, r_p, phi_p_m, theta_d, phi_d, f, c, Ts )

theta_range = [0:180]' ; % degrees
phi_range = [-180:180]' ; % degrees

P = length(r_p) ;
[ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;

[ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m ) ;

% Searching for beta_kaiser and ring_weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modified Kaiser Window
%-------------------------------------------------------------------------------------------

% num_iterations = 51 ; % for gradient descent
Delta_b_theta = zeros(num_iterations,length(f)) ;
Delta_b_phi = zeros(num_iterations,length(f)) ;
D = zeros(num_iterations,length(f)) ;
W = zeros(num_iterations,length(f)) ;

% initialize 
% --------------------------------------------------------------
beta_prime = init_beta_prime ;
beta = cell(1, P) ;
for p = 1 : P
    beta{p} = 100.^( beta_prime{p} ) ;
end

w = init_w ;

[ h ] = Modified_Bidir_Kaiser_CCA( M_active, beta, w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
% Metrics
F_low = -1000 ; F_high = -1000 ; % Hz (immaterial)
[ D(1,:) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;

[ b_theta, b_phi ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
Delta_b_theta(1,:) = b_theta - theta_BW ;
Delta_b_phi(1,:) = b_phi - phi_BW ;

[ W(1,:) ] = WNG_CCA(h, d) ;

% What parameter to optimize ?
% --------------------------------------------------------------------------------------------
theta_active = ( Delta_b_theta(1,:) < 0 ) + 0 ;
phi_active = ( Delta_b_phi(1,:) < 0 ) + 0 ;
D_active = not( theta_active + phi_active ) + 0 ; 

% gradient descent search
% --------------------------------------------------------------   
for iteration = 2 : num_iterations
    
    prev_beta_prime = beta_prime ;
    prev_beta = beta ;
    prev_w = w ; 
    prev_h = h ;
    
    prev_theta_active = theta_active ;
    prev_phi_active = phi_active ;
    prev_D_active = D_active ;

    
    % gradient of kaiser window parameter
    % --------------------------------------------------------------------------------------------
    gradient_beta_prime = cell(1,P) ;
    for p = 1:P
        
        if r_p(p) == 0
            gradient_beta_prime{p} = zeros(1, length(f) ) ;
            continue ;
        end

        % backward direction
        %------------------------------------------------------------
        tmp_beta_prime = beta_prime ;
        tmp_beta = beta ;
        tmp_beta_prime{p} = tmp_beta_prime{p} - gradient_stepsize/2 ; 
        %tmp_beta_prime{p} = tmp_beta_prime{p} .* ( tmp_beta_prime{p} > 0 ) ;
        tmp_beta{p} = 100.^(tmp_beta_prime{p}) ; 

        [ tmp_h ] = Modified_Bidir_Kaiser_CCA( M_active, tmp_beta, w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
        % Metrics
        F_low = -1000 ; F_high = -1000 ; % Hz (immaterial)
        [ tmp_D ] = DFanalytical_CCA( d, tmp_h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
        tmp_D = -log10(tmp_D) ;

        [ tmp_b_theta, tmp_b_phi ] = BW_CCA(tmp_h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
        tmp_Delta_b_theta = (tmp_b_theta - theta_BW)/90 ;
        tmp_Delta_b_phi = (tmp_b_phi - phi_BW)/90 ;

        %[ tmp_W ] = WNG_CCA(tmp_h, d) ;

        backward_parameter = theta_active.*(-tmp_Delta_b_theta) + phi_active.*(-tmp_Delta_b_phi) + D_active.*tmp_D ;


        % forward direction
        %------------------------------------------------------------
        tmp_beta_prime = beta_prime ;
        tmp_beta = beta ;
        tmp_beta_prime{p} = tmp_beta_prime{p} + gradient_stepsize/2 ; 
        tmp_beta_prime{p}( find(tmp_beta_prime{p} > 1) ) = 1 ;
        tmp_beta{p} = 100.^(tmp_beta_prime{p}) ; 

        [ tmp_h ] = Modified_Bidir_Kaiser_CCA( M_active, tmp_beta, w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
        % Metrics
        F_low = -1000 ; F_high = -1000 ; % Hz (immaterial)
        [ tmp_D ] = DFanalytical_CCA( d, tmp_h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
        tmp_D = -log10(tmp_D) ;

        [ tmp_b_theta, tmp_b_phi ] = BW_CCA(tmp_h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
        tmp_Delta_b_theta = (tmp_b_theta - theta_BW)/90 ;
        tmp_Delta_b_phi = (tmp_b_phi - phi_BW)/90 ;

        %[ tmp_W ] = WNG_CCA(tmp_h, d) ;

        forward_parameter = theta_active.*(-tmp_Delta_b_theta) + phi_active.*(-tmp_Delta_b_phi) + D_active.*tmp_D ;

        gradient_beta_prime{p} = (forward_parameter - backward_parameter)./gradient_stepsize ;

    end


    % gradient of ring weight parameter
    % --------------------------------------------------------------------------------------------
    gradient_w = cell(1, P) ;
    for p = 1:P

        % backward direction
        %------------------------------------------------------------
        tmp_w = w ;
        tmp_w{p} = tmp_w{p} - gradient_stepsize/2 ;  
        tmp_w{p} = tmp_w{p} .* ( tmp_w{p} > 0 ) ;

        [ tmp_h ] = Modified_Bidir_Kaiser_CCA( M_active, beta, tmp_w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
        % Metrics
        F_low = -1000 ; F_high = -1000 ; % Hz
        [ tmp_D ] = DFanalytical_CCA( d, tmp_h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
        tmp_D = -log10(tmp_D) ;

        [ tmp_b_theta, tmp_b_phi ] = BW_CCA(tmp_h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
        tmp_Delta_b_theta = (tmp_b_theta - theta_BW)/90 ;
        tmp_Delta_b_phi = (tmp_b_phi - phi_BW)/90 ;

        %[ tmp_W ] = WNG_CCA(tmp_h, d) ;

        backward_parameter = theta_active.*(-tmp_Delta_b_theta) + phi_active.*(-tmp_Delta_b_phi) + D_active.*tmp_D ;
        

        % forward direction
        %------------------------------------------------------------
        tmp_w = w ;
        tmp_w{p} = tmp_w{p} + gradient_stepsize/2 ; 

        [ tmp_h ] = Modified_Bidir_Kaiser_CCA( M_active, beta, tmp_w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
        % Metrics
        F_low = -1000 ; F_high = -1000 ; % Hz (immaterial)
        [ tmp_D ] = DFanalytical_CCA( d, tmp_h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
        tmp_D = -log10(tmp_D) ;

        [ tmp_b_theta, tmp_b_phi ] = BW_CCA(tmp_h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
        tmp_Delta_b_theta = (tmp_b_theta - theta_BW)/90 ;
        tmp_Delta_b_phi = (tmp_b_phi - phi_BW)/90 ;

        %[ tmp_W ] = WNG_CCA(tmp_h, d) ;

        forward_parameter = theta_active.*(-tmp_Delta_b_theta) + phi_active.*(-tmp_Delta_b_phi) + D_active.*tmp_D ;

        gradient_w{p} = (forward_parameter - backward_parameter)./gradient_stepsize ;

    end


    % Modify parameters each iteration
    % --------------------------------------------------------------
    regularization = 0.05 ;
    regularization = 1 - regularization ;
    for p = 1:P
%         value = regularization*beta_prime{p} - mu_beta_prime .* gradient_beta_prime{p} ;
%         sign_deviation = sign( value - beta_prime{p} ) ;
%         control_deviation = abs( value - beta_prime{p} ) > 0.2*beta_prime{p} ;
%         indices_deviation = find( control_deviation ) ;
%         tmp_value = beta_prime{p}.*(1 + 0.2*sign_deviation) ;
%         value( indices_deviation ) = tmp_value( indices_deviation ) ;
%         beta_prime{p} = value ;
        
        beta_prime{p} = regularization*beta_prime{p} - mu_beta_prime .* gradient_beta_prime{p} ;
        %beta_prime{p}( find(beta_prime{p} < 0) ) = -inf ;
        beta_prime{p}( find(beta_prime{p} > 1) ) = 1 ;
        beta{p} = 100.^( beta_prime{p} ) ;
    end

    total_w = zeros(1, length(f)) ;
    for p = 1:P
%         value = regularization*w{p} - mu_w .* gradient_w{p} ;
%         sign_deviation = sign( value - w{p} ) ;
%         control_deviation = abs( value - w{p} ) > 0.2*w{p} ;
%         indices_deviation = find( control_deviation ) ;
%         tmp_value = w{p}.*(1 + 0.2*sign_deviation) ;
%         value( indices_deviation ) = tmp_value( indices_deviation ) ;
%         w{p} = value ;        
        
        w{p} = regularization*w{p} - mu_w .* gradient_w{p} ;
        w{p}( find(w{p} < 0) ) = 0 ;
        total_w = total_w + w{p} ;
    end

    for p = 1:length(d)
        w{p} = w{p} ./ total_w ;
    end


    % Metrics after each iteration
    % --------------------------------------------------------------
    [ h ] = Modified_Bidir_Kaiser_CCA( M_active, beta, w, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
    % Metrics
    F_low = -1000 ; F_high = -1000 ; % Hz (immaterial)
    [ D(iteration,:) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;

    [ b_theta, b_phi ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
    Delta_b_theta(iteration,:) = b_theta - theta_BW ;
    Delta_b_phi(iteration,:) = b_phi - phi_BW ;

    [ W(iteration,:) ] = WNG_CCA(h, d) ;
    
   
    % What parameter to optimize in the next iteration ?
    % --------------------------------------------------------------------------------------------
    
    theta_active = ( Delta_b_theta(iteration,:) < 0 ) + 0 ;
    phi_active = ( Delta_b_phi(iteration,:) < 0 ) + 0 ;
    D_active = not( theta_active + phi_active ) + 0 ; 
    

%     % Increase Learning rate and Gradient stepsize when there is no change
%     % (except in Broadside scenario)
%     % --------------------------------------------------------------------------------------------
%     if and( iteration < 3, theta_d > 0 )
%         prev_increase = zeros(1, length(f) ) ;
%         check = zeros(1, length(f) ) ;
%         
%     elseif and( iteration >= 3, theta_d > 0 )
%         
%         check1 = ( sum( abs( diff( Delta_b_theta(iteration-2:iteration,:) ) ) ) < 1 ) ;
%         check2 = (sum( abs( diff( Delta_b_phi(iteration-2:iteration,:) ) ) ) < 1 ) ;
%         check = or(check1, check2) ;
%         check = check.*(not(prev_increase)) ;
%         
%         mu_beta_prime( find(check) ) = max(   min( 2.5*mu_beta_prime( find(check) ), 1.0 ) ,  0.15 ) ;
%         mu_w( find(check) ) = max(   min( 2.5*mu_w( find(check) ), 1.0 ) , 0.15 ) ;
%         gradient_stepsize( find(check) ) = max(   min( 2.5*gradient_stepsize( find(check) ), 1.0 ) , 0.25 ) ;
%         
%         prev_increase = check + 0 ;
%     end  


    % Increase Learning rate and Gradient stepsize when there is no change
    % (except in Broadside scenario)
    % --------------------------------------------------------------------------------------------
    if iteration < 3 
        prev_increase = zeros(1, length(f) ) ;
        check = zeros(1, length(f) ) ;
        
    elseif and( iteration >= 3, theta_d > 0 )
        
        check1 = ( sum( abs( diff( Delta_b_theta(iteration-2:iteration,:) ) ) ) < 1 ) ;
        check2 = (sum( abs( diff( Delta_b_phi(iteration-2:iteration,:) ) ) ) < 1 ) ;
        check = or(check1, check2) ;
        check = check.*(not(prev_increase)) ;
        
        mu_beta_prime( find(check) ) = max(   min( 2.5*mu_beta_prime( find(check) ), 1.0 ) ,  0.15 ) ;
        mu_w( find(check) ) = max(   min( 2.5*mu_w( find(check) ), 1.0 ) , 0.15 ) ;
        gradient_stepsize( find(check) ) = max(   min( 2.5*gradient_stepsize( find(check) ), 1.0 ) , 0.25 ) ;
        
        prev_increase = check + 0 ;
        
    elseif and( iteration >= 3, theta_d == 0 )
        
        check1 = ( sum( abs( diff( Delta_b_theta(iteration-2:iteration,:) ) ) ) < 1 ) ;
        check = check1 ;
        check = check.*(not(prev_increase)) ;
        
        mu_beta_prime( find(check) ) = max(   min( 2.5*mu_beta_prime( find(check) ), 1.0 ) ,  0.15 ) ;
        mu_w( find(check) ) = max(   min( 2.5*mu_w( find(check) ), 1.0 ) , 0.15 ) ;
        gradient_stepsize( find(check) ) = max(   min( 2.5*gradient_stepsize( find(check) ), 1.0 ) , 0.25 ) ;
        
        prev_increase = check + 0 ;
    end 
    
    
    % Decrease Learning rate and Gradient stepsize when parameter fluctuates between BW and DF
    % --------------------------------------------------------------------------------------------
    indices_freqs = ( abs(D_active - prev_D_active) > 0 ) + 0 ;
    
    % dont consider cases which require increasing rate 
    indices_freqs( find(check) ) = 0 ; 

    check1 = ( abs( diff( Delta_b_theta(iteration-1:iteration,:) ) ) >= 20 ) ;
    check2 = ( abs( diff( Delta_b_phi(iteration-1:iteration,:) ) ) >= 20 ) ;
    check = check1 + check2 ;
    indices_freqs = or( indices_freqs, check ) + 0 ;
    
    mu_beta_prime = mu_beta_prime .* (1 - indices_freqs*0.5) ; 
    mu_w = mu_w .* (1 - indices_freqs*0.5) ; 
    gradient_stepsize = gradient_stepsize .* (1 - indices_freqs*0.5) ; 
   
    theta_active( find(indices_freqs) ) = prev_theta_active( find(indices_freqs) ) ;
    phi_active( find(indices_freqs) ) = prev_phi_active( find(indices_freqs) ) ;
    D_active( find(indices_freqs) ) = prev_D_active( find(indices_freqs) ) ;
    
    % special cases 
    if theta_d > 0
        check = or( Delta_b_theta(iteration-1,:) == (360 - theta_BW) , Delta_b_phi(iteration-1,:) == (360 - phi_BW) ) + 0 ;
    elseif theta_d == 0
        check = ( Delta_b_theta(iteration-1,:) == (360 - theta_BW) ) + 0 ;
    end
    
    check = check.*(indices_freqs) ;
    indices_freqs( find(check) ) = 0 ;
    theta_active( find(check) ) = ( Delta_b_theta(iteration,find(check)) < 0 ) + 0 ;
    phi_active( find(check) ) = ( Delta_b_phi(iteration,find(check)) < 0 ) + 0 ;
    D_active( find(check) ) = not(  theta_active( find(check) ) + phi_active( find(check) )  ) + 0 ;

    
   % Revert back parameters 
   % --------------------------------------------------------------
    for p = 1:P
        beta_prime{p}( find(indices_freqs) ) = prev_beta_prime{p}( find(indices_freqs) ) ;
        beta{p}( find(indices_freqs) ) = prev_beta{p}( find(indices_freqs) ) ;
        w{p}( find(indices_freqs) ) = prev_w{p}( find(indices_freqs) ) ; 
        h{p}( : , find(indices_freqs) ) = prev_h{p}( : , find(indices_freqs) ) ;
    end

    Delta_b_theta(iteration, find(indices_freqs) ) = Delta_b_theta(iteration-1, find(indices_freqs) ) ;
    Delta_b_phi(iteration, find(indices_freqs) ) = Delta_b_phi(iteration-1, find(indices_freqs) ) ;
    D(iteration, find(indices_freqs) ) = D(iteration-1, find(indices_freqs) ) ;
    W(iteration, find(indices_freqs) ) = W(iteration-1, find(indices_freqs) ) ;
    
end
        
end


