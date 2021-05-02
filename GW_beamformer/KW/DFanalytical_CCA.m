%%% Find the Directivity Factor (DF) for each ring of Concentric Circular Array (UCCA)
%%% f : digital frequency vector, ranging betweeen (-0.5 0.5]
%%% c : speed of plane wave propagation (340 m/s)
%%% Fs : sampling period of the signals (8 kHz)
%%% d : Steering vector of the SOI at each frequency, for all concentric rings
%%% h : filters at different frequencies, for all concentric rings
%%% Gamma_distance : Matrix of Sensor-Distances

function [ D, D_bb, D_bb_low, D_bb_high, D_bb_ratio ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high )

    P = length(d) ;
    M_all = zeros(P,1) ;
    for p = 1 : P
        M_all(p) = size(d{p}, 1) ;
    end
    M_tot = sum(M_all) ;
    
    h_all = zeros(M_tot, length(f)) ;
    d_all = h_all ;

    FS = 1 / Ts ;
    
    % For each ring
    for p = 1 : P      
        
        % re-organize data
        if p == 1
            beg_sensor = 1 ;
        else
            beg_sensor = sum(M_all(1:p-1)) + 1 ;
        end
        end_sensor = sum(M_all(1:p)) ;

        d_all( beg_sensor : end_sensor , : ) = d{p} ;
        h_all( beg_sensor : end_sensor , : ) = h{p} ;

    end
    
    D = zeros( 1, length(f) ) ;
    num_bb = D ; den_bb = D ;
    for idx_f = 1 : length(f)

        d_f = d_all(:, idx_f) ;
        h_f = h_all(:, idx_f) ;

        Gamma_f = sinc( 2*f(idx_f)*(FS/c)*Gamma_distance ) ; % in matlab sinc(x) = sin(pi * x) / (pi * x)

        num = abs( h_f' * d_f )^2 ;
        den = abs( h_f' * Gamma_f * h_f ) ;
        
        num_bb(idx_f) = num ;
        den_bb(idx_f) = den ;
        
        % Narrowband Directivity 
        D(idx_f) = num / den ;
            
    end   
    
    % Broadband Directivity 
    D_bb = nanmean(num_bb) / nanmean(den_bb) ;
    
    % Broadband Directivity at low frequency
    idx_f_low = find( f*FS <= F_low ) ;
    D_bb_low = nanmean( num_bb(idx_f_low) ) / nanmean( den_bb(idx_f_low) ) ;
    
    % Broadband Directivity at high frequency
    idx_f_high = find( f*FS >= F_high ) ;
    D_bb_high = nanmean( num_bb(idx_f_high) ) / nanmean( den_bb(idx_f_high) ) ;
    
    % Broadband Directivity Ratio
    D_bb_ratio = D_bb_low / D_bb ;
    
end


