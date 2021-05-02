%%% BeamWidth for CCA
%%% theta_d, phi_d : DOA of SOI in degrees
%%% theta_d : [0, 90], phi_d : [0, 180]
%%% power_level_diff : 6 dB

function [ b_theta, b_phi, b_theta_bb, b_phi_bb, B_theta, B_phi ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff)

P = length(h) ; % number of rings, including the central sensor
M_all = zeros(P,1) ;
for p = 1 : P
        M_all(p) = size(h{p}, 1) ;
end
M_tot = sum(M_all) ;

num_f = length(f) ;

% re-organize beamformers
h_all = zeros(M_tot, num_f) ;
for p = 1 : length(h) % = P or P+1

    if p == 1
        beg_sensor = 1 ;
    else
        beg_sensor = sum(M_all(1:p-1)) + 1 ;
    end
    end_sensor = sum(M_all(1:p)) ;

    h_all( beg_sensor : end_sensor , : ) = h{p} ;

end

% Elevation
%-----------------------------------------------------------------------------------------
theta_range = [-90:90]' ; % Elevation angle of incidence of SOI
B_theta = zeros( length(theta_range), num_f ) ;  
for idx_theta_range = 1 : length(theta_range)

    theta = theta_range( idx_theta_range ) ;

    % evaluate steering vector for DOA (changing Elevation, Azimuth same as DOA)
    [ d ] = d_CCA( r_p, phi_p_m, theta, phi_d, f, c, Ts ) ;

    % re-organize steering vector
    d_all = zeros(M_tot, num_f) ;
    for p = 1 : P 

        if p == 1
            beg_sensor = 1 ;
        else
            beg_sensor = sum(M_all(1:p-1)) + 1 ;
        end
        end_sensor = sum(M_all(1:p)) ;

        d_all( beg_sensor : end_sensor , : ) = d{p} ;

    end

    B_theta( idx_theta_range, : ) = nansum( conj(h_all).* d_all , 1) ;

end
    
B_theta_dB = 10*log10( abs(B_theta).^2 ) ;



% Azimuth
%-----------------------------------------------------------------------------------------
phi_range = [-180:180]' ; % Azimuth angle of incidence of SOI
B_phi = zeros( length(phi_range), num_f ) ;  
for idx_phi_range = 1 : length(phi_range)

    phi = phi_range( idx_phi_range ) ;

    % evaluate steering vector for DOA (Elevation same as DOA, changing Azimuth)
    [ d ] = d_CCA( r_p, phi_p_m, theta_d, phi, f, c, Ts ) ;

    % re-organize steering vector
    d_all = zeros(M_tot, num_f) ;
    for p = 1 : P

        if p == 1
            beg_sensor = 1 ;
        else
            beg_sensor = sum(M_all(1:p-1)) + 1 ;
        end
        end_sensor = sum(M_all(1:p)) ;

        d_all( beg_sensor : end_sensor , : ) = d{p} ;

    end

    B_phi( idx_phi_range, : ) = nansum( conj(h_all).* d_all , 1) ;

end
    
B_phi_dB = 10*log10( abs(B_phi).^2 ) ;


% Calculate Beamwidths
%-----------------------------------------------------------------------------------------

% beampattern at DOA
[~, idx_theta_DOA] = min( abs( theta_range - theta_d ) ) ;
[~, idx_phi_DOA] = min( abs( phi_range - phi_d ) ) ;

B_theta_d = B_theta_dB(idx_theta_DOA, :) ;
B_phi_d = B_phi_dB(idx_phi_DOA, :) ;

B_theta_cutoff = B_theta_d - power_level_diff ;
B_phi_cutoff = B_phi_d - power_level_diff ;

b_theta = zeros(1, num_f) ;
b_phi = zeros(1, num_f) ;


% Elevation (theta_d : [0, 90])
%--------------------------------------------------------------------------
for idx_f = 1 : num_f
    
%     figure(1);
%     polarplot( theta_range*pi/180, 10.^((B_theta_dB(:,idx_f))/10) ) ; 
%     title( '$| \mathcal{B}(f,\theta,\phi_{ \mathrm{d} } ) |^2$' ) ;
%     rticks([0.25,1]) ; rticklabels({'-6 dB', '0 dB'}) ;
%     thetaticks( sort(unique([-90,0,theta_d,90]')) ) ; thetaticklabels( num2str(wrapTo180( sort(unique([-90,0,theta_d,90]')) )) ) ;
%     b=gca;
%     set(b, 'ThetaLim', [-90,90], 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top') ;
%     set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
%     a=findobj(gcf); % get the handles associated with the current figure
%     alllines=findall(a,'Type','line');
%     alltext=findall(a,'Type','text');
%     set(alllines,'Linewidth',2, 'MarkerSize', 10);
%     set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    tmp = find( B_theta_dB(:,idx_f) <= B_theta_cutoff(:,idx_f) ) ;
    if isempty(tmp)
        b_theta(idx_f) = 360 ;
        continue ;
    end
    tmp = theta_range(tmp) ;
    
    tmp_minus = tmp( find(tmp < theta_d ) ) ;
    % isempty(tmp_minus) % condition SHOULD never occur for theta_d > 0
    % special misdirected case - bad beamformer or precision error
    if isempty(tmp_minus) 
        tmp_minus = -180 ;
    else
        tmp_minus = max( tmp_minus ) ;
    end
    
%     % special misdirected case - bad beamformer
%     if isempty(tmp_minus) 
%         tmp_minus = theta_d ;
%     end
    
    tmp_plus = tmp( find(tmp > theta_d ) ) ;
    if isempty(tmp_plus)
        tmp_plus = 90 ;
        b_theta(idx_f) = 2 * abs(tmp_plus - tmp_minus) ;
        continue ;
    else
        tmp_plus = min( tmp_plus ) ;
    end

    b_theta(idx_f) = tmp_plus - tmp_minus ;
    
end


% Azimuth
%--------------------------------------------------------------------------      
for idx_f = 1 : num_f
    
%     figure(1);
%     polarplot( phi_range*pi/180, 10.^((B_phi_dB(:,idx_f))/10) ) ; 
%     title( '$| \mathcal{B}(f,\theta_{ \mathrm{d} },\phi ) |^2$' ) ;
%     rticks([0.25,1]) ; rticklabels({'-6 dB', '0 dB'}) ;
%     thetaticks( sort(unique([0,theta_d,90,180,270]')) ) ; thetaticklabels( num2str(wrapTo180( sort(unique([0,theta_d,90,180,270]')) )) ) ;
%     b=gca;
%     set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
%     a=findobj(gcf); % get the handles associated with the current figure
%     alllines=findall(a,'Type','line');
%     alltext=findall(a,'Type','text');
%     set(alllines,'Linewidth',2, 'MarkerSize', 10);
%     set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
    
    tmp = find( B_phi_dB(:,idx_f) <= B_phi_cutoff(:,idx_f) ) ;
    if isempty(tmp)
        b_phi(idx_f) = 360 ;
        continue ;
    end
    tmp = phi_range(tmp) ;
    
    tmp_minus = tmp( find(tmp < phi_d ) ) ;
    % isempty(tmp_minus) % condition SHOULD never occur for phi_d > 0
    % may occur due to precision error
    tmp_minus = max( tmp_minus ) ;
    
    if isempty(tmp_minus) 
        tmp_plus = tmp( find(tmp > phi_d ) ) ;
        tmp_plus = min( tmp_plus ) ;
        b_phi(idx_f) = 2 * abs( phi_d - tmp_plus ) ;
    else
        b_phi(idx_f) = 2 * ( phi_d - tmp_minus ) ;
    end

end


f = linspace(0,0.5, num_f)' ;
Delta_f = mean( abs( diff(f) ) ) ;

b_theta_bb = nanmean( b_theta ) ;
b_phi_bb = nanmean( b_phi ) ; 

end
