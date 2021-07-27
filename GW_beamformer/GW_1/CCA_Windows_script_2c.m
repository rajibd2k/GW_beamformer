% % Modified Gaussian Window beamformer for Concentric Circular Array (CCA)
% % Using Gradient Descent Algorithm 
% % Comparing GW-1 vs. KW vs. Method-1 vs. Method-2
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear all ; clc ; close all ;
% 
% design_name = {'GW-1', 'KW', 'Method-1', 'Method-2'} ;
% 
% c = 340 ; Ts = 1/16000 ; FS = 1/Ts ; 
% 
% % % choose frequency
% f = [0 : FS/256 : FS/2]' ; % Hz
% f = f/FS ;
% 
% theta_d = 45 ; %% Elevation - DOA of the SOI between [0,90]
% phi_d = 45 ; % Azimuth - DOA of the SOI between [0,90]
%  
% design = load( 'CCA_design' ) ;
% %design = load( 'CCA_design_Nyquist' ) ;
% 
% r_p = design.r_p ;
% phi_p_m = design.phi_p_m ;
% 
% clear design ;
% 
% active_rings = zeros(1,length(phi_p_m)) ;
% for p = 1 : length(phi_p_m) 
%     active_rings(p) = not( isempty(phi_p_m{p}) ) + 0 ; 
% end
% r_p = r_p( find(active_rings) ) ;
% phi_p_m = phi_p_m( find(active_rings) ) ;
% P = length(r_p) ;
% 
% [ Gamma_distance ] = GammaDistance_CCA( r_p, phi_p_m ) ;
% [ d ] = d_CCA( r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
% 
% % Modified Gaussian / Kaiser Window
% %-------------------------------------------------------------------------------------------
% M_active = 100 ;
% init_beta_prime = cell(1,P) ;
% for p = 1 : P
%     if r_p(p) == 0
%         init_beta_prime{p} = nan*ones(1,length(f)) ;
%     else
%         init_beta_prime{p} = 0.5*ones(1,length(f)) ; % initialize to 1 / 0 (all sensors active initially)
%     end
% end
% 
% init_w = cell(1,P) ;
% for p = 1 : P
%     init_w{p} = ( 1/P ) * ones(1, length(f)) ;
% end
% 
% num_iterations = 41 ;
% gradient_stepsize = 0.2*ones(1,length(f)) ;
% mu_beta_prime = 0.1*ones(1,length(f)) ; mu_w = 0.1*ones(1,length(f)) ;
% theta_BW = 40 ; 
% phi_BW = 40 ;
% power_level_diff = 6 ;
% 
% Delta_b_theta = zeros( num_iterations, length(f), length(design_name) ) ;
% Delta_b_phi = zeros( num_iterations, length(f), length(design_name) ) ;
% D = zeros( num_iterations, length(f), length(design_name) ) ;
% W = zeros( num_iterations, length(f), length(design_name) ) ;
% B = cell(1, length(design_name) ) ;
%     
% for idx_design_name = 1 : length(design_name)
%     
%     if strcmp(design_name{idx_design_name}, 'GW-1')
% 
%         [ h, ~, ~, ~, ~, Delta_b_theta(:,:,idx_design_name), Delta_b_phi(:,:,idx_design_name), D(:,:,idx_design_name), W(:,:,idx_design_name)] = GradDes_Gaussian_CCA( M_active, init_beta_prime, init_w, num_iterations, gradient_stepsize, mu_beta_prime, mu_w, theta_BW, phi_BW, power_level_diff, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
%         [theta_range, phi_range, B{idx_design_name}] = BP_CCA(h, r_p, phi_p_m, f, c, Ts) ;
% 
%     elseif strcmp(design_name{idx_design_name}, 'KW')
%         cd('..'); cd('KW') ;
%         [ h, ~, ~, ~, ~, Delta_b_theta(:,:,idx_design_name), Delta_b_phi(:,:,idx_design_name), D(:,:,idx_design_name), W(:,:,idx_design_name)] = GradDes_Kaiser_CCA( M_active, init_beta_prime, init_w, num_iterations, gradient_stepsize, mu_beta_prime, mu_w, theta_BW, phi_BW, power_level_diff, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
%         cd('..') ; cd('GW_1') ;
%         [theta_range, phi_range, B{idx_design_name}] = BP_CCA(h, r_p, phi_p_m, f, c, Ts) ;
% 
%         
%     elseif strcmp(design_name{idx_design_name}, 'Method-1')
%         
%         addpath('../Method_1') ;
%         bessel_epsilon = 0.01 ;     
%         [ h ] = Method_1( phi_BW, power_level_diff, bessel_epsilon, r_p, phi_p_m, theta_d, phi_d, f, c, Ts ) ;
%         
%         F_low = -1000 ; F_high = -1000 ;
%         [ D(end,:,idx_design_name) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
% 
%         [ b_theta, b_phi ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
%         Delta_b_theta(end,:,idx_design_name) = b_theta - theta_BW ;
%         Delta_b_phi(end,:,idx_design_name) = b_phi - phi_BW ;
% 
%         [ W(end,:,idx_design_name) ] = WNG_CCA(h, d) ;
% 
%         [theta_range, phi_range, B{idx_design_name}] = BP_CCA(h, r_p, phi_p_m, f, c, Ts) ;
%         
%         
%     elseif strcmp(design_name{idx_design_name}, 'Method-2')
%         addpath('../Method_2') ;
%         N = 12 ;
%         [ h ] = Method_2( phi_BW, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, N ) ;
%         
%         F_low = -1000 ; F_high = -1000 ;
%         [ D(end,:,idx_design_name) ] = DFanalytical_CCA( d, h, Gamma_distance, f, c, Ts, F_low, F_high ) ;
% 
%         [ b_theta, b_phi ] = BW_CCA(h, r_p, phi_p_m, theta_d, phi_d, f, c, Ts, power_level_diff) ;
%         Delta_b_theta(end,:,idx_design_name) = b_theta - theta_BW ;
%         Delta_b_phi(end,:,idx_design_name) = b_phi - phi_BW ;
% 
%         [ W(end,:,idx_design_name) ] = WNG_CCA(h, d) ;
% 
%         [theta_range, phi_range, B{idx_design_name}] = BP_CCA(h, r_p, phi_p_m, f, c, Ts) ;
%         
%     end
% 
% end
% 
% Delta_b_theta = reshape( Delta_b_theta(end,:,:), length(f), length(design_name) ) ;
% Delta_b_phi = reshape( Delta_b_phi(end,:,:), length(f), length(design_name) ) ;
% D = reshape( D(end,:,:), length(f), length(design_name) ) ;
% W = reshape( W(end,:,:), length(f), length(design_name) ) ;
% 
% save(['CCA_varying_beamformer'], 'Delta_b_theta', 'Delta_b_phi', 'D', 'W', 'B', 'theta_range', 'phi_range' ) ;
% %save(['CCA_Nyquist_varying_beamformer'], 'Delta_b_theta', 'Delta_b_phi', 'D', 'W', 'B', 'theta_range', 'phi_range' ) ;
% 
% exit ;
% 
% return ;

% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear ; clc ; %close all ;

design_name = {'CCA'; 'CCA-I'} ;

theta_d = 45 ; %% Elevation - DOA of the SOI between [0,90]
phi_d = 45 ; % Azimuth - DOA of the SOI between [0,90]

theta_BW = 40 ; 
phi_BW = 40 ;

beamformer_name = {'GW-1', 'KW', 'Method-1', 'Method-2'} ;

% % choose frequency
Ts = 1/16000 ; FS = 1/Ts ; % Hz
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

choice = 'normal' ;
if strcmp( choice , 'logscale' )
    frequencies = log2(f*FS) / log2(2) ; %log2
    fig_ticks = frequencies([1,4+1, 16+1, 32+1, 64+1, end]) ;
    fig_labels = 2.^( fig_ticks ) / 1000 ;
elseif strcmp( choice , 'normal' )
    frequencies = f*FS ;
    fig_ticks = linspace(0, FS/2, 5)' ;
    fig_labels = fig_ticks/1000 ;
end

figure();
for idx_design = 1 : length(design_name)

    if strcmp( design_name{idx_design} , 'CCA' )
        loadname = 'CCA_Nyquist_varying_beamformer' ;
    else
        loadname = 'CCA_varying_beamformer' ;
    end
    
    load(loadname) ;

    subplot(2,4,1+4*(idx_design-1)) ; 
    values = Delta_b_theta + theta_BW ;
    values = movmean(movmean(movmedian(values,10), 10), 10) ;
    plot(frequencies, values ) ; 
    if idx_design == 1
        title('elevation-beamwidth') ;
    else
        xlabel('$f$ (kHz)') ; 
    end
    ylabel('degrees') ; %axis('tight') ;
    xlim([ min(f*FS)+10, max(f*FS)]) ; 
    ylim([ 0, 360]) ; 
    yticks( sort(  unique( [ 0, theta_BW, round( linspace(min(ylim), max(ylim), 4) ) ] )  ) ) ;
    b=gca;
    set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    subplot(2,4,2+4*(idx_design-1)) ; 
    values = Delta_b_phi + theta_BW ;
    values = movmean(movmean(movmedian(values,10), 10), 10) ;
    plot(frequencies, values ) ; 
        if idx_design == 1
        title('azimuth-beamwidth') ;
    else
        xlabel('$f$ (kHz)') ; 
    end
    ylabel('degrees') ; %axis('tight') ; 
    xlim([ min(f*FS)+10, max(f*FS)]) ; 
    ylim([ 0, 360]) ; 
    yticks( sort(  unique( [ 0, phi_BW, round( linspace(min(ylim), max(ylim), 4) ) ] )  ) ) ;
    b=gca;
    set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    subplot(2,4,3+4*(idx_design-1)) ; 
    values = 10*log10( D ) ;
    values = movmean(movmean(movmedian(values,10), 10), 10) ;
    plot(frequencies, values ) ; 
    if idx_design == 1
        title('DF') ;
    else
        xlabel('$f$ (kHz)') ; 
    end
    ylabel('dB') ; %axis('tight') ;
    xlim([ min(f*FS)+10, max(f*FS)]) ; 
    ylim([ 0, 20]) ; 
    yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
    b=gca;
    set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    subplot(2,4,4+4*(idx_design-1)) ; 
    values = 10*log10( W ) ;
    values = movmean(movmean(movmedian(values,10), 10), 10) ;
    plot(frequencies, values ) ; 
    if idx_design == 1
        title('WNG') ;
    else
        xlabel('$f$ (kHz)') ; 
    end
    ylabel('dB') ; %axis('tight') ;
    xlim([ min(f*FS)+10, max(f*FS)]) ; 
    ylim([ 0, 22]) ; 
    yticks( round( linspace(min(ylim), max(ylim), 5), 1 ) ) ;
    hleg = legend(beamformer_name); 
    title(hleg, [design_name{idx_design},'~$ \theta_{\mathrm{d}} = \phi_{\mathrm{d}} = 45^o$ ~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$~Beamformer'],'Interpreter','Latex');
    b=gca;
    set (b, 'XTick', fig_ticks); set (b, 'XTickLabel', fig_labels );
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

end

return ;




clear ; clc ; %close all ;

loadname = 'CCA_varying_beamformer' ;
%loadname = 'CCA_Nyquist_varying_beamformer' ;

theta_d = 45 ; %% Elevation - DOA of the SOI between [0,90]
phi_d = 45 ; % Azimuth - DOA of the SOI between [0,90]

theta_BW = 40 ; 
phi_BW = 40 ;

load(loadname) ;

design_name = {'GW-1', 'KW', 'Method-1', 'Method-2'} ;

% % choose frequency
Ts = 1/16000 ; FS = 1/Ts ; % Hz
f = [0 : FS/256 : FS/2]' ; % Hz
f = f/FS ;

% Beampattern at 1 kHz
%---------------------------------------------------------------------------------------------

linestyles = {'-', '--', ':', '-.'} ;
[~,freq_index] = min(abs(f*FS - 1000)) ; %%
freq = f(freq_index)*FS/1000 ;
for idx_design = 1 : length(design_name)
    
    values = B{idx_design}(:,:,freq_index) ;
    values = abs( values ) ;
    values = values / max( max(values) ) ;
    values = reshape( values, length(theta_range), length(phi_range) ) ;
    values = values.^2 ;

    % 2-D Elevation Beampattern
    idx_phi_d = find(phi_range == theta_d) ;
    idx_theta_d = find(theta_range == theta_d) ;
    
    figure(2);
    subplot(1, 2, 1) ; 
    polarplot( theta_range*pi/180, values(:,idx_phi_d)/max(values(:,idx_phi_d)), linestyles{idx_design} ) ; hold on ;
    title( '$| \mathcal{B}(f,\theta,\phi_{ \mathrm{d} } ) |^2$' ) ;
    rticks([0.25,0.5,1]) ; rticklabels({'-6 dB', '-3 dB', '0 dB'}) ;
    %thetalim([0, 180]) ; thetaticks([0,45,90,180]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180]')) ) ;
    thetaticks([0,45,90,180,270]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180,270]')) ) ;
    b=gca;
    set(b, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top') ;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
    
    subplot(1, 2, 2) ; 
    polarplot( phi_range*pi/180, values(idx_theta_d,:)'/max(values(idx_theta_d,:)'), linestyles{idx_design} ) ; hold on ;
    title( '$| \mathcal{B}(f,\theta_{ \mathrm{d} },\phi) |^2$' ) ;
    rticks([0.25,0.5,1]) ; rticklabels({'-6 dB', '-3 dB', '0 dB'}) ;
    thetaticks([0,45,90,180,270]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180,270]')) ) ;
    hleg = legend(design_name); 
    title(hleg, ['$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = $',  num2str(phi_d),'$^o$', '~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$,~$fF_s =$', num2str(freq), '~kHz~Beamformer'],'Interpreter','Latex');
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');
    
end


% Beampattern at 6 kHz
%---------------------------------------------------------------------------------------------

linestyles = {'-', '--', ':', '-.'} ;
[~,freq_index] = min(abs(f*FS - 6000)) ; %%
freq = f(freq_index)*FS/1000 ;
for idx_design = 1 : length(design_name)

    values = B{idx_design}(:,:,freq_index) ;
    values = abs( values ) ;
    values = values / max( max(values) ) ;
    values = reshape( values, length(theta_range), length(phi_range) ) ;
    values = values.^2 ;

    % 2-D Azimuthal Beampattern
    idx_phi_d = find(phi_range == theta_d) ;
    idx_theta_d = find(theta_range == theta_d) ;

    figure(3);
    subplot(1, 2, 1) ; 
    polarplot( theta_range*pi/180, values(:,idx_phi_d)/max(values(:,idx_phi_d)), linestyles{idx_design} ) ; hold on ;
    title( '$| \mathcal{B}(f,\theta,\phi_{ \mathrm{d} } ) |^2$' ) ;
    rticks([0.25,0.5,1]) ; rticklabels({'-6 dB', '-3 dB', '0 dB'}) ;
    %thetalim([0, 180]) ; thetaticks([0,45,90,180]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180]')) ) ;
    thetaticks([0,45,90,180,270]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180,270]')) ) ;
    b=gca;
    set(b, 'ThetaDir', 'clockwise', 'ThetaZeroLocation', 'top') ;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

    subplot(1, 2, 2) ; 
    polarplot( phi_range*pi/180, values(idx_theta_d,:)'/max(values(idx_theta_d,:)'), linestyles{idx_design} ) ; hold on ;
    title( '$| \mathcal{B}(f,\theta_{ \mathrm{d} },\phi) |^2$' ) ;
    rticks([0.25,0.5,1]) ; rticklabels({'-6 dB', '-3 dB', '0 dB'}) ;
    thetaticks([0,45,90,180,270]') ; thetaticklabels( num2str(wrapTo180([0,45,90,180,270]')) ) ;
    hleg = legend(design_name); 
    title(hleg, ['$\theta_{\mathrm{d}} = \phi_{\mathrm{d}} = $',  num2str(phi_d),'$^o$', '~$\theta_{\mathrm{BW}} = \phi_{\mathrm{BW}} = 40^o$,~$fF_s =$', num2str(freq), '~kHz~Beamformer'],'Interpreter','Latex');
    b=gca;
    set(b,'FontName','Times New Roman','FontWeight','Bold','LineWidth',0.5,'FontSize',16);
    a=findobj(gcf); % get the handles associated with the current figure
    alllines=findall(a,'Type','line');
    alltext=findall(a,'Type','text');
    set(alllines,'Linewidth',2, 'MarkerSize', 10);
    set(alltext,'FontName','Times New Roman','FontWeight','Normal','FontSize',20,'Interpreter','Latex');

end

