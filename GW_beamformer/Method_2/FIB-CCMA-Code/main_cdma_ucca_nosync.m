
% % get some changeable parameters
% Filter_Para.c           = c;
% Filter_Para.P           = P;
% Filter_Para.theta_s     = theta_s;
% Filter_Para.b_vect      = b_vect;
% Filter_Para.RP_vect     = RP_vect;
% Filter_Para.NP_vect     = NP_vect;
% Filter_Para.MP_vect     = MP_vect;
% Filter_Para.phi0_vect   = phi0_vect;
% Filter_Para.theta_d     = theta_d;

j                       = sqrt(-1);
% compute the performance measures in every frequency points
% set some save buffers
% % 1) compute beampattern

BP_matrix           = zeros(length(theta),length(f_vect));
for f_index = 1 : length(f_vect)
    f                   = f_vect(f_index);
    Filter_Para.f       = f;
    h                   = comput_filter_ucca_nosync(Filter_Para);
    omega               = 2*pi*f;

    for index = 1 : size(theta)
        % compute the first ring
        Mp              = MP_vect(1);
        phi0            = phi0_vect(1);
        Mp_vector       = (0 : Mp-1)';
        phi_vector      = Mp_vector*2*pi/Mp + phi0;
        r               = RP_vect(1);
        varpi           = omega*r/c;
        df_theta        = exp(j*varpi*cos(theta(index)-phi_vector));
        
        if(P>1)
        % compute the second and other rings
        for p = 2:P
            Mp          = MP_vect(p);
            phi0        = phi0_vect(p);
            Mp_vector   = (0 : Mp-1)';
            phi_vector  = Mp_vector*2*pi/Mp + phi0;
            r           = RP_vect(p);
            varpi       = omega*r/c;
            df_theta_p  = exp(j*varpi*cos(theta(index)-phi_vector));
            df_theta    = [df_theta; df_theta_p];
        end
        end
        
        BP(index)       = abs(h'*df_theta);
    end
    %BP                  = BP/max(BP);
    BP                  = BP(:);
    BP_matrix(:,f_index) = BP;
end

% % 2) compute DF/WNG

G_wn                    = zeros(length(f_vect),1);
G_dn                    = zeros(length(f_vect),1);
for f_index = 1 : length(f_vect2)
    f                   = f_vect2(f_index);
    Filter_Para.f       = f;
    h                   = comput_filter_ucca_nosync(Filter_Para);
    omega               = 2*pi*f;

    % compute the first ring
        Mp              = MP_vect(1);
        phi0            = phi0_vect(1);
        Mp_vector       = (0 : Mp-1)';
        phi_vector      = Mp_vector*2*pi/Mp + phi0;
        r               = RP_vect(1);
        varpi           = omega*r/c;
        df_theta        = exp(j*varpi*cos(theta_s-phi_vector));
        
        if(P>1)
        % compute the second and other rings
        for p = 2:P
            Mp          = MP_vect(p);
            phi0        = phi0_vect(p);
            Mp_vector   = (0 : Mp-1)';
            phi_vector  = Mp_vector*2*pi/Mp + phi0;
            r           = RP_vect(p);
            varpi       = omega*r/c;
            df_theta_p  = exp(j*varpi*cos(theta_s-phi_vector));
            df_theta    = [df_theta; df_theta_p];
        end
        end
        BP_s       = abs(h'*df_theta);
    
    % % compute the diffuse noise pseudo-coherence matrix
     delta_matrix       = Compute_delta(RP_vect, MP_vect);
     Rv                 = sinc(2*f*delta_matrix/c);
     G_wn(f_index)      = BP_s/(h'*h + eps);
     G_dn(f_index)      = BP_s/(h'*Rv*h + eps);
end
