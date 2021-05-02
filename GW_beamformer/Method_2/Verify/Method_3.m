%%% Constant Beamwidth Beamforming for Uniform Concentric Circular Array (UCCA)
%%% "On the Design of Robust Steerable Frequency-Invariant Beampatterns with Concentric Circular Microphone Arrays"
%%% Optimize Azimuth only

%%% h : final beamformer, for all frequencies, for all the concentric rings
%%% M_all : Number of sensors, in each ring (including the central sensor)
%%% r_1 : radius of first ring (in cm)
%%% theta_range : [0:180]' (degrees)
%%% phi_range : [-180:180]' (degrees)

%%% Delta_b_phi : absolute beamwidth-Azimuth deviation (degree) of final beamformer, for all frequencies
%%% D : directivity of final beamformer, for all frequencies, for all iterations

%%% M_1 : number of sensors in first ring 
%%% central_sensor : whether sensor exists in center of UCCA (y/n)
%%% P : number of rings
%%% Delta_r : separation between rings
%%% theta_d : DOA-Elevation of SOI (in degrees)
%%% phi_d : DOA-Azimuth of SOI (in degrees)
%%% f : normalized frequencies [0, 0.5]
%%% c : speed of sound propagation (340 m/s)
%%% Ts : sampling period (= 1/Fs, Fs = 8000 Hz)
%%% N : order of Freq-invariant desired filter / bessel order

%%% phi_BW : desired beamwidth-Azimuth (in degrees)

function [ h, M_all, theta_range, phi_range] = Method_3( phi_BW, M_1, central_sensor, P, Delta_r, theta_d, phi_d, f, c, Ts, N )

theta_range = [0:180]' ; % degrees
phi_range = [-180:180]' ; % degrees

n = [-N:N]' ; % N is order of Freq-invariant desired filter / bessel order

% J
J = (1i).^(-n) ;
J = diag(J) ;

% Gamma
Gamma = exp(-1i*n*(phi_d*pi/180) ) ;
Gamma = diag(Gamma) ;

% Spatial filter
b_2N = ones(3,1) / 3 ;

% Calculate bar_Psi_Matrix at each frequency of each ring
[ bar_Psi, M_all ] = bar_Psi_Matrix( M_1, central_sensor, P, Delta_r, theta_d, f, c, Ts, N) ;

h_mat = zeros(sum(M_all), length(f)) ;
for idx_f = 1 : length(f)
    
    bar_Psi_f = bar_Psi(:,:,idx_f) ;
    h_f = bar_Psi_f' * pinv( bar_Psi_f * bar_Psi_f' ) * conj(J) * conj(Gamma) * b_2N ; 
    h_mat(:,idx_f) = h_f ;
    
end

h_all = sum( abs(h_mat), 1 ) ;
% normalizing 
h_mat = h_mat ./ (   ones( size(h_mat,1), 1 ) * h_all ) ;

% re-organizing the beamformers for each ring
h = cell(1, length(M_all)) ;
for p = 1 : length(h)
    
    M = [ (sum(M_all(1:(p-1))) + 1) : sum(M_all(1:p))]' ;
    h{p} = h_mat(M,:) ;
end

end


