%%% Calculate White Noise Gain (WNG) for CCA
%%% h : filters at different frequencies, for all concentric rings
%%% d : steering vectors at different frequencies, for all concentric rings

function [W, W_bb] = WNG_CCA(h, d)
    
P = size(d,2) ;
num_f = size( d{1}, 2 ) ;

M_all = zeros(P,1) ;
for p = 1 : P
    M_all(p) = size(d{p}, 1) ;
end
M_tot = sum(M_all) ;

% re-organize steering vectors and beamformers
d_all = zeros(M_tot, num_f) ;
h_all = zeros(M_tot, num_f) ;
for p = 1 : P

    if p == 1
        beg_sensor = 1 ;
    else
        beg_sensor = sum(M_all(1:p-1)) + 1 ;
    end
    end_sensor = sum(M_all(1:p)) ;

    d_all( beg_sensor : end_sensor , : ) = d{p} ;
    h_all( beg_sensor : end_sensor , : ) = h{p} ;

end

tmp_num = abs( sum( conj(h_all) .* d_all , 1) ).^2 ;
tmp_den = sum( conj(h_all) .* h_all , 1) ;
W = abs( tmp_num ./ tmp_den ) ;
W_bb = abs( nanmean(tmp_num) / nanmean(tmp_den) ) ;
    
end

