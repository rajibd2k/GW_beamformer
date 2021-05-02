%%% Delay-and-Sum (DS) beamformer for Concentric Circular Array (UCCA)
%%% d : steering vector of the SOI

function [ h ] = DS_CCA( d )
    
P = size(d,2) ;
h = cell(1,P) ;

M_all = zeros(P,1) ;
for p = 1 : P
    M_all(p) = size(d{p}, 1) ;
end
M_tot = sum(M_all) ;

for p = 1 : P
    h{p} = d{p} / M_tot ;
end
    
end


