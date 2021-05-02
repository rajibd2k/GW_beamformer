function h_filter = comput_filter_ucca_nosync(Filter_Para)
% compute the fixed CDMA filter with my Jacob's version

f           = Filter_Para.f;
c           = Filter_Para.c;
P           = Filter_Para.P;
theta_s     = Filter_Para.theta_s;
b_vect      = Filter_Para.b_vect;
RP_vect     = Filter_Para.RP_vect;
NP_vect     = Filter_Para.NP_vect;
MP_vect     = Filter_Para.MP_vect;
phi0_vect   = Filter_Para.phi0_vect;

theta_d     = Filter_Para.theta_d; 

j           =  sqrt(-1);
omega       = 2*pi*f;

N           = NP_vect(1);
M_sum       = sum(MP_vect);

B_matrix    = zeros(2*N+1,M_sum);
Jn_vect     = zeros(2*N+1,1);
Rs_vect     = zeros(2*N+1,1);

APN         = zeros(P,2*N+1);
for p = 1:P
    Np      = NP_vect(p);
    APN(p,N+1-Np:N+1+Np) = 1;
end

for n = -N:N
    Mp                  = MP_vect(1);
    phi0                = phi0_vect(1);
    Mp_vector           = (0 : Mp-1)';
    phi_vector          = Mp_vector*2*pi/Mp + phi0;
    r                   = RP_vect(1);
    varpi               = omega*r/c * sin(theta_d);
    Phi_vect            = APN(1,n+N+1)*besselj(n,varpi)*exp(-j*n*phi_vector);
    
    if(P>1)
    for p = 2:P
        Mp              = MP_vect(p);
        phi0            = phi0_vect(p);
        Mp_vector       = (0 : Mp-1)';
        phi_vector      = Mp_vector*2*pi/Mp + phi0;
        r               = RP_vect(p);
        varpi           = omega*r/c * sin(theta_d);
        AA              = APN(p,n+N+1)*besselj(n,varpi)*exp(-j*n*phi_vector);
        Phi_vect        = [Phi_vect;AA];
    end
    end
    
    index               = n+N+1;
    B_matrix(index,:)   = Phi_vect';
    Jn_vect(index)      = 1/(j^n);
    Rs_vect(index)      = exp(-j*n*theta_s);
end

J_matrix                = diag(Jn_vect);
Rs_matrix               = diag(Rs_vect);

h_filter = B_matrix'*inv(B_matrix*B_matrix')*conj(J_matrix)*conj(Rs_matrix)*b_vect;

end

