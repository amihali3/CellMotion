%% solve the PDE for A and M Molecules...  First solve assuming dA/dn =0,
%% later possibly add the CIM code to this... could also just do this by
%% setting everything outside the cell equal to zero when calculating the
%% derivatives if I send in u.

function [A,M ] = amPDE( A, M, dx, dt )
    DA      =   .0764;
    DM      =   .382;
    a       =   .084;
    b       =   1.146;
    e       =   .107;
    c       =   .0764;
    Ca       =   .001;
    Axx     =   ( A( 2:end-1, 3:end) - 2*A(2:end-1, 2:end-1)+ A(2:end-1,1:end-2) )/(dx^2);
    Ayy     =   ( A( 3:end,2:end-1 ) - 2*A(2:end-1,2:end-1) + A(1:end-2,2:end-1) )/(dx^2);
    Mxx     =   ( M( 2:end-1, 3:end) - 2*M(2:end-1, 2:end-1)+ M(2:end-1,1:end-2) )/(dx^2);
    Myy     =   ( M( 3:end,2:end-1 ) - 2*M(2:end-1,2:end-1) + M(1:end-2,2:end-1) )/(dx^2);
    
    A(2:end-1,2:end-1 )  = A(2:end-1,2:end-1 ) + dt * ...
                           ( DA * (Axx + Ayy ) -e*A(2:end-1, 2:end-1) -b*M(2:end-1, 2:end-1).*A(2:end-1, 2:end-1).^2 + Ca );
                       
    M(2:end-1,2:end-1 )  = M(2:end-1,2:end-1 ) + dt * ...
                           ( DM * (Mxx + Myy ) -c*M(2:end-1, 2:end-1) -  a *M(2:end-1, 2:end-1).*A(2:end-1, 2:end-1).^2 + Ca );
    
end