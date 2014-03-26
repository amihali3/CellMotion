%%CIM1_2d

function main()
close all;
a = 0; b = 1; c = 0; d = 1;
h = .25;

u = zeros( (b-a)/h + 1, (d-c)/h + 1 );
u(:,1)          = 0; % bottom
u(:,end)        = 0:h:(b-a); %top
u(1,2:end-1)    = 0; % left
u(end,2:end-1)  = (h:h:(d-c-h)); %right;


eGrid  = zeros( (b-a)/(h/2) + 3,(d-c)/(h/2) +3 );% twice as refined as the grid we solve for u
ex     =    a-h/2:h/2:b+h/2;
ey     =    c-h/2:h/2:d+h/2;

for cnt = 1:size(ey,2)
       eGrid(:, cnt ) = 0*(ex + ey(cnt) ) + ex;
end

n = (b -a)/h -1 ;
A   = -(1/h^2) *aMat(n, eGrid);   
f   = fMat( a, b, c, d, h );
f   = reshape( f(2:end-1, 2:end -1 ), numel(f(2:end-1, 2:end -1 )),1);
%%adj for known first and last row
f(1:n)  = f(1:n ) +(1/h^2)*u(2:end-1, 1) .*eGrid( 4:2:end-3,3);  %check
f( end-n+1:end ) = f( end-n+1:end ) + (1/h^2) * u(2:end-1,end) .* eGrid( 4:2:end-3,end-2); %check
%%adj for known first and last col
for cnt2 = 1:n
   ind1 = (cnt2-1) * n + 1;
   ind2 = ind1 + n - 1;
   f(ind1)    = f(ind1) + (1/h^2)* u(1, cnt2+1 ) * eGrid(3,(cnt2+1)*2);
   f(ind2)    = f(ind2) + (1/h^2)* u(end,cnt2+1) * eGrid(end-2,(cnt2+1)*2); 
end

uTemp = reshape( A\f, n,n );
u(2:end-1, 2:end-1) = uTemp;

figure(1);mesh(u);
figure(2);mesh( eGrid)

end


function A = aMat( n, eGrid )  
    A   = -(diag( reshape( ( eGrid( 3:2:end-4, 4:2:end-3 ) + eGrid( 5:2:end-2, 4:2:end-3) +...
            eGrid( 4:2:end-3, 3:2:end-4 ) + eGrid( 4:2:end-3, 5:2:end-2 )), n^2, 1 ) )); %check
    D2  =  diag( reshape( [eGrid( 5:2:end-4,4:2:end-3);zeros(1,n)] , n*n,1 ) );
    
    D3  =  diag( reshape( eGrid(4:2:end-3, 5:2:end-4)  , n * (n-1),1) );    
    
    A(1:end-1,2:end) = A(1:end-1,2:end)+ D2(1:end-1, 1:end-1 ); 
    A(2:end,1:end-1) = A(2:end,1:end-1)+ D2(1:end-1, 1:end-1 );
    
    A(1:end-n, n+1:end ) = A(1:end-n, n+1:end ) + D3;
    A(n+1:end, 1:end-n ) = A(n+1:end, 1:end-n ) + D3;
end

function f = fMat( a, b, c, d, h );
    f   =   zeros( (b-a)/h +1, (c-d)/h +1 );
    x   =   a:h:b;
    y   =   c:h:d;    
    for cnt = 1:size( y, 2);
        f(:, cnt ) = -y(cnt)+ 0*(x + y(cnt ));
    end
end