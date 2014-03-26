%%two dimensional Delta Function ( first order Smereka )

function dApprox = smeeka01( u, dx );
a   =   -1; b =1; c= -1; d = 1;
%C   =   1;
%dx  =   .02;
%dt  =   dx/2;
%xV  = (a-dx:dx:b+dx)';
%y%V  = (c-dx:dx:d+dx)';
%fT  =   (b-a)/dx*dt;
%ddDim =     2;
%epsilon = 1e-12;
%%-----------------------------
%u       = zeros( (b-a)/dx +1,(d-c)/dx +1 );
dApprox = zeros( (b-a)/dx +1,(d-c)/dx +1 );
%for cnt = 1:(d-c)/dx +1  % circle
 %   u( :,cnt ) = (xV(2:end-1).^2 + yV(cnt+1)^2).^.5 - .5;  %dist
%    u( :,cnt ) = (xV(2:end-1).^2 + yV(cnt+1)^2) - .5^2;
%end
u = padarray(u, [1 1], 0 );
u(:,1) = 2*u(:,2) - u(:,3 ) ; u(:,end)= 2*u(:,end-1) -u(:,end-2);
u(1,:) = 2*u(2,:) - u(3, :);  u(end,:)= 2*u(end-1,:) - u(end-2, : );


uxC     =   ( u( 3:end  , 2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
uyC     =   ( u( 2:end-1, 3:end  ) - u(2:end-1, 1:end-2) )/(2*dx);
uxB     =   ( u( 2:end-1, 2:end-1) - u(1:end-2, 2:end-1) )/(dx);
uyB     =   ( u( 2:end-1, 2:end-1) - u(2:end-1, 1:end-2) )/(dx);
uxF     =   ( u( 3:end  , 2:end-1) - u(2:end-1, 2:end-1) )/(dx);
uyF     =   ( u( 2:end-1, 3:end  ) - u(2:end-1, 2:end-1) )/(dx);
nGrad   =   (uxB.^2 +uyB.^2 + eps ).^.5;


%when an x change crosses boundary
deltX   = find( (u( 3:end, 2:end-1) .* u(2:end-1, 2:end-1)) <= 0 );
deltXn  = find( (u( 2:end-1, 2:end-1) .* u(1:end-2, 2:end-1)) <= 0 );
%when a y change crosses boundary
deltY   = find( (u( 2:end-1, 3:end) .* u(2:end-1, 2:end-1)) <= 0 );
deltYn  = find( (u( 2:end-1, 2:end-1) .* u(2:end-1, 1:end-2)) <= 0 );
ut = u(2:end-1, 2:end-1);

dApprox( deltX ) = abs( ut( deltX+1 ).* uxC(deltX )./(uxF(deltX).*nGrad(deltX ))/(dx^2 ) )';
dApprox( deltY ) = dApprox( deltY ) + abs( ut( deltY+size(ut,1) ).* uyC(deltY )./(uyF(deltY).*nGrad(deltY ))/(dx^2 ) );

dApprox( deltXn ) = dApprox( deltXn ) + abs( ut( deltXn-1 ).* uyC(deltXn )./(uyB(deltXn).*nGrad(deltXn ))/(dx^2 ) );
dApprox( deltYn ) = dApprox( deltYn ) + abs( ut( deltYn-size(ut,1)).* uyC(deltYn )./(uyB(deltYn).*nGrad(deltYn ))/(dx^2 ) ); 
                
%figure(1);mesh( dApprox );
%title( [ 'dx =', num2str( dx ), ': ' num2str( sum(sum(dApprox ) )*dx*dx ) ]);
% figure(2);contour( xV, yV,u', [0 0], 'k' );