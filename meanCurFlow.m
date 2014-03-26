%% Mean curvature flow: %see class 12 notes
a   =   -1; b =1; c= -1; d = 1;
C   =   1;
dx  =   .02;
dt  =   dx^2/8;
xV  = (a:dx:b)';
yV  = (c:dx:d)';
fT  =  .04;
pTime = fT/4;
crs =   [ 'm', 'r', 'b', 'c', 'k', 'y', 'm', 'y' ];
lgnd =  { 't = 0' };
%%-----------------------------
u0  = zeros( (b-a)/dx +1,(d-c)/dx +1 );

for cnt = 1:(d-c)/dx +1
   u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .5 );
%    u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .5^2;
end

 u = padarray(u0, [1 1], 0 );

% u(:,1) = u(:,2); u(:,end)= u(:,end-1);
%  u(1,:) = u(2,:); u(end,:)= u(end-1,:);
u(2:end-1,1) = 2*u(2:end-1,2)- u(2:end-1,3); u(2:end-1,end)= 2*u(2:end-1,end-1) - u(2:end-1,end-2);
u(1,2:end-1) = 2*u(2,2:end-1) - u(3,2:end-1); u(end,2:end-1)= 2*u(end-1,2:end-1) - u(end-2,2:end-1);
u(end, end) = 2*u(end-1, end-1) - u(end-2,end-2 );
u(1, end) = 2*u(2, end-1) - u(3,end-2 );
u(end,1) = 2*u(end-1, 2) - u(end-2,3 );
u(1, 1) = 2*u(2, 2) - u(3,3 );


figure(2); contour( xV, yV,u0', [0 0], 'k' );grid; hold on;

for t = dt:dt:fT
    ux      =   ( u( 3:end,  2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
    uy      =   ( u( 2:end-1,3:end)   - u(2:end-1, 1:end-2) )/(2*dx);
    uxx     =   ( u( 3:end,  2:end-1 )- 2*u(2:end-1,2:end-1) + u(1:end-2, 2:end-1) )/(dx^2);
    uyy     =   ( u( 2:end-1,3:end)   - 2*u(2:end-1,2:end-1) + u(2:end-1,1:end-2) )/(dx^2);
    uxy     =   (u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2))/(4*dx^2);
    cTerm   =   (ux.^2 .* uxx) + 2*( ux.*uy.*uxy ) + (uy.^2 .* uyy );
    gradM2  =   ux.^2 + uy.^2 + 1e-15; 
    toAdd   =   (uxx + uyy - cTerm./gradM2);
    u(2:end-1, 2:end-1)     =...
                u(2:end-1, 2:end-1) +C * dt*toAdd;  %THIS WAS POST HW FIX
    
    u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
    u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);
%     u(:,1) = u(:,2); u(:,end)= u(:,end-1);
%     u(1,:) = u(2,:); u(end,:)= u(end-1,:);
u(end, end) = 2*u(end-1, end-1) - u(end-2,end-2 );
u(1, end) = 2*u(2, end-1) - u(3,end-2 );
u(end,1) = 2*u(end-1, 2) - u(end-2,3 );
u(1, 1) = 2*u(2, 2) - u(3,3 );
  
    if (mod( t, pTime) == 0 )
        contour( xV, yV,u(2:end-1, 2:end-1)', [0 0], crs( round(t/pTime ))  );
        lgnd = [ lgnd; ['t = ',num2str( t ) ]];
    end

end
legend( lgnd );title( 'Zero Level Sets' );
hold off;
u = u(2:end-1, 2:end-1);

figure(1); mesh(xV, yV,u0');
title( 'u(x,0)' );xlabel( 'x-axis');ylabel( 'y-axis');

figure(3); mesh(xV, yV,u');
title( [ 'u(x, ', num2str( fT ), ' )' ] );xlabel( 'x-axis');ylabel( 'y-axis');