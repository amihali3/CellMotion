 %% (6/1/11, part 1 )
function main()
    a   =   -1; b =1; c= -1; d = 1;
    numSweeps = 1;
    C   =   1;
    dx  =   .05;
    dt  =   dx^2/4;
    xV  = (a:dx:b)';
    yV  = (c:dx:d)';
    fT  =  .01;
    pTime = fT/4;
    crs =   [ 'm', 'r', 'b', 'c', 'k', 'y', 'm', 'y' ];
    lgnd =  { 't = 0' };
    %%-----------------------------
    u0  = zeros( (b-a)/dx +1,(d-c)/dx +1 );
    for cnt = 1:(d-c)/dx +1
%        u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .5 );
       u0( :,cnt ) = xV.^2 +yV(cnt)^2 - .5^2 ;
    end
    u = padarray(u0, [1 1], 0 );
%  u(:,1) = u(:,2); u(:,end)= u(:,end-1);
%  u(1,:) = u(2,:); u(end,:)= u(end-1,:);
    u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
    u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);

    figure(2); contour( xV, yV,u0', [0 0], 'k' );grid; hold on;

    for t = dt:dt:fT
        ux      =   ( u( 2:end-1, 3:end) - u(2:end-1, 1:end-2) )/(2*dx);
        uy      =   ( u(3:end, 2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
        uxx     =   ( u( 2:end-1, 3:end) - 2*u(2:end-1, 2:end-1) + u(2:end-1, 1:end-2) )/(dx^2);
        uyy     =   ( u(3:end,2:end-1) - 2*u(2:end-1,2:end-1) + u(1:end-2,2:end-1) )/(dx^2);
        uxy     =   ( u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2) )/(4*dx^2);
        cTerm   =   ux.^2 .* uxx + 2*ux.*uy.*uxy + uy.^2.*uyy;
        gradM2  =   ux.^2 + uy.^2 + eps; 
        u(2:end-1, 2:end-1)     =...
                    u(2:end-1, 2:end-1) +C * dt*(uxx + uyy -cTerm./gradM2);
            
        u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
        u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);
        u =     sweepDist( u, dx, numSweeps );
       
        if (mod( t, pTime) == 0 )
            contour( xV, yV,u(2:end-1, 2:end-1)', [0 0], crs( round(t/pTime))  );
            lgnd = [ lgnd; ['t = ',num2str( t ) ]];
        end
    end
    u = u(2:end-1, 2:end-1);
    
      %%use secant method to estimate zero ( valid when symettric
    secCurve =  u( :, (size( u,2)-1)/2 + 1 );
    posToNeg = find( ( secCurve(1:end-1) > 0 ) & ( secCurve(2:end)<0 ) );
    uP  = secCurve( posToNeg );
    uN  = secCurve( posToNeg + 1 );
    xP  = xV( posToNeg );
    xN  = xV( posToNeg + 1 );
    xRad = xP - (uP)* (xN-xP )/(uN-uP);
   
    legend( lgnd );title( ['Integrated R = ',num2str( sqrt( .5^2 - fT*2) )] );
    xlabel( ['Final Radius from secant method: ', num2str( abs(xRad) ) ] );
    hold off;
    
    
    
    figure(1); mesh(xV, yV,u0');
    title( 'u(x,0)' );xlabel( 'x-axis');ylabel( 'y-axis');
    figure(3); mesh(xV, yV,u');
    title( [ 'u(x, ', num2str( fT ), ' )' ] );xlabel( 'x-axis');ylabel( 'y-axis');
end