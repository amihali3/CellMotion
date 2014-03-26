% This should reach an equilibrium when forces cancel out
% involves mean curvature and chem eq from notes #14
velocityS = 1;
a   =   -1.4; b =1.4; c= -1.4; d = 1.4;
dx  =   .02;
dt  =   dx^2/8;
xV  = (a:dx:b)';
yV  = (c:dx:d)';
fT  =  .0002;
pTime = fT/4;
crs =   [ 'm', 'r', 'b', 'c', 'k', 'y', 'm', 'y' ];
lgnd =  { 't = 0' };
epsilon = dx;
coordAtom = [0, 0 ];

%%-----------------------------
u0  = zeros( (b-a)/dx +1,(d-c)/dx +1 );
r0  =   u0;

for cnt = 1:(d-c)/dx +1
%    u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .5 );  % .5 diamond
%    u0( :,cnt ) = max( abs(xV)-1.5, abs(yV(cnt)) - 1.5 );  %1.5 diamond
 u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .8^2;
%     if(cnt< (b-a)/dx/2 )
%         u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .5 );
%     else
%         u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .5^2;
%     end
r0( :,cnt ) =  sqrt(( xV-coordAtom(1) ).^2 + (yV(cnt) - coordAtom(2))^2 );
     
end
u = padarray(u0, [1 1], 0 );
u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);

rI12m6    =   min((1./r0).^12 - (1./r0).^6, 150 );


figure(2); contour( xV, yV,u0', [0 0], 'k' );grid; hold on;

for t = dt:dt:fT
   
      
    
   %% MEAN CURVATURE SET-UP 
    ux      =   ( u( 3:end,  2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
    uy      =   ( u( 2:end-1,3:end)   - u(2:end-1, 1:end-2) )/(2*dx);
    uxx     =   ( u( 3:end,  2:end-1 )- 2*u(2:end-1,2:end-1) + u(1:end-2, 2:end-1) )/(dx^2);
    uyy     =   ( u( 2:end-1,3:end)   - 2*u(2:end-1,2:end-1) + u(2:end-1,1:end-2) )/(dx^2);
    uxy     =   (u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2))/(4*dx^2);
    cTerm   =   (ux.^2 .* uxx) + 2*( ux.*uy.*uxy ) + (uy.^2 .* uyy );
    gradM2  =   ux.^2 + uy.^2 + eps; 
    mcA     =   (uxx + uyy - cTerm./gradM2);
    
    %%Amount due to radius from the Atom
    bdX     = (u(2:end-1,2:end-1)-u(1:end-2,2:end-1))/dx;
    fdX     = (u(3:end  ,2:end-1)-u(2:end-1,2:end-1))/dx;
    bdY     = (u(2:end-1,2:end-1)-u(2:end-1, 1:end-2  ))/dx;
    fdY     = (u(2:end-1,3:end)-u(2:end-1, 2:end-1))/dx;

    adX     =   (bdX>=0 & fdX>=0).* bdX + ...
                (bdX<0  & fdX<0 ).* fdX + ...
                (bdX>0  & fdX<0 ).* max( abs(bdX), abs(fdX))+...
                (bdX<0  & fdX>0 ).* bdX *0;

    adY     =   (bdY>=0 & fdY>=0).* bdY + ...
                (bdY<0  & fdY<0 ).* fdY + ...
                (bdY>0  & fdY<0 ).* max( abs(bdY), abs(fdY))+...
                (bdY<0  & fdY>0 ).* bdY *0;
    
    if velocityS == 1
        d0      =    sweepDist( u,dx, 2 );
        rI12m6Swept = sweepVelocity( rI12m6, d0(2:end-1,2:end-1), dx,5);
    else 
        rI12m6Swept = rI12m6;
    end
            
    rAmount =   -1*sqrt(adX.^2+adY.^2) .* (rI12m6Swept );
    
    
    %%now step in time
    u(2:end-1, 2:end-1) =   u(2:end-1, 2:end-1) + dt * (mcA +rAmount);  
    
    u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
    u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);

    if (mod( t, pTime) == 0 )
        contour( xV, yV,u(2:end-1, 2:end-1)', [0 0], crs( round(t/pTime ))  );
        lgnd = [ lgnd; ['t = ',num2str( t ) ]];
    end

end
u = u(2:end-1, 2:end-1);
%%use secant method to estimate zero ( valid when symettric
% secCurve =  u( :, (size( u,2)-1)/2 + 1 );
% posToNeg = find( ( secCurve(1:end-1) > 0 ) & ( secCurve(2:end)<0 ) );
% uP  = secCurve( posToNeg );
% uN  = secCurve( posToNeg + 1 );
% xP  = xV( posToNeg );
% xN  = xV( posToNeg + 1 );
% xRad = xP - (uP)* (xN-xP )/(uN-uP);

% legend( lgnd );title( ['Integrated R = ',num2str( sqrt( .5^2 - fT*2) )] );
% xlabel( ['Final Radius from secant method: ', num2str( abs(xRad) ) ] );
hold off;


figure(1); mesh(xV, yV,u0');
title( 'u(x,0)' );xlabel( 'x-axis');ylabel( 'y-axis');

figure(3); mesh(xV, yV,u');
title( [ 'u(x, ', num2str( fT ), ' )' ] );xlabel( 'x-axis');ylabel( 'y-axis');