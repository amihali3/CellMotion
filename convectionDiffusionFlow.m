%  convection-diffusion flow
close all;
a   =   -1; b =1; c= -1; d = 1;
C   =   1;
dx  =   .02;
dt  =   dx^2/8;
xV  = (a:dx:b)';
yV  = (c:dx:d)';
fT  =    .009; %.02;
%%-----------------------------
u0  = zeros( (b-a)/dx +1,(d-c)/dx +1 );
for cnt = 1:(d-c)/dx +1
   u0( :,cnt ) = min((xV-.5).^2 + yV(cnt)^2 - .25^2, ...
                min( (xV+.5).^2 + yV(cnt)^2 - .25^2,...
                max(abs(xV)-.3 ,abs(yV(cnt))-.1 )));
end
u = padarray(u0, [1 1], 0 );
u(:,1) = u(:,2); u(:,end)= u(:,end-1);
u(1,:) = u(2,:); u(end,:)= u(end-1,:);

for t = dt:dt:fT
    %%first do mean curvature
    ux      =   ( u( 2:end-1, 3:end) - u(2:end-1, 1:end-2) )/(2*dx);
    uy      =   ( u(3:end, 2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
    uxx     =   ( u( 2:end-1, 3:end) - 2*u(2:end-1, 2:end-1) + u(2:end-1, 1:end-2) )/(dx^2);
    uyy     =   ( u(3:end,2:end-1) - 2*u(2:end-1,2:end-1) + u(1:end-2,2:end-1) )/(dx^2);
    uxy     =   (u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2))/(4*dx^2);
    cTerm   =   ux.^2 .* uxx + 2*ux.*uy.*uxy + uy.^2.*uyy;
    gradM2  =   ux.^2 + uy.^2 + eps;   
    %now upwind differencing for convection part
    bwX     =   ( u( 2:end-1, 2:end-1) - u(2:end-1, 1:end-2) )/dx;
    fwX     =   ( u( 2:end-1, 3:end)   - u(2:end-1, 2:end-1) )/dx;
    bwY     =   ( u( 2:end-1, 2:end-1) - u(1:end-2, 2:end-1) )/dx;
    fwY     =   ( u( 3:end,   2:end-1 )- u(2:end-1, 2:end-1) )/dx;
    
    uxUW= ( bwX > 0 & fwX >= 0 ).*bwX +...  % use BW diff
          ( bwX <= 0 & fwX <= 0 ) .*fwX +...  % use FW diff
          ( bwX <= 0 & fwX > 0 ) .*fwX * 0 +...  % use zero, both in wrong dir
          ( bwX > 0 & fwX < 0 ).*max(abs(fwX),bwX);%both good,choose larger
                   
    
    uyUW= ( bwY > 0 & fwY >= 0 ).*bwY +...  % use BW diff
          ( bwY <= 0 & fwY <= 0 ) .*fwY +...  % use FW diff
          ( bwY <= 0 & fwY > 0 ) .*fwY * 0 +...  % use zero, both in wrong dir
          ( bwY > 0 & fwY < 0 ).* max(abs(fwY),bwY);               
    
    %add in mean curvature
    u(2:end-1, 2:end-1)     =...
                u(2:end-1, 2:end-1) +C * dt*(uxx + uyy -cTerm./gradM2) ;
    
%     add in convection
    u(2:end-1, 2:end-1)     =...
                u(2:end-1, 2:end-1) - dt*(uxUW.^2 + uyUW.^2).^.5;     
            
    u(:,1) = u(:,2); u(:,end)= u(:,end-1);
    u(1,:) = u(2,:); u(end,:)= u(end-1,:);
end

u = u(2:end-1, 2:end-1);

figure(2); contour( xV, yV,u0', [0 0], 'k' );grid; hold on;
contour( xV, yV,u', [0 0], 'r' );
xlabel( 'x-axis');ylabel( 'y-axis');
legend( {'Initial Zero Level Set', 'Final Zero LS' } );
title( ['Final time: ' num2str(fT) ] );
hold off; 