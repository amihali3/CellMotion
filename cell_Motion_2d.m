%% cell motion 2-d (Mean curvature flow, %see class 12 notes)
function main()
a   =   -1; b =1; c= -1; d = 1;
C   =   1;
dx  =   .02;
dt  =   dx^2/8;
xV  = (a:dx:b)';
yV  = (c:dx:d)';
alp = 500;
beta= 500;
numSteps = 50;
fT  =  numSteps * dt;;
pTime = fT/4;
crs =   [ 'm', 'r', 'b', 'k', 'c', 'y', 'm', 'y' ];
lgnd =  { 't = 0' };
%%-----------------------------
u0  = zeros( (b-a)/dx +1,(d-c)/dx +1 );
A   = u0;
M   = u0;
for cnt = 1:(d-c)/dx +1
%    u0( :,cnt ) = max( abs(xV)-.495501, abs(yV(cnt)) - .495501 )+ .005;
   u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .5^2;
   A( :, cnt ) = cnt/((d-c)/dx +1);
   M( :, cnt ) = 1.1;%1-cnt/((d-c)/dx +1);
end
figure(2);subplot(3,2, 4); mesh( xV, yV, A.*( u0 <=0 ) );xlabel('x-axis');title( 'A at t0' )
u = padarray(u0, [1 1], 0 );
u(:,1) = u(:,2)+ u(:,2) -u(:,3); u(:,end)= u(:,end-1)+u(:,end-1)-u(:,end-2);
u(1,:) = u(2,:)+ u(2,:)-u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);
d0 = 0*u;

A(:,1) = A(:,2); A(:,end)= A(:,end-1);A(1,:) = A(2,:); A(end,:)= A(end-1,:);
M(:,1) = M(:,2); M(:,end)= M(:,end-1);M(1,:) = M(2,:); M(end,:)= M(end-1,:);

 subplot( 3,2, 1 );contour( xV, yV,u0, [0 0], 'k' );grid; hold on;

for t = dt:dt:fT
    ux      =   ( u( 2:end-1, 3:end) - u(2:end-1, 1:end-2) )/(2*dx);
    uy      =   ( u(3:end, 2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
    uxx     =   ( u( 2:end-1, 3:end) - 2*u(2:end-1, 2:end-1) + u(2:end-1, 1:end-2) )/(dx^2);
    uyy     =   ( u(3:end,2:end-1) - 2*u(2:end-1,2:end-1) + u(1:end-2,2:end-1) )/(dx^2);
    uxy     =   (u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2))/(4*dx^2);
    cTerm   =   ux.^2 .* uxx + 2*ux.*uy.*uxy + uy.^2.*uyy;
    gradM2  =   ux.^2 + uy.^2 + 1e-9; 
   
    ut      = C * (uxx + uyy - cTerm./gradM2) + (beta*A - alp* M).*sqrt(gradM2 );
   
    boundary     = smeeka01( u , dx );
    
    totExpansion = sum( sum( ut./sqrt( gradM2 ) .* boundary ) );
    boundarySize = sum( sum( boundary ) );
%     ut/sqrt(gradM2) is what we actually adjust, hence the mult through by sqrt(gradM2)
    ut = ut - sqrt(gradM2).*totExpansion/boundarySize;
   
    u(2:end-1, 2:end-1)     =   u(2:end-1, 2:end-1) + dt* ut;
    u(:,1) = u(:,2); u(:,end)= u(:,end-1);
    u(1,:) = u(2,:); u(end,:)= u(end-1,:);
  
    %Extend A and M into the new cell shape before solving PDE, maybe
    %instead of just extending by the normal to the surface we could use
    %the direction of the velocity on the surface.
    colInds = min(u ) < 5*dx;
    rowInds = min(u' ) < 5*dx;
    tempU   =   u(colInds, rowInds );
    
    d0(colInds, rowInds )      =   sweepDistBoundary( tempU,dx, 2 );
    
    d0(~colInds, : )= 10*dx;
    d0(:, ~rowInds )= 10*dx;
    
    Ap      =   sweepVelocityBoundary( A, d0(2:end-1,2:end-1), dx,1);
    Mp      =   sweepVelocityBoundary( M, d0(2:end-1,2:end-1), dx,1);
    A(u(2:end-1,2:end-1)>0)  =   Ap( u(2:end-1,2:end-1)>0 );
    M(u(2:end-1,2:end-1)>0)  =   Mp( u(2:end-1,2:end-1)>0 );
    
    [A,M ] = amPDE( A, M, dx, dt );

    %now implement CIM for AM here....
   
    if (mod( t, pTime) == 0 )
        contour( xV, yV,u(2:end-1, 2:end-1), [0 0], crs( round(t/pTime ))  );
        lgnd = [ lgnd; ['t = ',num2str( t ) ]];
    end

end
legend( lgnd );title( 'Zero Level Set/Cell Boundary' );xlabel('x-axis' );
hold off;
u = u(2:end-1, 2:end-1);

subplot(3,2, 3); mesh( xV, yV, A.*( u <=0 ) );xlabel('x-axis');title( 'A ( actin bundles, W in paper )' )
subplot(3,2, 5); mesh( xV, yV, M.*( u <=0 ) );xlabel('x-axis');title( 'M ( actin cross linked network (v)' );


% figure(1); mesh(xV, yV,u0');
% title( 'u(x,0)' );xlabel( 'x-axis');ylabel( 'y-axis');
% 
% figure(3); mesh(xV, yV,u');
% title( [ 'u(x, ', num2str( fT ), ' )' ] );xlabel( 'x-axis');ylabel( 'y-axis');



end

