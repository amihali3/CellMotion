% smooth curve... try to keep curve area the same while smoothing curve


function smoothCurve()
    a   =   -1.2; b =1.2; c= a; d = b;
    dx  =   .02;
    dt  =   dx^2/8;
    fT  =   .1;
    xV  = (a:dx:b)';yV  = (c:dx:d)';
    pTime = fT/4;
    epsilon = dx;
    crs =   [ 'm', 'r', 'b', 'c', 'k', 'y', 'm', 'y' ];
    lgnd =  { 't = 0' };
    for cnt = 1:(d-c)/dx +1
        u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .51 );  % .5 diamond
%    u0( :,cnt ) = max( abs(xV)-1.5, abs(yV(cnt)) - 1.5 );  %1.5 diamond
%  u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .5^2;
%     if(cnt< (b-a)/dx/2 )
%         u0( :,cnt ) = max( abs(xV)-.5, abs(yV(cnt)) - .5 );
%     else
%         u0( :,cnt ) = xV.^2 + yV(cnt)^2 - .5^2;
%     end     
    end
    u = padarray(u0, [1 1], 0 );
    
    figure(4); contour( xV, yV,u0', [0 0], 'k' );grid;hold on;
    
    for t = dt:dt:fT
        u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
        u(1,:) = 2*u(2,:) - u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);

        delU = calcDelta(u, dx);

          %% MEAN CURVATURE SET-UP 
        ux      =   ( u( 3:end,  2:end-1) - u(1:end-2, 2:end-1) )/(2*dx);
        uy      =   ( u( 2:end-1,3:end)   - u(2:end-1, 1:end-2) )/(2*dx);
        uxx     =   ( u( 3:end,  2:end-1 )- 2*u(2:end-1,2:end-1) + u(1:end-2, 2:end-1) )/(dx^2);
        uyy     =   ( u( 2:end-1,3:end)   - 2*u(2:end-1,2:end-1) + u(2:end-1,1:end-2) )/(dx^2);
        uxy     =   (u(3:end,3:end) - u(3:end,1:end-2) - u(1:end-2,3:end) + u(1:end-2,1:end-2))/(4*dx^2);
        cTerm   =   (ux.^2 .* uxx) + 2*( ux.*uy.*uxy ) + (uy.^2 .* uyy );
        gradM2  =   ux.^2 + uy.^2 + epsilon; 
        mCur    =   (uxx + uyy - cTerm./gradM2);

        %%find total div in order to get Havg
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

        gradUM  =   sqrt(adX.^2+adY.^2) ; 


        Havg    =   sum(sum(mCur.*delU))/(sum(sum(gradUM.*delU))+epsilon)
        
        if (mod( t, pTime) == 0 )
            contour( xV, yV,u(2:end-1, 2:end-1)', [0 0], crs( round(t/pTime ))  );
            lgnd = [ lgnd; ['t = ',num2str( t ) ]];
        end
        
        %%now step in time
        u(2:end-1, 2:end-1) =   u(2:end-1, 2:end-1) + dt * (mCur - Havg*gradUM ); 

    end
    hold off;
    u = u(2:end-1,2:end-1);
        
    figure(1);  mesh(xV, yV,delU');
    figure(2);  mesh(xV, yV, u');;title( 'uF');
    figure(3);  mesh(xV, yV, u0');title( 'u0');
    
end



function delU = calcDelta( u, dx )
    delU = zeros( size( u) );
    method = 2;
    switch method
        case    1
        %%neg to pos
        delU(1:end-1,:) =delU(1:end-1,:)+ (u(1:end-1,:)<0 & u(2:end,:)>0) .*  (u(2:end,:  ))./max(u(2:end,:)-u(1:end-1,:), eps );   
        delU(2:end,:  ) =delU(2:end,:  )+ (u(1:end-1,:)<0 & u(2:end,:)>0) .* -(u(1:end-1,:))./max(u(2:end,:)-u(1:end-1,:), eps );
    %     %pos to neg
        delU(1:end-1,:) =delU(1:end-1,:)+ (u(1:end-1,:)>0 & u(2:end,:)<0) .*  (u(2:end,:  ))./min(u(2:end,:)-u(1:end-1,:),-eps );   
        delU(2:end,:  ) =delU(2:end,:  )+ (u(1:end-1,:)>0 & u(2:end,:)<0) .* -(u(1:end-1,:))./min(u(2:end,:)-u(1:end-1,:),-eps );
    %     %%neg to pos
        delU(:,1:end-1) =delU(:,1:end-1)+ (u(:,1:end-1)<0 & u(:,2:end)>0) .*  (u(:,2:end  ))./max(u(:,2:end)-u(:,1:end-1),eps );   
        delU(:,2:end  ) =delU(:,2:end  )+ (u(:,1:end-1)<0 & u(:,2:end)>0) .* -(u(:,1:end-1))./max(u(:,2:end)-u(:,1:end-1),eps );
        %pos to neg
        delU(:,1:end-1) =delU(:,1:end-1)+ (u(:,1:end-1)>0 & u(:,2:end)<0) .*  (u(:,2:end  ))./min(u(:,2:end)-u(:,1:end-1),-eps );   
        delU(:,2:end  ) =delU(:,2:end  )+ (u(:,1:end-1)>0 & u(:,2:end)<0) .* -(u(:,1:end-1))./min(u(:,2:end)-u(:,1:end-1),-eps );
        delU( find( u == 0 ) ) = 1;
        
        case 2
        delU = 1- abs(u)./( abs(u) + dx^2 );
        
        
    end
    delU = delU(2:end-1,2:end-1);
    
end