%% HW 6 5/9
tic
a   =   -1; b =1; c= -1; d = 1;e=-1; f= 1;
dx  =   .02;
dt  =   dx/2;
xV  = (a:dx:b)';
yV  = (c:dx:d)';
zV  = (e:dx:f);

fT  =   .1;

epsS =  dx;

%%----------------------------------------

u0  = zeros( (b-a)/dx +1,(d-c)/dx +1, (f-e)/dx +1 );
dAct = u0;

for zcnt = 1:(f-e)/dx +1
    for ycnt = 1:(d-c)/dx +1
        dAct( :,ycnt, zcnt ) = (xV.^2 + yV(ycnt)^2 + zV(zcnt)^2).^.5 - .5;
        u0( :,ycnt, zcnt )   = (xV.^2 + yV(ycnt)^2 + zV(zcnt)^2) - .5^2;
    end
end
u = padarray(u0, [1 1 1], 0 );
u( :, :, 1) = u(:,:,2);u( :, :, end) = u(:,:,end-1);
u( :, 1, :) = u(:,2,:);u( :, end, :) = u(:,end-1, :);
u( 1, :, :) = u(2,:,:);u( end, :, :) = u(end-1,:, :);

sgnU = u./( u.^2 + epsS ).^.5;

%I can have u transformed into the dist function because the only
% think I still need from u is sgnU which I won't update
for t = dt:dt:fT
   
    bdX     =   ( u( 2:end-1, 2:end-1, 2:end-1) - u(1:end-2, 2:end-1, 2:end-1) )/dx;
    fdX     =   ( u( 3:end  , 2:end-1, 2:end-1) - u(2:end-1, 2:end-1, 2:end-1) )/dx;
    bdY     =   ( u( 2:end-1, 2:end-1, 2:end-1) - u(2:end-1, 1:end-2, 2:end-1) )/dx;
    fdY     =   ( u( 2:end-1, 3:end  , 2:end-1) - u(2:end-1, 2:end-1, 2:end-1) )/dx;
    bdZ     =   ( u( 2:end-1, 2:end-1, 2:end-1) - u(2:end-1, 2:end-1, 1:end-2) )/dx;
    fdZ     =   ( u( 2:end-1, 2:end-1, 3:end  ) - u(2:end-1, 2:end-1, 2:end-1) )/dx;
   
    dXPsgn     =    (bdX>=0 & fdX>=0).* bdX + ...
                    (bdX<0  & fdX<0 ).* fdX + ...
                    (bdX>0  & fdX<0 ).* max( abs(bdX), abs(fdX))+...
                    (bdX<0  & fdX>0 ).* bdX *0;
                
    dYPsgn     =    (bdY>=0 & fdY>=0).* bdY + ...
                    (bdY<0  & fdY<0 ).* fdY + ...
                    (bdY>0  & fdY<0 ).* max( abs(bdY), abs(fdY))+...
                    (bdY<0  & fdY>0 ).* bdY *0;
    
    dZPsgn     =    (bdZ>=0 & fdZ>=0).* bdZ + ...
                    (bdZ<0  & fdZ<0 ).* fdZ + ...
                    (bdZ>0  & fdZ<0 ).* max( abs(bdZ), abs(fdZ))+...
                    (bdZ<0  & fdZ>0 ).* bdZ *0;    
    
    dXNsgn     =    (bdX<=0 & fdX<=0).* bdX + ...
                    (bdX>0  & fdX>0 ).* fdX + ...
                    (bdX<0  & fdX>0 ).* max( abs(bdX), abs(fdX))+...
                    (bdX>0  & fdX<0 ).* bdX *0;
                
    dYNsgn     =    (bdY<=0 & fdY<=0).* bdY + ...
                    (bdY>0  & fdY>0 ).* fdY + ...
                    (bdY<0  & fdY>0 ).* max( abs(bdY), abs(fdY))+...
                    (bdY>0  & fdY<0 ).* bdY *0;
    
    dZNsgn     =    (bdZ<=0 & fdZ<=0).* bdZ + ...
                    (bdZ>0  & fdZ>0 ).* fdZ + ...
                    (bdZ<0  & fdZ>0 ).* max( abs(bdZ), abs(fdZ))+...
                    (bdZ>0  & fdZ<0 ).* bdZ *0;
                        

    u(2:end-1,2:end-1,2:end-1 ) =  u(2:end-1,2:end-1, 2:end-1 ) - dt * ...
                 ( (  ( sgnU(2:end-1, 2:end-1, 2:end-1 ) >0 ).* (dXPsgn.^2 + dYPsgn.^2 + dZPsgn.^2 ) + ... 
                      ( sgnU(2:end-1, 2:end-1, 2:end-1 ) <0 ).* (dXNsgn.^2 + dYNsgn.^2 + dZNsgn.^2 ) ... %% BW diff when positive
                  ).^.5 -1 ) .* sgnU(2:end-1, 2:end-1, 2:end-1 );  
                
                
    u( :, :, 1) = u(:,:,2);u( :, :, end) = u(:,:,end-1);
    u( :, 1, :) = u(:,2,:);u( :, end, :) = u(:,end-1, :);
    u( 1, :, :) = u(2,:,:);u( end, :, :) = u(end-1,:, :);

end
toc

u       = u(2:end-1, 2:end-1, 2:end-1 );
sgnU    = sgnU(2:end-1, 2:end-1, 2:end-1 );
isosurface(xV, yV, zV, u, 0);grid;
figure(2)
isosurface(xV, yV, zV, sgnU, 0);grid;


figure(3); contour( xV, yV,squeeze( u0(:,:,round((f-e)/dx/2)+1)'), [0 0] );hold on;
 contour( xV, yV,squeeze( u(:,:,round((f-e)/dx/2)+1)'), [0 0] ,'r');grid;
legend( {'Initial Zero level Set'; 'Final Zero Level Set' } ); title( [ 'Cross-Section of Zero Level Sets, Final time = ', num2str(fT) ]);
xlabel( 'x-axis' ); ylabel( 'y-axis' );hold off;

%%part b, test the accurracy of a few pts
fprintf('Calculate d value  = %f\n', u(    round((b-a)/dx/4)+1, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Actual d value     = %f\n', dAct( round((b-a)/dx/4)+1, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Calculate d value  = %f\n', u(    round((b-a)/dx/4)+2, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Actual d value     = %f\n', dAct( round((b-a)/dx/4)+2, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Calculate d value  = %f\n', u(    round((b-a)/dx/4)-1, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Actual d value     = %f\n', dAct( round((b-a)/dx/4)-1, round((d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Calculate d value  = %f\n',u(    round((2-1/sqrt(2))*(b-a)/dx/4),round((2-1/sqrt(2))*(d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );
fprintf('Actual d value     = %f\n', dAct(round((2-1/sqrt(2))*(b-a)/dx/4),round((2-1/sqrt(2))*(d-c)/dx/2)+1, round((f-e)/dx/2)+1 ) );

