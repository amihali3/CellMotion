%% Need to give this guy a padded u
function  u = sweepDistBoundary( u, dx, numSweeps)
    dy = dx;% assume symettry at this point
      tic;
    for cnt = 1: numSweeps
        %right ( not really r,l,u,d... just covering necessary sweeps )
        iFirst =    1;
        iStart =    size( u,1);  iEnd   =   1; icntr = -1;
        jStart =    size( u,2);  jEnd   =   1; jcntr = -1;      
        u   =   directionSweep( u, dx, dy, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );      
       
        %Down
        iFirst =    0;
        iStart =    1;  iEnd   =   size( u,1); icntr = 1;
        jStart =    1;  jEnd   =   size( u,2); jcntr = 1;           
        u   =   directionSweep( u, dx, dy, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
           
        %left
        iFirst =    1;
        iStart =    1;   iEnd   =   size( u,1 ); icntr = 1;
        jStart =    size( u,2 ) ;  jEnd   =   1; jcntr = -1;
        u   =   directionSweep( u, dx, dy, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
      
        %Up
        iFirst =    0;
        iStart =    size( u,1);   iEnd   =   1; icntr = -1;
        jStart =    1;  jEnd   =   size( u,2 ); jcntr = 1;
        u   =   directionSweep( u, dx, dy, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
    end
    toc;
end

function u = directionSweep( u, dx, dy, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr )
    dxHold = dx;
    dyHold = dy;
    noChangeThres   =  dx^2;
    if iFirst
        for i = iStart:icntr:iEnd
            for j = jStart:jcntr:jEnd
                dx  = dxHold;
                dy  = dyHold;
              
                if( u(i, j ) > noChangeThres )
                    if( i == 1)
                        xminN = u( 2, j );
                    elseif ( i== size(u,1))
                        xminN = u( i -1,j);
                    elseif( (u(i-1,j ) < 0) || (u( i+1, j ) < 0) )
                        iN      = i + 1 - 2*((u(i-1,j ) < 0) & (u( i+1, j ) < 0) & (u(i-1,j)>u(i+1,j)))...
                                    -2*(u(i+1,j) > 0 );
                        dx      = (dx * -1* u(i,j))/( u(iN,j) - u(i,j));
                        xminN = 0;
%                         dx    = (dx * -1* u(i,j))/(min( u(i-1,j), u(i+1,j ) ) - u(i,j));
                    else
                        xminN  =   min( u(i-1, j), u(i+1,j) );
                    end
                  
                    if( j == 1)
                        yminN = u( i, 2 );
                    elseif ( j== size(u,2))
                        yminN = u( i,j-1);
                    elseif( (u(i,j-1 ) < 0) || (u( i, j+1 ) < 0) )
                        jN      = j + 1 - 2*((u(i,j-1 ) < 0) & (u( i, j+1 ) < 0) & (u(i,j-1)>u(i,j+1)))...
                                    -2*(u(i,j+1) > 0 );
                        dy      = (dy * -1* u(i,j))/( u(i,jN)  - u(i,j));
                        yminN = 0;
%                         dy    = (dy * -1* u(i,j))/(min( u(i,j-1), u(i,j+1 ) ) - u(i,j));
                    else
                        yminN  =   min( u(i, j-1), u(i,j+1) );
                    end
                    a   =   dx^-2+dy^-2;
                    b   =   -2*xminN*dx^-2 - 2*yminN*dy^-2;
                    c   =   xminN^2*dx^-2+yminN^2*dy^-2 - 1;
                    
%                     uTemp  =  max(roots( [dx^-2+dy^-2, -2*xminN*dx^-2 - 2*yminN*dy^-2, xminN^2*dx^-2+yminN^2*dy^-2 - 1] ));
                    uTemp  =  (-b + sqrt(b^2 - 4*a*c))/(2*a);
                    
                    %    if( (xminN == 0) || (yminN == 0) )
%                         uTemp  = (xminN == 0 )*dx + (yminN == 0 ) * dy + (-dx - dy +uTemp) *(xminN == 0 && yminN == 0 );
%                     end
                  
                    if (uTemp > min( xminN, yminN ) && isreal( uTemp ) )
                        u( i, j ) = uTemp;
                    else
                        vMin = xminN * ( yminN >= xminN ) + yminN * ( yminN < xminN );
                        dv   = dx * ( yminN >= xminN ) + dy * ( yminN < xminN );
                        a = dv^-2;
                        b = -2*vMin*dv^-2 ;
                        c = vMin^2*dv^-2 - 1;
                        uTemp  = (-b + sqrt(b^2 - 4*a*c))/(2*a);
%                         uTemp  =  max(roots( [dv^-2, -2*vMin*dv^-2 , vMin^2*dv^-2 - 1] ));
%                         if (uTemp > min( xminN, yminN ) && isreal( uTemp ) )
                            u( i, j ) = uTemp;
%                         else
%                             fprintf('TROUBLE\n');
%                         end
                    end            

                elseif( u(i,j ) < -noChangeThres )
                  
                    if( i == 1)
                        xminN = -u( 2, j );
                    elseif ( i== size(u,1))
                        xminN = -u( i-1,j);
                    elseif( (u(i-1,j ) > 0) || (u( i+1, j ) > 0) )
                        iN      = i + 1 - 2*((u(i-1,j ) > 0) & (u( i+1, j ) > 0) & (u(i-1,j)<u(i+1,j)))...
                                    -2*(u(i+1,j) < 0 );
                        dx      = (dx * -1* u(i,j))/(u(iN,j) - u(i,j));
                        xminN = 0;
%                         dx    = (dx * -1* u(i,j))/(max( u(i-1,j), u(i+1,j ) ) - u(i,j));
                    else
                        xminN  =   min( -u(i-1, j), -u(i+1,j) );
                    end
                  
                    if( j == 1)
                        yminN = -u( i, 2 );
                    elseif ( j== size(u,2))
                        yminN = -u(i,j - 1);
                    elseif( (u(i,j-1 ) > 0) || (u( i, j+1 ) > 0) )
                        jN      = j + 1 - 2*((u(i,j-1 ) > 0) & (u( i, j+1 ) > 0) & (u(i,j-1)<u(i,j+1)))...
                                    -2*(u(i,j+1) < 0 );                   
                        dy      = (dy * -1* u(i,j))/( u(i,jN) - u(i,j));
                        yminN = 0;
%                         dy    = (dy * -1* u(i,j))/(max( u(i,j-1), u(i,j+1 ) ) - u(i,j));
                    else
                        yminN  =   min( -u(i, j-1), -u(i,j+1) );
                    end
                    a   =   dx^-2+dy^-2;
                    b   =   -2*xminN*dx^-2 - 2*yminN*dy^-2;
                    c   =   xminN^2*dx^-2+yminN^2*dy^-2 - 1;
                    
%                     uTemp  =  -1*max(roots( [dx^-2+dy^-2, -2*xminN*dx^-2 - 2*yminN*dy^-2, xminN^2*dx^-2+yminN^2*dy^-2 - 1] ));
                    uTemp  =  -1*(-b + sqrt(b^2 - 4*a*c))/(2*a);
%                     if( (xminN == 0) || (yminN == 0) )
%                         uTemp  = (xminN == 0 )* -dx + (yminN == 0 ) * -dy + (dx + dy + uTemp) *(xminN == 0 && yminN == 0 );
%                     end
                 
                    if (uTemp < min( -xminN, -yminN ) && isreal( uTemp ) )
                        u( i, j ) = uTemp;
                    else
                        vMin    = xminN * ( yminN >= xminN ) + yminN * ( yminN < xminN );
                        dv      = dx * ( yminN >= xminN ) + dy * ( yminN < xminN );
                        a = dv^-2;
                        b = -2*vMin*dv^-2 ;
                        c = vMin^2*dv^-2 - 1;
                        uTemp  =  -1*(-b + sqrt(b^2 - 4*a*c))/(2*a);
%                         uTemp   =  -1* max(roots( [dv^-2, -2*vMin*dv^-2 , vMin^2*dv^-2 - 1] ));
%                         if (uTemp < min( -xminN, -yminN ) && isreal( uTemp ) )
                            u( i, j ) = uTemp;
%                         else
%                             fprintf('TROUBLE\n');
%                         end
                    end
                end        
            end
        end
    else
        for j = jStart:jcntr:jEnd
            for i = iStart:icntr:iEnd
                dx  = dxHold;
                dy  = dyHold;
              
                if( u(i, j ) > noChangeThres )
                    if( i == 1)
                        xminN = u( 2, j );
                    elseif ( i== size(u,1))
                        xminN = u( i -1,j);
                    elseif( (u(i-1,j ) < 0) || (u( i+1, j ) < 0) )
                        iN      = i + 1 - 2*((u(i-1,j ) < 0) & (u( i+1, j ) < 0) & (u(i-1,j)>u(i+1,j)))...
                                    -2*(u(i+1,j) > 0 );
                        dx      = (dx * -1* u(i,j))/( u(iN,j) - u(i,j));
                        xminN = 0;
%                         dx    = (dx * -1* u(i,j))/(min( u(i-1,j), u(i+1,j ) ) - u(i,j));
                    else
                        xminN  =   min( u(i-1, j), u(i+1,j) );
                    end
                  
                    if( j == 1)
                        yminN = u( i, 2 );
                    elseif ( j== size(u,2))
                        yminN = u( i,j-1);
                    elseif( (u(i,j-1 ) < 0) || (u( i, j+1 ) < 0) )
                        jN      = j + 1 - 2*((u(i,j-1 ) < 0) & (u( i, j+1 ) < 0) & (u(i,j-1)>u(i,j+1)))...
                                    -2*(u(i,j+1) > 0 );
                        dy      = (dy * -1* u(i,j))/( u(i,jN)  - u(i,j));
                        yminN = 0;
%                         dy    = (dy * -1* u(i,j))/(min( u(i,j-1), u(i,j+1 ) ) - u(i,j));
                    else
                        yminN  =   min( u(i, j-1), u(i,j+1) );
                    end
                    a   =   dx^-2+dy^-2;
                    b   =   -2*xminN*dx^-2 - 2*yminN*dy^-2;
                    c   =   xminN^2*dx^-2+yminN^2*dy^-2 - 1;
                    
%                     uTemp  =  max(roots( [dx^-2+dy^-2, -2*xminN*dx^-2 - 2*yminN*dy^-2, xminN^2*dx^-2+yminN^2*dy^-2 - 1] ));
                    uTemp  =  (-b + sqrt(b^2 - 4*a*c))/(2*a);
                    
                    %    if( (xminN == 0) || (yminN == 0) )
%                         uTemp  = (xminN == 0 )*dx + (yminN == 0 ) * dy + (-dx - dy +uTemp) *(xminN == 0 && yminN == 0 );
%                     end
                  
                    if (uTemp > min( xminN, yminN ) && isreal( uTemp ) )
                        u( i, j ) = uTemp;
                    else
                        vMin = xminN * ( yminN >= xminN ) + yminN * ( yminN < xminN );
                        dv   = dx * ( yminN >= xminN ) + dy * ( yminN < xminN );
                        a = dv^-2;
                        b = -2*vMin*dv^-2 ;
                        c = vMin^2*dv^-2 - 1;
                        uTemp  =  (-b + sqrt(b^2 - 4*a*c))/(2*a);
%                         uTemp  =  max(roots( [dv^-2, -2*vMin*dv^-2 , vMin^2*dv^-2 - 1] ));
%                         if (uTemp > min( xminN, yminN ) && isreal( uTemp ) )
                            u( i, j ) = uTemp;
%                         else
%                             fprintf('TROUBLE\n');
%                         end
                    end            

                elseif( u(i,j ) < -noChangeThres )
                 
                    if( i == 1)
                        xminN = -u( 2, j );
                    elseif ( i== size(u,1))
                        xminN = -u( i-1,j);
                    elseif( (u(i-1,j ) > 0) || (u( i+1, j ) > 0) )
                        iN      = i + 1 - 2*((u(i-1,j ) > 0) & (u( i+1, j ) > 0) & (u(i-1,j)<u(i+1,j)))...
                                    -2*(u(i+1,j) < 0 );
                        dx      = (dx * -1* u(i,j))/(u(iN,j) - u(i,j));
                        xminN = 0;
%                         dx    = (dx * -1* u(i,j))/(max( u(i-1,j), u(i+1,j ) ) - u(i,j));
                    else
                        xminN  =   min( -u(i-1, j), -u(i+1,j) );
                    end
                  
                    if( j == 1)
                        yminN = -u( i, 2 );
                    elseif ( j== size(u,2))
                        yminN = -u(i,j - 1);
                    elseif( (u(i,j-1 ) > 0) || (u( i, j+1 ) > 0) )
                        jN      = j + 1 - 2*((u(i,j-1 ) > 0) & (u( i, j+1 ) > 0) & (u(i,j-1)<u(i,j+1)))...
                                    -2*(u(i,j+1) < 0 );                   
                        dy      = (dy * -1* u(i,j))/( u(i,jN) - u(i,j));
                        yminN = 0;
%                         dy    = (dy * -1* u(i,j))/(max( u(i,j-1), u(i,j+1 ) ) - u(i,j));
                    else
                        yminN  =   min( -u(i, j-1), -u(i,j+1) );
                    end
                    a   =   dx^-2+dy^-2;
                    b   =   -2*xminN*dx^-2 - 2*yminN*dy^-2;
                    c   =   xminN^2*dx^-2+yminN^2*dy^-2 - 1;
                    
%                     uTemp  =  -1*max(roots( [dx^-2+dy^-2, -2*xminN*dx^-2 - 2*yminN*dy^-2, xminN^2*dx^-2+yminN^2*dy^-2 - 1] ));
                    uTemp  =  -1*(-b + sqrt(b^2 - 4*a*c))/(2*a);
%                     if( (xminN == 0) || (yminN == 0) )
%                         uTemp  = (xminN == 0 )* -dx + (yminN == 0 ) * -dy + (dx + dy + uTemp) *(xminN == 0 && yminN == 0 );
%                     end
                 
                    if (uTemp < min( -xminN, -yminN ) && isreal( uTemp ) )
                        u( i, j ) = uTemp;
                    else
                        vMin    = xminN * ( yminN >= xminN ) + yminN * ( yminN < xminN );
                        dv      = dx * ( yminN >= xminN ) + dy * ( yminN < xminN );
                        a = dv^-2;
                        b = -2*vMin*dv^-2 ;
                        c = vMin^2*dv^-2 - 1;
                        uTemp  =  -1*(-b + sqrt(b^2 - 4*a*c))/(2*a);
%                         uTemp   =  -1* max(roots( [dv^-2, -2*vMin*dv^-2 , vMin^2*dv^-2 - 1] ));
%                         if (uTemp < min( -xminN, -yminN ) && isreal( uTemp ) )
                            u( i, j ) = uTemp;
%                         else
%                             fprintf('TROUBLE\n');
%                         end
                    end
                end        
            end
        end
             
    end
%     u(:,1) = 2*u(:,2)- u(:,3); u(:,end)= 2*u(:,end-1) - u(:,end-2);
%     u(1,:) = 2*u(2,:)- u(3,:); u(end,:)= 2*u(end-1,:) - u(end-2,:);
end