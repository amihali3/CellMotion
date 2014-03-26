%% Need to give this guy a padded u
function  vExt = sweepVelocityBoundary( v, d, dx, numSweeps)
epsilon = 1e-15;
% tic
    dy = dx;% assume symettry at this point    
    for cnt = 1: numSweeps        
        iFirst =    1;
        iStart =    size( v,2);  iEnd   =   1; icntr = -1;
        jStart =    size( v,1);  jEnd   =   1; jcntr = -1;   
        v   =   vSweep(  d, v, dx, dy ,epsilon, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
        
        iFirst =   0;
        iStart =    1;  iEnd   =   size( v,2); icntr = 1;
        jStart =    1;  jEnd   =   size( v,2); jcntr = 1;   
        v   =   vSweep(  d, v, dx, dy ,epsilon, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
        
        
        iFirst =    1;
        iStart =    1;  iEnd   =   size( v,2); icntr = 1;
        jStart =    size( v,1);  jEnd   =  1; jcntr = -1;   
        v   =   vSweep(  d, v, dx, dy ,epsilon, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );
        
        iFirst =    0;
        iStart =    size( v,2);  iEnd   =   1; icntr = -1;
        jStart =    1;  jEnd   =   size( v,2); jcntr = 1;   
        v   =   vSweep(  d, v, dx, dy ,epsilon, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr );

    end
    vExt = v;
% toc
end 
function v = vSweep( d, v, dx, dy, epsilon, iFirst, iStart, iEnd, icntr, jStart, jEnd, jcntr ) 

    dxHold = dx;
    dyHold = dy;
    if iFirst
        for i = iStart:icntr:iEnd
            for j = jStart:jcntr:jEnd
                dx  = dxHold;
                dy  = dyHold;                        
                if( d(i, j ) > 0 && (d(i,j) < 5*dx) )  %POSITIVE POINTS
                    if( i == 1)
                        ddx = d(1,j) - d(2,j);
                        vxN = v(2,j);
                    elseif( i== size(d,2) )
                        ddx = d(end,j) - d(end-1, j );
                        vxN = v(end-1, j );
                    elseif( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                        iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                            -2*(d(i+1,j) > 0 );

                        dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                        vtemp   = v(iN,j );
                        vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                        ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)

                    else
                        ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                        vxN     = v(i-1,j)* (d(i-1,j)<=d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                    end
                    if( j== 1)
                        ddy = d(i, 1 ) - d( i, 2 );
                        vyN = v(i, 2 );
                    elseif( j == size( d,2 ))
                        ddy = d(i,end) - d(i, end-1 );
                        vyN = v(i,end-1 );
                    elseif((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                        jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                            -2*(d(i,j+1) > 0 );

                        dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                        vtemp   = v(i,jN );
                        vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                        ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)

                    else
                        ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                        vyN     =  v(i,j-1)* (d(i,j-1)<=d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                    end           
                v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2); 
                                  
                end
            end
        end
    
    else  %then j is first
        for j = jStart:jcntr:jEnd
            for i = iStart:icntr:iEnd
                dx  = dxHold;
                dy  = dyHold;                        
                if( d(i, j ) > 0 && (d(i,j) < 5*dx) )  %POSITIVE POINTS
                    if( i == 1)
                        ddx = d(1,j) - d(2,j);
                        vxN = v(2,j);
                    elseif( i== size(d,2) )
                        ddx = d(end,j) - d(end-1, j );
                        vxN = v(end-1, j );
                    elseif( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                        iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                            -2*(d(i+1,j) > 0 );

                        dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                        vtemp   = v(iN,j );
                        vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                        ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)

                    else
                        ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                        vxN     = v(i-1,j)* (d(i-1,j)<=d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                    end
                    if( j== 1)
                        ddy = d(i, 1 ) - d( i, 2 );
                        vyN = v(i, 2 );
                    elseif( j == size( d,2 ))
                        ddy = d(i,end) - d(i, end-1 );
                        vyN = v(i,end-1 );
                    elseif((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                        jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                            -2*(d(i,j+1) > 0 );

                        dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                        vtemp   = v(i,jN );
                        vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                        ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)

                    else
                        ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                        vyN     =  v(i,j-1)* (d(i,j-1)<=d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                    end           

                v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2);    
                end
            end
        
        
        
        end

    end
end
