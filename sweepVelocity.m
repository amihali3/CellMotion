%% Need to give this guy a padded u
function  vExt = sweepVelocity( v, d, dx, numSweeps)
epsilon = 1e-15;
    dy = dx;% assume symettry at this point    
    for cnt = 1: numSweeps
        v   =   sweepVR(  d, v, dx, dy ,epsilon);
        
        v(:,1) = 2*v(:,2)- v(:,3); v(:,end)= 2*v(:,end-1) - v(:,end-2);
        v(1,:) = 2*v(2,:) - v(3,:); v(end,:)= 2*v(end-1,:) - v(end-2,:);
        
        v   =   sweepVDn( d, v, dx, dy,epsilon );
        v(:,1) = 2*v(:,2)- v(:,3); v(:,end)= 2*v(:,end-1) - v(:,end-2);
        v(1,:) = 2*v(2,:) - v(3,:); v(end,:)= 2*v(end-1,:) - v(end-2,:);
        v   =   sweepVL(  d, v, dx, dy,epsilon );
        
        v(:,1) = 2*v(:,2)- v(:,3); v(:,end)= 2*v(:,end-1) - v(:,end-2);
        v(1,:) = 2*v(2,:) - v(3,:); v(end,:)= 2*v(end-1,:) - v(end-2,:);
        
        v   =   sweepVUp( d, v, dx, dy,epsilon );
        v(:,1) = 2*v(:,2)- v(:,3); v(:,end)= 2*v(:,end-1) - v(:,end-2);
        v(1,:) = 2*v(2,:) - v(3,:); v(end,:)= 2*v(end-1,:) - v(end-2,:);
    end
    vExt = v;
end 
function v = sweepVUp( d, v, dx, dy, epsilon )
    dxHold = dx;
    dyHold = dy;
    for j = 2:size( v, 2)-1
        for i = 2:size( v, 1)-1
            dx  = dxHold;
            dy  = dyHold;                        
            if( d(i, j ) > 0 )  %POSITIVE POINTS
                if( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                        -2*(d(i+1,j) > 0 );
                
                    dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)<d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                end
                if((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                        -2*(d(i,j+1) > 0 );
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                    vyN     =  v(i,j-1)* (d(i,j-1)<d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                end           
                
            elseif( d(i,j ) < 0 )  %NEGATIVE PTs
                if( (d(i-1,j ) > 0) | (d( i+1, j ) > 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) > 0) & (d( i+1, j ) > 0) & (d(i-1,j)<d(i+1,j)))...
                        -2*(d(i+1,j) < 0 );
                    
                    dx      = (dx * -1* d(i,j))/(d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = -dx -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - max(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)>d(i+1,j)) + v(i+1,j)* (d(i+1,j)>d(i-1,j));
                end
                if((d(i,j-1 ) > 0) | (d( i, j+1) > 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) > 0) & (d( i, j+1 ) > 0) & (d(i,j-1)<d(i,j+1)))...
                        -2*(d(i,j+1) < 0 );                   
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN) - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = -dy -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - max( d(i,j-1), d(i, j+1) );
                    vyN      = v(i,j-1)* (d(i,j-1)>d(i,j+1)) + v(i,j+1)* (d(i,j+1)>d(i,j-1));
                end
                
            end
            v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2);   
        end
    end
end


function v = sweepVDn( d, v, dx, dy, epsilon )
    dxHold = dx;
    dyHold = dy;
    for j = size( v, 2)-1:-1:2
        for i = 2:size( v, 1)-1
            dx  = dxHold;
            dy  = dyHold;                        
            if( d(i, j ) > 0 )  %POSITIVE POINTS
                if( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                        -2*(d(i+1,j) > 0 );
                
                    dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)<d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                end
                if((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                        -2*(d(i,j+1) > 0 );
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                    vyN     =  v(i,j-1)* (d(i,j-1)<d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                end           
                
            elseif( d(i,j ) < 0 )  %NEGATIVE PTs
                if( (d(i-1,j ) > 0) | (d( i+1, j ) > 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) > 0) & (d( i+1, j ) > 0) & (d(i-1,j)<d(i+1,j)))...
                        -2*(d(i+1,j) < 0 );
                    
                    dx      = (dx * -1* d(i,j))/(d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = -dx -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - max(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)>d(i+1,j)) + v(i+1,j)* (d(i+1,j)>d(i-1,j));
                end
                if((d(i,j-1 ) > 0) | (d( i, j+1) > 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) > 0) & (d( i, j+1 ) > 0) & (d(i,j-1)<d(i,j+1)))...
                        -2*(d(i,j+1) < 0 );                   
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN) - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = -dy -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - max( d(i,j-1), d(i, j+1) );
                    vyN      = v(i,j-1)* (d(i,j-1)>d(i,j+1)) + v(i,j+1)* (d(i,j+1)>d(i,j-1));
                end
                
            end
            v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2);          
        end
    end
end

function v = sweepVR( d,v, dx, dy, epsilon )
    dxHold = dx;
    dyHold = dy;
    for i = size( v, 2)-1:-1:2
        for j = 2:size( v, 1)-1
            dx  = dxHold;
            dy  = dyHold;                        
            if( d(i, j ) > 0 )  %POSITIVE POINTS
                if( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                        -2*(d(i+1,j) > 0 );
                
                    dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)<d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                end
                if((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                        -2*(d(i,j+1) > 0 );
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                    vyN     =  v(i,j-1)* (d(i,j-1)<d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                end           
                
            elseif( d(i,j ) < 0 )  %NEGATIVE PTs
                if( (d(i-1,j ) > 0) | (d( i+1, j ) > 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) > 0) & (d( i+1, j ) > 0) & (d(i-1,j)<d(i+1,j)))...
                        -2*(d(i+1,j) < 0 );
                    
                    dx      = (dx * -1* d(i,j))/(d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = -dx -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - max(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)>d(i+1,j)) + v(i+1,j)* (d(i+1,j)>d(i-1,j));
                end
                if((d(i,j-1 ) > 0) | (d( i, j+1) > 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) > 0) & (d( i, j+1 ) > 0) & (d(i,j-1)<d(i,j+1)))...
                        -2*(d(i,j+1) < 0 );                   
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN) - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = -dy -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - max( d(i,j-1), d(i, j+1) );
                    vyN      = v(i,j-1)* (d(i,j-1)>d(i,j+1)) + v(i,j+1)* (d(i,j+1)>d(i,j-1));
                end
                
            end
            v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2);          
        end
    end
end

function v = sweepVL( d, v, dx, dy, epsilon )
    dxHold = dx;
    dyHold = dy;
    for i = size( v, 1)-1:-1:2
        for j = (size( v, 2)-1):-1:2
            dx  = dxHold;
            dy  = dyHold;                        
            if( d(i, j ) > 0 )  %POSITIVE POINTS
                if( (d(i-1,j ) < 0) | (d( i+1, j ) < 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) < 0) & (d( i+1, j ) < 0) & (d(i-1,j)>d(i+1,j)))...
                        -2*(d(i+1,j) > 0 );
                
                    dx      = (dx * -1* d(i,j))/( d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = dx +epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - min(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)<d(i+1,j)) + v(i+1,j)* (d(i+1,j)<d(i-1,j));
                end
                if((d(i,j-1 ) < 0) | (d( i, j+1) < 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) < 0) & (d( i, j+1 ) < 0) & (d(i,j-1)>d(i,j+1)))...
                        -2*(d(i,j+1) > 0 );
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN)  - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = dy + epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - min( d(i,j-1), d(i, j+1) );
                    vyN     =  v(i,j-1)* (d(i,j-1)<d(i,j+1)) + v(i,j+1)* (d(i,j+1)<d(i,j-1));
                end           
                
            elseif( d(i,j ) < 0 )  %NEGATIVE PTs
                if( (d(i-1,j ) > 0) | (d( i+1, j ) > 0) )
                    iN      = i + 1 - 2*((d(i-1,j ) > 0) & (d( i+1, j ) > 0) & (d(i-1,j)<d(i+1,j)))...
                        -2*(d(i+1,j) < 0 );
                    
                    dx      = (dx * -1* d(i,j))/(d(iN,j) - d(i,j));
                    vtemp   = v(iN,j );
                    vxN     =  v(i,j) + dx/dxHold*(vtemp - (v(i,j)));
                    ddx     = -dx -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddx     = d(i,j) - max(d(i-1,j), d(i+1,j) );
                    vxN     = v(i-1,j)* (d(i-1,j)>d(i+1,j)) + v(i+1,j)* (d(i+1,j)>d(i-1,j));
                end
                if((d(i,j-1 ) > 0) | (d( i, j+1) > 0) )
                    jN      = j + 1 - 2*((d(i,j-1 ) > 0) & (d( i, j+1 ) > 0) & (d(i,j-1)<d(i,j+1)))...
                        -2*(d(i,j+1) < 0 );                   
                    
                    dy      = (dy * -1* d(i,j))/( d(i,jN) - d(i,j));
                    vtemp   = v(i,jN );
                    vyN     =  v(i,j) + dy/dyHold*(vtemp - (v(i,j)));
                    ddy     = -dy -epsilon; % ddx is dij - (di-1j or di+1j)
                    
                else
                    ddy     = d(i,j) - max( d(i,j-1), d(i, j+1) );
                    vyN      = v(i,j-1)* (d(i,j-1)>d(i,j+1)) + v(i,j+1)* (d(i,j+1)>d(i,j-1));
                end
                
            end
            v(i,j)  =   (vxN*ddx/dx^2 + vyN*ddy/dy^2)/(ddx/dx^2+ddy/dy^2);     
        end
    end
end