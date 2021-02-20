function qt = curvisolve_scalar(qt,t,dT,cfl_des,grid,params)
%CURVISOLVE_SCALAR.   Sample advection solver with curvilinear coordinates
    
    % grid parameters
    XI = grid.XI; ET = grid.ET;
    mxi = grid.mxi; met = grid.met;
    dxi = grid.dxi; det = grid.det;

    % Transformation derivatives:
    dx_dxi = grid.dx_dxi;
    dx_det = grid.dx_det;
    dy_dxi = grid.dy_dxi;
    dy_det = grid.dy_det;

    % Advection velocities
    alphax = max( max( abs( 2*pi*XI.*cos(ET) ) ) );
    alphay = max( max( abs( 2*pi*XI.*sin(ET) ) ) );    
    alpha  = max(alphax,alphay);
    
    % Flux interface values
    Fp = zeros(mxi+1,met);
    Gp = zeros(mxi,met+1);

    % Fluxes (on the curvilinear mesh)
    % NOTE: These are qt values that are passed in NOT q values.  You need to convert these
    % in whatever flux you're using here!!
    %
    %    ft = @(Q,xi,et)( 2.*pi*Q*( dy_det(xi,et)*xi*sin(et) - dx_det(xi,et)*(-xi*cos(et)) ) );
    ft = @(Q,xi,et)( 0 );
    % gt = @(Q,xi,et)( 2.*pi*Q*(-dy_dxi(xi,et)*xi*sin(et) + dx_dxi(xi,et)*(-xi*cos(et)) ) );
    gt = @(Q,xi,et)( -2.*pi*Q );
    
    % ---- Time stepping loop ---- %
    Tend = t+dT;
    while( t < Tend )

        % TODO - choose the time step
        dt = min([Tend-t, cfl_des/sqrt( (alphax/dxi)^2 + (alphay/det)^2 ) ] );

        % Fluxes in the xi-direction
        for j=1:met
            
            for i=2:mxi
                Ql = qt(i-1,j); Qr = qt(i,j);
                xie = 0.5*(XI(i-1,j)+XI(i,j));
                ete = 0.5*(ET(i-1,j)+ET(i,j));
                Fp(i,j) = 0.5*( ft(Ql,xie,ete) + ft(Qr,xie,ete) ) - 0.5*alphax*(Qr-Ql);
            end

            % Boundary conditions:
            Ql = qt(1,j); Qr = qt(1,j);
            xie = XI(1,j)-0.5*dxi;  ete = ET(1,j);
            Fp(1,j) = 0.5*( ft(Ql,xie,ete) + ft(Qr,xie,ete) ) - 0.5*alphax*(Qr-Ql);

            Ql = qt(mxi,j);  Qr = Ql;
            xie = XI(mxi,j)+0.5*dxi;  ete = ET(mxi,j);            
            Fp(mxi+1,j) = 0.5*( ft(Ql,xie,ete) + ft(Qr,xie,ete) ) - 0.5*alphax*(Qr-Ql);
            
        end

        % Fluxes in the eta-direction
        for i=1:mxi
            
            for j=2:met
                Ql = qt(i,j-1); Qr = qt(i,j);
                xie = 0.5*(XI(i,j-1)+XI(i,j));
                ete = 0.5*(ET(i,j-1)+ET(i,j));
                Gp(i,j) = 0.5*( gt(Ql,xie,ete) + gt(Qr,xie,ete) ) - 0.5*alphay*(Qr-Ql);
            end

            % Boundary conditions:
            Ql = qt(i,met); Qr = qt(i,1);
            xie = XI(i,1);  ete = ET(i,1)-0.5*det;
            Gp(i,1) = 0.5*( gt(Ql,xie,ete) + gt(Qr,xie,ete) ) - 0.5*alphay*(Qr-Ql);
            Gp(i,met+1) = Gp(i,1);
            
        end

        % Update the solution
        for i=1:mxi
            for j=1:met
                qt(i,j) = qt(i,j) - (dt/dxi)*(Fp(i+1,j)-Fp(i,j)) - (dt/det)*(Gp(i,j+1)-Gp(i,j));
            end
        end
        t = t+dt;
        
    end
    
end
