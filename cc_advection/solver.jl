function output_soln(t,qt,nf,outputdir)
#OUTPUT_SOLN.  Output the solution to file

    if( ~isdir(outputdir) )
        mkdir( outputdir )
    end

    if( nf==0 )
        cp("params.ini",string(outputdir,"/params.ini"),force=true)        
        open( string(outputdir, "/parameters_ini_filename"), "w" ) do io
            write(io,"params.ini\n")
        end
    end
    
    # Create an output file - save the time
    fname = string(outputdir,"/q",lpad(nf,4,"0"),".dat")
    open(fname, "w") do io

        # Save the time
        write(io, string(t,"\n") )

        # Save each element of qt
        for i=1:length(qt)
            write(io, string(qt[i],"\n") )
        end

    end

    return 0
    
end

function cc_advec_solve(t,qt,Tend,dxvec,uvec)
#CC_ADVEC_SOLVE - single time step on the constant coefficient advection equation:

    # Grid parameters
    dx = dxvec[1]
    dy = dxvec[2]

    u1 = uvec[1]
    u2 = uvec[2]
    
#   # Grab the number of elements
#   if( length( size(qt) ) > 2 )
#       (mx,my,meqn) = size(qt)
#   else
#       (mx,my) = size(qt)
#       meqn = 1
#   end
    
    # Storage for the fluxes
    Fp = zeros(mx+1,my)
    Gp = zeros(mx,my+1)

    
    # TOOD - REPLACE WITH ACTUAL PARAMETERS HERE
    alphax = abs(u1)
    alphay = abs(u2)
    while( t < Tend )

        dt = minimum([Tend-t, cfl_des/sqrt( (alphax/dx)^2 + (alphay/dy)^2 ) ] );
#        println("cfl_des",cfl_des)
        
        # Fluxes in the xi-direction
        Threads.@threads for j=1:my
#        for j=1:my

            for i=2:mx
                Ql = q[i-1,j]; Qr = q[i,j]
#               xie = 0.5*(XI[i-1,j]+XI[i,j])
#               ete = 0.5*(ET[i-1,j]+ET[i,j])
#                Fp[i,j] = 0.5*( ft[Ql,xie,ete] + ft[Qr,xie,ete] ) - 0.5*alphax*(Qr-Ql)
                Fp[i,j] = 0.5*u1*( Ql + Qr ) - 0.5*alphax*(Qr-Ql)
                
            end

            # Boundary conditions:
            Ql = q[mx,j]; Qr = q[1,j]
#            xie = XI[1,j]-0.5*dxi;  ete = ET[1,j]
            Fp[1,j] = 0.5*u1*( Ql + Qr ) - 0.5*alphax*(Qr-Ql)
            Fp[mx+1,j] = Fp[1,j]

        end

        # Fluxes in the eta-direction
        Threads.@threads for i=1:mx
#        for i=1:mx

            for j=2:my
                Ql = q[i,j-1]; Qr = q[i,j]
#               xie = 0.5*(XI[i,j-1]+XI[i,j])
#               ete = 0.5*(ET[i,j-1]+ET[i,j])
                Gp[i,j] = 0.5*u2*( Ql + Qr ) - 0.5*alphay*(Qr-Ql)
            end

            # Boundary conditions:
            Ql = q[i,my]; Qr = q[i,1]
#            xie = XI[i,1];  ete = ET[i,1]-0.5*det()
            Gp[i,1] = 0.5*u2*( Ql + Qr ) - 0.5*alphay*(Qr-Ql)
            Gp[i,my+1] = Gp[i,1]

        end
        
        # Update the solution
        Threads.@threads for i=1:mx
#        for i=1:mx
            for j=1:my
                q[i,j] = q[i,j] - (dt/dx)*(Fp[i+1,j]-Fp[i,j]) - (dt/dy)*(Gp[i,j+1]-Gp[i,j]);
            end
        end
        
       t = t+dt
       #@printf("t = %f; Tend = %f\n", t, Tend)
   end
    
    return q
    
end
