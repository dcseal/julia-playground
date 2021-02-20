"""
Sample 2D advection solver on curvilinear coordinates

See: https://docs.julialang.org/en/v1/manual/documentation/
- I need to properly document all the routines I call in this code
"""

using ConfParser
using Printf
include("solver.jl")

#function read_params_ini()
#: https://github.com/JuliaIO/ConfParser.jl

# Parse the input parameters
conf = ConfParse("params.ini")
parse_conf!(conf)

# get and store config parameters
const mx = parse( Int, retrieve(conf, "grid", "mx" ) )
const my = parse( Int, retrieve(conf, "grid", "my" ) )

const xlow  = parse( Float64, retrieve(conf, "grid", "xlow" ) )
const xhigh = parse( Float64, retrieve(conf, "grid", "xhigh" ) )    
const ylow  = parse( Float64, retrieve(conf, "grid", "ylow" ) )
const yhigh = parse( Float64, retrieve(conf, "grid", "yhigh" ) )

const Tf = parse( Float64, retrieve(conf, "default", "Tfinal" ) )
const cfl_des = parse( Float64, retrieve(conf, "default", "cfl_des" ) )
const nout = parse( Float64, retrieve(conf, "default", "nout" ) )

const outputdir = string(retrieve(conf, "default", "outputdir" ), ".", string(mx), ".", string(my) )

# ---- Derived Parameters ---- #
const dx = (xhigh-xlow)/mx
const dy = (yhigh-ylow)/my
const dxvec = [dx dy]

# Grid cell centers
const x  = range(xlow+0.5*dx, xhigh-0.5*dx, length = mx )
const y  = range(ylow+0.5*dy, yhigh-0.5*dy, length = my )

# define our curvilinear coordinates [or computational graph] values
#const mxi = mx
#const met = my

# Computational grid
#xi_min = 0; xi_max = 1
#et_min = 0; et_max = 1

# derived parameters
#const dxi = (xi_max-xi_min)/mxi
#const det = (et_max-et_min)/met

#=
# Grid cell centers
const xi  = range(xi_min+0.5*dxi, xi_max-0.5*dxi, length = mxi )
const et  = range(et_min+0.5*det, et_max-0.5*det, length = met )

# Grid edges
const xievec  = range(xi_min, xi_max, length = mxi+1 )
const etevec  = range(et_min, et_max, length = met+1 )

# Physical grid (defiend through a mapping)

# perturbation values
const ex = 0.01
const ey = 0.02
const ax = 2
const ay = 4

# define our physical grid()
xval = (xi,eta) -> xi  + ex*sin(2*pi*ax*eta)
yval = (xi,eta) -> eta + ey*sin(2*pi*ay*xi )

# TODO - fix this
jacob = (xi,eta) -> 1.0

# Transformation derivatives:
dx_dxi = (xi,et) -> ( 1.0*ones(size(xi)) )
dx_det = (xi,et) -> ( (2.0*pi*ax)*ex*cos(2.0*pi*ax*et) )
dy_dxi = (xi,et) -> ( (2.0*pi*ay)*ey*cos(2.0*pi*ay*xi) )
dy_det = (xi,et) -> ( 1.0*ones(size(xi)) )

# define our local Jacobian value Jt & qt [qt representing q tilda]
qt = q
for i=1:mxi
    for j=1:met
        qt[i,j] = q[i,j] * jacob( xi[i], et[j] )
    end
end

dxivec = [dxi,det]
=#

# problem specific parameters; in the physical domain
u1 = 2.0;      u2 = 1.0
uvec = [u1 u2]

# Initial conditions
q = zeros(Float64,mx,my)

# Smooth (periodic) initial condition
# qinit = .1 + sin(4*pi.*XX) + cos(2*pi.*YY)

# Cosine bump function
xc = 0.4
yc = 0.5
for i = 1:mx
    for j = 1:my
        #        rc = sqrt( (xval(xi[i],et[j]) - xc )^2  + (yval(xi[i],et[j])-yc)^2 )
        rc = sqrt( (x[i]-xc)^2 + (y[j]-yc)^2 )
        if ( rc < 0.3 )
            q[i,j] = cos(5/3*pi*rc)^6
        end
    end
end

# ---- Main time stepping loop ---- #
t = 0
dT = Tf/nout
output_soln(t,q,0,outputdir)  # Initial conditions
for nf = 1:nout

    global q = cc_advec_solve(t,q,t+dT,dxvec,uvec)
    global t  = t+dT     # See: https://docs.julialang.org/en/v1/manual/variables-and-scoping/
    # q  = qt.*Jt

    # Save the current frame
    output_soln(t,q,Int(nf),outputdir)
    @printf("---- Finished frame number %03d of %03d ----\n", nf, nout )
    
end
