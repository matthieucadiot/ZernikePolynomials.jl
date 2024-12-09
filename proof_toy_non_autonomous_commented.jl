using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, MATLAB, JLD2

include("C:/Users/Matthieu/Desktop/code/julia/Navier_Stokes/list_functions.jl")

# Newton's method for solving a nonlinear system related to a PDE
# Inputs:
#   U: Initial guess for the solution
#   MMT, iMMT: Forward and inverse transformation matrices
#   N: Number of modes/resolution parameter
#   k: Jacobi parameter
#   m: Wave number/mode parameter
#   precision: Convergence tolerance
function Newton(U, MMT, iMMT, N, k, m, precision)
    # Compute inverse Laplacian operator with Dirichlet boundary conditions
    Δ⁻¹ = mid.(iLaplacian_Dirichlet(N, N, m))
    
    # Initialize operators for mode conversion and reduction
    R_m = 1.0*I[1:N+1,1:N+1]  # Reduction operator
    C_m = 1.0*I[1:N+1,1:N+1]  # Conversion operator
    
    # Build composite operators through recursive application
    for j=0:m-1
        R_m = mid.(R_minus(k,2*m-j,N,N,0))*R_m
        C_m = mid.(conversion_m(k,m+j,N,N,0))*C_m
    end
    
    # Compute nonlinear term (U^2 in transformed space)
    U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))
    
    # Compute F(U) = U + Δ⁻¹*R_m*U^2
    F = U + Δ⁻¹*R_m*U2
    
    # Compute Jacobian DF = I + 2Δ⁻¹*R_m*diag(U)*C_m
    DF = I + 2*Δ⁻¹*R_m*mid.(multiplication_0(C_m*U,N,N,MMT,iMMT))*C_m
    nf = norm(F)

    # Newton iteration loop
    k = 0
    while (nf>precision)&&(k<30) 
        U = U - DF\F  # Newton step
        # Recompute F and DF at new point
        U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))
        F = U + Δ⁻¹*R_m*U2
        DF = I + 2*Δ⁻¹*R_m*mid.(multiplication_0(C_m*U,N,N,MMT,iMMT))*C_m
        nf = norm(F)
        k+=1
    end
    return U,nf
end

# Compute F(U) for the toy problem 
function F_toy(U, MMT,iMMT,N,k,m)
   # Compute r^m operator and conversion operator
    R_m = 1.0*I[1:N+1,1:N+1]
    C_m = 1.0*I[1:N+1,1:N+1]
    for j=0:m-1
        R_m = mid.(R_minus(k,2*m-j,N,N))*R_m
        C_m = mid.(conversion_m(k,m+j,N,N))*C_m
    end
    # computation of the Dirichlet inverse laplacian
    Δ⁻¹ = mid.(iLaplacian_Dirichlet(N,N,m)) 
    U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))
    return U + Δ⁻¹*R_m*U2
end

# Main function for computer-assisted proof using radii polynomials
# This implements a rigorous verification of solutions using interval arithmetic
function proof_R_to_the_m(U,N,k,m,precis,eps_N)
    # Set precision for interval arithmetic
    setprecision(precis)
    
    # Get transformation matrices
    MMT, iMMT = passage_matrix(2*N+m+1,k,2*m,precis,eps_N)

    # Convert inputs to interval arithmetic
    U = interval.(U)
    Ubig = interval.(big.(U))
    U2 = [U;interval.(zeros(N+m+1))]
    U2big = [Ubig;interval.(big.(zeros(N+m+1)))]

    # Build operators in interval arithmetic
    R_m = interval.(big.(1.0*I[1:2*N+m+2,1:2*N+m+2]))
    C_m = interval.(big.(1.0*I[1:2*N+m+2,1:2*N+m+2]))
    for j=0:m-1
        R_m = R_minus(k,2*m-j,2*N+m+1,2*N+m+1,1)*R_m
        C_m = conversion_m(k,m+j,2*N+m+1,2*N+m+1,1)*C_m
    end

    # Compute inverse Laplacian and derivative operators
    Δ⁻¹ = iLaplacian_Dirichlet(N,N,m) 
    DG = interval(2)*multiplication_0(C_m*U2,2*N+m+1,2*N+m+1,interval.(Float64.(inf.(MMT),RoundDown),Float64.(sup.(MMT),RoundUp)),interval.(Float64.(inf.(iMMT),RoundDown),Float64.(sup.(iMMT),RoundUp)))*C_m
    DG = interval.(Float64.(inf.(DG),RoundDown),Float64.(sup.(DG),RoundUp))

    # Compute inverse of I + Δ⁻¹*R_m*DG
    A = interval.(inv(I + mid.(Δ⁻¹)*mid.(R_m[1:N+1,1:N+1])*mid.(DG[1:N+1,1:N+1])))

    # Compute Y0 bound (residual evaluation)
    UU = iMMT*((MMT*C_m*U2big).*(MMT*C_m*U2big))
    Y0  = norm(A*(U+ iLaplacian_Dirichlet(2*N+m+1,N,m)*(R_m*UU)),1) 
    Y01 =  norm((iLaplacian_Dirichlet(2*N+m+1,2*N+m+1,m)*R_m*UU)[N+2:end],1)
    Y0 = Y0 + Y01

    display("Y0 bound")
    display(Y0)

    # Compute Z1 bound (linear part)
    UU = interval.(Float64.(inf.(UU),RoundDown),Float64.(sup.(UU),RoundUp))
    Z1 = maximum([opnorm(LinearOperator(interval.(1.0*I[1:N+1,1:N+1]) - A*( interval.(1.0*I[1:N+1,1:N+1]) +iLaplacian_Dirichlet(N+1,N,m)*R_m[1:N+2,1:N+2]*DG[1:N+2,1:N+1])),1)  opnorm(LinearOperator(A*iLaplacian_Dirichlet(N+1,N,m)*R_m[1:N+2,1:N+2]*DG[1:N+2,N+2:end]),1)]) + interval(2)/((interval(2)*interval(N)+1+interval(abs(m)))^2)*norm(U,1)

    display("Z1 bound")
    display(Z1)

    # Compute Z2 bound (nonlinear part)
    Z2 = interval(2)*(opnorm(LinearOperator(A*iLaplacian_Dirichlet(N+1,N,m)*R_m[1:N+2,1:N+2]),1) + interval(2)/((interval(2)*interval(N)+1+interval(abs(m)))^2))

    display("Z2 bound")
    display(Z2)

    # Compute radii polynomial coefficients
    α = interval(sup(interval(1)/interval(2) + Y0*Z2/(interval(1)-Z1)^2))
    r = interval(2)*α*Y0/(interval(1)-Z1)

    # Check verification conditions
    if inf(1- Z1)>0
        if sup(Z2*r+Z1) < 1
            if inf(1/2*Z2*r^2 - (1-Z1)*r + Y0) < 0
                display("The computer-assisted proof was successful")
                display("Radius of contraction :")
                display(r)
            else
                display("failure: discriminant is negative")
            end
        else
            display("r is too big")
        end
    else
        display("failure: 1-z > 0")
    end
end

# The rest of the code contains test cases for different values of m (wave numbers)
# Each test case:
# 1. Sets up the problem parameters (N, k, m)
# 2. Creates transformation matrices
# 3. Initializes a guess solution
# 4. Runs Newton iteration to find a numerical solution
# 5. Attempts to prove existence of a true solution near the numerical one





################# Proofs and data for toy problems ########################


########### m=0 data #####################################################

N = 36
k=0
m=0
MMT, iMMT = passage_matrix(N,k,m,128,1e-30)
MMT = mid.(MMT)
iMMT = mid.(iMMT)

display("MMT and iMMT created")

U = zeros(N+1)
U[1] = 3.4
U[2] = -3.4
U[3] = 0.5
U[4] = -0.5
U,nf = Newton(20*U, MMT,iMMT,N,k,m,1e-30)

display(norm(U,1))

# r = 0:0.001:1
# theta = 2π*r
# Z = zeros(length(r),length(theta))

# for i = 1:length(theta)
#     Z[:,i] = evaluation(U,k,m,N,r)*cos(m*theta[i])
# end

# mat"
#     [T,R] = meshgrid($theta,$r);
#     [X,Y,Z] = pol2cart(T,R,$Z);
#     h = surf(X,Y,Z);
#     set(h,'LineStyle','none')"

    precis = 256 ; eps_N = 1e-30
    proof_R_to_the_m(U,N,k,m,precis,eps_N)  

########### m=1 data #####################################################
N = 36
k=0
m=1
MMT, iMMT = passage_matrix(N,k,2*m,128,1e-16)
MMT = mid.(MMT)
iMMT = mid.(iMMT)
# MMT = Float64.(MMT) ; iMMT = Float64.(iMMT)

display("MMT and iMMT created")
U = zeros(N+1)
U[1] = 0.89
U[2] = -0.89
U[3] = 0.24
U[4] = -0.24
U,nf = Newton(20*U, MMT,iMMT,N,k,m,1e-20)

display(norm(U,1))

# r = 0:0.001:1
# theta = 2π*r
# Z = zeros(length(r),length(theta))

# for i = 1:length(theta)
#     Z[:,i] = evaluation(U,k,m,N,r)*cos(m*theta[i])
# end


# mat"
#     [T,R] = meshgrid($theta,$r);
#     [X,Y,Z] = pol2cart(T,R,$Z);
#     h = surf(X,Y,Z);
#     set(h,'LineStyle','none')"


precis = 256 ; eps_N = 1e-30
proof_R_to_the_m(U,N,k,m,precis,eps_N)  

########################################################################


########### m=2 data #####################################################
N = 36
k=0
m=2
MMT, iMMT = passage_matrix(N,k,2*m,128,1e-16)
MMT = mid.(MMT);
#  MMT = Float64.(MMT)
iMMT = mid.(iMMT); 
# iMMT = Float64.(iMMT)

display("MMT and iMMT created")
U = zeros(N+1)
U[1] = 0.94
U[2] = -0.94
U[3] = 0.36
U[4] = -0.36
U,nf = Newton(50*U, MMT,iMMT,N,k,m,1e-20)

display(norm(U,1))

# r = 0:0.001:1
# theta = 2π*r
# Z = zeros(length(r),length(theta))

# for i = 1:length(theta)
#     Z[:,i] = evaluation(U,k,m,N,r)*cos(m*theta[i])
# end


# mat"
#     [T,R] = meshgrid($theta,$r);
#     [X,Y,Z] = pol2cart(T,R,$Z);
#     h = surf(X,Y,Z);
#     set(h,'LineStyle','none')"


precis = 256 ; eps_N = 1e-30 ; 
proof_R_to_the_m(U,N,k,m,precis,eps_N)  

########################## m=20 ################################################

# For the case m=20, more coefficients are needed and consequently the proof might take a longer time 


# k=0 ; m=20
# eps_N = 1e-30
# precis = 256
#  U = load("U_20.jld2","U")
# N = 75 ; V = zeros(N+1) ; V[1:size(U)[1]] = U

# MMT, iMMT = passage_matrix(N,k,2*m,precis,eps_N)
# MMT = mid.(MMT); MMT = Float64.(MMT)
# iMMT = mid.(iMMT); iMMT = Float64.(iMMT)

# V,nf = Newton(V, MMT,iMMT,N,k,m,1e-15)

# proof_R_to_the_m(V,N,k,m,precis,eps_N) 