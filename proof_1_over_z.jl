using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, Plots, MATLAB

# need to include the code list_functions
include("path_to_the_repository/list_functions.jl")

function Newton(U, MMT, iMMT, N, precision)
    # Implementation of Newton's method for solving the nonlinear equation
    # Parameters:
    # U: Initial guess
    # MMT, iMMT: Transformation matrices for spectral method
    # N: Truncation parameter
    # precision: Desired accuracy
    
    # Initialize operators
    Δ⁻¹ = mid.(iLaplacian_Dirichlet(N,N,1))  # Inverse Laplacian operator
    iR = mid.(iR_plus(N,N))                   # Auxiliary operator R⁺
    C_m = mid.(conversion_m(0,1,N,N,0))       # Conversion operator
    
    # Compute initial residual
    U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))     # Nonlinear term
    F = U + Δ⁻¹*iR*U2                         # Function evaluation
    DF = I + 2*Δ⁻¹*iR*mid.(multiplication_0(C_m*U,N,N,MMT,iMMT))*C_m  # Derivative
    nf = norm(F)

    # Newton iteration
    k = 0
    while (nf>precision)&&(k<20) 
        U = U - DF\F  # Newton step
        # Recompute residual and derivative
        U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))
        F = U + Δ⁻¹*iR*U2
        DF = I + 2*Δ⁻¹*iR*mid.(multiplication_0(C_m*U,N,N,MMT,iMMT))*C_m
        nf = norm(F)
        k+=1
    end
    return U
end

function F_toy(U, MMT, iMMT, N)
    # Evaluates the nonlinear operator F
    # This is the main equation we're trying to solve
    
    iR = mid.(iR_plus(N,N))
    C_m = mid.(conversion_m(0,1,N,N))
    Δ⁻¹ = mid.(iLaplacian_Dirichlet(N,N,1)) 
    
    U2 = iMMT*((MMT*C_m*U).*(MMT*C_m*U))  # Quadratic term
    return U + Δ⁻¹*iR*U2
end

function proof_one_over_R(U,N,k,m)
    # Main function implementing the computer-assisted proof
    # Uses interval arithmetic to prove existence of solution
    
    setprecision(164)  # Set high precision for interval computations
    
    # Initialize transformation matrices
    MMT, iMMT = passage_matrix(2*N+1,k,2*m,164,1e-30)
    
    # Convert input to interval arithmetic
    U = interval.(U)
    Ubig = interval.(big.(U))
    U2 = [U;interval.(zeros(N+1))]
    U2big = [Ubig;interval.(big.(zeros(N+1)))]
    
    # Initialize operators
    iR1 = iR_plus(2*N+1,2N+1)
    iR0 = iR_plus_0(2*N+1,2*N+1)
    C_m = conversion_m(0,1,2N+1,2N+1,1)
    
    # Compute bounds for the verification
    
    # Y0 bound - measures the residual
    Δ⁻¹ = iLaplacian_Dirichlet(N,N,m) 
    DG = interval(2)*multiplication_0(C_m*U2,2*N+1,2*N+1,
          interval.(Float64.(inf.(MMT),RoundDown),Float64.(sup.(MMT),RoundUp)),
          interval.(Float64.(inf.(iMMT),RoundDown),Float64.(sup.(iMMT),RoundUp)))*C_m
    DG = interval.(Float64.(inf.(DG),RoundDown),Float64.(sup.(DG),RoundUp))
    
    # Compute approximate inverse and its norm
    A = interval.(inv(I + mid.(Δ⁻¹)*mid.(iR1[1:N+1,1:N+1])*mid.(DG[1:N+1,1:N+1])))
    norm_A = opnorm(LinearOperator(A),1)
     
    # Compute Y0 bound
    UU = iMMT*((MMT*C_m*U2big).*(MMT*C_m*U2big))
    Y0  = norm(A*(U+ iLaplacian_Dirichlet(2*N+1,N,m)*(iR1*UU)),1) + 
          norm((iLaplacian_Dirichlet(2*N+1,2*N+1,m)*iR1*UU)[N+2:end],1)

          display("Y0 bound")
          display(Y0)      
    
    # Compute Z1 bound - measures how close A is to inverse
    UU = interval.(Float64.(inf.(UU),RoundDown),Float64.(sup.(UU),RoundUp))
    MMT0, iMMT0 = passage_matrix(2*N+1,k,1,128,1e-30)
    DG_R_inv = interval(2)*multiplication_0(conversion_m(0,0,2N+1,2N+1,0)*iR0*U2,
                2*N+1,2*N+1,interval.(Float64.(inf.(MMT0),RoundDown),
                Float64.(sup.(MMT0),RoundUp)),interval.(Float64.(inf.(iMMT0),
                RoundDown),Float64.(sup.(iMMT0),RoundUp)))
    
    Z1 = maximum([opnorm(LinearOperator(interval.(1.0*I[1:N+1,1:N+1])-
         A*(interval.(1.0*I[1:N+1,1:N+1])+iLaplacian_Dirichlet(N+1,N,m)*
         iR1[1:N+2,:]*DG[:,1:N+1])),1) opnorm(LinearOperator((A*
         iLaplacian_Dirichlet(N+1,N,m)*DG_R_inv[1:N+2,N+2:end])),1)]) + 
         norm(iR0*U2,1)/(interval(2)*(interval(N)^2))

         display("Z1 bound")
         display(Z1)     
    
    # Compute Z2 bound - related to tail estimates
    Z2 = interval(2)*opnorm(LinearOperator(A*iLaplacian_Dirichlet(N,N,m)*
         iR1[1:N+1,1:N+1]),1) + interval(16)/interval(15)*norm_A/(interval(N+3)) + 
         interval(1)/(interval(2*N))

         display("Z2 bound")
         display(Z2)     
    
    # Final verification step using radii polynomial approach
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

# Main script execution
N = 36   # Truncation parameter
k = 0    # Polynomial parameter
m = 1    # Polynomial parameter



r = 0:0.001:1
theta = 2π*r
Z = zeros(length(r),length(theta))

for i = 1:length(theta)
    Z[:,i] = evaluation(U,k,m,N,r)*cos(m*theta[i])
end


mat"
    [T,R] = meshgrid($theta,$r);
    [X,Y,Z] = pol2cart(T,R,$Z);
    h = surf(X,Y,Z);
    set(h,'LineStyle','none')"


# Initialize transformation matrices
MMT, iMMT = passage_matrix(N,k,2*m,80,1e-16)
MMT = mid.(MMT) ; MMT= Float64.(MMT)
iMMT = mid.(iMMT); iMMT= Float64.(iMMT)

# Execute the computer-assisted proof
U = load("U_1_over_z.jld2","U") # loading the approximate solution
proof_one_over_R(U,N,k,m)  
