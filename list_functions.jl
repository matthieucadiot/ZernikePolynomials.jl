



function eig_multiplication(N,k,m,id)
    if id ==0
        n = interval.(big.(0:N))
    else 
        n = (0:N)
    end
    γ = 2*(n.+1).*(n.+( k + m +1)).*(2*n.+(k+m))
    

    c = 2*(n.+k).*(n.+m).*(2*n.+(k+m+2))
    b = -(k.^2-m.^2)*(2*n.+(k+m+1))
    a = (2*n.+(k+m+1)).*(2*n.+(m+k+2)).*(2*n.+(k+m))
    if (k==0)&&(m==0)
        a[1] = 1
    end

    σ = b./a
    
    θ = sqrt.(γ[1:end-1].*a[1:end-1]).*sqrt.(c[2:end]./a[2:end])./a[1:end-1]

    if (k==0)&&(m==0)
        θ[1] = sqrt(c[2]/a[2])
    end


    return Tridiagonal(θ,σ,θ)
end


function normalization_D(k,N)
    n=0:N
    return Diagonal(interval.(binomial.(n.+k,n)))
end


function newton_eig(B,U,eps_N)
    
    V = coefficients(component(U,2))
    V0 = V 
    v = component(U,1)[1]
    norm_f = 1
    
    k=0
    while (norm_f > eps_N)&&(k<20)
        F = [dot(V-V0,V0);B*V-v*V]
        DF = [0  V0';-V (B-v*I)]
        U = U - DF\F 
        V = U[2:end] 
        v = U[1]
        norm_f = norm(F)
        k= k+1
    end
    return U
end







function proof_eigen(B,U)

    pspace = space(U)

    V = component(U,2) 
    v = component(U,1)[1]

    T = LinearOperator(pspace,pspace, [[0   (coefficients(V))'] ; [-coefficients(V) (B - v*I)]] )
   

# # # # ################ Z BOUND ######################################################

    A = inv(mid.(T)) 
    A = interval.(A)
    norm_A =opnorm(A,1)
    
    
    Z = opnorm(I - A*T,1)
################# Y BOUND ######################################################

    Y = norm(coefficients(A)*[0;B*coefficients(V)-v*coefficients(V)],1)
    
################## Z2 BOUND ######################################################

    Z2 = norm_A
    
################## Verification ###################################################

x = 1;
if 1- sup(Z)>0
    if 1-sup(2*Y*Z2) > 0
      rmin = sup((1-Z - sqrt((1-Z)^2-2*Y*Z2))/Z2)
      rmax = inf((1-Z + sqrt((1-Z)^2-2*Y*Z2))/Z2)
      if rmin<rmax
        
      else
        display("rmin>=rmax")
        x=0
      end
    else
      display("failure: discriminant is negative")
      x=0
    end
  else
      display("failure: 1-z > 0")
      x=0
  end

  return [x==1,rmin]

end


function computation_W(k,m,N,id)
    if id==0
        n = 0:N
        if k==0
            n = interval.(big.(0:N))
            interval(big(2^(m+1)))./(interval(big(2))*n.+interval(big(m+1)))
        else
            return interval(big(2^(k+m+1)))./((interval(big(2))*interval.(big.(0:N))) .+interval(big(k+m+1))).*interval.(factorial.(big.(n.+k))).*interval.(factorial.(big.(n.+m)))./interval.((factorial.(big.(n.+(k+m)))).*interval.(factorial.(big.(n)))) 
        end
    else
        n = 0:N
        if k==0
            n = interval.((0:N))
            interval(2^(m+1))./(interval(2)*n.+interval(m+1))
        else
            return interval(2)^(k+m+1)./(interval(2)*n .+interval(k+m+1)).*factorial.(n.+k).*factorial.(n.+m)./(factorial.(n.+(k+m)).*factorial.(n)) 
        end
    end
end


function proof_full_eigen(E,B,N,pspace,eps_N)
    M = interval.(big.(ones(N+1,N+1)))
    d = interval.(big.(ones(N+1)))
    for m = 1:N+1
    
        v = E.values[m] 
        V = E.vectors[:,m] 
        U = Sequence(pspace,(vec([big(v);big.(V)]))) 
        U = newton_eig(mid.(B),mid.(U),eps_N)
        U = interval.(U)
        
        P,rmin = proof_eigen(B,U)
        
        if P != 1
            break
        end
        M[:,m] = interval.(inf.(coefficients(component(U,2))).-rmin,sup.(coefficients(component(U,2))).+rmin)
        d[m] = interval(inf(component(U,1)[1])-rmin,sup(component(U,1)[1])+rmin)
    end
    return M,d
end





###########################################################################################################################


function passage_matrix(N,k,m,precision_proof,eps_N)

    setprecision(precision_proof)
    pspace = ParameterSpace()×ParameterSpace()^(N+1) 
    B = mid.(eig_multiplication(N,k,m,1))

    E = eigen(Hermitian(Matrix(B)))
    n = interval.(big.((1:2*N)))

    γ = interval(big(2))*(n.+interval(big(1))).*(n.+interval(big(k + m +1))).*(interval(big(2))*n.+interval(big(k+m)))
    c = interval(big(2))*(n.+interval(big(k))).*(n.+interval(big(m))).*(interval(big(2))*n.+interval(big(k+m+2)))./γ
    b = -(interval(big(k)).^2-interval(big(m)).^2)*(interval(big(2))*n.+interval(big(k+m+1)))./γ
    a = (interval(big(2))*n.+interval(big(k+m+1))).*(interval(big(2))*n.+interval(big(m+k+2))).*(interval(big(2))*n.+interval(big(k+m)))./γ
    a = [interval(big(1.0));a]
    b = [interval(big(0.0));b]
    c = [interval(big(0.0));c]

    

    mbig = interval(big(m))
    kbig = interval(big(k))
    B = Matrix(eig_multiplication(N,kbig,mbig,0))
    eigenvectors,eigenvalues = proof_full_eigen(E,B,N,pspace,eps_N)
   

    MMT = interval.(big.(ones(N+1,N+1))) 
    MMT[:,1] = interval.(big.(ones(N+1)))
    MMT[:,2] = (kbig+interval(big(1))).+interval(big(0.5))*(kbig+mbig+interval(big(2)))*(eigenvalues.-interval(big(1)))
    for i = 3:N+1
        MMT[:,i] = (a[i-1]*eigenvalues.-b[i-1]).*MMT[:,i-1] - c[i-1]*MMT[:,i-2]
     end


    # MMT = interval.(big.(ones(N+1,N+1))) 
    # A1 = interval.(big.(zeros(N+1,N+1))) 
    # for i = 1:N+1
    #     eig0 = eigenvalues[i]
    #     A1[1,1] = interval(big(1)) ; A1[2,2] = interval(big(1))
    #     for j = 3:N+1
    #         A1[j,j-2] = c[j-1]
    #         A1[j,j-1] = b[j-1]-a[j-1]*eig0
    #         A1[j,j] = interval(big(1))
    #     end
    #     b1 = interval.(big.(zeros(N+1)))
    #     b1[1] = interval(big(1))
    #     b1[2] = (kbig+interval(big(1))) .+ interval(big(0.5))*(kbig+mbig +interval(big(2)))*(eig0-interval(big(1)))

    #     MMT[i,:] = solve_linear(A1,b1,precision_proof)
    
    # end





    # computation of the wights for the quadrature formula
    ω = interval.(big.(zeros((N+1))))
    for j=0:N
        ω[j+1] = eigenvectors[1,j+1]^2
    end
    ω = interval(big(2^(k+m+1)))*interval(factorial(big(k)))*interval(factorial(big(m)))/interval(factorial(big(k+m+1)))*ω

    # computation of the weights for the orthogonal basis (inner product weights)
    W =  computation_W(k,m,N,0)


    iMMT = (interval(big(1))./W).*(MMT').*(ω')
    # MMT = interval.(Float64.(inf.(MMT),RoundDown),Float64.(sup.(MMT),RoundUp) )
    # iMMT = interval.(Float64.(inf.(iMMT),RoundDown),Float64.(sup.(iMMT),RoundUp) )

    return MMT,iMMT
end


function solve_linear(M,b,precis)
    setprecision(precis)

    x = interval.(big.(mid.(M)\mid.(b))) 
    Minv = interval.((inv((mid.(M)))))
    N = size(b)[1]
    Id = interval.(big.(1.0*I[1:N,1:N]))

    Z1 = opnorm(LinearOperator(Id - M*Minv),Inf)


    if inf(interval(1)-Z1) >0
        Y0 = norm(Minv*(M*x-b),Inf)
        rmin = big.(sup(Y0/inf(1-Z1)))
        
        return interval.(inf.(x).-rmin,sup.(x).+rmin)
    else
        return NaN
    end
end



function iLaplacian_Dirichlet(M,N,m)
    n1 = interval.(0:N)
    n2 = interval.(0:N-1)
    n3 = interval.(1:N)
    an = -interval.(ones(N+1))./(interval(2)*(interval(2)*n1 .+ (interval(m)+interval(2))).*(interval(2)*n1 .+ interval(m)))
    bn = interval.(ones(N))./(interval(4)*(interval(2)*n2.+interval(m+1)).*(interval(2)*n2.+interval(m+2)))
    an[1] = -interval(1)/interval(4*(m+1)*(m+2))
    bn[1] = interval(1)/interval(4*(m+1)*(m+2))
    
    cn = interval.(ones(N))./(interval(4)*(interval(2)*n3.+interval(m)).*(interval(2)*n3.+interval(m+1)))

    L1 = Tridiagonal(bn,an,cn)
        if M == N+1
            en = interval.(zeros(N+1))
            en[end] = interval(1)/(interval(4)*(interval(2)*(interval(N+1))+interval(m))*(interval(2)*(interval(N+1))+interval(m+1)))
            L = [L1 en]
        elseif M > N+1
            en = interval.(zeros(N+1))
            en[end] = interval(1)/(interval(4)*(interval(2)*(interval(N+1))+interval(m))*(interval(2)*(interval(N+1))+interval(m+1)))
            L = [L1 en interval.(zeros(N+1,M-N-1))]
        else
            L = L1
        end
    return L
end



function Laplacian_plus2(M,N,k,m)

    n = interval.(1:N)
    return Tridiagonal(interval.(zeros(N)),interval.(zeros(N+1)),4*(n .+ m).*(n .+ (k+m+1)))

end



function iLaplacian(M,N,k,m,α)

    # need M >= N !!!!!!!!!!
    if N>M
        display("M needs to be greater that N")
        return 1
    end

    n1 = interval.(0:N)
    n2 = interval.(0:N-1)
    n3 = interval.(1:N)
    bn = ones(N+1)./(2*(2*n1 .+ (m+k+2)).*(2*n1 .+ (m+k)))
    an = ones(N)./(4*(2*n2.+(m+k+1)).*(2*n2.+(m+k+2)))
    cn = ones(N)./(4*(2*n3.+(m+k)).*(2*n3.+(m+k+1)))

    bn[1] = 0
    an[1] = 1/(4*(m+1)*(k+m+2))
    an[2] = (k+m+2)/(4*(k+m+4)*(k+m+3)*(m+2))
    bn[2] = 1/(4*(k+m+2)*(k+m+3))+1/(4*(k+m+3)*(k+m+4))
    cn[1] = 0

    L1 = Tridiagonal(an,-bn,cn)
        if M == N+1
            en = interval.(zeros(N+1))
            en[end] = 1/(4*(2*(interval(N)+1)+m+k)*(2*(interval(N)+1)+m+k+1))
            L = [L1 en]
        elseif M > N+1
            en = interval.(zeros(N+1))
            en[end] = 1/(4*(2*(interval(N)+1)+m+k)*(2*(interval(N)+1)+m+k+1))
            L = [L1 en zeros(N+1,M-N-1)]
        else
            L = L1
        end

        n1 = interval.(0:M)
    
        bn = ones(M+1)./(2*(2*n1 .+ (m+k+2)).*(2*n1 .+ (m+k)))
        an = ones(M+1)./(4*(2*n1.+(m+k+1)).*(2*n1.+(m+k+2)))
        cn = ones(M+1)./(4*(2*n1.+(m+k)).*(2*n1.+(m+k+1)))
    
        cn[1] = 0
        bn[1] = 0
        an[1] = 1/(4*(m+1)*(k+m+2))
        an[2] = (k+m+2)/(4*(k+m+4)*(k+m+3)*(m+2))
        bn[2] = 1/(4*(k+m+2)*(k+m+3))+1/(4*(k+m+3)*(k+m+4))
        cn[2] = 0
        
    γ = interval.(zeros(M+1)) 
    γ[1] =   α[2]*an[1] - α[1]*bn[1]
    for j=1:M
        γ[j+1] = α[j+2]*an[j+1] - α[j+1]*bn[j+1] + α[j]*cn[j+1]
    end
    L = Matrix(L)
    L[1,:] = -γ'
    return L
end







function multiplication_0(U,N1,N2,MMT,iMMT)
    N = size(U)[1]-1
    DG = interval.(zeros(N1+1,N2+1))
    if N >= N2
        U = U[1:N2+1]
    else
        U = [U ; interval.(zeros(N2-N))]
    end
    
    for j = 0:N1
        DG[:,j+1] = iMMT*((MMT*U).*(MMT*(interval.(1.0*I[1:N2+1,j+1]))))
    end
    # we can try building I as a matrix and vectorize this computation to be a matrix multiplication only
    return DG
end





function conversion_k(k,m,M,N)
    
    n1 = interval.(0:N)
    n2 = interval.(1:N)
    CN = Tridiagonal(interval.(zeros(N)),(n1.+interval(k+m+1))./(interval(2)*n1.+interval(k+m+1)),-(n2.+interval(m))./(interval(2)*n2.+interval(k+m+1)))
    if N==M
        return CN
    elseif M == N+1
        return [CN [interval.(zeros(N)); -interval(N+1+m)/(interval(2)*interval(N+1)+interval(k+m+1))]]
    else
        return [CN [interval.(zeros(N)); -interval(N+1+m)/(interval(2)*interval(N+1)+interval(k+m+1))] interval.(zeros(N+1,M-N-1))] 
    end
end


function sigma(N,N0)
   
    b = 1.0*ones(N-N0)
    
    return Tridiagonal(b,1.0*zeros(N-N0+1),1.0*zeros(N-N0))
end



function conversion_m(k,m,M,N,id)
    

    if id==0
    n1 = interval.(0:N)
    n2 = interval.(1:N)
    CN = Tridiagonal(interval.(zeros(N)),(n1.+interval(k+m+1))./(interval(2)*n1.+interval(k+m+1)),(n2.+interval(k))./(interval(2)*n2.+interval(k+m+1)))
    if N==M
        return CN
    elseif M == N+1
        return [CN [interval.(zeros(N)); interval(N+1+k)/(interval(2)*interval(N+1)+interval(k+m+1))]]
    else
        return [CN [interval.(zeros(N)); interval(N+1+k)/(interval(2)*interval(N+1)+interval(k+m+1))] interval.(zeros(N+1,M-N-1))] 
    end



    else
        n1 = interval.(big.(0:N))
        n2 = interval.(big.(1:N))
        CN = Tridiagonal(interval.(big.(zeros(N))),(n1.+interval(big(k+m+1)))./(interval(big(2))*n1.+interval(big(k+m+1))),(n2.+interval(big(k)))./(interval(big(2))*n2.+interval(big(k+m+1))))
        if N==M
            return CN
        elseif M == N+1
            return [CN [interval.(big.(zeros(N))); interval(big(N+1+k))/(interval(big(2))*interval(big(N+1))+interval(big(k+m+1)))]]
        else
            return [CN [interval.(big.(zeros(N))); interval(big(N+1+k))/(interval(big(2))*interval(big(N+1))+interval(big(k+m+1)))] interval.(big.(zeros(N+1,M-N-1)))] 
        end


    end
end




function R_plus(k,m,M,N)
    
    n1 = interval.(0:N)
    n2 = interval.(1:N)
    RN = Tridiagonal(interval.(zeros(N)),(n1.+(k+m+1))./(2*n1.+(k+m+1)),(n2.+k)./(2*n2.+(k+m+1)))
    if N==M
        return RN
    elseif M == N+1
        return [RN [zeros(N); (N+1+k)/(2*(N+1)+k+m+1)]]
    else
        return [RN [zeros(N); (N+1+k)/(2*(N+1)+k+m+1)] zeros(N+1,M-N-1)] 
    end
end


function R_minus(k,m,M,N,id)
    if id ==0
        n1 = interval.(0:N)
        n2 = interval.(0:N-1)
        RN = Tridiagonal((n2.+interval(1))./(interval(2)*n2.+interval(k+m+1)),(n1.+ interval(m))./(interval(2)*n1.+interval(k+m+1)),interval.(zeros(N)))

    else
        n1 = interval.(big.(0:N))
        n2 = interval.(big.(0:N-1))
        RN = Tridiagonal((n2.+interval(big(1)))./(interval(big(2))*n2.+interval(big(k+m+1))),(n1.+ interval(big(m)))./(interval(big(2))*n1.+interval(big(k+m+1))),interval.(big.(zeros(N))))

    end
    # if N==M
    #     return RN
    # elseif M == N+1
    #     return [RN [zeros(N); (N+1+k)/(2*(N+1)+k+m+1)]]
    # else
    #     return [RN [zeros(N); (N+1+k)/(2*(N+1)+k+m+1)] zeros(N+1,M-N-1)] 
    # end
    return RN
end



function iR_plus(M,N)
    # compute the inverse of R_plus for the case k=0 and m=1
    R = interval.(zeros(N+1,M+1))
    for i=0:N
        for j=0:M
            if j >= i
                R[i+1,j+1] = interval(2)*interval((-1)^(i+j))*(interval(i)+interval(1))^2/((interval(j)+interval(1))*(interval(j)+interval(2)))
            end
        end
    end
    return R
end



function iR_plus_0(M,N)
    # compute the inverse of R_plus for the case k=0 and m=1
    R = interval.(zeros(N+1,M+1))
    for i=0:N
        for j=0:M
            if j >= i
                R[i+1,j+1] = interval((-1)^(i+j))*(interval(2)*interval(i)+interval(1))/((interval(j)+interval(1)))
            end
        end
    end
    return R
end





function evaluation(U,k,m,N,r)

    n = 1:2*N

    γ = 2*(n.+1).*(n.+( k + m +1)).*(2*n.+(k+m))
    c = 2*(n.+k).*(n.+m).*(2*n.+(k+m+2))./γ
    b = -(k.^2-m.^2)*(2*n.+(k+m+1))./γ
    a = (2*n.+(k+m+1)).*(2*n.+(m+k+2)).*(2*n.+(k+m))./γ
    a = [1;a]
    b = [0;b]
    c = [0;c]
    R = length(r)
    f = zeros(R)
    for i = 1:R
        A1 = zeros(N+1,N+1) 
        A1[1,1] = 1 ; A1[2,2] = 1
            for j = 3:N+1
                A1[j,j-2] = c[j-1]
                A1[j,j-1] = b[j-1]-a[j-1]*(2*r[i]^2-1)
                A1[j,j] = 1
            end
        b1 = zeros(N+1)
        b1[1] = 1
        b1[2] = (k+1) + 1/2*(k+m +2)*(2*r[i]^2-2)
        f[i] = sum((A1\b1).*U)    
    end
    return  (r.^m).*f
end




