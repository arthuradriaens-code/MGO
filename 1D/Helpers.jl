using LinearAlgebra

function GramSchmidt(A)
    N = size(A)[1]; 
    V = zeros(N,N); 
    V[1,:] = A[1,:]; 
    for i in 2:N
        V[i,:] = A[i,:]
        for j in 1:(i-1)
            V[i,:] -= dot(V[j,:],A[i,:])*V[j,:]/dot(V[j,:],V[j,:]);
        end
    end
    V
end
