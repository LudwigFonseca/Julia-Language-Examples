function retro_substituicao(A, b)
    n = size(A)[2]
    x = zeros(n)
    x[n] = b[n]/A[n,n]
    for i = n-1:-1:1
        temp = 0
        for k = i+1:n
            temp += A[i,k]*x[k]
        end
        x[i] = (b[i]-temp)/A[i,i]
    end
    return x
end


function substituicao(A, b)
    n = size(A)[2]
    x = zeros(n)
    x[1] = b[1]/A[1,1]
    for i = 2:n
        temp = 0
        for k = 1:i-1
            temp += A[i,k]*x[k]
        end
        x[i] = (b[i]-temp)/A[i,i]
    end
    return x
end


function eliminacao(A)
    row, col = size(A)
    for j in 1:row-1
        for i in j+1:row
            m = A[i,j]/A[j,j]
            A[i,:] = A[i,:]-m*A[j,:]
        end
    end
    return A
end

function fatoracao_lu(A)
    row, col = size(A)
    L = eye(row)

    for j in 1:row-1
        for i in j+1:row
            L[i,j] = A[i,j]/A[j,j]
            A[i,:] = A[i,:]-L[i,j]*A[j,:]
        end
    end

    U = A[:, 1:end-1]
    return L, A
end

function resolve_fatlu(A)
    for i in 6:8
        L, U = fatoracao_lu(A[:, [1:5; i]])
        println("A fatoração LU estará coreta se L*U for igual a A ")
        println(L*U, A)

        y = substituicao(L, b)
        x = retro_substituicao(U, y)
        println("Valor do y é $y o valor do x é $x")
    end
end
