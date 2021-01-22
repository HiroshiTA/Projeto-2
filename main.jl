using LinearAlgebra
# 10x + 3y - 2z = 57
# 2x + 8y -1z = 20
# x + y + 5z = -4
#A = [10 3 -2; 2 8 -1; 1 1 5]
#B = [57, 20, -4]
#C = [0, 0, 0]

using LinearAlgebra

function Método_Jacobi(A, B, C, max_iter = 100, ϵ= 1e-3) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)  
    i = 1
    j = 1
    k = 1     #número de iterações
    v = zeros(0)    #vetor que recebe os x1, x2, ..., xn.
    erro = 1
    print(ϵ)
    while (k <= max_iter) || (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)
        while i <= m
            a = A[i,j]                #elementos da Diagonal Principal
            b = B[j]
            o = (b/a) 
            E = A[i,:]
            E[i] = 0
            n = (dot(E,C))/a     
            x = o - n
            push!(v, x)
            i = i + 1
            j = j + 1
            println(a)
        end
        println(v)
              
        #Para achar o maior dos x dentro do vetor v."
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        #Para achar a maior distancia entre x^(k) - x^(k-1)"
        distancia = abs.(v - C)
        maiord = distancia[1] 
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiord          #Erro Relativo
        C = v                           #Atualização da matriz para a proxima iterção: x^(k-1) = x^(k)
        k = k + 1
    end
    return v
end


function Método_Gauss_Seidel(A, B, C, max_iter = 100, E = 1e-3) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)
    k = 1    #número de iterações 
    v = zeros(0)
    erroR = E
    while (k <= max_iter) || (erroR < E)
        i = 1
        j = 1
        v = zeros(0)   #vetor que recebe os x1, x2, ..., xn.
        while i <= m
            a = A[i,j]        #elementos da Diagonal Principal
            b = B[j]
            o = (b/a)     
            E = A[i,:]
            E[i] = 0
            p = (dot(E,C))/a     
            x = o - p
            push!(v,x)
            C[i] = x
            i = i + 1
            j = j + 1
        end
        
        # Acharemos o maior elemento do Vetor de v e a maior distância entre o chute inicial de x e os valores achados.
        maiorx = abs(v[1])               
        distancia = abs.(v - C)
        maiord = distancia[1] 
        
        #Para o maior x
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        
        # Para a maior distância
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        erroR = maiord / maiorx             #Erro Relativo
        C = v
        k = k + 1
    end
    return v
end


