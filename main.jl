# Dado o sistema: Ax = B

using LinearAlgebra

function metodo_jacobi(A, B, C, max_iter = 100, ϵ= 1e-5) # C: vetor do chute inicial.
    m,n = size(A)  
    k = 1                           #número de iterações
    v = zeros(0)                   
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        v = zeros(0)                #vetor que recebe os x1, x2, ..., xn, en cada iteração.
        for i = 1:m
            a = A[i,i]
            b = B[i]
            o = (b/a) 
            E = A[i,:]              # E: Matriz A com os elementos da diagonal valendo 0.
            E[i] = 0
            p = (dot(E,C))/a     
            x = o - p
            push!(v, x)
        end
              
        #Para achar o maior dos x dentro do vetor v
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        
        #A maior distancia entre as soluções dos x".
        distancia = abs.(v - C)
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx          #Erro Relativo
        C = v                           #Atualização do vetor C.
        k = k + 1
    end
    return "Números de Iterações = $k", "Vetor solução = $v"
end


function metodo_gauss_seidel(A, B, C, max_iter = 100, ϵ = 1e-5)
    m,n = size(A)
    k = 1
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        v = zeros(0)
        C2 = copy(C)
        for i = 1:m
            a = A[i,i]
            b = B[i]
            o = (b/a)     
            E = A[i,:]
            E[i] = 0
            p = (dot(E,C))/a     
            x = o - p
            push!(v,x)
            C[i] = x                #Atualização de cada x_i
        end
        
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        
        distancia = abs.(v - C2)
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx
        C = v
        k = k + 1
    end
    return "Números de Iterações = $k", "Vetor solução = $C"
end


function metodo_SOR(A, B, C, ω = 2, max_iter = 100, ϵ = 1e-5)        # w : coeficiente de relaxamento.
    m,n = size(A)
    k = 1
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        v = zeros(0)
        C2 = copy(C)
        for i = 1:m
            a = A[i,i]
            b = B[i]     
            E = A[i,:]
            R = b - dot(E,C)                #Resíduo
            x = C[i] + (ω / a) * R
            push!(v,x)
            C[i] = x
        end
    
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        
        distancia = abs.(v - C2)      
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx          #Erro Relativo
        C = v
        k = k + 1
    end
    return "Números de Iterações = $k", "Vetor solução = $C"
end


# Verificação se a solução do sistema converge ou não. Esta função é apenas uma condição suficiente no método.

# metodo = 1 ou 2    
    # 1 = Jacobi
    # 2 = Gauss-Seidel

function converge(A, metodo)    
    soma = 0
    r = 0
    m, n = size(A)
    z = ones(m)
    
    if (metodo == 1) || (metodo == 2)
    #Critério das linhas
        s = 0
        p = 0
        for a = 1:m
            for b = 1:n
                if b == a
                    r = abs(A[a, a])
                else
                    s = s + abs(A[a, b])
                end
                b = b + 1
            end
            if (s/r) >= p
                p = s/r
            end
            a = a + 1
        end
        if p < 1
            soma = soma + 1
            println("Pelo critério das linhas a matriz converge.")
        else
            println("Pelo critério das linhas a matriz não converge.")
        end
    end
    
    if (metodo == 2)
    #Critério de sassenfeld
        s = 0
        p = 0
               
        for a = 1:m
            s = 0
            for b = 1:n
                if b == a
                    r = abs(A[a, b])
                else
                    s = s + z[b]* abs(A[a, b])
                end
            end
            z[a]= s/r
            if (s/r) >= p
                p = s/r
            end
        end
        if p < 1
            soma = soma + 1    
            println("Pelo critério de Sassenfeld a matriz converge.")
        else
            println("Pelo critério de Sassenfeld a matriz não converge.")
        end
    end
    if soma != 0
        println("Conclusão: Converge.")
    else
        println("Conclusão: Pode ser que não convirja.")
    end
end
