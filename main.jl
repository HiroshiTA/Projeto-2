# Dado o sistema: Ax = B

using LinearAlgebra

function metodo_jacobi(A, B, C, max_iter = 100, ϵ= 1e-5) # C: vetor do chute inicial.
    m,n = size(A)  
    k = 1                           #número de iterações
    v = zeros(0)                   
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)                #vetor que recebe os x1, x2, ..., xn, en cada iteração.
        while i <= m
            a = A[i,j]
            b = B[j]
            o = (b/a) 
            E = A[i,:]              # E: Matriz A com os elementos da diagonal valendo 0.
            E[i] = 0
            n = (dot(E,C))/a     
            x = o - n
            push!(v, x)
            i = i + 1
            j = j + 1
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
    println("Números de Iterações = $k")
    println("Vetor solução = $v")
end


function metodo_gauss_seidel(A, B, C, max_iter = 100, ϵ = 1e-5)
    m,n = size(A)
    k = 1
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)
        C2 = copy(C)
        while i <= m
            a = A[i,j]
            b = B[j]
            o = (b/a)     
            E = A[i,:]
            E[i] = 0
            p = (dot(E,C))/a     
            x = o - p
            push!(v,x)
            C[i] = x                #Atualização de cada x_i
            i = i + 1
            j = j + 1
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
    println("Números de Iterações = $k")
    println("Vetor solução = $v")
end


function metodo_SOR(A, B, C, ω = 1.25, max_iter = 100, ϵ = 1e-5)        # w : coeficiente de relaxamento.
    m,n = size(A)
    k = 1
    erro = 1.0
    while (k < max_iter) && (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)
        C2 = copy(C)
        while i <= m
            a = A[i,j]
            b = B[j]     
            E = A[i,:]
            R = b - dot(E,C)                #Resíduo
            x = C[i] + (ω / a) * R
            push!(v,x)
            C[i] = x
            i = i + 1
            j = j + 1
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
    println("Números de Iterações = $k")
    println("Vetor solução = $v")
end



function converge(A, metodo)
    
    # 1 = Jacobi
    # 2 = Gauss-Seidel
    
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
