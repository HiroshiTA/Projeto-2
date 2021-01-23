using LinearAlgebra
# 10x + 3y - 2z = 57
# 2x + 8y -1z = 20
# x + y + 5z = -4
#A = [10 3 -2; 2 8 -1; 1 1 5]
#B = [57, 20, -4]
#C = [0, 0, 0]

function Método_Jacobi(A, B, C, max_iter = 100, ϵ= 1e-3) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)  
    i = 1
    j = 1
    k = 1     #número de iterações
    v = zeros(0)    #vetor que recebe os x1, x2, ..., xn.
    erro = 1
    while (k <= max_iter) && (erro > ϵ)
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
        end
        #println(v)
              
        #Para achar o maior dos x dentro do vetor v
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        #println("maior x = $maiorx")
        
        #A maior distancia entre as soluções dos x".
        distancia = abs.(v - C)
        #println("|x^(k) - x^(k-1)| = $distancia")
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx          #Erro Relativo
        #println("Erro Relativo = $erro")
        C = v                           #Atualização da nova matriz: x^k-1 = x^k.
        #println("x^k = $C")
        k = k + 1
    end
    return v
end


function Método_Gauss_Seidel(A, B, C, max_iter = 100, ϵ = 1e-3) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)
    k = 1    #número de iterações 
    v = zeros(0)
    erro = 1.0
    while (k <= max_iter) && (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)   #vetor que recebe os x1, x2, ..., xn.
        C2 = copy(C)
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
            println(a)
        end
        println(v)
        
        #Para achar o maior dos x dentro do vetor v
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        println("maior x = $maiorx")
        println(C2)
        
        #A maior distancia entre as soluções dos x".
        
        distancia = abs.(v - C2)
        println("|x^(k) - x^(k-1)| = $distancia")
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx          #Erro Relativo
        println("Erro Relativo = $erro")
        C = v                           #Atualização da nova matriz: x^k-1 = x^k.
        println("x^k = $C")
        k = k + 1
        println("")
    end
    return v
end


function Método_SOR(A, B, C, max_iter = 100, ω = 1.25, ϵ = 1e-4) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)
    k = 1    #número de iterações 
    v = zeros(0)
    erro = 1.0
    while (k <= max_iter) && (erro > ϵ)
        i = 1
        j = 1
        v = zeros(0)   #vetor que recebe os x1, x2, ..., xn.
        C2 = copy(C)
        while i <= m
            a = A[i,j]        #elementos da Diagonal Principal
            b = B[j]     
            E = A[i,:]
            R = b - dot(E,C) 
            x = C[i] + (ω / a) * R
            push!(v,x)
            C[i] = x
            i = i + 1
            j = j + 1
            println(R)
        end
        
        #Para achar o maior dos x dentro do vetor v
        maiorx = abs(v[1])
        for h = 2:length(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        println("maior x = $maiorx")
        
        #A maior distancia entre as soluções dos x".
        
        distancia = abs.(v - C2)      
        maiord = distancia[1]         
        for l = 2:length(distancia)
            if maiord < abs(distancia[l])
                maiord = abs(distancia[l])
            end
        end
        
        erro = maiord / maiorx          #Erro Relativo
        println("Erro Relativo = $erro")
        C = v                           #Atualização da nova matriz: x^k-1 = x^k.
        println("x^k = $C")
        k = k + 1
        println("")
    end
    return v
end



function converge(A, metodo)
    
    # 1 = Jacobi
    # 2 = Gauss-Seidel
    # 3 = SOR
    
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
            println("Pelo critério das linhas a matriz converge")
        else
            println("Pelo critério das linhas a matriz não converge")
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
            println("Pelo critério de Sassenfeld a matriz converge")
        else
            println("Pelo critério de Sassenfeld a matriz não converge")
        end
    end
    if soma != 0
        println("Conclusão: Converge")
    else
        println("Conclusão: Pode ser que não converga")
    end
end



