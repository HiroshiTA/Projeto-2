using LinearAlgebra
# 10x + 3y - 2z = 57
# 2x + 8y -1z = 20
# x + y + 5z = -4
# chute inicial: C = (x,y,z)=(0,0,0)
# x = 57/10 + (- 3y + 2z)/10
# y = 20/8 + (- 2x + z)/8
# z = -4/5 + (- x - y)/5
# Ax = B
# A = 100 3 - 2 
#     2  8 -1 
#     1  1  5
# B = [ 57  20  -4]
# Recebe a matriz A(coeficientes),vetor B(soluções) e o vetor C (chute inicial).Vetor X = [x1 ...xn] x1 
# (primeiro elemento da matriz X) é igual ao primeiro elemento da matriz
# B (b11) dividido pelo coeficiente de x1 (a11) mais o somatório do oposto dos coeficientes (a1j) multiplicados pelo vetor C
# - exceto o primeiro elemento x1, tudo isso dividido pelo coeficiente de a11.

a = A[1,1] = 10
b = B[1] = 57
m = (b/a) = 57/10
A[1,1]=0
n = (dot(A[1,:],C))/a
x1 = m - n

function(A::Matrix, B::Vector, C = zeros(0), max_iter = 100, E = 1e-3)
m,n = size(A)
i = 1
j = 1
C = zeros(m)
v = zeros(m)
D = zeros(m,n)
a = 1 
for k = 1:max_iter
    while i <= m
        a = A[i,j]
        b = B[j]
        m = (b/a)
        A[i,j] = 0
        n = (dot(A[i,:],C))/a
        x = m - n
        push!(x,v)
        i = i + 1
        j = j + 1
    end
    D[:,a] = v
end
return v
#reverse.(v)

