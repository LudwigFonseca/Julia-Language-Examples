#Método de Euler
function Euler(df, f0, ti, tf, n)
    f = Array{Float64}(n)
    h = (tf - ti) / n
    f[1] = f0 + h*df(f0, ti)
    tk = ti + h
    for k = 1:(n-1)
        f[k + 1] = f[k] + h*df(f[k], tk)
        tk = tk + h
    end
    return vcat(f0, f)
end

#Método de Euler Modificado ou RUnge Kutta 2ª Orde
function Runge_Kutta_2(df, f0, ti, tf, n)
    f = Array{Float64}(n)
    h = (tf - ti) / n
    k1 = df(f0, ti)
    k2 = df(f0 + k1*h/2, ti + h/2)
    f[1] = f0 + k2*h
    tk = ti + h
    for k = 1:(n-1)
        k1 = df(f[k], tk)
        k2 = df(f[k] + k1*h/2, tk + h/2)
        f[k+1] = f[k] + k2*h
        tk = tk + h
    end
    return vcat(f0, f)
end

#Método de RUnge Kutta 3ª Ordem
function Runge_Kutta_3(df, f0, ti, tf, n)
    f = Array{Float64}(n)
    h = (tf - ti) / n
    k1 = h*df(f0, ti)
    k2 = h*df(f0 + k1/2, ti + h/2)
    k3 = h*df(f0 + k2*3/4, ti + 3*h/4)
    f[1] = f0 + (2/9)*k1 + (1/3)*k2 + (4/9)*k3
    tk = ti + h
    for k = 1:(n-1)
        k1 = h*df(f[k], tk)
        k2 = h*df(f[k] + k1/2, tk + h/2)
        k3 = h*df(f[k] + k2*3/4, tk + h/2)
        f[k + 1] = f[k] + (2/9)*k1 + (1/3)*k2 + (4/9)*k3
        tk = tk + h
    end
    return vcat(f0, f)
end

#Método de RUnge Kutta 4ª Ordem
function Runge_Kutta_4(df, f0, ti, tf, n)
    f = Array{Float64}(n)
    h = (tf - ti) / n
    k1 = h*df(f0, ti)
    k2 = h*df(f0 + k1/2, ti + h/2)
    k3 = h*df(f0 + k2/2, ti + h/2)
    k4 = h*df(f0 + k3, ti + h)
    f[1] = f0 + (k1 + 2*(k2 + k3) + k4)*(1/6)
    tk = ti + h
    for k = 1:(n-1)
        k1 = h*df(f[k], tk)
        k2 = h*df(f[k] + k1/2, tk + h/2)
        k3 = h*df(f[k] + k2/2, tk + h/2)
        k4 = h*df(f[k] + k3, tk + h)
        f[k + 1] = f[k] + (k1 + 2*(k2 + k3) + k4)*(1/6)
        tk = tk + h
    end
    return vcat(f0, f)
end

using Plots

#Equação Diferencial Ordinária
df(y, t) = -t/y
f0 = 20
ti = 0
tf = 16
n = 4

#Comparação dos Métodos Com a Resposta Analítica
y = collect(ti:(tf-ti)/n:tf)
y1(y) = 0.1*exp.(2*y)

y2 = Euler(df, f0, ti, tf, n)
y3 = Runge_Kutta_2(df, f0, ti, tf, n)
y4 = Runge_Kutta_3(df, f0, ti, tf, n)
y5 = Runge_Kutta_4(df, f0, ti, tf, n)
plot(y, y1,  label="Analítico")
plot!(y, y2, label="Euler")
plot!(y, y3, label="Runge Kutta 2")
plot!(y, y4, line=(:dot,3), label="Runge Kutta 3")
plot!(y, y5, line=(:dot,4), label="Runge Kutta 4")

#Análise do Erro Relativo
erro1 = abs(y1(y) - y2)
erro2 = abs(y1(y) - y3)
erro3 = abs(y1(y) - y4)
erro4 = abs(y1(y) - y5)
plot(y, erro1, label="Erro Euler")
plot!(y, erro2, label="Erro Runge Kutta 2")
plot!(y, erro3, line=(:dot,3), label="Erro Runge Kutta 3")
plot!(y, erro4, line=(:dot,4), label="Ero Runge Kutta 4")
