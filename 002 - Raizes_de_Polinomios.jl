# medotdo da bisseccao
function met_bis(a,b,f,tol=1.0e-2)
    if a>=b
        error("entre com a < b")
    elseif !(f(a)*f(b) < 0)
        error("entre com a,b tal que f(a)*f(b) < 0")
    end
    itermax = 100
    iter= 0
    xk= (a+b)/2
    while !(abs(f(xk))< tol) && iter < itermax
        if f(xk)*f(a) < 0
            b = xk
        elseif f(xk)*f(b) < 0
            a = xk
        end
        xk=(a+b)/2
        iter += 1 #iter = iter +1
    end
    println("$iter  é o numero de interações")
    if abs(f(xk)) < tol
        println("$xk é zero de f")
    end
    return xk
end

# Método de Newton
function newton(xk , f , df , tol = 1.0e-5)
    itermax = 100;
    iter= 0;
    while !(abs(f(xk))< tol) && iter < itermax
        xk = xk - (f(xk)/df(xk))
        iter += 1
        println(iter ,":", f(xk))
    end
    println("$iter  é o numero de interações")
    if abs(f(xk)) < tol
        println("$xk é zero de f")
    end
    return xk
end

#Método de Newton Secante
function newton_secante(xk, xk1, f, tol = 1.0e-5)
    itermax = 100;
    iter= 0;
    while !(abs(f(xk))< tol) && iter < itermax
        mk = (f(xk)-f(xk1))/(xk-xk1)
        xk1=xk
        xk = xk - (f(xk)/mk)
        iter += 1 #iter = iter +1
        println(iter ,":", f(xk))
    end
    println("$iter  é o numero de interações")
    if abs(f(xk)) < tol
        println("$xk é zero de f")
    end
    return xk
end
