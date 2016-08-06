function getparameters(filename::AbstractString)
    f = open(filename,"r")
    dict = Dict{AbstractString, AbstractString}()
    for ln in eachline(f)
        m = match(r"^\s*(?P<key>\w+)\s+(?P<value>[\w+-.]+)", ln)
        if m != nothing
            dict[m[:key]] = m[:value]
        end
    end
    close(f)
    return dict
end

function print_parameters(parameters)
    println("Parameters")
    for (key, value) in parameters
        println("$key ==> $value")
    end
end


function density(array, x, y)
    return sum(array[x,y,:])
end

function density(array)
    return sum(array,3)
end


function velocity(array, x, y, rho)
    return sum(array[x,y,:] * c)/rho
end

function equilibrium(rho, u, c, w)
    cu = cell(size(c,2))
    for i = 1:size(c,2)
        cu[i] = 3 * dot(c[:,i],u)
    end
    usqr = 3/2 * dot(u,u)
    feq = cell(size(c,2))
    for i = 1:size(c,2)
        feq[i] = w[i]*rho(1+cu[i]+0.5*cu[i]^2 - usqr)
    end
end

function stream!(src, dsc, c)
    for i = 1:size(c,2)
        dsc[:,:,i] = circshift(src[:,:,i], -c[:,i])
    end
end

function collision!(array, c, omega, w)
    rho = density(array)
    for x = 1:size(array,1)
        for y = 1:size(array,2)
            u = velocity(array, x, y, rho[x,y])
            feq = equilibrium(rho[x,y], u, c, w)
            array[x,y,:] = array[x,y,:] - omega * (array[x,y,:] - feq[:])
        end
    end
end

c = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 1 1 1 0 0]
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
q = size(c,2)

parameters = getparameters("test")
print_parameters(parameters)
omega = parse(Float64, parameters["omega"])

println("c: $c")

x = parse(Int, parameters["xlength"])
y = parse(Int, parameters["ylength"])

src = Array(Float64, x, y, size(c,2))
dsc = Array(Float64, x, y, size(c,2))


while true
    stream!(src, dsc, c)
    collision!(dsc, c, omega, w)
    tmp = src
    src = dsc
    dsc = tmp
end
