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
    return sum(array[:,x,y])
end

function density(array)
    return sum(array,1)
end


function velocity(array, x, y, rho)
    return (c * array[:,x,y])/rho
end

function equilibrium(rho, u)
    cu = Array(Float64, q)
    for i = 1:q
        cu[i] = 3 * dot(c[:,i],u)
    end
    usqr = 3/2 * dot(u,u)
    feq = Array(Float64, q)
    for i = 1:q
        feq[i] = w[i]*rho*(1+cu[i]+0.5*cu[i]^2 - usqr)
    end
    return feq
end

function stream!(src, dsc)
    for i = 1:q
        dsc[i,:,:] = circshift(src[i,:,:], -c[:,i])
    end
end

function collision!(array, c, omega, w)
    rho = density(array)
    rho = reshape(rho, x_length, y_length)
    for x = 1:x_length
        for y = 1:y_length
            u = velocity(array, x, y, rho[x,y])
            feq = equilibrium(rho[x,y], u)
            array[:,x,y] = array[:,x,y] - omega * (array[:,x,y] - feq[:])
        end
    end
end


function boundary!(array)

end


c = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 1 1 1 0 0]
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
q = size(c,2)

parameters = getparameters("test")
print_parameters(parameters)
omega = parse(Float64, parameters["omega"])

println("c: $c")
println("q: $q")

x_length = parse(Int, parameters["xlength"])
y_length = parse(Int, parameters["ylength"])

src = Array(Float64, q, x_length, y_length)
dsc = Array(Float64, q, x_length, y_length)

# Init
feq = equilibrium(1.0,[0.0,0.0])
for x = 1:x_length
    for y = 1:y_length
        src[:,x,y] = feq[:]
    end
end


for i = 1:100
    println("Iteration: $i")
    stream!(src, dsc)
    collision!(dsc, c, omega, w)
    dsc, src = src, dsc
end
