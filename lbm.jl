using WriteVTK, Images


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
    println("Parameters:")
    for (key, value) in parameters
        println("\t$key ==> $value")
    end
end


function density(array, x, y)
    return sum(array[:,x,y])
end

function density(array)
    return squeeze(sum(array,1),1)
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
    for x = 2:x_length+1, y = 2:y_length+1
    #for i = 1:q
        #dsc[i,:,:] = circshift(src[i,:,:], -c[:,i])
        dsc[2,x,y] = src[2,x-1,y]
        dsc[3,x,y] = src[3,x,y-1]
        dsc[4,x,y] = src[4,x+1,y]
        dsc[5,x,y] = src[5,x,y+1]
        dsc[6,x,y] = src[6,x-1,y-1]
        dsc[7,x,y] = src[7,x+1,y-1]
        dsc[8,x,y] = src[8,x+1,y+1]
        dsc[9,x,y] = src[9,x-1,y+1]
    end
end

function collision!(array, c, omega, w)
    #rho = density(array)
    for x = 2:x_length+1, y  = 2:y_length+1
        rho = sum(array[:,x,y])
        u = velocity(array, x, y, rho)
        feq = equilibrium(rho, u)
        array[:,x,y] = array[:,x,y] - omega * (array[:,x,y] - feq[:])
    end
end


function boundary!(array) #so far only periodic
    array[:,1,:] = array[:,x_length+1,:]
    array[:,x_length+2,:] = array[:,2,:]

    array[:,:,1] = array[:,:,y_length+1]
    array[:,:,y_length+2] = array[:,:,2]
end

function findfirstcolumn(A, v)
    index = findfirst(A[1,:],v[1])
    found = false
    while index != 0 && found == false
        found = true
        for i = 2:size(v)[1]
            if A[i,index] != v[i]
                found = false
                break
            end
        end
        if found == true
            return index
        end
        index = findnext(A[1,:], v[1], index+1)
    end
    return 0
end

function obstacles!(f, obstacles)
    for x in 1:size(obstacles,1), y in 1:size(obstacles,2) 
        if obstacles[x,y] == true
            for i = 2:q
                f[i,x+1,y+1] = f[noslip[i], x+c[1,i]+1, y+c[2,i]+1]
            end
        end
    end
end

function visualization(array, obstacles, iteration)
    vtk_file = vtk_grid("sim_$iteration", size(obstacles,1), size(obstacles,2), size(obstacles,3))
    rho = density(array[:,2:size(obstacles,1)+1, 2:size(obstacles,2)+1])
    u = Array(Float64, 2, size(obstacles,1), size(obstacles,2))
    for x in 1:size(obstacles,1), y in 1:size(obstacles,2)
        u[1,x,y], u[2,x,y] = velocity(array, x+1, y+1, rho[x,y])
    end
    vtk_point_data(vtk_file, rho, "Pressure")
    vtk_point_data(vtk_file, u, "Velocity")
    vtk_save(vtk_file)
end


c = [0 1 0 -1 0 1 -1 -1 1; 0 0 1 0 -1 1 1 -1 -1]
w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
q = size(c,2)

noslip = Array(Int, q)
for i=1:q
    noslip[i] = findfirstcolumn(c, -c[:,i])
end
@show noslip

parameters = getparameters(ARGS[1])
#img = raw(load(ARGS[2]))
img = load(ARGS[2])
obstacles = convert(Array, img)'
obstacles = !map(Bool, obstacles)
print_parameters(parameters)
omega = parse(Float64, parameters["omega"])

println("c: $c")
println("q: $q")

x_length = size(obstacles, 1)
y_length = size(obstacles, 2)

println("x_length = $x_length")
println("y_length = $y_length")

src = Array(Float64, q, x_length+2, y_length+2)
dsc = similar(src)

# Init
feq = equilibrium(1.0,[0.1,0.0])
for x = 1:x_length+2
    for y = 1:y_length+2
        src[:,x,y] = feq[:]
        dsc[:,x,y] = feq[:]
    end
end
visualization(dsc, obstacles, 0)

# Time loop
for i = 1:30
    @show i
    obstacles!(src,obstacles)
    stream!(src, dsc)
    collision!(dsc, c, omega, w)
    #boundary!(src)
    visualization(dsc, obstacles, i)
    dsc, src = src, dsc
end
