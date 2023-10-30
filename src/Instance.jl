"""
    Struct of Lot Sizing Problem with uncertain demand
"""
mutable struct Instance
    T::Int                  # Planning horizon
    f::Float64              # Fixed ordering cost
    c::Float64              # Unit ordering cost
    h::Float64              # Unit holding cost per period
    b::Float64              # Unit backlogging cost per period
    D::Vector{Float64}      # Historical demand vector
    beta::Vector{Float64}   #
    theta::Vector{Float64}  #

    write::Function

    function Instance(T,f,c,h,b,D,beta,theta)
        this = new(T,f,c,h,b,D,beta,theta)
        this.write = function(name::String)
            if !isdir(name)
                mkdir(name)
            end
            f = open(string(name,"/costs.txt"),"w")
                print(f,string(T," "))
                print(f,string(this.f," "))
                print(f,string(c," "))
                print(f,string(h," "))
                println(f,string(b," ")) # Next line after printing b
                println(f,beta)
                println(f,theta)
            close(f)
            df = DataFrame(d = D)
            CSV.write(string(name,"/demands.csv"),df)
        end
        return this
    end

    function Instance(name::String;path="Instances/",evaluate=false)
		D = Vector{Float64}(DataFrame(CSV.File(string(path,name,"/demands.csv"))).d)
		f = open(string(path,name,"/costs.txt"),"r")
			c = parse.(Float64,Base.split(readline(f)))
            beta = parse.(Float64,Base.split(readline(f)[2:end-1]))
            theta = parse.(Float64,Base.split(readline(f)[2:end-1]))
		close(f)
        return Instance(c[1],c[2],c[3],c[4],c[5],D,beta,theta)
	end
end

function evaluate_sol(u::Vector{Float64},inst::Instance,D_test::Matrix{Float64},X,B,H,C)
    N = length(D_test[1,:])
    # X = Matrix{Float64}(undef,inst.T,N)
    # B = zeros(Float64,N)
    # H = zeros(Float64,N)
    # C = zeros(Float64,N)

    for i in 1:N
        C[i] = B[i] = H[i] = 0
        for t in 1:inst.T
            X[t,i] = (sum(u[k] - D_test[k,i] for k in 1:t))
            C[i] += (u[t] > 0)*inst.f + inst.c*u[t] + max(inst.h*X[t,i],-inst.b*X[t,i])
            if X[t,i] < 0
                B[i] += -X[t,i]
            elseif X[t,i] > 0
                H[i] += X[t,i]
            end
        end
        B[i] = B[i]/sum(D_test[:,i])*100
    end
    return X, C, H, B
end


function evaluate_adjustable_sol(u,x,inst::Instance,D_test::Matrix{Float64},X,B,H,C)
    N = length(D_test[1,:])
    # X = Matrix{Float64}(undef,inst.T,N)
    # B = zeros(Float64,N)
    # H = zeros(Float64,N)
    # C = zeros(Float64,N)
    for i in 1:N
        C[i] = B[i] = H[i] = 0
        for t in 1:inst.T
            X[t,i] = (sum(value.(u[k+1,1]) + sum(value.(u[k+1,j+1])*D_test[j,i] for j in 1:k) - D_test[k,i] for k in 1:t))
            C[i] += x[t]*inst.f + inst.c*(value.(u[t+1,1]) + sum(value.(u[t+1,j+1])*D_test[j,i] for j in 1:t)) + max(inst.h*X[t,i],-inst.b*X[t,i])
            if X[t,i] < 0
                B[i] += -X[t,i]
            elseif X[t,i] > 0
                H[i] += X[t,i]
            end
        end
        B[i] = B[i]/sum(D_test[:,i])*100
    end
    return X, C, H, B
end
