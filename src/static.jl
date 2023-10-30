################################################################
# Author : Benoit Loger
# Modified : 30/10/2023
#
# Main file for static models
#
################################################################
using JuMP, CPLEX, Distributions, CSV, DataFrames, LinearAlgebra, Statistics
include("Instance.jl")
include("environment.jl")
include("algorithms.jl")
include("Models/Static/Bertsimas.jl")
include("Models/Static/ARMA.jl")


INSTANCE_PATH = "../Instances/"
SOLUTION_PATH = "../Solutions/"
CUT_PATH = "../Cuts/"

"""
    test()

"""
function test()
    i = 1
    # Generate a test_instance
    generate_instance()

    # Set solver attributes
    set_CPLEX_solver()
    set_solver_attribute("CPXPARAM_Threads",2)

    # Load the instance
    inst = Instance("$(i)",path=INSTANCE_PATH)

    # # Generate and solve budget based model
    # LSP = Bertsimas(inst,PARAM)
    # cutting_plane(LSP,verbose=true)
    # solution = LSP.sol
    #
    # # Save the solution
    # if !isdir("$(SOLUTION_PATH)$(i)")
    #     mkdir("$(SOLUTION_PATH)$(i)")
    # end
    # f = open("$(SOLUTION_PATH)/Bertsimas_$(i).txt","w")
    #     println(f,solution)
    # close(f)

    # sol = get_classic_sol(i)

    # Generate and solve budget based model
    println(inst.D)

    LSP = LSPARMA(inst,inst.beta,inst.theta,10.0)
    cutting_plane(LSP,verbose=true)
    solution = LSP.sol

    # Save the solution
    if !isdir("$(SOLUTION_PATH)$(i)")
        mkdir("$(SOLUTION_PATH)$(i)")
    end
    f = open("$(SOLUTION_PATH)/ARMA_$(i).txt","w")
        println(f,solution)
    close(f)
end

"""
    get_classic_sol(i::Int)

    Return the solution obtain with the Budget based uncertainty set if it exist
"""
function get_classic_sol(i::Int;dir=SOLUTION_PATH)
    AR1_file = open("$(dir)Bertsimas_$(i).txt","r")
    while true
        try
            line_array = parse.(Float64,split(readline(AR1_file)[2:end-1],","))
            return line_array[1:end]
        catch
            break
        end
    end
end

"""
    get_arma_sol(i::Int)

    Return the solution obtain with the Budget based uncertainty set if it exist
"""
function get_arma_sol(i::Int;dir=SOLUTION_PATH)
    AR1_file = open("$(dir)ARMA_$(i).txt","r")
    while true
        try
            line_array = parse.(Float64,split(readline(AR1_file)[2:end-1],","))
            return line_array[1:end]
        catch
            break
        end
    end
end


#========================================================#
#                TO DELETE LATER                         #
#========================================================#
"""
    generate_instance()

    Generate a test instance (used to test functions)
"""
function generate_instance()
    i = 1
    d = demand_ar1(100.,0.5,10.0,N=100)
    inst = Instance(14,0.0,0.1,0.3,8*0.3,d,[0.5],[0.0])
    inst.write("$(INSTANCE_PATH)$(i)")
end


"""
    demand_ar1(mu::Float64, beta1::Float64, Std::Float64; N::Int=1000, d0::Float64=100.)

Generate a vector of N demand following an AR(1) process d_t = beta0 + beta1 d_{t-1} + epsilon_t.
Keyword argument d0 is the starting point of the vector (i.e d0 -> [d_1,d_2,d_3,...,d_N]).
"""
function demand_ar1(mu::Float64, beta1::Float64, Std::Float64; N::Int=1000, d0::Float64=-1.0,mode="")
    beta0 = mu - mu*beta1
    b = sqrt(12*Std^2)/2
    a = -sqrt(12*Std^2)/2
    if d0 != -1
        if mode != "uni"
            d = [max(0,beta0 + beta1*d0 + rand(Normal(0,Std)))]
        else

            d = [max(0,beta0 + beta1*d0 + rand(Uniform(a,b)))]
        end
    else
        if mode != "uni"
            d = [max(0,beta0 + beta1*mu + rand(Normal(0,Std)))]
        else
            d = [max(0,beta0 + beta1*mu + rand(Uniform(a,b)))]
        end
    end
    for t in 1:N-1
        if mode != "uni"
            push!(d,max(0,beta0 + beta1*d[t] + rand(Normal(0,Std))))
        else
            push!(d,max(0,beta0 + beta1*d[t] + rand(Uniform(a,b))))
        end
    end
    return d
end
