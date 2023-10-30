################################################################
# Author : Benoit Loger
# Modified : 05/04/2023
#
# This file contains models for
# our study about the Robust Lot Sizing Problem
#
################################################################
"""
    Master problem structure
"""
mutable struct LSPARMA
    inst::Instance
    Gamma::Float64
    beta::Vector{Float64}
    theta::Vector{Float64}
    Std::Float64
    UB::Float64
    LB::Float64
	sol::Vector{Float64}

    # Master problem and decision variables
    master::Model
    x::Vector{VariableRef}
    u::Vector{VariableRef}
    y::Vector{Vector{VariableRef}}
    z::VariableRef

    # Sub problem and decision variables
    sub::Model
    d::Vector{VariableRef}
    LSPARMA(inst::Instance,beta::Vector{Float64},theta::Vector{Float64},Std::Float64) = new(inst,2.0,beta,theta,Std)

end

function build_master_problem(LSP::LSPARMA)
    LSP.master = Model(SOLVER)
    D = [mean(LSP.inst.D) for t in 1:LSP.inst.T]
    LSP.u = @variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0)     # (1)
    LSP.x = @variable(LSP.master, [i in 1:LSP.inst.T], Bin)     # (1)
    LSP.y = Vector{VariableRef}[]
    push!(LSP.y,@variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0))
    LSP.z = @variable(LSP.master,lower_bound = 0.0)
    @objective(LSP.master, Min, LSP.z )    # (2)
	@constraint(LSP.master,con,LSP.z >= 0)
    @constraint(LSP.master,LSP.z >= sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[1][t]  for t in 1:LSP.inst.T))    # (3)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[1][t] >= LSP.inst.h * (sum(LSP.u[k] - D[k] for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[1][t] >= -LSP.inst.b * (sum(LSP.u[k] - D[k] for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.u[t] <= LSP.x[t]*sum(D))


    for att in OPTIMIZER_ATRIBUTES
        set_optimizer_attribute(LSP.master, att[1], att[2])
    end
    if SILENT_MASTER
        set_silent(LSP.master)
    end
    LSP.UB = Inf
end

function solve_master_problem(LSP::LSPARMA)
    optimize!(LSP.master)
    # Stop if master problem is infeasible
    if termination_status(LSP.master) == INFEASIBLE
        error("Master problem infeasible, cutting plane generation stoped !")
    end
    LSP.LB = objective_value(LSP.master)
	LSP.sol = value.(LSP.u)
end
#
function build_sub_problem(LSP::LSPARMA)
    # TODO : Changer après implémentation de l'apprentissage
    beta = LSP.inst.beta
    theta = LSP.inst.theta
    beta0 = 50.0
    Epsilon = zeros(Float64,length(LSP.inst.D))
    Std = LSP.Std
    P, Q = length(beta), length(theta) # Les instances sont telle que p == q (quitte à ajouter des 0.0)

    LSP.sub = Model(SOLVER)

    @variable(LSP.sub,I[t in 1:LSP.inst.T] >= 0) # Holding cost at each period
    @variable(LSP.sub,B[t in 1:LSP.inst.T] >= 0) # Backlogging cost at each period

    @variable(LSP.sub,p[t in 1:LSP.inst.T], Bin)
    @variable(LSP.sub,r[t in 1:LSP.inst.T], Bin)
    LSP.d = @variable(LSP.sub,[t in 1:LSP.inst.T],lower_bound=0.0)

    @variable(LSP.sub, -3.0*Std <= epsilon[t in 1:LSP.inst.T] <= 3.0*Std)

    @objective(LSP.sub,Max,sum(LSP.inst.f*value.(LSP.x[t]) + LSP.inst.c*value.(LSP.u[t]) + I[t] + B[t] for t in 1:LSP.inst.T))

    @constraint(LSP.sub,[k in 1:LSP.inst.T], I[k] >= LSP.inst.h * (sum(value.(LSP.u[i]) - LSP.d[i] for i in 1:k)))
    @constraint(LSP.sub,[k in 1:LSP.inst.T], I[k] <= LSP.inst.h * (sum(value.(LSP.u[i]) - LSP.d[i] for i in 1:k)) + 100000*(1-p[k]))
    @constraint(LSP.sub,[k in 1:LSP.inst.T], I[k] <= 100000*(p[k]))

    @constraint(LSP.sub,[k in 1:LSP.inst.T], B[k] >= -LSP.inst.b * (sum(value.(LSP.u[i]) - LSP.d[i] for i in 1:k)))
    @constraint(LSP.sub,[k in 1:LSP.inst.T], B[k] <= -LSP.inst.b * (sum(value.(LSP.u[i]) - LSP.d[i] for i in 1:k)) + 100000*(p[k]))
    @constraint(LSP.sub,[k in 1:LSP.inst.T], B[k] <= 100000*(1-p[k]))


    # On distingue 3 cas :
    # t = 1
        # Dépend uniquement des valeurs historiques et de epsilon[1]
        @constraint(LSP.sub, LSP.d[1] == beta0 + sum(beta[k]*LSP.inst.D[end-(k-1)] for k in 1:P) + epsilon[1])

    # t > 1 && <= max(p,q)
        # Dépend des valeus historiques, des valeurs de d précédentes et de epsilon[t-q:t]
        if P > 1
            @constraint(LSP.sub,[t in 2:P], LSP.d[t] == beta0 + sum(beta[k]*d[t-k] for k in 1:t-1) + sum(beta[k]*LSP.inst.D[end-(k-t)] for k in t:P) +
            sum(theta[k]*epsilon[t-k] for k in 1:t-1) + sum(theta[k]*Epsilon[end-(k-t)] for k in t:P) + epsilon[t])
        end
    # t > max(p,q)
        # Dépend uniquement des valeurs de d précédentes et des epsilon[t-q:t]
        @constraint(LSP.sub,[t in P+1:LSP.inst.T], LSP.d[t] == beta0 + sum(beta[k]*LSP.d[t-k] for k in 1:P) + sum(theta[k]*epsilon[t-k] for k in 1:Q) + epsilon[t])

    # Budget of uncertainty
    @constraint(LSP.sub,[t in 1:LSP.inst.T],sum(epsilon[k] for k in 1:t) <= 2*Std*sqrt(t))
    @constraint(LSP.sub,[t in 1:LSP.inst.T],sum(epsilon[k] for k in 1:t) >=  -2*Std*sqrt(t))

    for att in OPTIMIZER_ATRIBUTES
        set_optimizer_attribute(LSP.sub, att[1], att[2])
    end
    if SILENT_SUB
        set_silent(LSP.sub)
    end

end

function update_sub_problem(LSP::LSPARMA)
    build_sub_problem(LSP)
end

function solve_sub_problem(LSP::LSPARMA)
    optimize!(LSP.sub)
    if termination_status(LSP.sub) == INFEASIBLE
        error("Sub problem infeasible, cutting plane generation stoped at iteration $iteration !")
    end
    LSP.UB = min(LSP.UB, objective_value(LSP.sub))
end

function is_optimal(LSP::LSPARMA)
    return (LSP.UB - LSP.LB)/LSP.LB <= 0.0001
end

function update_master_problem(LSP::LSPARMA)
    iteration = length(LSP.y)

    push!(LSP.y,@variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0))
    @constraint(LSP.master,LSP.z >= sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[iteration+1][t]  for t in 1:LSP.inst.T))    # (3)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= LSP.inst.h * (sum(LSP.u[k] - value.(LSP.d[k]) for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= -LSP.inst.b * (sum(LSP.u[k] - value.(LSP.d[k]) for k in 1:t)))
    x = [LSP.sol[t] > 0 ? 1 : 0 for t in 1:LSP.inst.T]
	set_start_value.(LSP.x, x)
	u = LSP.sol
	set_start_value.(LSP.u, u)
end

function push_cut(LSP::LSPARMA,d::Vector{Float64})
    iteration = length(LSP.y)
    push!(LSP.y,@variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0))
    @constraint(LSP.master,LSP.z >= sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[iteration+1][t]  for t in 1:LSP.inst.T))    # (3)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= LSP.inst.h * (sum(LSP.u[k] - d[k] for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= -LSP.inst.b * (sum(LSP.u[k] - d[k] for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.u[t] <= LSP.x[t]*sum(d))
end
