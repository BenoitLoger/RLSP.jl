################################################################
# Author : Benoit Loger
# Modified : 05/04/2023
#
# This file contains models for
# our study about the Robust Lot Sizing Problem
#
################################################################
"""
    Robust Lot Sizing Model with uncertain demands and budget based uncertainty set
"""
mutable struct Bertsimas
    inst::Instance
    Gamma::Vector{Float64}
    epsilon::Float64 # Optimality tolerance
    UB::Float64
    LB::Float64
	M::Float64
    sol::Vector{Float64}
    tol::Float64 # Optimality tolerance

    # Master problem and decision variables
    master::Model
    x::Vector{VariableRef}
    u::Vector{VariableRef}
    y::Vector{Vector{VariableRef}}
    z::VariableRef

    # Sub problem and decision variables
    sub::Model
    d::Vector{VariableRef}
    Bertsimas(inst::Instance,Gamma::Vector{Float64}) = new(inst,Gamma,0.001)

end

function build_master_problem(LSP::Bertsimas)
    # Get the nominal case
    D = [mean(LSP.inst.D) for t in 1:LSP.inst.T]
    DD = [2.0*std(LSP.inst.D) for t in 1:LSP.inst.T]
	LSP.M = sum(D .+ DD)
    # Problem initialization
    LSP.master = Model(SOLVER)
    # Problem variables
    LSP.u = @variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0)       # u[t] : quantity ordered at period t
    LSP.x = @variable(LSP.master, [i in 1:LSP.inst.T], Bin)                     # x[t] : 1 if a order is placed in period t

    LSP.y = Vector{VariableRef}[]                                               # y[t] : holding/backlogging cost of period t
    LSP.z = @variable(LSP.master,lower_bound = 0.0)                             # z : worst cost among scenarios (epigraph formulation)
    # Problem constraints : general
    @objective(LSP.master, Min, LSP.z )                                                                                     # (1)
    @constraint(LSP.master,BigM[t in 1:LSP.inst.T], LSP.u[t] <= LSP.x[t]*LSP.M)                                    		    # (2)
    # Problem constraints : scenario related
    push!(LSP.y,@variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0))                                              # (3)
    @constraint(LSP.master,LSP.z >= sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[1][t]  for t in 1:LSP.inst.T))  # (4)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[1][t] >= LSP.inst.h * (sum(LSP.u[k] - D[k] for k in 1:t)))            # (5)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[1][t] >= -LSP.inst.b * (sum(LSP.u[k] - D[k] for k in 1:t)))           # (6)
    @constraint(LSP.master,LSP.x[1] == 1)

    # Constraint description
    # (1) : Minimizing the worst cost
    # (2) : Big M constraint linking u[t] and x[t]
    # (3) : Add y[t] variables for the nominal scenario
    # (4) : Epigraph formulation of Min Max model
    # (5) : Holding cost constraint
    # (6) : Backlogging cost constraint

    # Set Optimizer attributes
    for att in OPTIMIZER_ATRIBUTES
        set_optimizer_attribute(LSP.master, att[1], att[2])
    end
    if SILENT_MASTER
        set_silent(LSP.master)
    end
    LSP.UB = Inf
end

function solve_master_problem(LSP::Bertsimas)
    optimize!(LSP.master)
    # Stop if master problem is infeasible
    if termination_status(LSP.master) == INFEASIBLE
        error("Master problem infeasible, cutting plane generation stoped !")
    end
    LSP.LB = objective_value(LSP.master)
    LSP.sol = value.(LSP.u)
end
#
function build_sub_problem(LSP::Bertsimas)
    D = [mean(LSP.inst.D) for t in 1:LSP.inst.T]
    DD = [2.0*std(LSP.inst.D) for t in 1:LSP.inst.T]

    LSP.sub = Model(SOLVER)

    @variable(LSP.sub,I[t in 1:LSP.inst.T] >= 0) # Holding cost at each period
    @variable(LSP.sub,B[t in 1:LSP.inst.T] >= 0) # Backlogging cost at each period

    @variable(LSP.sub,p[t in 1:LSP.inst.T], Bin)
    LSP.d = @variable(LSP.sub,[t in 1:LSP.inst.T],lower_bound=0.0)
    @variable(LSP.sub,0 <= z_pos[t in 1:LSP.inst.T] <= 1)
    @variable(LSP.sub,0 <= z_neg[t in 1:LSP.inst.T] <= 1)

    @objective(LSP.sub,Max,sum(LSP.inst.f*value.(LSP.x[t]) + LSP.inst.c*value.(LSP.u[t]) + I[t] + B[t] for t in 1:LSP.inst.T))

    @constraint(LSP.sub,[t in 1:LSP.inst.T], I[t] >= LSP.inst.h * (sum(value.(LSP.u[k]) - LSP.d[k] for k in 1:t)))
    @constraint(LSP.sub,[t in 1:LSP.inst.T], I[t] <= LSP.inst.h * (sum(value.(LSP.u[k]) - LSP.d[k] for k in 1:t)) + 10000*(1-p[t]))
    @constraint(LSP.sub,[t in 1:LSP.inst.T], I[t] <= 10000*(p[t]))

    @constraint(LSP.sub,[t in 1:LSP.inst.T], B[t] >= -LSP.inst.b * (sum(value.(LSP.u[k]) - LSP.d[k] for k in 1:t)))
    @constraint(LSP.sub,[t in 1:LSP.inst.T], B[t] <= -LSP.inst.b * (sum(value.(LSP.u[k]) - LSP.d[k] for k in 1:t)) + 10000*(p[t]))
    @constraint(LSP.sub,[t in 1:LSP.inst.T], B[t] <= 10000*(1-p[t]))

    @constraint(LSP.sub,[t in 1:LSP.inst.T], LSP.d[t] == D[t]+DD[t]*(z_pos[t]-z_neg[t]))
    @constraint(LSP.sub,[t in 1:LSP.inst.T],sum(z_pos[k]+z_neg[k] for k in 1:t) <= LSP.Gamma[t]) # Budget for Demands

    for att in OPTIMIZER_ATRIBUTES
        set_optimizer_attribute(LSP.sub, att[1], att[2])
    end
    if SILENT_SUB
        set_silent(LSP.sub)
    end

end

function update_sub_problem(LSP::Bertsimas)
    build_sub_problem(LSP)
end

function solve_sub_problem(LSP::Bertsimas)
    optimize!(LSP.sub)
    if termination_status(LSP.sub) == INFEASIBLE
        error("Sub problem infeasible, cutting plane generation stoped at iteration !")
    end
    LSP.UB = min(LSP.UB, objective_value(LSP.sub))
end

function is_optimal(LSP::Bertsimas)
    return (LSP.UB - LSP.LB)/LSP.LB <= LSP.tol
end

function update_master_problem(LSP::Bertsimas)
    iteration = length(LSP.y)
    # Add the y variables corresponding to the curent iteration (Collumn generation)
    push!(LSP.y,@variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0))

    # (Constraint generation)
    @constraint(LSP.master,LSP.z >= sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[iteration+1][t]  for t in 1:LSP.inst.T))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= LSP.inst.h * (sum(LSP.u[k] - value.(LSP.d[k]) for k in 1:t)))
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.y[iteration+1][t] >= -LSP.inst.b * (sum(LSP.u[k] - value.(LSP.d[k]) for k in 1:t)))

    # Update the BigM constraints
	LSP.M = max(LSP.M,sum(value.(LSP.d)))
	delete.(LSP.master,LSP.master[:BigM])
	LSP.master[:BigM] = @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.u[t] <= LSP.x[t]*LSP.M)
end

###################################################################################################
