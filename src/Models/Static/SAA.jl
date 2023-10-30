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
mutable struct SAA
    inst::Instance
    D::Matrix{Float64}

    # Master problem and decision variables
    master::Model
    x::Vector{VariableRef}
    u::Vector{VariableRef}
    y::Matrix{VariableRef}

    # Sub problem and decision variables
    SAA(inst::Instance,D::Matrix{Float64}) = new(inst,D)
end

function build_master_problem(LSP::SAA)
    n,N = size(LSP.D)
    # Problem initialization
    LSP.master = Model(SOLVER)
    # Problem variables
    LSP.u = @variable(LSP.master, [i in 1:LSP.inst.T], lower_bound = 0.0)           # u[t] : quantity ordered at period t
    LSP.x = @variable(LSP.master, [i in 1:LSP.inst.T], Bin)                         # x[t] : 1 if a order is placed in period t
    LSP.y = @variable(LSP.master, [t in 1:LSP.inst.T,i in 1:N], lower_bound = 0.0)  # x[t] : 1 if a order is placed in period t

    # Problem constraints : general
    @objective(LSP.master, Min, sum(1/N * sum(LSP.inst.f*LSP.x[t] + LSP.inst.c * LSP.u[t] + LSP.y[t,i] for t in 1:LSP.inst.T) for i in 1:N))                                                                                     # (1)
    @constraint(LSP.master,[t in 1:LSP.inst.T], LSP.u[t] <= LSP.x[t]*10000)                                      # (2)
    # Problem constraints : scenario related
    @constraint(LSP.master,[t in 1:LSP.inst.T,i in 1:N], LSP.y[t,i] >= LSP.inst.h * (sum(LSP.u[k] - LSP.D[k,i] for k in 1:t)))            # (5)
    @constraint(LSP.master,[t in 1:LSP.inst.T,i in 1:N], LSP.y[t,i] >= -LSP.inst.b * (sum(LSP.u[k] - LSP.D[k,i] for k in 1:t)))           # (6)

    # Set Optimizer attributes
    for att in OPTIMIZER_ATRIBUTES
        set_optimizer_attribute(LSP.master, att[1], att[2])
    end
    println(SILENT_MASTER)
    if SILENT_MASTER
        set_silent(LSP.master)
    end
end

function solve_master_problem(LSP::SAA)
    optimize!(LSP.master)
end

###################################################################################################
