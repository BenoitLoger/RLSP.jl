"""
    Defining the solver used
"""
SOLVER = CPLEX.Optimizer
function set_CPLEX_solver()
    try
        @eval using CPLEX
        global SOLVER = CPLEX.Optimizer
    catch
        @warn("Failed to import CPLEX using $SOLVER instead")
    end
end

function set_GUROBI_solver()
    try
        @eval using Gurobi
        global SOLVER = Gurobi.Optimizer
    catch
        @warn("Failed to import Gurobi using $SOLVER instead")
    end
end
"""
    Add a solver attribute
"""
function set_solver_attribute(name::String,value)
    push!(OPTIMIZER_ATRIBUTES,(name,value))
end
OPTIMIZER_ATRIBUTES = []
"""
    verbose mode of solvers for master and sub problems
"""
SILENT_MASTER = true
SILENT_SUB = true

PARAM = [0.98,1.56,2.11,2.63,3.17,3.68,4.23,4.69,5.21,5.68,6.17,6.63,7.12,7.59,7.99,8.46,8.9,9.41,9.93,10.39]
