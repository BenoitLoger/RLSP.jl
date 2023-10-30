################################################################
# Author : Benoit Loger
# Modified : 30/10/2023
#
# This file implement the cutting plane algorithm used to solve the RULSP
#
################################################################
"""

"""
function cutting_plane(problem;verbose::Bool=true,log_file::String="",cut_file_d::String="",cut_file_lt::String="",solution_file::String="",time_limit=600)
    # Initializing optimization models
    build_master_problem(problem)
    iteration = 0
    optimal = false
    total_time = 0

    # Open log and cut file if necessary
    if log_file != ""
        log = open(log_file,"w")
        	println(log,"Iteration,Time,LB,UB,Gap")
		close(log)
    end
    if cut_file_d != ""
        cut_d = open(cut_file_d,"w")
        close(cut_d)
    end
    if cut_file_lt != ""
        cut_lt = open(cut_file_lt,"w")
        close(cut_lt)
    end
    solution = []
    # Repeat until an Optimal Robust solution is found
    while !optimal && total_time <= time_limit
        iteration_time = (@timed begin
            # Solving master problem
            solve_master_problem(problem)
            solution = value.(problem.u)

            # update_sub_problem and solve it
            update_sub_problem(problem)
            solve_sub_problem(problem)

            # Check optimality of current solution and update master if necessary
            optimal = is_optimal(problem)
            if !optimal && total_time <= time_limit # Do not modify if we reach the time limit
                update_master_problem(problem)
            end
			if cut_file_d != ""
				cut_d = open(cut_file_d,"a")
				println(cut_d,value.(problem.d))
				close(cut_d)
			end
			if cut_file_lt != ""
				cut_lt = open(cut_file_lt,"a")
				println(cut_lt,round.(value.(problem.lt)))
				close(cut_lt)
			end

            iteration += 1
        end)[2]
        total_time += iteration_time
        if verbose
            println("| Iteration $iteration | $(round(total_time,digits=2)) | $(round(problem.LB,digits=2)) | $(round(problem.UB,digits=2)) | $(min(100.0,round((problem.UB-problem.LB)/problem.LB*100,digits=2)))")
        end
        if log_file != ""
			log = open(log_file,"a")
            	println(log,"$iteration,$(round(total_time,digits=2)),$(round(problem.LB,digits=2)),$(round(problem.UB,digits=2)),$(min(100.0,round((problem.UB-problem.LB)/problem.LB*100,digits=2)))")
			close(log)
        end

    end
    return total_time, solution
end
