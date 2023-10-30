using StatsPlots, CSV, DataFrames
include("Instance.jl")
include("../data_factory.jl")

####################################################################
#                      TESTING FUNCTIONS                           #
####################################################################

function plot_guillaume(dir::String,set::Vector{Int};N::Int=10000,file="../../Plots/CVAR.pdf")
	beta = 0.3
	conf = DataFrame(CSV.File("../../Instances/$(dir)/conf.csv"))
	inst = Instance("../../Instances/$(dir)/1",path="")

	# Initialize Vectors and Matrices
	D = Matrix{Float64}(undef,inst.T,N)
	X = Matrix{Float64}(undef,inst.T,N)
	C = Matrix{Float64}(undef,inst.T,N)
	H = Vector{Float64}(undef,N)
	B = Vector{Float64}(undef,N)

	df = DataFrame(alpha=String[],method=String[],CVar=Float64[])
	plt = plot(ylabel="Relative superquantile \$\\rho_.(\\alpha) \$",xlabel="Value of \$\\alpha\$")
	for alpha in [0.0,0.25,0.5,0.75,0.90]
		for i in set
			inst = Instance("../../Instances/$(dir)/$i",path="")
			for j in 1:N
				D[:,j] = demand_ar1(conf.mu[i],conf.beta[i],conf.Std[i],N=inst.T,d0=conf.d_0[i])
			end
			println(inst.f)
			# Plot solution SAA
			C_SAA = Float64[]
			for j in 1:10
				SAA_sol = []
				try
					SAA_sol = get_SAA_sol(i,dir=dir,N=100,j=j)
				catch
					SAA_sol = get_SAA_sol(i,dir=dir,N=200,j=j)
				end
				X, C, H, B = evaluate_sol_seperate(SAA_sol,inst,D, X, C, H, B)
				C_SAA = vcat(C_SAA,[sum(C[:,i]) for i in 1:N])
			end
			ref_saa = tvar(C_SAA,alpha)
			# push!(df,[string(alpha),"SAA",tvar(C_SAA,alpha)])

			# Plot solution Classic
			Classic_sol = get_Classic_sol(i,2.0,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(Classic_sol,inst,D, X, C, H, B)
			Costs_rob = [sum(C[:,i]) for i in 1:N]
			# push!(df,[string(alpha),"\$\\mathcal{U}^{\\Gamma}_d\$",tvar(Costs,alpha)])
			ref_rob = tvar(Costs_rob,alpha)

			# Plot solution AR1
			AR1_sol = get_AR1_sol(i,2.0,2.0,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
			Costs = [sum(C[:,i]) for i in 1:N]
			ref_ar1 = tvar(Costs,alpha)
			push!(df,[string(alpha),"\$\\rho_{SAA}(\\alpha)\$",(ref_saa-ref_ar1)/ref_ar1])
			push!(df,[string(alpha),"\$\\rho_{\\Gamma}(\\alpha)\$",(ref_rob-ref_ar1)/ref_ar1])
		end

	end
	plt = groupedboxplot!(df.alpha,df.CVar,groups=df.method,dpi=300,legend=:outertopright)
	savefig(plt,file)
end

function plot_q(dir::String,set::Vector{Int};N::Int=10000)
	beta = 0.3
	conf = DataFrame(CSV.File("../../Instances/$(dir)/conf.csv"))
	inst = Instance("../../Instances/$(dir)/1",path="")

	# Initialize Vectors and Matrices
	D = Matrix{Float64}(undef,inst.T,N)
	X = Matrix{Float64}(undef,inst.T,N)
	C = Matrix{Float64}(undef,inst.T,N)
	H = Vector{Float64}(undef,N)
	B = Vector{Float64}(undef,N)

	df = DataFrame(alpha=String[],method=String[],CVar=Float64[])
	plt = plot(ylabel="Relative \$CVar_{\\alpha}\$",xlabel="Value of \$\\alpha\$")
	for alpha in [0.5,0.75,0.90,0.95]
		for i in set
			inst = Instance("../../Instances/$(dir)/$i",path="")
			for j in 1:N
				D[:,j] = demand_ar1(conf.mu[i],conf.beta[i],conf.Std[i],N=inst.T,d0=conf.d_0[i])
			end
			println(inst.f)
			# Plot solution SAA
			C_SAA = Float64[]
			for j in 1:1
				SAA_sol = get_SAA_sol(i,dir=dir,N=200,j=j)
				X, C, H, B = evaluate_sol_seperate(SAA_sol,inst,D, X, C, H, B)
				C_SAA = vcat(C_SAA,[sum(C[:,i]) for i in 1:N])
			end
			ref_saa = quantile(C_SAA,alpha)
			# push!(df,[string(alpha),"SAA",tvar(C_SAA,alpha)])

			# Plot solution Classic
			Classic_sol = get_Classic_sol(i,2.0,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(Classic_sol,inst,D, X, C, H, B)
			Costs_rob = [sum(C[:,i]) for i in 1:N]
			# push!(df,[string(alpha),"\$\\mathcal{U}^{\\Gamma}_d\$",tvar(Costs,alpha)])
			ref_rob = quantile(Costs_rob,alpha)

			# Plot solution AR1
			AR1_sol = get_AR1_sol_v4(i,2.0,2.0,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
			Costs = [sum(C[:,i]) for i in 1:N]
			ref_ar1 = quantile(Costs,alpha)
			push!(df,[string(alpha),"SAA",ref_saa])
			push!(df,[string(alpha),"Robust",ref_rob])
			push!(df,[string(alpha),"AR1",ref_ar1])

		end

	end
	plt = groupedboxplot!(df.alpha,df.CVar,groups=df.method,dpi=300,legend=:outertopright,ylim=(-1:1))
	# savefig(plt,"CVAR.pdf")
end


function evaluate(dir,set;N=10000)
	conf = DataFrame(CSV.File("../../Instances/$(dir)/conf.csv"))
	i = set[1]
    # Initialize Vectors and Matrices
    D = Matrix{Float64}(undef,conf.T[i],N)
    X = Matrix{Float64}(undef,conf.T[i],N)
    C = Matrix{Float64}(undef,conf.T[i],N)
    H = Vector{Float64}(undef,N)
    B = Vector{Float64}(undef,N)


	df = DataFrame(i=Int[],method=String[],gamma=Float64[],Avg=Float64[],median=Float64[],Std=Float64[],Q95=Float64[],Q99=Float64[],beta=Float64[],sigma=Float64[],s=Float64[],N=Float64[],d0=Float64[],value=Float64[])
	for i in set
		inst = Instance("../../Instances/$(dir)/$i",path="")
		for j in 1:N
	        D[:,j] = demand_ar1(conf.mu[i],conf.beta[i],conf.Std[i],N=inst.T,d0=conf.d_0[i])
	    end
		# SAA
		C_SAA = Float64[]
		for j in 1:10
		SAA_sol = get_SAA_sol(i,dir=dir,N=100,j=j)
		X, C, H, B = evaluate_sol_seperate(SAA_sol,inst,D, X, C, H, B)
		C_SAA = vcat(C_SAA,[sum(C[:,i]) for i in 1:N])
		end
		push!(df,[i,"SAA",0.0,mean(C_SAA),median(C_SAA),std(C_SAA),quantile(C_SAA,0.95),quantile(C_SAA,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],0.0])

		# CLASSIC ROBUST
		# for Gamma in 0.0:inst.T
			Classic_sol = get_Classic_sol(i,2.0,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(Classic_sol,inst,D, X, C, H, B)
			Costs = [sum(C[:,i]) for i in 1:N]
			push!(df,[i,"Rob",2.0,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],Classic_sol[end-1]])
		# end

		# AR(1) Robust
		for Gamma in collect(2.0)
			println(i," ", Gamma)
			AR1_sol = get_AR1_sol(i,2.0,Gamma,2.0,dir=dir)
			X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
			Costs = [sum(C[:,i]) for i in 1:N]
			push!(df,[i,"AR1",Gamma,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],AR1_sol[end-1]])
		end

		# for Gamma in collect(2.0)
		# 	println(i," ", Gamma)
		# 	AR1_sol = get_AR1_sol_v3(i,2.0,Gamma,2.0,dir=dir)
		# 	X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
		# 	Costs = [sum(C[:,i]) for i in 1:N]
		# 	push!(df,[i,"AR1_diff",Gamma,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],AR1_sol[end-1]])
		# end
		#
		# for Gamma in collect(2.0)
		# 	println(i," ", Gamma)
		# 	AR1_sol = get_AR1_sol_v4(i,2.0,Gamma,2.0,dir=dir)
		# 	X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
		# 	Costs = [sum(C[:,i]) for i in 1:N]
		# 	push!(df,[i,"AR1_free",Gamma,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],AR1_sol[end-1]])
		# end
		#
		# for Gamma in collect(2.0)
		# 	println(i," ", Gamma)
		# 	AR1_sol = get_AR1_sol_v5(i,2.0,Gamma,2.0,dir=dir)
		# 	X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
		# 	Costs = [sum(C[:,i]) for i in 1:N]
		# 	push!(df,[i,"AR1_Var",Gamma,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99),conf.beta[i],conf.Std[i],inst.f,length(inst.D),inst.D[end],AR1_sol[end-1]])
		# end
	end
	return df
end

# function evaluate(dir,i;N=10000)
# 	conf = DataFrame(CSV.File("../../Instances/$(dir)/conf.csv"))
# 	inst = Instance("../../Instances/$(dir)/$i",path="")
#
#     # Initialize Vectors and Matrices
# 	println(conf)
# 	println(conf.T[i])
#     D = Matrix{Float64}(undef,inst.T,N)
#     X = Matrix{Float64}(undef,inst.T,N)
#     C = Matrix{Float64}(undef,inst.T,N)
#     H = Vector{Float64}(undef,N)
#     B = Vector{Float64}(undef,N)
#
# 	inst = Instance("../../Instances/$(dir)/$i",path="")
# 	for j in 1:N
# 		D[:,j] = demand_ar1(conf.mu[i],conf.beta[i],conf.Std[i],N=inst.T,d0=conf.d_0[i])
# 	end
#
# 	df = DataFrame(i=Int[],method=String[],gamma=Float64[],Avg=Float64[],median=Float64[],Std=Float64[],Q95=Float64[],Q99=Float64[])
# 	f = open("solutions.txt","r")
# 	j=1.0
# 	while !eof(f)
# 		line_array = parse.(Float64,split(readline(f)[2:end-1],","))
# 		AR1_sol = line_array[1:inst.T]
# 		X, C, H, B = evaluate_sol_seperate(AR1_sol,inst,D, X, C, H, B)
# 		Costs = [sum(C[:,i]) for i in 1:N]
# 		push!(df,[i,"AR1",j,mean(Costs),median(Costs),std(Costs),quantile(Costs,0.95),quantile(Costs,0.99)])
# 		j += 1
# 	end
# 	return df
# end

function summary(set)
	df = DataFrame(CSV.File("../../results/Results.csv"))
	res = DataFrame(i=Int[],method=String[],gamma=Float64[],Avg=Float64[],median=Float64[],Std=Float64[],Q95=Float64[],Q99=Float64[],beta=Float64[],sigma=Float64[],s=Float64[],b=Float64[])
	for i in set
		# SAA
		sub = filter(row -> row.i == i && row.method == "SAA",df)
		push!(res,sub[1,1:12])
		println(sub)
		# Classic
		#If we need the best average solution
		# sub = filter(row -> row.i == i && row.method == "Rob",df)
		# index = findmin(sub.Avg)[2]
		# push!(res,sub[index,:])
		sub = filter(row -> row.i == i && row.method == "Rob",df)
		index = findmin(sub.Q99)[2]
		push!(res,sub[index,1:12])

		# AR1
		#If we need the best average solution
		# sub = filter(row -> row.i == i && row.method == "AR1",df)
		# index = findmin(sub.Avg)[2]
		# push!(res,sub[index,:])
		sub = filter(row -> row.i == i && row.method == "AR1",df)
		index = findmin(sub.Q99)[2]
		push!(res,sub[index,1:12])

	end
	return res
end

function table()
	sep = " & "
	line = "\\\\"
	por = "\\%"
	df = DataFrame(CSV.File("../../results/Results_2.csv"))
	sub = filter(row -> row.N == 100 && row.sigma == 10.0,df)
	betas = unique(sub.beta)
	println(betas)
		for s in unique(sub.s)
		select = filter(row -> row.s == s,sub)
			for beta in unique(select.beta)
			select2 = filter(row -> row.beta == beta,select)
			SAA = filter(row -> row.method == "SAA",select2)
			Rob = filter(row -> row.method == "Rob",select2)
			AR1 = filter(row -> row.method == "AR1",select2)

			saa_avg, saa_q95, saa_q99 = mean(SAA.Avg), mean(SAA.Q95), mean(SAA.Q99)
			rob_avg, rob_q95, rob_q99 = mean(Rob.Avg), mean(Rob.Q95), mean(Rob.Q99)
			ar1_avg, ar1_q95, ar1_q99 = mean(AR1.Avg), mean(AR1.Q95), mean(AR1.Q99)
			println(round(select2.s[1]),sep,round(select2.beta[1],digits=2),sep,round(saa_avg,digits=2),sep,
			round((rob_avg-saa_avg)/saa_avg*100,digits=2),por,sep,round((ar1_avg-saa_avg)/saa_avg*100,digits=2),por,sep,
			round(saa_q95,digits=2),sep,round((rob_q95-saa_q95)/saa_q95*100,digits=2),por,sep,round((ar1_q95-saa_q95)/saa_q95*100,digits=2),por,sep,
			round(saa_q99,digits=2),sep,round((rob_q99-saa_q99)/saa_q99*100,digits=2),por,sep,round((ar1_q99-saa_q99)/saa_q99*100,digits=2),por,line)
		end
	end
	return sub
end

####################################################################
#                      USEFULL FUNCTIONs                           #
####################################################################
"""
    get_AR1_sol(i::Int,zeta::Float64,gamma::Float64)

    Return the solution obtain with the AR1 based uncertainty set with a given value of parameter
"""
function get_AR1_sol(i::Int,zeta::Float64,gamma::Float64,alpha::Float64;dir="test")
    AR1_file = open("../../Solutions/$(dir)/AR1_$(i).txt","r")
    while true
        try
            line_array = parse.(Float64,split(readline(AR1_file)[2:end-1],","))
            if line_array[end] == gamma
                return line_array[1:end]
            end
        catch
            break
        end
    end
end


"""
    get_Classic_sol(i::Int,gamma::Float64)

    Return the solution obtain with the AR1 based uncertainty set with a given value of parameter
"""
function get_Classic_sol(i::Int,gamma::Float64,alpha::Float64;dir="test")
    AR1_file = open("../../Solutions/$(dir)/classic_$(i).txt","r")
    while true
        try
            line_array = parse.(Float64,split(readline(AR1_file)[2:end-1],","))
            if line_array[end] == gamma
                return line_array[1:end]
            end
        catch
            break
        end
    end
end

"""
    get_SAA_sol(i::Int,zeta::Float64,gamma::Float64)

    Return the solution obtain with the SAA based model
"""
function get_SAA_sol(i::Int;dir="test",solutions_path="../../Solutions",N::Int=500,j::Int=1)
    AR1_file = open("$(solutions_path)/$(dir)/SAA_$(i)_$(N)_$(j).txt","r")
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
    get_SAA_sol(i::Int,zeta::Float64,gamma::Float64)

    Return the solution obtain with the SAA based model
"""
function get_Normal_sol(i::Int;dir="test",solutions_path="../../Solutions",N::Int=500,j::Int=1)
    AR1_file = open("$(solutions_path)/$(dir)/Normal_$(i)_$(N)_$(j).txt","r")
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
    get_SAA_sol(i::Int,zeta::Float64,gamma::Float64)

    Return the solution obtain with the SAA based model
"""
function get_RAAA_sol(i::Int;dir="test")
    AR1_file = open("../../Solutions/$(dir)/RAAA_$(i).txt","r")
    while true
        try
            line_array = parse.(Float64,split(readline(AR1_file)[2:end-1],","))
            return line_array[1:end]
        catch
            break
        end
    end
end
