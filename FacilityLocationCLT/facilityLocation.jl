###
# Facility Location Example in Julia
###
# Develops a facility location Example
# Emphasis on usng data to construct sets and tune

using JuMP
using Distributions
using Gurobi

type FacilityData
    fixedCost::Float64
    capCost::Float64
    revenue::Float64
    shipCost::Float64
    maxCap::Float64
    facLocs::Array{Float64, 2}
    demandLocs::Array{Float64, 2}
end

#Random locations for both facilities and locations on 10 x 10 square
function FacilityData(numFacs::Int64, numDemands::Int64; 
                        square_size=10., fixedCost=120., capCost=2., revenue=10., shipCost=1., maxCap=1e3)
    facLocs = rand(Uniform(0, square_size), (2, numFacs))

    #place half of demands in little clusters at top and bottom
    demandLocs = rand(Uniform(-2, 2), (2, int(numDemands/4))) + 2
    demandLocs2 = rand(Uniform(-2, 2), (2, int(numDemands/4))) + 8
    demandLocs = hcat(demandLocs, demandLocs2)
    demandLocs = hcat(demandLocs, rand(Uniform(0, square_size), (2, numDemands - size(demandLocs, 2))))
    return FacilityData(fixedCost, capCost, revenue, shipCost, maxCap, facLocs, demandLocs)
end

numFacs(data::FacilityData) = size(data.facLocs, 2)
numDemands(data::FacilityData) = size(data.demandLocs, 2)
dist(f, d, data::FacilityData) = norm(data.facLocs[:, f] - data.demandLocs[:, d])

type FacilityModel
    m
    z
    x
    caps
    costComp
    shipComp
end
FacilityModel(m, z, x, caps, costComp) = FacilityModel(m, z, x, caps, costComp, AffExpr())

function createProblem(data::FacilityData)
    m = Model(solver=GurobiSolver(OutputFlag=false))
    @defVar(m, z[1:numFacs(data)], Bin)
    @defVar(m, 0<= x[1:numFacs(data), 1:numDemands(data)] <= 1)
    @defVar(m, 0 <= cap[1:numFacs(data)] <= data.maxCap)
    @defVar(m, costComp >= 0)
    for d = 1:numDemands(data)
        @addConstraint(m, sum{x[f, d], f = 1:numFacs(data)} <= 1)
    end
    for f = 1:numFacs(data)
        @addConstraint(m, cap[f] <= data.maxCap * z[f])
    end
    @addConstraint(m, costComp == data.fixedCost * sum{z[f], f=1:numFacs(data)} 
                                    + data.capCost * sum{cap[f], f=1:numFacs(data)} )
    return FacilityModel(m, z, x, cap, costComp)
end

#Solve facility location for a fixed value of demand
function addDemandConstr!(facModel::FacilityModel, data::FacilityData, demands::Array{Float64, 1})
    if size(demands, 1) != numDemands(data)
        error("Demand sequence wrong length")
    end

    m = facModel.m; x = facModel.x
    cNumDemands = numDemands(data); cNumFacs = numFacs(data)
    @defVar(m, shipComp)
    @addConstraint(m, shipComp >= 
        sum{(data.shipCost * dist(f, d, data) - data.revenue) * demands[d] * x[f, d], 
                f=1:cNumFacs, d=1:cNumDemands})

    facModel.shipComp = shipComp
    for f = 1:cNumFacs
        @addConstraint(m, sum{demands[d] * x[f, d], d=1:cNumDemands} <= facModel.caps[f])
    end
    @setObjective(m, Min, facModel.costComp + shipComp)
end

addDemandConstr!(facModel::FacilityModel, data::FacilityData, avgDemand::Float64) = 
    addDemandConstr!(facModel, data, avgDemand * ones(numDemands(data)))

function solve(facModel::FacilityModel, data::FacilityData)
    m = facModel.m; status = JuMP.solve(m)
    println("Objective value: ", getObjectiveValue(m))
    println("Total Capacity: ", sum([getValue(facModel.caps)[f] for f=1:numFacs(data)]))
    # for f = 1:numFacs(data)
    #     if getValue(facModel.z[f]) > .5
    #         println("Facility: ", f, "\t", getValue(facModel.caps[f]))
    #     end
    # end

    #extract all the xs and then manipulate... seems like a bug
    xvals = getValue(facModel.x)
    # for d = 1:numDemands(data)
    #     println("Demand Fraction: ", d, "\t", sum(xvals[:, d]))
    # end
    return getObjectiveValue(m), getValue(facModel.shipComp), getValue(facModel.z), getValue(facModel.caps)
end

#Used to solve the second stage problem.  Not most efficient way.
function fixFacilities(facModel::FacilityModel, zs, caps)
    for f = 1:size(facModel.z, 1)
        setLower( facModel.z[f], zs[f] ); setUpper( facModel.z[f], zs[f] )
        setLower( facModel.caps[f], caps[f] ); setUpper( facModel.caps[f], caps[f] )
    end
end

type CLTSet
    # C' C = Sigma
    C::Array{Float64, 2} 
    lbounds::Array{Float64, 1}
    ubounds::Array{Float64, 1}
    sigmas::Array{Float64, 1}
    means::Array{Float64, 1}
    Gamma::Float64
end

#Assume data comes by columns
#VG Ask Iain about right way to do this idiomaticallly
function CLTSet(data_set, support_size, Gamma, bSmart)
    numDemands = size(data_set, 1)
    means = Float64[mean(data_set[i, :]) for i=1:numDemands]

    if bSmart
        C = chol(cov(transpose(data_set)))
        means = C\means  # solves mean_zeta = C * mean_u
        sigmas = ones(numDemands)
        lbounds = Float64[means[i] - support_size for i=1:numDemands]
        ubounds = Float64[means[i] + support_size for i=1:numDemands]
   else
        C = eye(numDemands)
        sigmas = Float64[sqrt(var( data_set[i, :])) for i=1:numDemands]
        lbounds = Float64[means[i] - support_size * sigmas[i] for i=1:numDemands]
        ubounds = Float64[means[i] + support_size * sigmas[i] for i=1:numDemands]
    end
    return CLTSet(C, lbounds, ubounds, sigmas, means, Gamma)
end

function addRobConstr!(m, xbar, rhs, Uset::CLTSet)
    numDemands = size(Uset.lbounds, 1)
    @defVar(m, theta1 >= 0)
    @defVar(m, theta2 >= 0)
    @defVar(m, p[1:numDemands] >= 0.)
    @defVar(m, q[1:numDemands] >= 0.)
    for d = 1:numDemands
        @addConstraint(m, q[d] - p[d] + (theta1 - theta2)/Uset.sigmas[d] == 
                            sum{Uset.C[i, d] * xbar[i], i = 1:numDemands}) 
    end
    @addConstraint(m, sum{Uset.ubounds[d] * q[d], d=1:numDemands} - 
                      sum{Uset.lbounds[d] * p[d], d=1:numDemands} + 
                      Uset.Gamma * sqrt(numDemands) * (theta1 + theta2) + 
                      sum([Uset.means[d] / Uset.sigmas[d] for d = 1:numDemands]) * (theta1 - theta2) <= rhs)
end

#Solve facility location for a CLT Set of demand
function addDemandConstr!(facModel::FacilityModel, data::FacilityData, Uset::CLTSet)
    cNumDemands = numDemands(data)
    m = facModel.m; x = facModel.x
    @defVar(m, shipComp)
    xbar = Array(AffExpr, cNumDemands)

    #long way...  Ask Iain
    for d = 1:cNumDemands
        xbar[d] = AffExpr()
        for f = 1:numFacs(data)
            xbar[d] += (data.shipCost * dist(f, d, data) - data.revenue) * x[f, d]
        end
    end
    addRobConstr!(m, xbar, shipComp, Uset)

    facModel.shipComp = shipComp
    for f = 1:numFacs(data)
        addRobConstr!(m, x[f, :], facModel.caps[f], Uset)
    end
    @setObjective(m, Min, facModel.costComp + shipComp)
end


### just playing for now
function genDemand(numDemands::Int64, numSims::Int64)
    #For now a normal market factor, plus exponential terms
    mkt = rand(Normal(10, 2), numSims)
    t = broadcast(+, mkt, rand(Exponential(3), numSims, numDemands) )
    return transpose(t)
end

srand(8675309)
const cNumDemands = 20; const cNumFacs = 10;
const data = FacilityData(cNumFacs, cNumDemands, revenue=12.)

#gen some insample data
const demand_in = genDemand(cNumDemands, 50)


#solve a nominal problem
facModelNom = createProblem(data)
d_hat = Array(Float64, cNumDemands)
for i = 1:cNumDemands
    d_hat[i] = mean(demand_in[i, :])
end
addDemandConstr!(facModelNom, data, d_hat)
println("Means")
obj, shipComp, zvals, capvals = solve(facModelNom, data)
total_caps = Array(Float64, 1)
total_caps[1] = sum([capvals[f] for f in 1:cNumFacs])

#create some out of sample demand realizations
const numSims = 100
const demand_out = genDemand(cNumDemands, numSims)

#Nominal stage2
profitsNom = zeros(numSims)
for iRun = 1:numSims
    facModel2 = createProblem(data)
    fixFacilities( facModel2, zvals, capvals )
    addDemandConstr!(facModel2, data, demand_out[:, iRun])
    profitsNom[iRun] = solve(facModel2, data)[1]
    profitsNom[iRun] += obj - shipComp
end

const width = 2; const Gamma = 2.;
## Similar thing for the naive CLT
facModelCLT1 = createProblem(data)
Uset1 = CLTSet(demand_in, width, Gamma, false)
addDemandConstr!(facModelCLT1, data, Uset1)
println("CLT1")
objCLT1, shipCompCLT1, zvalsCLT1, capvalsCLT1 = solve(facModelCLT1, data)
push!(total_caps, sum([capvalsCLT1[f] for f in 1:cNumFacs]))

#CLT1 Stage 2
profits1 = zeros(numSims)
for iRun = 1:numSims
    facModelCLT1_ = createProblem(data)
    fixFacilities( facModelCLT1_, zvalsCLT1, capvalsCLT1 )
    addDemandConstr!(facModelCLT1_, data, demand_out[:, iRun])
    profits1[iRun] = solve(facModelCLT1_, data)[1]
    profits1[iRun] += objCLT1 - shipCompCLT1
end

###  A Better CLT Model
profits2 = Array(Float64, (numSims, 11) )
gamma_grid = linspace(1.5, 4, 11)

for iGamma in 1:11
    facModelCLT2 = createProblem(data)
    Uset2 = CLTSet(demand_in, gamma_grid[iGamma], gamma_grid[iGamma], true)
    addDemandConstr!(facModelCLT2, data, Uset2)
    println("CLT2 ", iGamma)
    objCLT2, shipCompCLT2, zvalsCLT2, capvalsCLT2 = solve(facModelCLT2, data)
    push!(total_caps, sum([capvalsCLT2[f] for f in 1:cNumFacs]))

    #CLT2 Stage 2
    for iRun = 1:numSims
        facModelCLT2_ = createProblem(data)
        fixFacilities( facModelCLT2_, zvalsCLT2, capvalsCLT2 )
        addDemandConstr!(facModelCLT2_, data, demand_out[:, iRun])
        profits2[iRun, iGamma] = solve(facModelCLT2_, data)[1]
        profits2[iRun, iGamma] += objCLT2 - shipCompCLT2
    end
end

#write out some things
using DataFrames
writetable("demandLocs.csv", DataFrame(data.demandLocs))
writetable("facilityLocs.csv", DataFrame(data.facLocs))

all_profits = hcat(profitsNom, profits1, profits2)
writetable("profits.csv", DataFrame(all_profits))

writetable("capacities.csv", DataFrame(total_caps))


