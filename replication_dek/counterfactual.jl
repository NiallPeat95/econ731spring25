#
#   Using hat algebra to solve DEK
#
#   Nels Lind, 2/13/2025
#

cd("/Users/niallpeat/Documents/GitHub/econ731spring25/replication_dek/")
using Pkg; Pkg.add(["FileIO","DataFrames","Chain","Plots","Distributions","LinearAlgebra", "RData"])
using Pkg; Pkg.activate("."); Pkg.instantiate()
using FileIO, DataFrames, Chain, Plots
using Distributions, LinearAlgebra, RData

# location of data files on your computer
datadir = "/Users/niallpeat/Dropbox/International Trade/Data/"

# load WIOD
years = 2000:2014
df = DataFrame()
for year in years
    t0 = time()
    print("\nLoading WIOD data for $year...")
    file = datadir*"WIOD/WIOTS_in_R/WIOT$(year)_October16_ROW.RData"
    tdf = load(file)["wiot"]
    # tdf.year .= year
    append!(df,tdf)
    print("took $(round(time()-t0,digits=3)) seconds.")
end

# make the row number "RNr" variable an integer
df.RNr = string.(map(Int64,df.RNr))
df.Year = map(Int64,df.Year)
rename!(df,:RNr => :i)

# Within each year, WIOD is organized into inputs/factors (rows) and
#   outputs/final-uses (columns) by industries × origins + destination
#   factors (rows), and industries × destinations + final-uses × destinations (columns)
#
# When loading from the RData files, we don't get full information on the columns, but
#   inspection of the XLSB versions of the WIOD we can organize as follows.

# For any given column (industry/use × destination), the non-industry inputs are
#   the rows corresponding to:
totalNames = @chain begin
    df[df.Country .== "TOT",[:IndustryCode,:IndustryDescription,:i]]
    unique
    Dict(x.i => x.IndustryDescription for x in eachrow(_))
end
Mtotal = length(keys(totalNames))

# For any given column (industry/use × destination), the industry intermediate inputs are
#   the rows corresponding to:
industryNames = @chain begin
    df[df.Country .!== "TOT",[:IndustryCode,:IndustryDescription,:i]]
    unique
    Dict(x.i => x.IndustryDescription for x in eachrow(_))
end
Mindustries = length(keys(industryNames))

# The set of countries is
countries = filter(x->x != "TOT",unique(df.Country))
N = length(countries)

# For any given row (industry × origin + destination factor), a column is either
#   an industry within a country on the output side or a final use within a country

# The first 5 columns correspond to information about the rows/year
names(df[:,1:5])

# The remaining columns have variable names corresponding to an identifying integer
#   paired with a country
names(df[:,5+1:end])

# In these labels, the first 56 integers correspond to the 56 industry intermediate inputs,
#   while the remaining correspond to final uses. The intermediate input demands are organized
#   by country.
names(df[:,5+1:5+Mindustries])
names(df[:,5+Mindustries+1:5+2*Mindustries])
names(df[:,5+2*Mindustries+1:5+3*Mindustries])

# This continues across all N countries, with the final country being the rest-of-world aggregate.
names(df[:,5+(N-1)*Mindustries+1:5+N*Mindustries])

# From there, the columns are final uses, which have identifying integers of 57 to 61.
names(df[:,5+N*Mindustries+1:end])

# The XLSB version of the data describes these final use columns. The labels are:
finalUseNames = Dict(
    "57" => "Final consumption expenditure by households",
    "58" => "Final consumption expenditure by non-profit organisations serving households",
    "59" => "Final consumption expenditure by government",
    "60" => "Gross fixed capital formation",
    "61" => "Changes in inventories and valuables"
)

allNames = merge(totalNames,industryNames,finalUseNames)

# reshape from wide to long
df = stack(df,names(df)[6:end-1])

rename!(df,:Country => :o,:i => :input,:Year => :t)
df.d = map(x->x[1:3],df.variable)
df.output = map(x->x[4:end],df.variable)
df.inputDesc = map(x->allNames[x],df.input)
df.outputDesc = map(x->allNames[x],df.output)
select!(df,Not(:variable,:IndustryCode,:IndustryDescription))

# manufacturing
for x in eachrow(unique(df[:,[:input,:inputDesc]]))
    @show x.input,x.inputDesc
end
manufacturingNames = Dict("$i" => industryNames["$i"] for i=5:22)

# Y_n = W_n*L_n
# X_n = Y_n + D_n
# Xm_n = Ym_n + Dm_n
# Xm_n = α*X_n + (1-β)*Ym_n
# Ym_n = ∑_d Π_nd * Xm_d

#  Dm_n = Xm_n - Ym_n
#       = Xm_n - ∑_d Π_nd * Xm_d
#       = (1 - Π_nn)*Xm_n - ∑_d≠n Π_nd * Xm_d
#       = ∑_o≠n Π_on)*Xm_n - ∑_d≠n Π_nd * Xm_d

#
# GrossOutput_n = ∑_ijd M_injd + ∑_id F_ind
# Expenditure_n = ∑_ijo M_iojn + ∑_io F_ion
# IntermediateUse_n = ∑_ioj M_iojn
# FinalExpenditure_n = ∑_io F_ion
# Deficit_n = ∑_ijk≠n (M_injk - M_ikjn) + ∑_ik≠n (F_ink - F_ikn)
# ValueAdded_n = GrossOutput_n - IntermediateUse_n
#              = ∑_ijk (M_injk - M_ikjn) + ∑_id F_ind 
#              = ∑_ijk≠n (M_injk - M_ikjn) + ∑_id≠n F_ind + ∑_i F_inn
#              = ∑_ijk≠n (M_injk - M_ikjn) + ∑_id≠n F_ind + X_n - ∑_io≠n F_ion
#              = Exports_n - Imports_n + FinalExpenditure_n
#              = Deficit_n + FinalExpenditure_n

dfCountry = @chain begin
    filter(x->in(x.input,keys(industryNames)),df)
    groupby([:o,:t])
    combine(:value => sum => :GrossOutput)
    rename(:o=>:n)
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(industryNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :Expenditure)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(industryNames)) && !in(x.output,keys(finalUseNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :IntermediateUse)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(industryNames)) && in(x.output,keys(finalUseNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :FinalExpenditure)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry.ValueAdded = dfCountry.GrossOutput - dfCountry.IntermediateUse
dfCountry.Deficit = dfCountry.ValueAdded - dfCountry.FinalExpenditure

extrema(dfCountry.Deficit)

# GrossOutputManuf_n = ∑_jd ∑_i∈Manuf M_injd + ∑_d ∑_i∈Manuf F_ind
# ExpenditureManuf_n = ∑_jo ∑_i∈Manuf M_iojn  + ∑_o ∑_i∈Manuf F_ion
# IntermediateUseManuf_n = ∑_oj ∑_i∈Manuf M_iojn
# FinalExpenditureManuf_n = ∑_k ∑_i∈Manuf F_ikn
# DeficitManuf_n = ∑_jk≠n ∑_i∈Manuf ( M_ikjn - M_injk ) + ∑_k≠n ∑_i∈Manuf ( F_ikn - F_ink )
#                = ∑_jk ∑_i∈Manuf ( M_ikjn - M_injk ) + ∑_k ∑_i∈Manuf ( F_ikn - F_ink )
#                = ∑_jk ∑_i∈Manuf M_ikjn  + ∑_k ∑_i∈Manuf F_ikn - GrossManufOutput_n
#                = ExpenditureManuf_n - GrossManufOutput_n

dfCountry = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)),df)
    groupby([:o,:t])
    combine(:value => sum => :GrossOutputManuf)
    rename(:o=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :ExpenditureManuf)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)) && !in(x.output,keys(finalUseNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :IntermediateUseManuf)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)) && in(x.output,keys(finalUseNames)),df)
    groupby([:d,:t])
    combine(:value => sum => :FinalExpenditureManuf)
    rename(:d=>:n)
    leftjoin(dfCountry,_,on=[:n,:t])
end
dfCountry.ValueAddedManuf = dfCountry.GrossOutputManuf - dfCountry.IntermediateUseManuf
dfCountry.DeficitManuf = dfCountry.ExpenditureManuf - dfCountry.GrossOutputManuf

# β = value added share of manufacturing
# α = manufactured good share of final expenditure

dfCountry.α = dfCountry.FinalExpenditureManuf ./ dfCountry.FinalExpenditure
dfCountry.β = dfCountry.ValueAddedManuf ./ dfCountry.GrossOutputManuf

mean(dfCountry.α)
mean(dfCountry.β)

histogram(dfCountry.α)
histogram(dfCountry.β)

dfTrade = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)),df)
    groupby([:o,:d,:t])
    combine(:value => sum => :TradeFlowManuf)
    leftjoin(_,rename(dfCountry[:,[:n,:t,:ExpenditureManuf]],:n=>:d),on=[:d,:t])
end
dfTrade.TradeShareManuf = dfTrade.TradeFlowManuf ./ dfTrade.ExpenditureManuf


include("DEK.jl")

α = .188
#α = dfCountry.α
β = .312
#β = dfCountry.β
θ = 8.28
t = 2004
countries = sort(unique(dfCountry.n))
Π = reshape(sort(filter(x->x.t==t,dfTrade),reverse([:o,:d])).TradeShareManuf,N,N)
Y = sort(filter(x->x.t==t,dfCountry),:n).ValueAdded
D = sort(filter(x->x.t==t,dfCountry),:n).Deficit
Ym = sort(filter(x->x.t==t,dfCountry),:n).GrossOutputManuf
Dm = sort(filter(x->x.t==t,dfCountry),:n).DeficitManuf
m = DEK(Π,Y,D,Dm,α,β,θ)

# DEK counterfactual of zeroing out current accounts assuming CA = - D
Dm′ = Dm - D
D′ = zeros(N)
T̂ = ones(N)
τ̂ = ones(N,N)
Ŵ = tâtonnment(m,T̂,τ̂,D′,Dm′,report=true,reportrate=1,λ=.001, maxit=1e6,tol=1e-7)
P̂ = prices(m,Ŵ,T̂,τ̂)
fsz = 8
psz = (600,500)
scatter(-D./Y,Ŵ./P̂,legend=false,
    alpha=0,series_annotations=Plots.series_annotations(countries,Plots.font(fsz,:black)),size=psz,fontfamily=:Times)