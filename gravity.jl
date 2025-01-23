
#
#   Exploring gravity relationship in trade data:
#       X_od ∝ Y_o * Y_d / Dist_od
#
#   Nels Lind, 1/21/2025
#

cd("/Users/nelslind/Dropbox/teaching/emory/2024-2025/Econ 731 Spring 2025/code/econ731spring25/")
using Pkg; Pkg.instantiate()
using FileIO, DataFrames, FixedEffectModels, RegressionTables, Chain, Plots, Statistics

# location of data files on your computer
datadir = "/Users/nelslind/Dropbox/data/"

# load trade data, three equivalent ways of doing so using standard vs "pipe" notation
# df = DataFrame(load(datadir*"tradedata/tradeDataSITC.csv"))
# df = load(datadir*"tradedata/tradeDataSITC.csv") |> DataFrame
df = datadir*"tradedata/tradeDataSITC.csv" |> load |> DataFrame

# penn world tables and gravity covariates
pwt = datadir*"PWT/pwt1001.dta" |> load |> DataFrame
dist = datadir*"CEPII/dist_cepii.dta" |> load |> DataFrame

# aggregate trade data to o-d-t level, two equivalent notations, standard vs "chain"
# aggdf = combine(groupby(df,[:o,:d,:t]),:value => sum => :value)
aggdf = @chain begin
    df
    groupby([:o,:d,:t])
    combine(:value => sum => :value)
end

# next lets merge in GDP and population data from PWT

# look at variables in pwt
println(names(pwt))

# merge first for exporters then for importers
vars = [:countrycode,:year,:rgdpe,:pop]
aggdf = leftjoin(aggdf,rename(pwt[:,vars],:countrycode => :o,:year => :t),on=[:o,:t])
rename!(aggdf,:rgdpe => :Yot,:pop => :Not)
aggdf = leftjoin(aggdf,rename(pwt[:,vars],:countrycode => :d,:year => :t),on=[:d,:t])
rename!(aggdf,:rgdpe => :Ydt,:pop => :Ndt)

# USA imports in 2005
sample = @chain begin
    filter(:t => ==(2005),aggdf)
    filter(:d => ==("USA"),_)
end
x = log.(sample.Yot) .- mean(skipmissing(log.(sample.Yot)))
y = log.(sample.value) .- mean(skipmissing(log.(sample.value)))
plt = scatter(x,y,legend=false)
plot!(plt,x,x)

# USA exports in 2005
sample = @chain begin
    filter(:t => ==(2005),aggdf)
    filter(:o => ==("USA"),_)
end
x = log.(sample.Ydt) .- mean(skipmissing(log.(sample.Ydt)))
y = log.(sample.value) .- mean(skipmissing(log.(sample.value)))
plt = scatter(x,y,legend=false)
plot!(plt,x,x)

# merge in CEPII
println(names(dist))
aggdf = leftjoin(aggdf,rename(dist[:,[:iso_o,:iso_d,:contig,:comlang_off,:comlang_ethno,:distw]],:iso_o => :o,:iso_d => :d),on=[:o,:d])

# USA exports in 2005
sample = @chain begin
    filter(:t => ==(2005),aggdf)
    filter(:o => ==("USA"),_)
end
plt = scatter(log.(sample.distw),log.(sample.value ./ sample.Ydt),legend=false)
x = log.(sample.distw) .- mean(skipmissing(log.(sample.distw)))
y = log.(sample.value ./ sample.Ydt) .- mean(skipmissing(log.(sample.value ./ sample.Ydt)))
plt = scatter(x,y,legend=false)
plot!(plt,x,-x)

# USA exports in 2005
sample = @chain begin
    filter(:t => ==(2005),aggdf)
    filter(:d => ==("USA"),_)
end
plt = scatter(log.(sample.distw),log.(sample.value ./ sample.Yot),legend=false)
x = log.(sample.distw) .- mean(skipmissing(log.(sample.distw)))
y = log.(sample.value ./ sample.Yot) .- mean(skipmissing(log.(sample.value ./ sample.Yot)))
plt = scatter(x,y,legend=false)
plot!(plt,x,-x)

# classic gravity regression
vcov = Vcov.cluster(:o,:d)
f = @formula log(value) ~ log(Yot) + log(Ydt) + log(distw)

# gravity regression consistent with CES model, use FE's to control for market access and multilateral resistance
result1 = reg(aggdf,f,vcov)
f = @formula log(value) ~ fe(o)&fe(t) + fe(d)&fe(t) + log(distw)
result2 = reg(aggdf,f,vcov)

regtable(result1,result2)