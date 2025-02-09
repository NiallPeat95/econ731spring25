#   Exploring gravity relationship in trade data:
#       X_od ∝ Y_o * Y_d / Dist_od
#
#   Niall Peat, 1/21/2025
#
cd("/Users/niallpeat/Desktop/Second Year Coursework/ECON731 - International Trade")

using Pkg; Pkg.instantiate()
using Pkg
Pkg.add("StatFiles")

import Pkg; Pkg.add("CSVFiles")

using CSVFiles,FileIO, DataFrames, FixedEffectModels, RegressionTables, Chain, Plots

# location of data files on your computer
datadir = "/Users/niallpeat/Desktop/Second Year Coursework/ECON731 - International Trade/data/"

# load trade data, three equivalent ways of doing so using standard vs "pipe" notation
# df = DataFrame(load(datadir*"tradedata/tradeDataSITC.csv"))
# df = load(datadir*"tradedata/tradeDataSITC.csv") |> DataFrame
df = datadir*"tradeDataSITC.csv" |> load |> DataFrame

# penn world tables and gravity covariates
pwt = datadir*"pwt1001.dta" |> load |> DataFrame
dist = datadir*"dist_cepii.dta" |> load |> DataFrame

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

#merg in dist_cepii
println(names(dist))
aggdf = leftjoin(aggdf,rename(dist[:,[:iso_o,:iso_d,:contig,:comlang_off,:comlang_ethno,:distw]],:iso_o => :o,:iso_d => :d),on=[:o,:d])

# USA exports in 2005
sample = @chain begin
    filter(:t => ==(2005),aggdf)
    filter(:d => ==("USA"),_)
end
plt = scatter(log.(sample.distw),log.(sample.value ./ sample.Yot),legend=false)
#plot!(plt,log.(sample.Yot),log.(sample.Yot))

#Classic gravity Regression
vcov = Vcov.cluster(:o,:d)
f = @formula log(value) ~ log(Yot) + log(Ydt) + log(distw)
result1 = reg(aggdf,f,vcov)

f = @formula log(value) ~ fe(o)&fe(t)+fe(d)&fe(t)+log(distw)
result2 = reg(aggdf,f,vcov)

regtable(result1,result2) 



