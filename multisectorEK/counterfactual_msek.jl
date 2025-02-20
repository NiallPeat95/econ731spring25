#
#   Using hat algebra to solve DEK
#
#   Nels Lind, 2/18/2025
#

username = "nelslind"
cd("/Users/$username/Dropbox/teaching/emory/2024-2025/Econ 731 Spring 2025/code/econ731spring25/")
using Pkg
Pkg.activate("."); Pkg.instantiate()
using FileIO, DataFrames, Chain, Plots
using Distributions, LinearAlgebra

# location of data files on your computer
datadir = "/Users/$username/Dropbox/data/"

# load WIOD
years = 2000:2014
df = DataFrame()
for year in years
    t0 = time()
    print("\nLoading WIOD data for $year...")
    file = datadir*"WIOD/WIOTS_in_R/WIOT$(year)_October16_ROW.RData"
    tdf = load(file)["wiot"]
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
Mtotal = length(keys(inputNames))

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

dfCountry.α = dfCountry.FinalExpenditureManuf ./ dfCountry.FinalExpenditure
dfCountry.β = dfCountry.ValueAddedManuf ./ dfCountry.GrossOutputManuf

mean(dfCountry.α)
mean(dfCountry.β)

dfTrade = @chain begin
    filter(x->in(x.input,keys(manufacturingNames)),df)
    groupby([:o,:d,:t])
    combine(:value => sum => :TradeFlowManuf)
    leftjoin(_,rename(dfCountry[:,[:n,:t,:ExpenditureManuf]],:n=>:d),on=[:d,:t])
end
dfTrade.TradeShareManuf = dfTrade.TradeFlowManuf ./ dfTrade.ExpenditureManuf

save(,dfCountry)







# # verify that TOT corresponds to output
# dfOutput = @chain begin
#     filter(x->in(x.output,keys(industryNames)) && x.inputDesc == "Output at basic prices",df)
#     select(_,[:output,:d,:t,:value])
#     rename(_,:output=>:industry,:d=>:n)
# end

# # verify that TOT corresponds to output
# @chain begin
#     filter(x->in(x.input,keys(industryNames)),df)
#     select(_,[:input,:o,:t,:TOT])
#     unique
#     rename(_,:input=>:industry,:o=>:n)
#     leftjoin(dfOutput,_,on=[:industry,:n,:t])
#     extrema(_.value .- _.TOT)
# end

# dfIntermediateConsumption = @chain begin
#     filter(x->in(x.output,keys(industryNames)) && x.inputDesc == "Total intermediate consumption",df)
#     select(_,[:output,:d,:t,:value])
#     rename(_,:output=>:industry,:d=>:n)
# end

# # aggregate 




# # verify that total intermediate consumption correspond to aggregating over rows
# @chain begin
#     filter(x->in(x.input,keys(industryNames)) && in(x.output,keys(industryNames)),df)
#     groupby([:output,:d,:t])
#     combine(:value => sum => :SUM)
#     rename(_,:output=>:industry,:d=>:n)
#     leftjoin(dfIntermediateConsumption,_,on=[:industry,:n,:t])
#     extrema(_.value .- _.SUM)
# end

# # final expenditure
# dfFinalExpenditure = @chain begin
#     filter(x->in(x.input,keys(industryNames)) && in(x.output,keys(finalUseNames)),df)
#     groupby([:input,:d,:t])
#     combine(:value => sum => :value)
#     rename(_,:input=>:industry,:d=>:n)
# end

# # what industries are non-traded?
# @chain begin
#     filter(x->in(x.input,keys(industryNames)),df)
#     groupby([:input,:inputDesc])
#     combine([:value,:o,:d] => ( (v,o,d) -> sum((o.!=d).*v )/sum(v)) => :importShare)
#     rename(_,:input=>:industry)
#     sort(_,[:importShare])
#     plot(1:56,_.importShare)
# end



# @chain begin
#     rename(dfIntermediateConsumption,:value=>:M)
#     leftjoin(dfGrossOutput,_,on=[:industry,:n,:t])
# end



# α is the share of final goods expenditure on intermediates


# df




# # aggregate final uses
# df.input = map(x->in(x,keys(industryNames)) ? "M" : x,df.input)
# df.output = map(x->in(x,keys(industryNames)) ? "M" : x,df.output)
# df.output = map(x->in(x,keys(finalUseNames)) ? "F" : x,df.output)
# df = @chain begin
#     df
#     groupby([:o,:d,:t,:input,:output])
#     combine(:value => sum => :value,:TOT => sum => :total)
# end

# unique(df.inputDesc)
# unique(df.outputDesc)

# df

# filter(x->x.o == "TOT" && x.d == "USA" && x.t == 2014,df)
# filter(x->x.o != "TOT" && x.d == "USA" && x.t == 2014,df)

# # there is no data on output or value added for final goods
# extrema(filter(x->x.inputDesc == "Value added at basic prices" && x.output == "F",df).value)
# extrema(filter(x->x.inputDesc == "Output at basic prices" && x.output == "F",df).value)

# # income is 


# tdf = @chain begin
#     filter(x->x.output == "F",df)
# end



# tdf = @chain begin
#     filter(x->x.output == "M",df)
#     groupby([:input,:d,:t])
#     combine(:value => sum => :value)
# end



# extrema(filter(x->x.output=="final" && x.o!="TOT",df).value)
# extrema(filter(x->x.output=="final" && x.o!="TOT",df).total)
# extrema(filter(x->x.output=="final" && x.o=="TOT",df).value)
# extrema(filter(x->x.output=="final" && x.o=="TOT",df).total)

# allNames["M"] = "Intermediate Goods"
# allNames["F"] = "Final Goods"
# df.inputDesc = map(x->allNames[x],df.input)
# df.outputDesc = map(x->allNames[x],df.output)


# dfOutput = @chain begin
#     filter(x->x.inputDesc=="Output at basic prices",df)
#     filter(x->x.output!="F",_)
#     extrema(_.value)
#     # groupby([:d,:t])
#     # combine(:value => sum => :value)
# end


# dfIntermediatesInFinal = @chain begin
#     filter(x->x.output=="final" && x.o !="TOT",df)
#     groupby([:d,:t])
#     combine(:value => sum => :value)
# end




# filter(x->x.IndustryCode == "VA",df)


# # what goods are non traded?
# df.trade = df.o != df.d
# @chain begin
#     df
#     filter(x->x.o!="TOT",_)
#     groupby(_,[:input,:trade])
#     collapse(:value => sum)
# end




# # aggregate industries
# df.in = map(x->in(x,namesIndustries.i) ? "0" : x,df.input)
# df.out = map(x->in(x,namesIndustries.i) ? "0" : x,df.output)
# unique(df.in)
# unique(df.out)
# df = @chain begin
#     df
#     groupby([:in,:out,:o,:d,:t])
#     combine(:value => sum => :value,:TOT => sum => :total)
# end

# inNames = Dict(x.i => x.IndustryDescription for x in eachrow(namesInputs))
# inNames["0"] = "Intermediate"

# outNames = Dict(
#     "0"  => "Production",
#     "57" => "Final consumption expenditure by households",
#     "58" => "Final consumption expenditure by non-profit organisations serving households",
#     "59" => "Final consumption expenditure by government",
#     "60" => "Gross fixed capital formation",
#     "61" => "Changes in inventories and valuables"
# )

# df.inDescription = map(x->inNames[x],df.in)
# df.outDescription = map(x->outNames[x],df.out)


# display(filter(x->x.in != "0" && x.t == 2010,df))
# display(filter(x->x.d == "USA" && x.t == 2010,df))

# unique(df[:,[:in,:inDescription]])
# unique(df[:,[:out,:outDescription]])


# # intermediate inputs (in=="0") aggregate across o to total intermeadiate consumption (in="65")
# @chain begin
#     df
#     filter(x->in(x.in,["0","65"]),_)
#     groupby([:in,:out,:d,:t])
#     combine(:value => sum => :value)
# end


# # everything adds up to output (not quite!)
# tdf = copy(df)
# tdf.in = map(x->x=="73",tdf.in)
# tdf.o = map(x->x=="TOT",tdf.o)
# @chain begin
#     tdf
#     groupby([:in,:out,:o,:d,:t])
#     combine(:value => sum => :value)
#     filter(x->x.d == "USA" && x.t == 2014,_)
# end

# # everything adds up to output
# tdf = copy(df)
# tdf.in = map(x->in(x,["0","65","70","73"]) ? x : "Other",tdf.in)
# tdf.inDesc = map(x->x=="Other" ? x : inNames[x],tdf.in)
# tdf.o = map(x->x=="TOT",tdf.o)
# @chain begin
#     tdf
#     groupby([:in,:inDesc,:out,:o,:d,:t])
#     combine(:value => sum => :value)
#     filter(x->x.out == "0" && x.d == "USA" && x.t == 2014,_)
# end


# extrema(filter(x->x.o == "TOT" && x.d == "USA" && x.t == 2014,df).total)
# extrema(filter(x->x.o == "TOT" && x.d == "USA" && x.t == 2014,df).total)


# # combine other into value added
# tdf = copy(df)
# # tdf.in = map(x->in(x,["0","65","70","73"]) ? x : "70",tdf.in)
# # tdf.inDescription = map(x->x=="Other" ? x : inNames[x],tdf.in)
# tdf.out = map(x->x=="0",tdf.out)
# tdf.outDescription = map(x->x ? "Intermediate" : "Final",tdf.out)
# tdf = @chain begin
#     tdf
#     groupby([:in,:inDescription,:out,:outDescription,:o,:d,:t])
#     combine(:value => sum => :value)
# end

# inNames

# dfY = select(filter(x->x.o=="TOT" && x.outDescription=="Intermediate" && x.in=="73",tdf),[:d,:t,:value])
# rename!(dfY,:value=>:Y)

# extrema(filter(x->x.o=="TOT" && x.outDescription=="Final" && x.in=="73",tdf).value)
# extrema(filter(x->x.o=="TOT" && x.outDescription=="Final" && x.in=="70",tdf).value)

# dfXf = select(filter(x->x.o!="TOT" && x.outDescription=="Final",tdf),[:o,:d,:t,:value])
# rename!(dfXf,:value=>:Xf)


# rename!(dfOutput,:value=>:Y)



# dfTotalFinal = filter(x->x.o=="TOT" && x.outDescription=="Final",tdf)
# dfTotalFinal.in = "X".*dfTotalFinal.in
# dfTotalFinal = unstack(dfTotalFinal,[:d,:t],:in,:value)
# extrema(dfTotalFinal.X73)



# # can drop 67 68 69
# namesInputs


# dfTotal.variable = map(x->Dict("65" => "M","70"=>"VA","73"=>"Y")[x.in]*"_"*x.outDescription,eachrow(dfTotal))
# dfTotal = unstack(dfTotal,[:d,:t],:variable,:value)





# dfTotal = filter(x->x.o=="TOT",tdf)
# dfTotal.variable = map(x->Dict("65" => "M","70"=>"VA","73"=>"Y")[x.in]*"_"*x.outDescription,eachrow(dfTotal))
# dfTotal = unstack(dfTotal,[:d,:t],:variable,:value)


# dfInt = filter(x->x.o!="TOT",tdf)
# dfInt = unstack(dfInt,[:o,:d,:t],:outDescription,:value)

# dfInt = leftjoin(dfInt,dfTotal[:,[:d,:t,:M_Intermediate,:M_Final]],on=[:d,:t])
# dfInt.share_intermediate = dfInt.Intermediate ./ dfInt.M_Intermediate
# dfInt.share_final = dfInt.Final ./ dfInt.M_Final
# select!(dfInt,[:o,:d,:t,:share_intermediate,:share_final])

# scatter(dfInt.share_intermediate,dfInt.share_final,legend=false)

# extrema(dfTotal.M_Intermediate .+ dfTotal.VA_Intermediate - dfTotal.Y_Intermediate)
# dfInt