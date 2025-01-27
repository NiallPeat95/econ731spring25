#
#   Creates a dataset on bilateral trade flows and tariffs. Uses HS92 level
#       trade flow data in the BACI data set from CEPII and SITC revision 2
#       level tariff data from WTD. Additionally, need a concordance between
#       HS92 and SITC2, which comes from WITS. All these raw data sources
#       should be in the "datadir."
#
#   Raw data sources:
#       1. BACI @ CEPII: https://www.cepii.fr/DATA_DOWNLOAD/baci/data/BACI_HS92_V202401b.zip
#       2. WTD: https://www.dropbox.com/sh/2qwxckjxzj7zdek/AACo4xeB__rwYpRPOi-CUo4ha?dl=0
#       3. HS0 to SITC2 Concordance: https://wits.worldbank.org/data/public/concordance/Concordance_H0_to_S2.zip
#
#   Creates three files in "datadir":
#       1. datadir*"tradedata/countryNames.csv"
#       2. datadir*"tradedata/tradeDataSITC.csv"
#       3. datadir*"tradedata/tradeDataSITC2Digit.csv"
#
#   Nels Lind, 1/21/2025
#

using FileIO, DataFrames, Chain, Statistics

datadir = "/Users/nelslind/Dropbox/data/"

# years to load
years = 1995:2010

country_codes = datadir*"CEPII/BACI/BACI_HS92_V202401b/country_codes_V202401b.csv" |> load |> DataFrame
hs_to_sitc = @chain begin
    datadir*"WITS/concordances/Concordance_H0_to_S2/JobID-10_Concordance_H0_to_S2.csv"
    load(_,colparsers=Dict(Symbol("HS 1988/92 Product Code") => String,Symbol("SITC Revision 2 Product Code") => String))
    DataFrame
    rename(Symbol("HS 1988/92 Product Code") => :hs,Symbol("SITC Revision 2 Product Code") => :sitc)
end

df = DataFrame()
cases = []
for year in years
    println("Loading WTD data for $year...")
    @time tariffs = @chain begin
        datadir*"WTD/disaggregate_tariffs_base/tariff$(year)_base_v1.dta"
        load
        DataFrame
        select([:sitc4,:eiso,:iiso,:year,:t1,:t1_pref])
        rename(:sitc4 => :s4,:eiso => :o, :iiso => :d, :year => :t)
    end
    tariffs.t = map(x->ismissing(x) ? missing : convert(Int64,x),tariffs.t)
    tariffs.t1_pref = map(x->ismissing(x) ? missing : convert(Float64,x),tariffs.t1_pref)

    println("Loading BACI trade flow data for $year...")
    @time flows = @chain begin
        load(datadir*"CEPII/BACI/BACI_HS92_V202401b/BACI_HS92_Y$(year)_V202401b.csv",
            colparsers=Dict(:k => String))
        DataFrame
        select([:t,:i,:j,:k,:v])
        leftjoin(_,rename(country_codes[:,Not(:country_iso2)],:country_code => :i,
            :country_name => :exporter,:country_iso3 => :o),on=[:i])
        leftjoin(_,rename(country_codes[:,Not(:country_iso2)],:country_code => :j,
            :country_name => :importer,:country_iso3 => :d),on=[:j])
        rename(:v => :value,:k => :hs)
        select(Not([:i,:j]))
        leftjoin(hs_to_sitc[:,[:hs,:sitc]],on=[:hs])
    end

    # 334 and 33452 exist in flows, while we have 3341, 3342, 3343, 3344, 3345 in tariffs
    #   1. reduce 3341, 3342, 3343, 3344 to 334 in tariffs
    #   2. sitc codes in flows that are missing correspond to hs codes that didn't match, classify as "XXXXX"
    #   3. truncate 5 digit codes in flows to 4 digit to match with 4 digit codes in tariffs
    tariffs.s = map(x->x ∈ ["3341","3342","3343","3343"] ? "334" : x,tariffs.s4)
    flows.sitc = map(x->ismissing(x) ? "XXXXX" : x, flows.sitc)
    flows.s = map(x->length(x)==5 ? x[1:4] : x,flows.sitc)

    tariffs_s_codes = unique(tariffs.s)
    flows_s_codes = unique(flows.s)
    cases = union(cases,setdiff(flows_s_codes,tariffs_s_codes))

    tariffs = @chain begin
        tariffs
        groupby([:s,:o,:d,:t])
        combine(:t1 => mean => :t1,:t1_pref => mean => :t1_pref)
    end
    append!(df,leftjoin(flows,tariffs,on=[:s,:o,:d,:t]))
end
# cases is the accumulated set of all codes in flows that aren't matched to a code in tariffs.
@show cases
df = @chain begin
    df
    groupby([:s,:o,:d,:t,:importer,:exporter])
    combine(:value => sum => :value,:t1 => mean => :t1,:t1_pref => mean => :t1_pref)
end

# about 20% missing for the two tariff measures
@show mean(ismissing.(df.value))
@show mean(ismissing.(df.t1))
@show mean(ismissing.(df.t1_pref))

df = sort(df[:,[:s,:o,:d,:t,:exporter,:importer,
    :value,:t1,:t1_pref]],reverse([:s,:o,:d,:t]))

dir = datadir*"tradedata/"
isdir(dir) || mkpath(dir)

# save space by splitting off country names, can merge back in as needed
countries = unique(vcat(rename(df[:,[:o,:exporter]],:o=>:n,:exporter=>:country),
    rename(df[:,[:d,:importer]],:d=>:n,:importer=>:country)))

save(dir*"countryNames.csv",countries)
save(dir*"tradeDataSITC.csv",df[:,Not([:importer,:exporter])])

# aggregate to two digit level

# merge in first year values
df = @chain begin
    filter(:t => ==(1995),df)
    select([:s,:o,:d,:value])
    rename(_,:value => :value95)
    leftjoin(df,_,on=[:s,:o,:d])
end

# reduce to two digit codes
df.s = map(x->x[1:2],df.s)

# aggregate by summing value, and taking the mean of tariffs with various weights
function weightedMean(w,x) 
    if all(ismissing.(w .* x))
        return missing
    else
        return sum(skipmissing(w .* x)) / sum(skipmissing(w .* ismissing.(x)))
    end
end
df = @chain begin
    df
    groupby([:s,:o,:d,:t,:exporter,:importer])
    combine(
        :value => sum => :value,
        :t1 => (x->mean(skipmissing(x))) => :t1_mean,
        [:value,:t1] => weightedMean => :t1_wmean,
        [:value95,:t1] => weightedMean => :t1_wmean95,
        :t1_pref => (x->mean(skipmissing(x))) => :t1_pref_mean,
        [:value,:t1_pref] => weightedMean => :t1_pref_wmean,
        [:value95,:t1_pref] => weightedMean => :t1_pref_wmean95
    )
end

save(dir*"tradeDataSITC2Digit.csv",df[:,Not([:importer,:exporter])])