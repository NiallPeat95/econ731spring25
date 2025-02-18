#
#   Using hat algebra to solve DEK
#
#   Nels Lind, 2/13/2025
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
    # tdf.year .= year
    append!(df,tdf)
    print("took $(round(time()-t0,digits=3)) seconds.")
end

colnames = names(df)[6:end-1]
countries = unique(map(x->x[1:3],colnames))
sectors = unique(map(x->x[4:end],colnames))

N = length(countries)



display(colnames[1:56])
display(colnames[1:56])

# number of commodities
M = 56

display(colnames[M+1:M+M])
display(colnames[1:N*M])
display(colnames[N*M+1:end])


# number of final uses
K = 5
display(colnames[N*M+1:N*M+10])
display(colnames[N*M+1:N*M+N*K])

dfInt = df[:,vcat(1:5,5+1:5+N*M)]


industries = unique(df[:,[:IndustryCode,:IndustryDescription]])

print(industries)










@time Ŵ = tâtonnment(m,D′,Dm′,report=true,reportrate=1)

