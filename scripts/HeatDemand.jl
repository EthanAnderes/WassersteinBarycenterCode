# ============================================
#  Barycenter example
# ============================================


# --------------------------------------------
#  Load Packages, set paths and load the data.
# --------------------------------------------

using DataFrames
using JuMP
using Clp
using PyCall
using PyPlot


# ---- Set directory paths to code and the paper
figurepath = joinpath(pwd(),"figures")
codepath   = joinpath(pwd(),"scripts")



# -----------  Load the data
# data from https://en.wikipedia.org/wiki/Climate_of_California
citydf = DataFrame(
	city = [
		:Bakersfield,
		:Eureka,
		:Fresno,
		:LosAngeles,
		:Sacramento,
		:SanBernardino,
		:SanFrancisco,
		:SanJose,
		:SouthLakeTahoe ],
	lat_y      = [35.3667,   40.8019,   36.7500,   34.0500,   38.5556,	 34.1000,   37.7833,	37.3382,	38.9400],
	lon_x      = [-119.0167, -124.1636, -119.7667, -118.2500, -121.4689, -117.3000, -122.4167,	-121.8863,	-119.9769],
	population = [363_630,   26_913,    509_924,   3_884_000, 479_686,	 2_088_000, 837_442,	998_537,	21_387],
	jan = [57, 56, 55, 71, 54, 69, 58, 61, 42],
	feb = [64, 57, 62, 72, 61, 71, 61, 64, 43],
	mar = [70, 57, 68, 75, 66, 75, 63, 67, 48],
	apr = [76, 58, 75, 78, 72, 80, 64, 73, 53],
	may = [84, 61, 84, 80, 80, 85, 66, 77, 63],
	jun = [92, 63, 92, 83, 87, 92, 68, 82, 72],
	jul = [98, 64, 99, 90, 92, 10, 68, 84, 80],
	aug = [96, 65, 97, 92, 92, 103,69, 84, 80],
	sep = [91, 64, 91, 91, 88, 100, 71, 81, 73],
	oct = [80, 62, 79, 84, 78, 86, 70, 76, 62],
	nov = [66, 58, 65, 78, 64, 75, 64, 65, 50],
	dec = [57, 55, 55, 71, 54, 68, 58, 61, 42],
)


# ---------------------------------------------------
# Transform temp data and package into a custom type
# ---------------------------------------------------
# This is basically a named array type.
# It is designed to contain the images in the columns (indexed by month)
immutable city_x_month_image <: AbstractMatrix{Float64}
	imagematrix::Array{Float64,2}
	city_nms ::Array{Symbol,1}
	month_nms::Array{Symbol,1}
end

# here is a constructor
# This allows one to
function city_x_month_image(df::DataFrame, months::Symbol...)
	nrw = length(df[:city])
	ncl = length(months)
	transform_normalize(tempu, pop)  = (z=pop.*((tempu.-72).^2); return z./sum(z))
	imageslocal = zeros(Float64, (nrw, ncl))
	for cl = 1:ncl
		imageslocal[:,cl] = transform_normalize(citydf[months[cl]], citydf[:population])
	end
	return city_x_month_image(imageslocal, df[:city], [months...])
end
city_x_month_image(df::DataFrame) = city_x_month_image(df::DataFrame, :jan,:feb,:mar,:apr,:may,:jun,:jul,:aug,:sep,:oct,:nov,:dec)

Base.size(im::city_x_month_image)   = size(im.imagematrix)
Base.length(im::city_x_month_image) = length(im.imagematrix)
Base.getindex(im::city_x_month_image, i::Symbol, j::Symbol) = im.imagematrix[find(im.city_nms .== i)[1], find(im.month_nms .== j)[1]]
Base.getindex(im::city_x_month_image, i::Symbol, j) = im.imagematrix[find(im.city_nms .== i)[1], j]
Base.getindex(im::city_x_month_image, i, j::Symbol) = im.imagematrix[i, find(im.month_nms .== j)[1]]
Base.getindex(im::city_x_month_image, i, j) = im.imagematrix[i, j]
months(im::city_x_month_image) = im.month_nms
cities(im::city_x_month_image) = im.city_nms

# make an instance of this type
images = city_x_month_image(citydf, :dec, :jan, :feb, :mar, :jun, :jul, :aug, :sep)
# now we can index as follows: images[:Sacramento, :jan]






# ---------------------------------------------------
#   Get barycenter support
# ---------------------------------------------------

# The following object iterates over the number of ways to choose 8 locations out of 20.
function get_support(citydf, images)
	ncities      = length(cities(images))
	nmonths      = length(months(images))
	ndividers    = ncities - 1          # 8
	nones        = nmonths + ndividers  # 20
	itercomb     = combinations(1:nones, ndividers)

	bary_support = zeros(Float64, (2, length(itercomb))) # columns give the support coordinates
	for (index, value) in enumerate(itercomb)
		tmp        = vcat(value[1], diff(sort(value))) .- 1 # this gives the counts for the first 8  partitions
		citycounts = vcat(tmp, nmonths - sum(tmp))
		for i = 1:ncities
			bary_support[:, index] .+= citycounts[i].*[citydf[i, :lon_x], citydf[i, :lat_y]]
		end
		bary_support[:, index] .*= (1.0 / nmonths)
	end
	return bary_support
end
@time bary_support = get_support(citydf, images)


# ------------------------------------------------------------
# Define and solve the Wasserstein Barycenter linear program
# -----------------------------------------------------------

# ------- initialize distance matrix and variable size
N = size(images, 2)        # number of images
n = size(images, 1)        # each image size
b = size(bary_support, 2)  # barycenter size
zk = (citydf[:lon_x] .- bary_support[1,:]).^2  + (citydf[:lat_y] .- bary_support[2,:]).^2

# ------- define the solver
m = Model(solver = ClpSolver(SolveType=5))

# ------ set variables and constraints
@defVar(m, p[1:b] >= 0 )
@defVar(m, mk[1:n, 1:b, 1:N] >= 0 )
for i = 1:n, k = 1:N
	@addConstraint(m, sum{mk[i, j, k], j = 1:b} == images[i, k])
end
for j = 1:b, k = 1:N
	@addConstraint(m, sum{mk[i, j, k], i = 1:n} == p[j])
end

# ------ set objective
@setObjective(m, Min, sum{zk[i, j] * mk[i, j, k], i=1:n, j=1:b, k=1:N})

# ------- solve it
@time status = solve(m)

# ----- get the solution and the support
pval  = getValue(p[:])
mkval = Array(Float64, n, b, N)
for kntr = 1:N
	mkval[:,:,kntr] = getValue(mk[:,:,kntr])
end
psupport = find(pval)




# ------------------------------------------------------------
# Generate figure1* (the individual month demands)
# ------------------------------------------------------------


@pyimport mpl_toolkits.basemap as basemap
Basemap =  basemap.Basemap
lat_0, lon_0 =37.0, -119.3
map = Basemap(
	projection="merc",
	lat_0=lat_0,
	lon_0=lon_0,
    resolution = "h",
    area_thresh=1000,
    llcrnrlon=lon_0 - 5.2 ,
    llcrnrlat=lat_0 - 5.2,
    urcrnrlon=lon_0 + 5.2 ,
    urcrnrlat=lat_0 + 5.2
    )

# downloaded from http://www.census.gov/cgi-bin/geo/shapefiles2010/main
shapename = "tl_2010_06_state10"
shapefile = joinpath(codepath, "$shapename/$shapename")


# ---- plot individual months
function monthplots(setmonth::Symbol)
	map[:readshapefile](shapefile, shapename, linewidth=1.5)
	map[:scatter](citydf[:lon_x], citydf[:lat_y], latlon = true, s=5000*images[:,find(months(images) .== setmonth)], c=[0,0,0], alpha=0.5)
	map[:scatter](citydf[:lon_x], citydf[:lat_y], latlon = true,
		s=50,
		linewidths=2,
		marker="x",
		c=[1,0,0],
		alpha=1)
	axis("off")
end

monthplots(:feb)
savefig(joinpath(figurepath, "figure1a.pdf"), dpi=300, bbox_inches="tight", transparent=true)


monthplots(:mar)
savefig(joinpath(figurepath, "figure1b.pdf"), dpi=300, bbox_inches="tight", transparent=true)


monthplots(:jun)
savefig(joinpath(figurepath, "figure1c.pdf"), dpi=300, bbox_inches="tight", transparent=true)


monthplots(:jul)
savefig(joinpath(figurepath, "figure1d.pdf"), dpi=300, bbox_inches="tight", transparent=true)



# -------------------------------------------------------------------
# Generate figure2* (the barycenter and the possible support points)
# -------------------------------------------------------------------

# ---- barycenter  plot
map[:readshapefile](shapefile, shapename, linewidth=1.5)
map[:scatter](bary_support[1,psupport][:], bary_support[2,psupport][:],  latlon = true,
	s=5000*pval[psupport],
	c=[0,0,0],
	alpha=0.5)
map[:scatter](citydf[:lon_x], citydf[:lat_y], latlon = true,
	s=50,
	linewidths=2,
	marker="x",
	c=[1,0,0],
	alpha=1)
axis("off")
savefig(joinpath(figurepath, "figure1e.pdf"), dpi=300, bbox_inches="tight", transparent=true)


# ---- barycenter support,
map[:readshapefile](shapefile, shapename, linewidth=1.5)
map[:scatter](bary_support[1,:][:], bary_support[2,:][:], latlon = true, s = 10,  marker=".", c=[1,1,1], alpha=0.2)
map[:scatter](citydf[:lon_x], citydf[:lat_y], latlon = true,
	s=50,
	linewidths=2,
	marker="x",
	c=[1,0,0],
	alpha=1)
axis("off")
savefig(joinpath(figurepath, "figure1f.pdf"), dpi=300, bbox_inches="tight", transparent=true)





# -------------------------------------------------------------------
# Generate figure3* (the transportation maps)
# -------------------------------------------------------------------

# barycntr_trnsprt_map:
#	- Each row corresponds to a barycenter support point
#   - Columns will correspond to month and contain the name of the servicing city.
barycntr_trnsprt_map = DataFrame(bary_supp_indx = psupport)

for iim = 1:length(months(images))   # iim indexes the months
	findcity = Array(Int64, length(psupport)) # initialize the container that  this should contain entries 1,2,...,9 which correspond to the rows of images.
	for kr = 1:length(psupport)
		_, findcity[kr] = findmax(mkval[:, psupport[kr], iim])
	end
	barycntr_trnsprt_map[months(images)[iim]] = cities(images)[findcity]
end
@assert barycntr_trnsprt_map[:bary_supp_indx] == psupport


# cit_srvic_fun:
#   - A function cit_srvic_map(:month, :city) => the barycenter support points which serve that city for that month
cit_srvic_fun(mnth, cit) =  (barycntr_trnsprt_map[cit .== barycntr_trnsprt_map[mnth], :bary_supp_indx])
for mnth in months(images)
	@assert length(psupport) == sum([length(cit_srvic_fun(mnth,city)) for city in cities(images)])
end

# ---- transportation into :sink_city on month :month
function makequiver(monthindex::Int)
	println(months(images)[monthindex])
	map[:readshapefile](shapefile, shapename, linewidth=1.5)
	for sink_city in [:SanFrancisco, :Sacramento, :LosAngeles, :SanBernardino]
		month      = months(images)[monthindex] # pick a month [:mar, :sep]
		b_sply_inx = cit_srvic_fun(month, sink_city)
		suppl_x, suppl_y  = map[:__call__](bary_support[1, b_sply_inx][:], bary_support[2, b_sply_inx][:])
		sinkl_x, sinkl_y  = map[:__call__](citydf[find(citydf[:city].==sink_city)[1], :lon_x], citydf[find(citydf[:city].==sink_city)[1],:lat_y])
		map[:scatter](suppl_x, suppl_y,
			s=5000*pval[b_sply_inx],
			c=[0,0,0],
			alpha=0.15)
		function smallt(u,v)
			nrm = √(u.^2 + v.^2)
			return 30_000 .* u ./ nrm, 30_000 .* v ./ nrm
		end
		sub1, sub2 = smallt(sinkl_x .- suppl_x, sinkl_y .- suppl_y)
		quiver(suppl_x, suppl_y, sinkl_x .- suppl_x .- sub1, sinkl_y .- suppl_y .- sub2,
			angles="xy",
			scale_units="xy",
			scale = 1,
			width = 0.004,
			headwidth = 4.9,
			headlength =  3.0,
			headaxislength =  3.5,
			minshaft =  2.5, #1.5
			color = [0,0,0],
			alpha=0.9)
		map[:scatter](sinkl_x, sinkl_y,
			s=50,
			linewidths=2,
			marker="x",
			c=[1,0,0],
			alpha=1)
		axis("off")
	end
end

makequiver(4)
savefig(joinpath(figurepath, "figure1g.pdf"), dpi=300, bbox_inches="tight", transparent=true)

makequiver(6)
savefig(joinpath(figurepath, "figure1h.pdf"), dpi=300, bbox_inches="tight", transparent=true)
