


init = list(
	length = 4, #length of region in which to generate particles
	width = 4,
	height = 1e-8,

	density = 0.1, #particle density; can also be function. (0,1].
	init_seed = 1 #set.seed() for particle distribution
)

rules = list(
	structure_seed = 1, #set.seed() for structure generation; preferred over init_seed
	dimension = 2, #links per node; does the model use links, triangles, or tetrahedra?

	link_range = 0.18,
	mobility = 0.04,
	contraction = 1.0,
	cohesion = 0,

	branching = 1
)
#deal with bias, genesis, thickness later
#thickness should go in the "visualization suite"

#you can also pass heatmaps into any of these with the format function(pos) [pos=c(x,y,z)]





parameter = function(rule) ifelse(is.function(rule),rule,function(x,y,z) rule)


point_colnames = c(
	'x', #x position of point
	'y',
	'z',
	'loose?' #does the point belong to a node yet?
)
node_colnames = c(
	'p1', #pid of first connected point
	'p2',
	'p3',
	'timeout',
	'generation', #starts at 1
	'parent', #nid of parent node
	'child1', #nid of first child node
	'child2',
	'child3',
	'active',
	'x', #average x position of points
	'y',
	'z'
)


#to enable 1D stuff, we want to let active go above one and decrease each time it branches


Points = function(init){
	set.seed(init$seed)
	n_points = 2048 * ceiling(init$length) * ceiling(init$width) * ceiling(init$height)
	pts = matrix(0,n_points,length(point_colnames)) #giving a matrix real colnames slows it down signficantly

	#distribute points randomly
	X = runif(n_points,0,init$length)
	Y = runif(n_points,0,init$width)
	Z = runif(n_points,0,init$height)

	#enforce density by thinning
	pid = 1
	for (p in 1:n_points){
			pos = c(X[p],Y[p],Z[p])
			if (parameter(init$density)(pos) > runif(1)) {
				pts[pid,] = c(pos,0)
				pid = pid + 1
			}
	}

	#return the remaining points
	pts[1:(pid-1),]
}

Nodes = function(rules){
	n_nodes = 65536 #dynamically updated; the size of this matrix increases by 65536 if it gets close to filling
	nodes = matrix(0,n_nodes,length(node_colnames))

	nodes
}

enchunk = function(pos) paste0(pos, collapse=',')
Chunk = function(pts,size){ #fuzz: include some points in multiple if they sit on the border, in case they drift
	size=1
	warning('chunking nonfunctional atm')
	chunks = list(size=size) #3D list
	for (p in 1:nrow(pts)){
		chunkpos = enchunk(floor(pts[p,1:3]/size))
		if (is.null(chunks[[chunkpos]])) chunks[[chunkpos]] = p
		else chunks[[chunkpos]] = c(chunks[[chunkpos]], p)
	}
	chunks
}

Sim = function(init,rules){
	sim = list(init=init,rules=rules)
	sim$points = Points(init) #uses init$init_seed for distribution
	sim$nodes = Nodes(rules)
	set.seed(rules$seed)

	sim$chunks = Chunk(sim$points, rules$link_range * 2)

	#first link
	sim$points = rbind(c(init$length/2,init$width/2,init$height/2,1),sim$points) #point in the center
	sim$nodes[1,] = c(rep(1,rules$dimension),rep(0,3-rules$dimension),1,1,0,0,0,0,1,sim$points[1,1:3]) #node connected to only that point
	sim$nid = 2 #first open spot in sim$nodes

	sim
}


grow = function(sim, nid){
	node = sim$nodes[nid,]

	# ---- Get Nearby Points ----
	pos = node[11:13]
	chunkpos = floor(pos / sim$chunks$size)
	subchunkpos = sign((pos / sim$chunks$size - chunkpos) - 0.5)

	#get all point pids from a 2x2x2 cube of chunks around the node
	close_pids = c(
		sim$chunks[[enchunk(chunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(0,0,1) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(0,1,0) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(0,1,1) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(1,0,0) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(1,0,1) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(1,1,0) * subchunkpos)]],
		sim$chunks[[enchunk(chunkpos + c(1,1,1) * subchunkpos)]]
	)

	# ---- Calculate Distance ----
	pts = sim$points[close_pids,1:3,drop=FALSE]
	diffs = t(t(pts) - pos) #casting hack; might be way slower
	dists = sqrt(rowSums((diffs)^2))
	dists[sim$points[close_pids,4] != 0] = dists[sim$points[close_pids,4] != 0] / sim$rules$cohesion #apply cohesion

	# ---- Create New Nodes ----
	closest_point = which.min(dists)
	if (length(closest_point) && dists[closest_point] < sim$rules$link_range) {#if a closest valid point exists to link to
		pid = close_pids[closest_point]
		sim$points[pid,4] = 1 #no longer free

		for (p in 1:sim$rules$dimension){ #add number of nodes equal to dimension
			sim$nodes[nid, 7 + p] = sim$nid #add child
			pts = c(pid, node[(1:3)[-p]]) #exclude one of the old points
			sim$nodes[sim$nid,] = c(pts, 10, node[5] + 1, nid, 0, 0, 0, 1, 0, 0, 0) #pos will get updated during contraction #ADD MOBILITY HERE
			sim$nid = sim$nid + 1
		}
	}

	sim
}

contract = function(sim){
	active_nids = which(sim$nodes[,4] > 0)
	n_active = length(active_nids)
	if (n_active == 0) return(sim) #exit if none active
	sim$nodes[active_nids, 4] = sim$nodes[active_nids, 4] - 1 #decrement timer

	#for each node: xyz xyz xyz xyz for points 1-4
	glob_points = matrix(NA,n_active,9)
	pts = list()
	for (p in 1:3) if (sum(sim$nodes[active_nids,p]) > 0) glob_points[, (p-1)*3 + 1:3] = sim$points[sim$nodes[active_nids,p],1:3]

	avg_pos = cbind(
		rowMeans(glob_points[,c(1,4,7),drop=FALSE],na.rm=TRUE),
		rowMeans(glob_points[,c(2,5,8),drop=FALSE],na.rm=TRUE),
		rowMeans(glob_points[,c(3,6,9),drop=FALSE],na.rm=TRUE)
	)

	sim$nodes[active_nids,11:13] = avg_pos #update average pos of nodes

	glob_avg = cbind(avg_pos,avg_pos,avg_pos) #duplicate to overlay point_glob

	glob_velocity = (glob_avg - glob_points) * 0.01 #CONTRACTION

	for (n in 1:nrow(glob_points)){
		for (p in 1:3){
			sim$points[sim$nodes[n,p],1:3] = sim$points[sim$nodes[n,p],1:3] + glob_velocity[n, (p-1)*3 + 1:3]
		}
	}

	sim
}


#need some renaming of active vs live
tick = function(sim){
	live_nids = which(sim$nodes[,10] == 1)
	n_live = length(live_nids)

	# ---- Grow Nodes ----
	to_grow = live_nids[0.05 > runif(n_live)] #feeling weird about this. would rather not be tick-based and instead always fire randomly, but then contraction speed is an issue...
	if (length(to_grow)) for (n in 1:length(to_grow)) sim = grow(sim, to_grow[n])
	sim$nodes[to_grow,10] = 0 #deactivate

	# ---- Contract Nodes ----
	sim = contract(sim)

	sim
}


view_network = function(sim){
	real_nids = which(sim$nodes[,1] > 0)
	n_real = length(real_nids)
	plot(NULL,xlim=c(0,sim$init$length),ylim=c(0,sim$init$width))
	for (n in 1:n_real){
		node = sim$nodes[n,]
		pnode = sim$nodes[node[6],]
		cnode = sim$node[node[8],]
		#lines(c(node[11],pnode[11]),c(node[12],pnode[12]))

		polygon(c(sim$points[node[1],1],sim$points[node[2],1],sim$points[node[3],1]),c(sim$points[node[1],2],sim$points[node[2],2],sim$points[node[3],1]))
		#polygon(c(sim$points[cnode[1],1],sim$points[node[1],1],sim$points[node[2],1]),c(sim$points[cnode[1],2],sim$points[node[1],2],sim$points[node[2],2]),col='black')
	}
}

sim=Sim(init,rules)
for(i in 1:1000){
	cat('\r',i,'     ')
	sim=tick(sim)
}
view_network(sim)

