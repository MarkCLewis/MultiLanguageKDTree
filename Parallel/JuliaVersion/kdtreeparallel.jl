module kdtreeparallel
# KD-Tree Sequential, single struct for nodes
# Audrey & Clarissa

# Constants
const MAX_PARTS = 7
const THETA = 0.3

mutable struct Body
    p::Array{Float64} #p[0] = x, p[1] = y...
    v::Array{Float64} #v[0] = vx, v[1] = vy...
    m::Float64
end


abstract type KDTree end

struct KDLeaf <: KDTree
    # for leaves
    particles::Array{Int64}
end

struct KDInternal <: KDTree
    # for internal nodes
    split_dim::Int64
    split_val::Float64
    m::Float64
    cm::Array{Float64}
    size::Float64
    left::KDTree
    right::KDTree
end

function build_tree(indices::Array{Int64}, start::Int64, ending::Int64, system::Array{Body})::KDTree
    np = ending - start
    #println(np, " ", start, " ", ending)
    if np <= MAX_PARTS
        node = KDLeaf(zeros(Int64, np))
        for i in 0:np-1
            node.particles[i+1] = indices[start + i]
        end
        node
    else
	#println("Build internal")
        minp = [1e100, 1e100, 1e100]
        maxp = [-1e100, -1e100, -1e100]
        m = 0.0
        cm = [0.0, 0.0, 0.0]
        for i in start:ending-1
            m += system[indices[i]].m
            cm += system[indices[i]].m * system[indices[i]].p
            minp = min.(minp, system[indices[i]].p)
            maxp = max.(maxp, system[indices[i]].p)
        end
        cm /= m
        split_dim = 1
        if maxp[2] - minp[2] > maxp[split_dim] - minp[split_dim]
            split_dim = 2
        end
        if maxp[3] - minp[3] > maxp[split_dim] - minp[split_dim]
            split_dim = 3
        end
        size = maxp[split_dim] - minp[split_dim]
        # partition time
        mid::Int64 = div((start + ending), 2) # int division mid bs
        s = start
        e = ending
        while s + 1 < e
            pivot = rand(s:e-1)
            swapTmp = indices[s]
            indices[s] = indices[pivot]
            indices[pivot] = swapTmp
            low = s+1
            high = e-1
            while low <= high
                if system[indices[low]].p[split_dim] < system[indices[s]].p[split_dim]
                    low += 1
                else
                    swapTmp2 = indices[low]
                    indices[low] = indices[high]
                    indices[high] = swapTmp2
                    high-= 1
                end
            end
            swapTmp3 = indices[s]
            indices[s] = indices[high]
            indices[high] = swapTmp3
            if high < mid
                s = high + 1
            elseif high > mid
                e = high
            else
                s = e
            end
        end
        split_val = system[indices[mid]].p[split_dim]
        #recursion on kids
        left = build_tree(indices, start, mid, system)
        right = build_tree(indices, mid, ending, system)
        
        node = KDInternal(split_dim, split_val, m, cm, size, left, right)
        
        node
    end
end

    function calc_pp_accel(system::Vector{Body}, i::Int64, j::Int64, acc::Vector{Float64})
        d = system[i].p - system[j].p
        dist = sqrt(sum(d.^2))
        magi = -system[j].m / (dist^3) 
        acc += d * magi
    end

    function accel_recur(cur_node::KDLeaf, p::Int64, system::Vector{Body}, acc::Vector{Float64})
        for i in cur_node.particles
            if i != p
                calc_pp_accel(system, p, i, acc)
            end
        end
    end
    function accel_recur(cur_node::KDInternal, p::Int64, system::Vector{Body}, acc::Vector{Float64})
        d = system[p].p - cur_node.cm
        dist_sqr = sum(d.^2)
        if cur_node.size * cur_node.size < THETA^2 * dist_sqr
            dist = sqrt(dist_sqr)
            magi = -cur_node.m / (dist_sqr * dist)
            acc += d * magi
        else
            accel_recur(cur_node.left, p, system, acc)
            accel_recur(cur_node.right, p, system, acc)
        end
    end

    function calc_accel(p::Int64, tree:: KDTree, system::Vector{Body}, acc::Vector{Float64})
        accel_recur(tree, p, system, acc)
    end

    function print_tree(step::Int64, tree::KDTree, system::Array{Body})
        function print_node(n::KDLeaf, file::IO)
            println(file, "L $(n.num_parts)")
                        for i in n.particles
                            println(i)
                            println(file, "$(system[p].p[1]) $(system[p].p[2]) $(system[p].p[3])")
                        end
        end

        function print_node(n::KDInternal, file::IO)
            println(file, "I $(n.split_dim) $(n.split_val)")
            print_node(n.left, file)
            print_node(n.right, file)
        end

        fname = "tree$step.txt"
        try
            open(fname, "w") do file
                print_node(tree, file)
            end
        catch ex
            println(ex)
        end
    end


    function simple_sim(system::Vector{Body}, dt::Float64, steps::Int64)
        nb::Int64 = length(system)
        acc = zeros(nb, 3)
        indices = collect(1:nb)
    
        for step in 1:steps
            tree = build_tree(indices, 1, nb+1, system)
            Threads.@threads for i in 1:nb
                calc_accel(i, tree, system, acc[i, :])
            end
            Threads.@threads for i in 1:nb
                system[i].v += dt * acc[i,:]
                system[i].p += dt * system[i].v
                acc[i, :] .= 0.0
            end
            # print_tree(step, tree, system)
        end
    end

    function circular_orbits(n::Int64)::Vector{Body}
        first = Body([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0)
        bods = [first]
        for i in 1:n
            d = .1 + (i * 5.0 / n)
            v = sqrt(1.0 / d)
            theta = rand(Float64)*2*pi
            x = d * cos(theta)
            y = d * sin(theta)
            vx = -v * sin(theta)
            vy = v * cos(theta)
            temp = Body([x, y, 0.0], [vx, vy, 0] , 1.0e-7)
            push!(bods, temp)
        end
        bods
        #print(bods)
    end
    # to run:
    # julia main.jl #particles #steps
    if !isinteractive()
        steps = parse(Int64, ARGS[1])
        n = parse(Int64, ARGS[2])
        dt = 1e-3
        system = circular_orbits(n)
        simple_sim(system, dt, steps)
    end

end
