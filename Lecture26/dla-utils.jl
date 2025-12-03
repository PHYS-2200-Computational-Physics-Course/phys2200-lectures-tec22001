"""
    i, j = random_circular_position(c, r)

Generate cartezian pixel coordinates of a random point on a circle
of radius r and center c.
"""
function random_circular_position(c, r)
    theta = 2 * pi * rand()
    sn, cs = sincos(theta)
    i = Int(floor(c[1] + r * cs))
    j = Int(floor(c[2] + r * sn))
    return i, j
end

"""
    xc, yc = center_of_mass(dish)

Calculate the position (in pixels) of the center of mass of the dendrite in the dish.
"""
function center_of_mass(dish)
    x = 0
    y = 0
    n = 0
    for index in findall(!iszero, dish)
        x += index[1]
        y += index[2]
        n += 1
    end
    return div(x, n), div(y, n)
end

"""
    xcn, ycn = update_center_of_mass(dish, cold, padd)

Re-calculate the position (in pixels) of the center of mass of the dendrite
in the dish after an additional molecule got attached.
"""
@inline function update_center_of_mass(dish, cold, padd)
    n = nzz(dish)
    x = cold[1] * n + padd[1]
    y = cold[2] * n + padd[2]
    return div(x, n+1), div(y, n+1)
end

"""
    rg = radius_of_gyration(dish, c)

Calculate the radius of gyration (in pixels) about point c (which is
expected to be pixel coordinates of the center of mass) as a measure
of spatial spread of the dish configuration.
"""
function radius_of_gyration(dish, c)
    r2 = 0.0
    n = 0
    for index in findall(!iszero, dish)
        r2 += (c[1] - index[1])^2 + (c[2] - index[2])^2
        n += 1
    end
    return sqrt(r2 / n)
end

# basevecs are determined by the geometry of the dish lattice. Those are directions
# to the neighbours as well as the directions of diffusion steps.

"""
    inew, jnew = step_random_direction(i, j, basevecs)

Move a particle with pixel coordinates (i, j) one step in a random direction.
Return the new coordinates of the particle.
"""
function step_random_direction(i, j, basevecs)
    direction = rand(1:length(basevecs))
    shift = basevecs[direction]
    return (i,j) .+ shift 
end

"""
    isfree = is_free(i, j, dish, basevecs)

Check if the particle (i, j) is 'touching' a fixed particle.
"""
function is_free(i, j, dish, k, basevecs)
    isfree = true
    for dir in basevecs
        if dish[(i,j) .+ dir] != 0
            # in contact with a fixed particle
            dish[i, j] = k
            isfree = false
            break
        end
    end
    return isfree
end
    
"""
    stepno, count = add_particle!(dish, basevecs, c, r, k)

Randomly place a free particle (labeled by index k) at the edge of
a circle of radius r and center, c. Let the particle diffuse until it get
attached to the growing structure. If the particle diffuses too far away
from c, discard it and inject a new particle instead.

Return the number of random-walk steps, stepno, and the number of attempts
to add a particle.
"""
function add_particle!(dish, basevecs, c, r, k)

    # Place a particle at the edge of a circle of radius r
    i, j = random_circular_position(c, r)

    # Check if the position (i, j) is not already occupied
    if dish[i, j] != 0
        error("A new particle $k is added to the position ($i, $j) that is already occupied.")
    end
    
    count = 1        # number of attempts to add a particle
    stepno = 0       # number of steps until attaching

    # Check if the particle is 'touching' a fixed one
    isfree = is_free(i, j, dish, k, basevecs)

    C = div.(size(dish) .- 1, 2)
    R2 = (C[1] - 10)^2
    
    # Random walk untill the particle is trapped or diffuses away
    while isfree

        stepno += 1
        i, j = step_random_direction(i, j, basevecs)

        # Check if the particle moved too far away from the center of the dish
        if (i - C[1])^2 + (j - C[2])^2 >= R2
            if stepno > 1
                # If the particle diffused away, replace it with a new random particle.
                # Reset stepno, increase count
                i, j = random_circular_position(c, r)
                # Check if the position (i, j) is not already occupied
                if dish[i, j] != 0
                    error("A new particle $k is added to the position ($i, $j) that is already occupied.")
                end
                stepno = 0
                count += 1
                continue
            else
                error("The cluster is too big for the current dish. Increase the dish size.")
            end
        end

        # Check if the particle is 'touching' a fixed one
        isfree = is_free(i, j, dish, k, basevecs)
    end

    return (i, j), stepno, count
end

"""
    dish, rads = dish_init(K, m)

K is the radius of the dish (in pixels), m is the number of molecules in the cluster,
dish(2K+1, 2k+1) stores the growing cluster, rads(m) is the storage for gyration radia
of the cluster.
"""
function dish_init(K, m)
    N = 2*K + 1      # dish size, in pixels, expected to be an odd number
    dish = spzeros(Int, N, N);
    dish[K, K] = 1   # seed in the center of the dish

    rads = zeros(m)  # preallocate storage for characteristic cluster sizes
    rads[1] = 0.0

    return dish, rads
end

"""
    show_cluster(dish)

Visualize the cluster `dish` obtained via the process of diffusion-limited aggregation.
The coloring of the clusterr corresponds to the "age" of aggregation.
"""
function show_cluster(dish)
    
    figure()
    imshow(dish, cmap="tab20") # cmap="cool", "viridis", "magma", ...
    # axis(false)
    m = nnz(dish) # The size of the cluster
    titl = @sprintf("%d particles", m)
    # title(titl, font="monospace")
    title(titl)

    return nothing
end

"""
    grow_cluster!(dish, basevecs, rads, m, animate, anistep)

Grow and optionally visualize a dendrite cluster. dish is the cluster storage, basevecs are
the allowed directions of the diffusion steps, rads(m) is the vector of gyration radia of 
the cluster, m is the number of molecules to add to the cluster, animate = (true|false) is
the visualization switch: if animate is true, every 'anistep' of the growth process is shown.
"""
function grow_cluster!(dish, basevecs, rads, m, animate, anistep)
    
    if animate
        fig = figure()
        img = imshow(dish)
        # axis(false)
    end

    K = div(size(dish)[1] - 1, 2)
    
    for i = 2:m
        c = center_of_mass(dish)
        rads[i] = radius_of_gyration(dish, c)
        r = Int64(round(rads[i]))
        R = 2*r + 15
        # Check if the agregate became too big for the dish
        if R >= K - 15
            break
        end
        add_particle!(dish, basevecs, c, R, i)
    
        if animate && i % anistep == 0
            img.set_array(dish)
            titl = @sprintf("%6d particles", i)
            # title(titl, font="monospace")
            title(titl)
            display(fig)
            sleep(0.01)
            IJulia.clear_output(true)
        end
    end

end
