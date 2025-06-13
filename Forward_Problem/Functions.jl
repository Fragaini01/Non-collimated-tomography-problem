function regrid_wall(wall , step, order = 1)
    sector_one_r = wall.r[2:6] 
    sector_one_z = wall.z[2:6]

    sector_two_r = reverse(wall.r[6:14])
    sector_two_z = reverse(wall.z[6:14])
    
    sector_three_r = wall.r[14:16]
    sector_three_z = wall.z[14:16]

    sector_four_r = wall.r[16:end-1]
    sector_four_z = wall.z[16:end-1]
    
    radii_1 = collect(sector_one_r[1]:step:sector_one_r[end])
    radii_2 = collect(sector_two_r[1]:step:sector_two_r[end])
    spline_3 = collect(minimum(sector_three_z):step:maximum(sector_three_z))

    spline_1 = Spline1D(sector_one_r, sector_one_z, k=order)(radii_1)
    spline_2 = Spline1D(sector_two_r, sector_two_z, k=order)(reverse(radii_2))
    radii_3 = collect(range(minimum(sector_three_r), stop=maximum(sector_three_r), length=length(spline_3)))

    radii = [radii_1; reverse(radii_2); reverse(radii_3) ; sector_four_r]
    zs = [spline_1; spline_2; reverse(spline_3) ; sector_four_z]

    return radii, zs

end

function line_two_points(p1, p2)
    p1 = [p1[1] , 0 , p1[3]]
    p2 = [p2[1] , 0 , p2[3]]
    # vector from p1 to p2
    v = (p2 - p1) / norm(p2 - p1)
    return v 
end

function angle_two_lines(m1, m2)
    angle = acos(round((dot(m1, m2) / (norm(m1) * norm(m2))), digits=8))
    cross_product = cross(m1, m2)
    return -sign(cross_product[2]) * angle
end

function target(wall_point, plasma_point, m2)
    m1 = line_two_points(wall_point, plasma_point)
    angle = angle_two_lines(m1, m2)
    return angle
end


function index_to_point(index, wall_r, wall_z)
    return [wall_r[index] , 0 , wall_z[index]]
end

function find_zero_bisection(voxel, detector_position, grid_r, grid_z, first_guess , verbose=true, tol=.1, max_iter=6000)
    if verbose
        println("Finding zero using bisection method...")
    end
    m_plasma_detector = line_two_points(detector_position, voxel)
    if verbose
        println("m_plasma_detector: ", m_plasma_detector)
        println("First guess: ", first_guess)
        println("low index: ", first_guess[1])
        println("high index: ", first_guess[2])
    end

    index_low = Int(first_guess[1])
    index_high = Int(first_guess[2])
    low = index_to_point(index_low, grid_r, grid_z)
    high = index_to_point(index_high, grid_r, grid_z)
    while(target(low , voxel, m_plasma_detector) * target(high , voxel, m_plasma_detector)) > 0
        if verbose
            println("Target function at low: ", target(low , voxel, m_plasma_detector))
            println("Target function at high: ", target(high , voxel, m_plasma_detector))
        end
        index_high = rand(1:length(grid_r))
        high = index_to_point(index_high, grid_r, grid_z)
    end

    for i in 1:max_iter
        mid_index = Int(floor(index_low + (index_high - index_low) / 2))
        if mid_index > length(grid_r)
            mid_index = mid_index - length(grid_r)
        elseif mid_index < 1
            mid_index = mid_index + length(grid_r)
        end
        mid = index_to_point(mid_index, grid_r, grid_z)
        if target(low , voxel, m_plasma_detector) * target(mid , voxel, m_plasma_detector) < 0
            index_high = mid_index
            high = mid
        else
            index_low = mid_index
            low = mid
        end

        if abs(target(low , voxel, m_plasma_detector)) < tol
            if verbose
                println("Bisection method converged after ", i, " iterations.")
            end
            return low, high, index_low, index_high
        end
        if abs(target(high , voxel, m_plasma_detector)) < tol
            if verbose
                println("Bisection method converged after ", i, " iterations.")
            end
            return high, low, index_high, index_low
        end

        if abs(index_high - index_low) < 2
            if verbose
                println("Bisection method converged after ", i, " iterations.")
            end
            index_high = rand(1:length(grid_r))
            return high, low, index_high, rand(1:length(grid_r))
        end

        if i % 1 == 0 && verbose
            println("Iteration: ", i, " Low: ", low, " High: ", high, "index_low: ", index_low, " index_high: ", index_high)
        end
    end


    return low, high, index_low, index_high
end

function a_ij(voxel, detector_position, radii, zs, first_wall_size, second_wall_size, cross_section_1, cross_section_2, first_guess)
    # R z plan
    points_in_the_wall = find_zero_bisection(voxel, detector_position, radii, zs, first_guess,false)
    linear_wall = line_two_points(points_in_the_wall[1], points_in_the_wall[2])
    line_detector_voxel = line_two_points(detector_position, voxel)
    angle_wall = angle_two_lines(line_detector_voxel, linear_wall)


    # Total distance
    x_1 = abs(first_wall_size / sin(angle_wall))
    x_2 = abs(second_wall_size / sin(angle_wall))
    x_voxel = voxel[1] * cos(voxel[2])
    y_voxel = voxel[1] * sin(voxel[2])
    x_detector = detector_position[1] * cos(detector_position[2])
    y_detector = detector_position[1] * sin(detector_position[2])

    A = sqrt((x_voxel - x_detector)^2 + (y_voxel - y_detector)^2)
    costheta = abs(-voxel[1]^2 + A^2 + detector_position[1]^2) / (2 * A * detector_position[1])
    costheta = abs(detector_position[1]^2 + voxel[1]^2 - A^2) / (2 * detector_position[1] * voxel[1])
    if costheta > 1.0 || costheta < 0
        println("Warning: cos(theta) out of bounds: ", costheta)
        costheta = clamp(costheta, -1.0, 1.0)
    end
    if -cross_section_1 * x_1 - cross_section_2 * x_2 > 0 
        println("Warning: Negative cross section: ", -cross_section_1 * x_1 - cross_section_2 * x_2)
    end
    if x_1 <0 || x_2 < 0
        println("Warning: Negative distance: x1: ", x_1, " x2: ", x_2)
    end
    
    d = sqrt((detector_position[1]- points_in_the_wall[1][1])^2 + (detector_position[3]- points_in_the_wall[1][3])^2)
    n = first_wall_size/ second_wall_size
    x_1 = d * n / (n+1)
    x_2 = d - x_1
    if x_1 < 0 || x_2 < 0
        println("Warning: Negative distance: x1: ", x_1, " x2: ", x_2)
    end
    total_attenuation = exp((-cross_section_1 * x_1 - cross_section_2 * x_2) * costheta)


    return total_attenuation, points_in_the_wall[3], points_in_the_wall[4], points_in_the_wall[1]
end

function plasma_gridder(dr, dphi, dz, wall_r, wall_z)
    # Create a grid of points in the plasma
    r = collect(minimum(wall_r):dr:maximum(wall_r)) # skip dr, will return on it later
    z = collect(minimum(wall_z):dz:maximum(wall_z))
    phi_grid = collect(range(-180,stop=180,length=360)).*(pi/180.0) # Toroidal angle in degrees
    phi_grid = phi_grid[1:end-1] # Remove phi=180 (since it is equivalent to phi=-180)
    total_grid = []
    grid = []
    for i in 1:length(r)
        for j in 1:length(z)
            push!(grid, [r[i], 0, z[j]])
        end
    end

    #grid = grid[[in_boundary(wall , point[1] , point[3]) for point in grid]]
    total_grid = []
    for phi in phi_grid
        for i in grid 
            push!(total_grid, [i[1], phi, i[3]])
        end
    end
    # Sort total_grid by the second column (phi)
    total_grid = total_grid[sortperm([point[2] for point in total_grid]), :]

    return total_grid
end

function plasma_gridder_henrik(nR , nz, wall)
    nphi = 360
    R_grid = collect(range(minimum(wall.r),stop=maximum(wall.r),length=nR)) # Major radius in meters
    phi_grid = collect(range(-180,stop=180,length=nphi)).*(pi/180.0) # Toroidal angle in degrees
    phi_grid = phi_grid[1:end-1] # Remove phi=180 (since it is equivalent to phi=-180)
    z_grid = collect(range(minimum(wall.z),stop=maximum(wall.z),length=nz)) # Vertical coordinate in meters
    grid = []

    RZ_grid = [[r, z] for r in R_grid, z in z_grid]
    RZ_grid = reshape(RZ_grid, :)

    RZ_grid = [pt for pt in RZ_grid if in_boundary(wall, pt[1], pt[2])]

    n = length(RZ_grid)
    m = length(phi_grid)
    Matrix = Array{Float64,2}(undef, n*m, 3)
    for (j, phi) in enumerate(phi_grid)
        for (i, pt) in enumerate(RZ_grid)
            Matrix[(j-1)*n + i, :] = [pt[1], phi, pt[2]]
        end
    end
    # Sort Matrix by the second column (phi)
    Matrix = Matrix[sortperm(Matrix[:,2]), :]
    return Matrix
end

function point_to_grid(point, grid_plasma) 
    point_plasma = findminimum(grid_plasma, x -> norm(x - point))
    index_plasma = findfirst(x -> grid_plasma[x] == point_plasma, grid_plasma)
    return index_plasma
end

function place_detector(wall_r, wall_z , first_wall_size, second_wall_size, index_wall)
    first_point_wall = [wall_r[index_wall], 0, wall_z[index_wall]]
    if first_point_wall[3] < -2
        println("Warning: First point wall is below -2 m: ", first_point_wall)
    end
    # find the angular coefficient of the wall in that point considering the point after and before
    m_wall = line_two_points([wall_r[index_wall - 1], 0, wall_z[index_wall - 1]], [wall_r[index_wall + 1], 0, wall_z[index_wall + 1]])
    m_ortogonal = [m_wall[3], 0, -m_wall[1]]
    # find line from the point in the wall with -1/m coefficient, evalueted in first_wall_size + second_wall_size
    # find the detector position
    dr = (first_wall_size + second_wall_size) * m_ortogonal[1]
    dz = (first_wall_size + second_wall_size) * m_ortogonal[3]
    r_detector = first_point_wall[1] + dr
    z_detector = first_point_wall[3] + dz
    
    return [r_detector, 0, z_detector]
end

function detectors_on_file(detectors, filename = "detectors.csv")
    filename = "../my_code/"*filename
    open(filename, "w") do io
        println(io, "r,phi,z")
        for det in detectors
            println(io, "$(det[1]),$(det[2]),$(det[3])")
        end
    end
end

function detectors_positions(wall_r, wall_z , first_wall_size, second_wall_size, number_of_detectors)
    detectors = []

    #Find the step 
    interesting_r = wall_r[wall_z .> -2]
    interesting_z = wall_z[wall_z .> -2]
    total_length = [sqrt((interesting_r[i] - interesting_r[i+1])^2 + (interesting_z[i] - interesting_z[i+1])^2) for i in 1:length(interesting_r)-1]
    total_length = sum(total_length)
    step = total_length / number_of_detectors

    i = 1
    while i < length(interesting_r) && length(detectors) < number_of_detectors
        direction_of_the_wall = line_two_points([interesting_r[i], 0, interesting_z[i]], [interesting_r[i+1], 0, interesting_z[i+1]])
        next_z = interesting_z[i] + step * direction_of_the_wall[3] / norm(direction_of_the_wall)
        next_r = interesting_r[i] + step * direction_of_the_wall[1] / norm(direction_of_the_wall)

        # Find the index that is closest to the values next_r and next_z
        distances = [sqrt((interesting_r[j] - next_r)^2 + (interesting_z[j] - next_z)^2) for j in i:min(i+2000,length(interesting_r)-1)]
        min_dist, min_idx = findmin(distances)
        next_index = min_idx + i 
        @assert next_index > i "Next index must be greater than current index"
        i = next_index
        detector = place_detector(interesting_r, interesting_z , first_wall_size, second_wall_size, i)
        push!(detectors, detector)
        
    end
    return detectors
end

function mask_grid_by_slice(plasma_grid, points_seen_detector, dr, dz) 
    mask = fill(0, length(plasma_grid))
    slice_size = Int(floor(length(plasma_grid) / 360))
    max_phi = maximum([point[2] for point in points_seen_detector])
    min_phi = minimum([point[2] for point in points_seen_detector])

    first_index_voxel = findfirst(x -> isapprox(x[2], min_phi; atol=1e-8), plasma_grid)
    last_index_voxel = findlast(x -> isapprox(x[2], max_phi; atol=1e-8), plasma_grid)
    
    interesting_voxels = plasma_grid[first_index_voxel:last_index_voxel]
    println("Number of voxels to check : " , length(interesting_voxels))
    max_dist = sqrt(dr^2 + dz^2) / 2

    voxel_before = [0, -99  , 0 ]

    slice = [x for x in points_seen_detector if isapprox(x[2], interesting_voxels[1][2]; atol=1e-8)]
    # First slices
    for voxel in interesting_voxels
        idx = findfirst(x -> x == voxel, interesting_voxels) + first_index_voxel - 1

        if !isapprox(voxel_before[2], voxel[2]; atol=1e-8)
            slice = [x for x in points_seen_detector if isapprox(x[2], voxel[2]; atol=1e-8)]
        end
        for x in slice
            if sqrt((x[1] - voxel[1])^2 + (x[3] - voxel[3])^2) < max_dist
                mask[idx] = 1
                break
            end
        end
        voxel_before = voxel
    end
    return mask
end

function mask_grid_interpolate(voxels_seen, plasma_grid_henrik, plasma_grid, regridded_wall)
    
    slice = [[p[1],p[2] , p[3]]  for p in plasma_grid]
    first_slice = plasma_grid[1:Int(length(plasma_grid) / 359)]


    slice_matrix = reduce(hcat, slice)'
    xi = slice_matrix
    points = (plasma_grid_henrik[:, 1], plasma_grid_henrik[:, 2], plasma_grid_henrik[:, 3])

    voxels_seen = sort(voxels_seen, by = x -> x[2])
    voxels_seen_matrix = reduce(hcat, voxels_seen)'

    minphi = minimum(voxels_seen_matrix[:, 2])
    maxphi = maximum(voxels_seen_matrix[:, 2])

    plasma_grid_henrik_used = plasma_grid_henrik[(plasma_grid_henrik[:, 2] .>= minphi) .& (plasma_grid_henrik[:, 2] .<= maxphi), :]
    # Trova l'indice del primo elemento di plasma_grid_henrik_used all'interno di plasma_grid_henrik
    first_point = plasma_grid_henrik_used[1, :]
    min_index = findfirst(row -> all(row .== first_point), eachrow(plasma_grid_henrik))
    tol = 1e-1

    indices_in_full = Int[]
    # Confronta le due matrici: plasma_grid_henrik_used e voxels_seen_matrix

    # Build a KDTree for fast nearest neighbor search
    tree = KDTree(voxels_seen_matrix[:, 1:3]')  # Only use the first two columns (x, y)
    inds = inrange(tree, plasma_grid_henrik_used[:, 1:3]', tol)

    for (i, found) in enumerate(inds)
        if !isempty(found)
            push!(indices_in_full, i)
        end
    end

    samples = zeros(Int, size(plasma_grid_henrik, 1))
    if !isempty(indices_in_full)
        #println("Indices in full: ", indices_in_full)
        samples[indices_in_full .+ (min_index .- 1)] .= 1
    end

    samples = collect(samples)


    # Convert points and xi to the format expected by scipy.interpolate.NearestNDInterpolator
    nearest_interp = interpolate.NearestNDInterpolator(points, samples)
    interpolated = nearest_interp(xi)

    # Trova gli indici validi (non NaN)
    valid_idx = .!isnan.(interpolated)

    # Arrotonda solo i valori validi
    interpolated_int = Array{Int}(undef, length(interpolated))
    interpolated_int[valid_idx] = round.(Int, interpolated[valid_idx])
    # Cambia il segno dei punti che non corrispondono a punti visti (cioè fuori dal wall)
    for i in 1:length(first_slice)
        if !in_boundary(regridded_wall, first_slice[i][1], first_slice[i][3])
            for j in 0:358
                interpolated_int[i + j * slice_size] = 0
            end
        end
    end


    # check the first column by comparing with the first column in voxel seen 

    
    return interpolated_int
end


function a_to_W(a , plasma_grid, detectors, out_file , verbose=true)
    if verbose
        println("Converting a to W...")
    end
    W = a'  # [nDetectors, nGridPoints]

    Ed_array = collect(1:length(detectors))  # O le coordinate reali dei detectors
    Ed_array_units = "-"  # O l'unità corretta

    slice_size = Int(length(plasma_grid) / 359)
    first_slice = plasma_grid[1:slice_size]  
    # Coordinate dei punti della griglia
    println(typeof(plasma_grid))
    println(size(plasma_grid))
    D1_array = unique([p[1] for p in first_slice])  # R
    D1_array_units = "m"
    D2_array = unique([p[3] for p in first_slice])  # z
    D2_array_units = "m"
    
    println(typeof(D1_array))
    println(typeof(D2_array))

    output = out_file * ".jld2"
    if verbose
        println("Saving W to file: ", output)
    end
    W = reshape(W , (length(detectors) , length(D1_array) , length(D2_array)))

    @save output W Ed_array Ed_array_units D1_array D1_array_units D2_array D2_array_units

end

function compute_measurements(intensities , in_file,  out_file, number_of_detectors = 99)

    # read W from file
    in_file = in_file * ".jld2"
    if isfile(in_file )
        @load (in_file ) W Ed_array Ed_array_units D1_array D1_array_units D2_array D2_array_units
    else
        println("File not found: ", in_file )
        return
    end
    if length(W) == 0
        println("W is empty, cannot compute measurements.")
        return
    end
    if length(intensities) == 0
        println("Intensities is empty, cannot compute measurements.")
        return
    end
    W = reshape(W, length(Ed_array), length(D1_array) * length(D2_array))
    
    columns_to_keep = round.(Int, range(1, length=number_of_detectors, stop=size(W,1)))
    W = W[columns_to_keep, :]
    Ed_array = Ed_array[columns_to_keep]
    
    println("W loaded from file: ", in_file)
    println("Number of detectors: ", size(W, 1))
    println("Number of grid points: ", size(W, 2))
    println("Number of intensities: ", length(intensities))
    
    # Compute the measurements
    measurements = W * intensities
    # Save the measurements to a file
    output_file = out_file * ".jld2"
    e = randn(size(measurements))
    e = e / norm(e)
    err = ones(length(measurements)) 
    S_units = "-"
    err_units = "-"
    Ed_array = collect(1:length(measurements)) 
    Ed_array_units = "-" 
    println(size(measurements))
    println(size(W))
    println(size(intensities))
    # Ed_array e Ed_array_units già caricati dal file W

    #Adding noise 
    measurements = measurements+0.01*norm(measurements)*e
    output_file = out_file * "_S.jld2"
    if number_of_detectors != 99
        output_file = out_file * "_S_$(number_of_detectors).jld2"
        base = split(in_file, r"\.(?=[^\.]*$)")[1]
        W = reshape(W , (number_of_detectors , length(D1_array) , length(D2_array)))

        @save base * "_$(number_of_detectors).jld2" W Ed_array Ed_array_units D1_array D1_array_units D2_array D2_array_units
        println("new W saved to ", base * "_$(number_of_detectors).jld2")
    end
    @save output_file S=measurements S_units=S_units err=err err_units=err_units Ed_array Ed_array_units
    println("Measurements saved to ", output_file)


    true_values_file = out_file * "_true_values.jld2"
    @save true_values_file values = intensities values_units = "-" D1_array D1_array_units D2_array D2_array_units 

end