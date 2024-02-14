
using JSON, PETScBinaryIO
import PETScBinaryIO: read_vec

export Read_Cell, Read_Vel, Read_DMDA_coords, Read_DMDA_data, get_dimensions_json, reshape_DMDA, ReadVec

function Read_Cell(step_dir::String; int_type=Int32, scalar_type=Float64)

    json_file  = step_dir*".dmda-pressure_dmda.json"  # center json (has dimensions)
    coord_file = step_dir*".dmda-cell.coords.vec"     # coordinates of P
    pres_file  = step_dir*".dmda-Xp"                  # values @ center
 
    # read coordinates of vertexes
    n, ndof = get_dimensions_json(json_file);
    coord_vec  = ReadVec(coord_file, int_type=int_type, scalar_type=scalar_type);  # read vec
    Coord = reshape_DMDA(coord_vec, n.+1, dim); # note that the size is onemore
    
    # Read data
    Data = Read_DMDA_data(pres_file, n, ndof; int_type=int_type, scalar_type=scalar_type);

    return Coord, Data
 end

 
function Read_Vel(step_dir::String; int_type=Int32, scalar_type=Float64)

   json_file  = step_dir*".dmda-velocity_dmda.json"  # velocity json (has dimensions)
   coord_file = step_dir*".dmda-velocity_dmda_coords.pbvec" # coordinates of V
   velocity_file = step_dir*".dmda-Xu" # values

   Coord, n = Read_DMDA_coords(json_file, coord_file; int_type=int_type, scalar_type=scalar_type);
   Vel = Read_DMDA_data(velocity_file, n, length(Coord); int_type=int_type, scalar_type=scalar_type);
    

    return Coord, Vel
end


"""
    Coord, n = Read_DMDA_coords(json_file::String, coord_file::String; int_type=Int32, scalar_type=Float64)

Reads a PETSc DMDA dataset as specified in `json_file` with coordinates file `coord_file`
"""
function Read_DMDA_coords(json_file::String, coord_file::String; int_type=Int32, scalar_type=Float64)
    
    # Retrieve names of files
    n, ndof = get_dimensions_json(json_file);

    # read coordinates
    coord_vec  = ReadVec(coord_file, int_type=int_type, scalar_type=scalar_type);  # read vec
    Coord = reshape_DMDA(coord_vec, n, dim);

    return Coord, n
end

"""

Reads a DMDA `file` with data of size `n` and `dof` degrees of freedom
"""
function Read_DMDA_data(file::String, n::Vector{Int}, dof::Int; int_type=Int32, scalar_type=Float64)
    
    # read velocity data
    vec  = ReadVec(file, int_type=int_type, scalar_type=scalar_type);  
    Data = reshape_DMDA(vec, n, dof);

    return Data
end



"""
returns the dimensions of the JSON file `json_file`
"""
function get_dimensions_json(json_file::String)
    json_data = JSON.parsefile(json_file) # read JSON file
    
    n = []
    ndof = 1
    if  any(contains.(keys(json_data),"DMDA"))
        dmda = json_data["DMDA"];   # info contained in JSON file
        
        # get some useful info
        dim  = dmda["dim"]
        ndof = dmda["ndof"]
        stencilWidth = dmda["stencilWidth"]

        # get the number of points in every direction
        n = [ dmda["directions"][i]["points"] for i = 1:dim ]   # compacter than a loop, same result

    elseif any(contains.(keys(json_data),"FVDA"))
        fvda = json_data["FVDA"];
        dim  = fvda["dim"]
        
    else

        error("unknown file format")
    end

    return n, ndof
end

"""
    A = reshape_DMDA(vec::Vector, n::Vector{Int}, dof::Int )

Reshape the vector `vec` into an vector `A` of length `dof` where each entry is an array of `size(n)`.   
A[1] can, for example, contain a 3D array with x-coordinates. In 3D, `A` would be a vector of `length(A)==3`
"""
function reshape_DMDA(vec::Vector, n::Vector{Int}, dof::Int )
    @assert length(vec) == prod(n)*dof

    A = fill(zeros(n...),dof)    #initialize arrays
    for i = 1:dof
        A[i] = reshape(vec[i:dof:end], n...);
    end
    return A
end


"""
    vec =  ReadVec(file::String; int_type=Int32, scalar_type=Float64)

This reads a PETSc vector from file
"""
function  ReadVec(file::String; int_type=Int32, scalar_type=Float64)
    io = open(file,"r")
    class_id = ntoh(read(io, int_type))
    @assert class_id == 1211214

    vec = PETScBinaryIO.read_prefix_vec(io, int_type, scalar_type)
    close(io)

    return vec
end


