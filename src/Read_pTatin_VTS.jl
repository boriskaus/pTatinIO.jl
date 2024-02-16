# This reads VTS files, created by ptatin3d-pyviztools, back into julia
using ReadVTK, GeophysicalModelGenerator, WriteVTK, LightXML

export pTatin_CartData, Compress_pTatin_Simulation, Compress_pTatin_VTS_File, Read_pTatin_VTS_File


"""
    Compress_pTatin_Simulation(pvd_filename::String; Dir=pwd())

This compresses all the files of a pTatin simulation.
"""
function Compress_pTatin_Simulation(pvd_filename::String; Dir=pwd())
    if pvd_filename[end-2:end]!="pvd"
        error("Filename should end on *.pvd; currently is $pvd_filename")
    end
    CurDir = pwd()
    cd(Dir)

    # Load the *.pvd simulation ---
    # Normally we would use PVDFile for that (from the ReadVTK package), but for some reason, 
    # the ptatin3d-pyviztools writes a too old version of this fileformat. 
    # It's a quite simple fileformat, so we simply mimic this here
    xml_file_contents = read(pvd_filename, String)
    
    # Open file and ensure that it is a valid VTK file
    xml_file = LightXML.parse_string(xml_file_contents)
    root = LightXML.root(xml_file)
    @assert LightXML.name(root) == "VTKFile"

    # Extract attributes (use `required=true` to fail fast & hard in case of unexpected content)
    file_type = attribute(root, "type", required = true)
   
    # Ensure matching file types
    if file_type != "Collection"
        error("Unsupported PVD file type: ", file_type)
    end

    version = VersionNumber(attribute(root, "version", required = true))
    
    pieces = root[file_type][1]["DataSet"]
    n_pieces = length(pieces)
    vtk_filenames = Vector{String}(undef, n_pieces)
    directories = Vector{String}(undef, n_pieces)
    timesteps = Vector{Float64}(undef, n_pieces)
    for i in 1:n_pieces
        file_dir = attribute(pieces[i], "file", required = true)
        vtk_filenames[i] = file_dir
        directories[i] = dirname(file_dir)
        timesteps[i] = parse(Float64, attribute(pieces[i], "timestep", required = true))
    end
    pvd = PVDFile(pvd_filename, file_type, vtk_filenames, directories, timesteps)
    # -----------------------------


    # Now loop through the files and a) compress each of them, b) add them to a new PVD file
    vtk_filenames_compressed =  Vector{WriteVTK.DatasetFile}(undef, n_pieces)
    pvd_filename_compressed = pvd_filename[1:end-4]*"_compressed.pvd"
    pvd_compressed = paraview_collection(pvd_filename_compressed)
    for i in 1:n_pieces
        vtk = Compress_pTatin_VTS_File(pvd.vtk_filenames[i]);
        pvd_compressed[pvd.timesteps[i]] = vtk
    end
    
    # and save the new PVD file
    vtk_save(pvd_compressed)

    cd(CurDir)

    return pvd_filename_compressed
end

"""
    cart_data = pTatin_CartData(FileName::String; Dir::String=pwd(), to_km=true)

This opens a pTatin `*.vts` file and imports it as GeophysicalModelGenerator `CartData` dataset.
Note that we interpolate vertex data (usually velocity) to cell centers. `to_km` indicates whether we scale 
the length from meters to kilometers. 
"""
function pTatin_CartData(FileName::String; Dir::String=pwd(), to_km=true)

    # read file
    coord, data_points, data_cell = Read_pTatin_VTS_File(FileName; DirName_base=Dir);

    # interpolate vertex -> cell data
    coord_cell = Vertex2Cell(coord)
    data_cell_interp = Vertex2Cell(data_points);

    # Create CartData struct
    data_cell_merged = merge(data_cell, data_cell_interp[1])  # merge the datasets into one
    if to_km
        scale_length = 1e-3;
    else
        scale_length = 1;
    end

    cart_data =  CartData(coord_cell[1]*scale_length,coord_cell[2]*scale_length,coord_cell[3]*scale_length,data_cell_merged);

    return  cart_data
end

"""
    coords, data_points, data_cell = Read_pTatin_VTS_File(FileName::String; DirName_base::String=pwd())

This reads a pTatin *.vts file and returns the coordinates `coords`, points data (at vertexes) and cell data.

"""
function Read_pTatin_VTS_File(FileName::String; DirName_base::String=pwd())
    CurDir = pwd();

    DirName, File = split_pathname(DirName_base, FileName)
    
    # read data from parallel rectilinear grid
    cd(DirName)
    vts = VTKFile(File)
    cd(CurDir)

    # Fields stored in this data file:    
    data_point = get_point_data(vts)
    data_cell = get_cell_data(vts)

    names_point = keys(data_point)
    names_cell = keys(data_cell)

    # Read coordinates
    coords = get_coordinates(vts)

    # Read all the point data (vectors are shaped accordingly)
    data_points = NamedTuple();
    for FieldName in names_point
        dat = ReadField_3D_VTS(vts, FieldName, isCell=false);
        data_points = merge(data_points,dat)               
    end

    # Read all the cell data
    data_cell = NamedTuple();
    for FieldName in names_cell
        dat = ReadField_3D_VTS(vts, FieldName, isCell=true);
        data_cell = merge(data_cell,dat)               
    end

    return coords, data_points, data_cell
end


# The FileName can contain a directory as well; deal with that here
function split_pathname(DirName_base::String, FileName::String)
   FullName = joinpath(DirName_base,FileName)
   id       = findlast("/", FullName)[1];
   
   DirName  = FullName[1:id-1]
   File     = FullName[id+1:end]
   return DirName, File
end

"""
    data::NamedTuple = ReadField_3D_VTS(vts, FieldName; isCell=true)

Reads the data with name `FieldName`  from the `vts` file structure and returns it as a NamedTuple
"""
function ReadField_3D_VTS(vts, FieldName; isCell=true)
   
    # Get the data
    if isCell
        data_f = get_cell_data(vts)
    else
        data_f = get_point_data(vts)
    end

    # Load data
    data_Field  = get_data_reshaped(data_f[FieldName], cell_data=isCell)
    if typeof(data_Field[1])==UInt8
        data_Field = Int64.(data_Field)
    end

    if length(size(data_Field))>3
        # we are dealing with vector fields, so shape it accordingly
        data_Tuple = (ntuple(i->data_Field[i,:,:,:], size(data_Field)[1]),)
    else
        data_Tuple  = (data_Field,)
    end

    # Extract name
    name = filter(x -> !isspace(x), FieldName)  # remove white spaces
    
    id   = findfirst("[", name)
    if !isnothing(id)
        name = name[1:id[1]-1]      # strip out "[" signs
    end

    # Create NamedTuple
    data_out = NamedTuple{(Symbol(name),)}(data_Tuple,);

    return data_out
end


"""
    vtk = Compress_pTatin_VTS_File(FileName::String; DirName_base::String=pwd())

Opens a pTatin `*.vts` file and saves it again but in compressed form.
"""
function Compress_pTatin_VTS_File(FileName::String; DirName_base::String=pwd())

    coord, data_point, data_cell = Read_pTatin_VTS_File(FileName; DirName_base=DirName_base);

    CurDir = pwd();
    DirName, File = split_pathname(DirName_base, FileName)
    cd(DirName)

    # new name
    FileName_compress = FileName[1:end-4]*"_compress"*FileName[end-3:end];

    vtk = vtk_grid(FileName_compress, coord)    # open file

    # add vertex data
    for (i,name) = enumerate(keys(data_point))
        vtk[String(name)] = data_point[i]
    end

    # add cell data
    for (i,name) = enumerate(keys(data_cell))
        vtk[String(name)] = data_cell[i]
    end

    # Save and close vtk file.
    output = vtk_save(vtk)

    cd(CurDir)

    return vtk
end

    
