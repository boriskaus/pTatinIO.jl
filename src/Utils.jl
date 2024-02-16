# 

"""
    data_cell::NamedTuple = Vertex2Cell(data_vertex::NamedTuple)

Interpolates vertex data to cell centers when both are given as NamedTuples
"""
function Vertex2Cell(data_vertex::NamedTuple)
    data_cell =  NamedTuple()
    N = length(data_vertex)
    names = keys(data_vertex)
    for i = 1:N
        data_cell_new = NamedTuple{(names[i],)}((Vertex2Cell(data_vertex[i]),))

        data_cell = (data_cell..., data_cell_new)
    end
    return data_cell
end


"""
    data_cell::NTuple = Vertex2Cell(data_vertex::NTuple{N,T})

Interpolates vertex data to cell centers when both are given as Tuples
"""
function Vertex2Cell(data_vertex::NTuple{N,T}) where {N,T<:Array}
    data_cell = ();
    for i = 1:N
        data_cell = (data_cell..., Vertex2Cell(data_vertex[i]))
    end
    return data_cell
end

"""
    data_cell = Vertex2Cell(data_vertex::Array{_T,3})

Interpolate 3D array from vertexes -> centers  
"""
function Vertex2Cell(data_vertex::Array{_T,3}) where _T<:Number

    N = size(data_vertex)
    data_cell = zeros(N .- 1)
    for I in CartesianIndices(data_cell)
        data_cell[I] = 0.125*(data_vertex[I + CartesianIndex(0,0,0)] + 
                              data_vertex[I + CartesianIndex(1,0,0)] +  
                              data_vertex[I + CartesianIndex(0,1,0)] +  
                              data_vertex[I + CartesianIndex(1,1,0)] +
                              data_vertex[I + CartesianIndex(0,0,1)] + 
                              data_vertex[I + CartesianIndex(1,0,1)] +  
                              data_vertex[I + CartesianIndex(0,1,1)] +  
                              data_vertex[I + CartesianIndex(1,1,1)]);
    end

    return data_cell
end
