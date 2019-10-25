import JLD2

function ego_save(ego::EGO, filename::AbstractString)
    JLD2.jldopen(filename, "w") do io
        io["ego/x"] = ego.gp.x
        io["ego/y"] = ego.gp.y
    end
    return length(ego)
end

function append!(ego::EGO, filename::AbstractString)
    JLD2.jldopen(filename, "r") do io
        x = io["ego/x"]
        y = io["ego/y"]
        append!(ego, x, y)
    end
end