export StatusTag

# Enum for differnt floe status
@enum StatusTag begin
    active = 1
    remove = 2
    fuse = 3
end

# Status struct
mutable struct Status
    tag::StatusTag
    fuse_idx::Vector{Int}
end

# Default status is active 
Status(; tag = active) = Status(tag, Vector{Int}())  # active floe