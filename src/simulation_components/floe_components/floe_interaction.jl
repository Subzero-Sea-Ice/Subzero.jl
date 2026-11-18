export InteractionFields, floeidx, xforce, yforce, xpoint, ypoint, torque, overlap 

# Enum to index into floe interactions field with more intuituve names
@enum InteractionFields begin
    floeidx = 1
    xforce = 2
    yforce = 3
    xpoint = 4
    ypoint = 5
    torque = 6
    overlap = 7
end

# Index into interactions field with InteractionFields enum objects
Base.to_index(s::InteractionFields) = Int(s)

# Create a range of interactions field columns with InteractionFields enum objects
Base.:(:)(a::InteractionFields, b::InteractionFields) = Int(a):Int(b)