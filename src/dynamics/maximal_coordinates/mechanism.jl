struct Mechanism
    bodies::Vector{Body}
    joints::Vector{Joint}
end

function Mechanism()
    Mechanism(Body[], Joint[])
end
