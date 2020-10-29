using JuMP, Cbc

struct Arc
    i::Int64
    j::Int64
    c::Int64
end

struct Instance
    n::Int64
    m::Int64
    s::Int64
    t::Int64
    arcs::Array{Arc, 1}
end

struct SolutionArc
    i::Int64
    j::Int64
    flow::Int64
end

struct Solution
    v::Int64
    arcs::Array{SolutionArc, 1}
end

function readInstance(instanceFilePath::String)
    n = 0
    m = 0
    s = 0
    t = 0
    arcs = Arc[]

    count = 0
    for word in eachline(instanceFilePath)
        if count == 0
            n = parse(Int64, word)
        elseif count == 1
            m = parse(Int64, word)
        elseif count == 2
            s = parse(Int64, word)
        elseif count == 3
            t = parse(Int64, word)
        else 
            wordArray = split(word, " ")
            push!(arcs, Arc(parse(Int64, wordArray[1]), parse(Int64, wordArray[2]), parse(Int64, wordArray[3])))
        
        end
        count = count + 1
    end

    push!(arcs, Arc(t, s, typemax(Int64)))
    return Instance(n, m, s, t, arcs)
end

function solve(instance::Instance)
    PFM_model = Model(Cbc.Optimizer)
    set_optimizer_attribute(PFM_model, "logLevel", 1)

    arc_vertices = []
    for arc in instance.arcs
        if isempty(arc_vertices)
            push!(arc_vertices, arc.i)
            push!(arc_vertices, arc.j)
        elseif !isassigned(arc_vertices, arc.i)
            push!(arc_vertices, arc.i)
        end
        if !isassigned(arc_vertices, arc.j)
            push!(arc_vertices, arc.j)
        end
    end

    arc_pairs = []
    for arc in instance.arcs
        push!(arc_pairs, (arc.i, arc.j))
    end

    @variable(PFM_model, 0 <= x[i=1:instance.n, j=1:instance.n;(i, j) in arc_pairs], integer = true)    
    for arc in instance.arcs
       set_upper_bound(x[arc.i, arc.j], arc.c) 
    end

    @objective(PFM_model, Min, -x[instance.t, instance.s])

    for i in arc_vertices
        inFlow = 0
        outFlow = 0
        if i in arc_vertices
            outFlow = sum(x[i, j] for j in arc_vertices if (i, j) in  arc_pairs)
        end
        if i in arc_vertices
            inFlow = sum(x[j, i] for j in arc_vertices if (j, i) in  arc_pairs)
        end
        @constraint(PFM_model, (outFlow - inFlow) == 0)
    end

    JuMP.optimize!(PFM_model)

    arcsSolution = []
    for arc in instance.arcs
        if value(x[arc.i, arc.j]) > 0
            auxArc = SolutionArc(arc.i, arc.j, value(x[arc.i, arc.j]))
            push!(arcsSolution, auxArc)
        end
    end

    return Solution(objective_value(PFM_model), arcsSolution)
end

instance = readInstance("instancias/instance1.txt")
z = solve(instance)
println(z)