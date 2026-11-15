#=
Polygon clipping on polygons stored in arrays (`GOG.ArrayPolygon`), so that it can run in GPU
kernels: the intersection of two polygons, as `intersect_polys` (`GO.intersection`), which
uses GeometryOps' Foster-Hormann clipping. Like the functions of `GeometryOpsGPU`, each
function gives the same result as the GeometryOps function of the same name (without the
`!`), bit for bit,
and the functions are tested against them in `test/test_device_clipping.jl`.

GeometryOps clips two polygons A and B in steps:
1. `_build_a_list`: a list of nodes for A: its vertices, with the points where B's edges
   intersect A's edges inserted in order along each edge.
2. `_build_b_list`: the same list for B, with each intersection node linked to the
   equivalent node in the other list (`neighbor`).
3. `_classify_crossing!`, `_flag_ent_exit!`, `_index_crossing_intrs!`: mark which
   intersections cross from the inside to the outside, and in which direction.
4. `_trace_polynodes`: walk the lists from crossing to crossing to form the overlap regions.

This file has steps 1 and 2 so far. GeometryOps keeps the lists in growing vectors of
`PolyNode`s; here they are fixed-size arrays (`_PolyNodes`), and a function returns -1 if
they are too small.

`intersect_polys(p1, p2, FT)` always clips in Float64, also for Float32 floes, because
`GO.intersection` doesn't pass the float type on to the clipping algorithm. So the clipping
type `T` is Float64 to match it, although the functions work with any `T`.
=#

# Where an intersection node is in an overlapping chain of edges, as GeometryOps'
# `start_chain`, `end_chain`, `not_endpoint` (`EndPointType`)
const _START_CHAIN = Int8(1)
const _END_CHAIN = Int8(2)
const _NOT_ENDPOINT = Int8(3)

#=
Lists of clipping nodes, as GeometryOps' `Vector{PolyNode{T}}`: node m of list k has the
fields of a `PolyNode` in column m of row k of these arrays:
- `point[k, m, :]`: the node's point (x, y)
- `inter[k, m]`: whether it's an intersection point
- `neighbor[k, m]`: for an intersection point, the index of the equivalent node in the
  other list (while building the lists, first the index of the other polygon's edge)
- `idx[k, m]`: for a crossing point, its index among the crossing points
- `ent_exit[k, m]`: whether the polygon enters the other polygon at this point
- `crossing[k, m]`: whether the boundaries cross here, rather than touch
- `endpoint[k, m]`: `_START_CHAIN`, `_END_CHAIN`, or `_NOT_ENDPOINT`
- `fracs[k, m, :]`: for an intersection point, the fractions along the edges of both
  polygons where it is
`k` is the list that this value refers to; make one for another list with
`_PolyNodes(nodes, k)`. It's isbits on the device, so kernels can create and pass it.
=#
struct _PolyNodes{PA<:AbstractArray{<:Any, 3}, BA<:AbstractMatrix{Bool},
        IA<:AbstractMatrix{<:Integer}, EA<:AbstractMatrix{Int8}}
    point::PA
    fracs::PA
    inter::BA
    ent_exit::BA
    crossing::BA
    neighbor::IA
    idx::IA
    endpoint::EA
    k::Int
end

Adapt.@adapt_structure _PolyNodes

# `nlists` lists of up to `capacity` nodes with points of type T, on the CPU (move them to a
# device with `adapt`)
_PolyNodes(::Type{T}, nlists, capacity) where T = _PolyNodes(
    zeros(T, nlists, capacity, 2), zeros(T, nlists, capacity, 2),
    zeros(Bool, nlists, capacity), zeros(Bool, nlists, capacity),
    zeros(Bool, nlists, capacity), zeros(Int32, nlists, capacity),
    zeros(Int32, nlists, capacity), fill(_NOT_ENDPOINT, nlists, capacity), 0,
)

# The same arrays, for list k
_PolyNodes(l::_PolyNodes, k) = _PolyNodes(l.point, l.fracs, l.inter, l.ent_exit, l.crossing,
    l.neighbor, l.idx, l.endpoint, Int(k))

_capacity(l::_PolyNodes) = size(l.inter, 2)

#= A node as a NamedTuple with the fields of GeometryOps' `PolyNode`, with its defaults, so
that the code can follow GeometryOps' `PolyNode{T}(; point, inter, ...)`. Copies with some
fields changed, as `PolyNode(node; crossing = true)`, are
`merge(node, (; crossing = true))`. =#
_poly_node(::Type{T}; point = (zero(T), zero(T)), inter = false, neighbor = 0, idx = 0,
    ent_exit = false, crossing = false, endpoint = _NOT_ENDPOINT,
    fracs = (zero(T), zero(T))) where T =
    (; point, inter, neighbor, idx, ent_exit, crossing, endpoint, fracs)

# Node m of the list
function _get_node(l::_PolyNodes, m)
    k = l.k
    return (;
        point = (l.point[k, m, 1], l.point[k, m, 2]),
        inter = l.inter[k, m],
        neighbor = Int(l.neighbor[k, m]),
        idx = Int(l.idx[k, m]),
        ent_exit = l.ent_exit[k, m],
        crossing = l.crossing[k, m],
        endpoint = l.endpoint[k, m],
        fracs = (l.fracs[k, m, 1], l.fracs[k, m, 2]),
    )
end

# Set node m of the list
function _set_node!(l::_PolyNodes, m, node)
    k = l.k
    l.point[k, m, 1], l.point[k, m, 2] = node.point
    l.inter[k, m] = node.inter
    l.neighbor[k, m] = node.neighbor
    l.idx[k, m] = node.idx
    l.ent_exit[k, m] = node.ent_exit
    l.crossing[k, m] = node.crossing
    l.endpoint[k, m] = node.endpoint
    l.fracs[k, m, 1], l.fracs[k, m, 2] = node.fracs
    return
end

#= A list of node indices, as GeometryOps' `a_idx_list`: entry m of list k is `idx[k, m]`.
`k` is the list that this value refers to; make one for another list with
`_IndexList(list, k)`. =#
struct _IndexList{IA<:AbstractMatrix{<:Integer}}
    idx::IA
    k::Int
end

Adapt.@adapt_structure _IndexList

# `nlists` lists of up to `capacity` indices, on the CPU
_IndexList(nlists, capacity) = _IndexList(zeros(Int32, nlists, capacity), 0)
_IndexList(l::_IndexList, k) = _IndexList(l.idx, Int(k))

_capacity(l::_IndexList) = size(l.idx, 2)
Base.getindex(l::_IndexList, m) = Int(l.idx[l.k, m])
Base.setindex!(l::_IndexList, x, m) = (l.idx[l.k, m] = x; l)

# Point j of the polygon, converted to T, as GeometryOps' `_tuple_point(p, T)`
function _tuple_point(::Type{T}, p::GOG.ArrayPolygon, j) where T
    x, y = GOG.getpoint(p, j)
    return (T(x), T(y))
end

#= Sort nodes `first:last` of the list by the fraction along their edge of the list's
polygon, as `sort!(nodes, by = x -> x.fracs[1])`: a stable insertion sort, so nodes with the
same fraction stay in order, as with Julia's (stable) default sort. =#
function _sort_nodes_by_frac!(l::_PolyNodes, first, last)
    for m in (first + 1):last
        node = _get_node(l, m)
        j = m - 1
        while j >= first && isless(node.fracs[1], l.fracs[l.k, j, 1])
            _set_node!(l, j + 1, _get_node(l, j))
            j -= 1
        end
        _set_node!(l, j + 1, node)
    end
    return
end

#= Sort entries `1:n` of the index list by `by(index)`, as `sort!(list, by = by)`: a stable
insertion sort, as Julia's default sort. =#
function _sort_indices!(l::_IndexList, n, by = identity)
    for m in 2:n
        x = l[m]
        key = by(x)
        j = m - 1
        while j >= 1 && isless(key, by(l[j]))
            l[j + 1] = l[j]
            j -= 1
        end
        l[j + 1] = x
    end
    return
end

#=
The node list of polygon A for clipping with polygon B, as `GO._build_a_list(alg, T, a, b;
exact = True())` with the default algorithm `FosterHormannClipping(Planar())`. For each edge
of A, in order: a node for its first point, then the points where B's edges intersect it,
sorted along the edge. A vertex of A that lies on an edge of B becomes an intersection node
instead. Edges of length zero are skipped, but still counted in the edge numbers that
`neighbor` holds.

Writes the nodes to `a` and the indices of the intersection nodes to `a_idx`. Returns the
number of nodes, the number of indices, and the number of intersection points that will be
new nodes in B's list (GeometryOps' `n_b_intrs`), or (-1, 0, 0) if `a` or `a_idx` is too
small.
=#
function _build_a_list!(a::_PolyNodes, a_idx::_IndexList, ::Type{T}, poly_a::GOG.ArrayPolygon,
        poly_b::GOG.ArrayPolygon) where T
    a_count = 0  # number of nodes in a
    n_idx = 0  # number of indices in a_idx
    n_b_intrs = 0
    prev_counter = 0  # index of the node of the current edge's first point
    for i in 1:(GOG.npoint(poly_a) - 1)
        a_pt1, a_pt2 = _tuple_point(T, poly_a, i), _tuple_point(T, poly_a, i + 1)
        a_pt1 == a_pt2 && continue
        a_count == _capacity(a) && return -1, 0, 0
        a_count += 1
        _set_node!(a, a_count, _poly_node(T; point = a_pt1))
        prev_counter = a_count
        for j in 1:(GOG.npoint(poly_b) - 1)
            b_pt1, b_pt2 = _tuple_point(T, poly_b, j), _tuple_point(T, poly_b, j + 1)
            b_pt1 == b_pt2 && continue
            line_orient, intr1, intr2 = GOG._intersection_point(T, (a_pt1, a_pt2), (b_pt1, b_pt2))
            line_orient == GOG.line_out && continue
            if line_orient == GOG.line_cross  # intersection point that isn't a vertex
                int_pt, fracs = intr1
                (a_count == _capacity(a) || n_idx == _capacity(a_idx)) && return -1, 0, 0
                a_count += 1
                n_b_intrs += 1
                _set_node!(a, a_count, _poly_node(T; point = int_pt, inter = true,
                    neighbor = j, crossing = true, fracs))
                n_idx += 1
                a_idx[n_idx] = a_count
            else
                (_, (α1, β1)) = intr1
                # whether a1 or b1 is an intersection node in a
                add_a1 = α1 == 0 && 0 ≤ β1 < 1
                a1_β = add_a1 ? β1 : zero(T)
                add_b1 = β1 == 0 && 0 < α1 < 1
                b1_α = add_b1 ? α1 : zero(T)
                # collinear, overlapping edges have a second intersection point
                if line_orient == GOG.line_over
                    (_, (α2, β2)) = intr2
                    if α2 == 0 && 0 ≤ β2 < 1
                        add_a1, a1_β = true, β2
                    end
                    if β2 == 0 && 0 < α2 < 1
                        add_b1, b1_α = true, α2
                    end
                end
                if add_a1
                    n_b_intrs += a1_β == 0 ? 0 : 1
                    _set_node!(a, prev_counter, _poly_node(T; point = a_pt1, inter = true,
                        neighbor = j, fracs = (zero(T), a1_β)))
                    n_idx == _capacity(a_idx) && return -1, 0, 0
                    n_idx += 1
                    a_idx[n_idx] = prev_counter
                end
                if add_b1
                    a_count == _capacity(a) && return -1, 0, 0
                    n_idx == _capacity(a_idx) && return -1, 0, 0
                    a_count += 1
                    _set_node!(a, a_count, _poly_node(T; point = b_pt1, inter = true,
                        neighbor = j, fracs = (b1_α, zero(T))))
                    n_idx += 1
                    a_idx[n_idx] = a_count
                end
            end
        end
        # order the edge's intersection points along the edge
        _sort_nodes_by_frac!(a, prev_counter + 1, a_count)
    end
    return a_count, n_idx, n_b_intrs
end

#=
The node list of polygon B, as `GO._build_b_list(alg, T, a_idx_list, a_list, n_b_intrs,
b)`, from A's list made by `_build_a_list!`: B's vertices, with the intersection nodes of
A's list inserted along B's edges. A vertex of B that is an intersection point becomes an
intersection node instead. Links each intersection node of A's list with its node in B's
list through `neighbor`, and puts `a_idx` back in increasing order.

Writes the nodes to `b`. Returns the number of nodes, or -1 if `b` is too small.
=#
function _build_b_list!(b::_PolyNodes, a::_PolyNodes, a_idx::_IndexList, n_idx,
        ::Type{T}, poly_b::GOG.ArrayPolygon) where T
    # order the intersection nodes by their place along B
    _sort_indices!(a_idx, n_idx, x -> a.neighbor[a.k, x] + a.fracs[a.k, x, 2])
    intr_curr = 1
    b_count = 0
    b_pt1 = (zero(T), zero(T))
    for i in 1:GOG.npoint(poly_b)
        b_pt2 = _tuple_point(T, poly_b, i)
        if i ≤ 1 || b_pt1 == b_pt2  # don't repeat points
            b_pt1 = b_pt2
            continue
        end
        b_count == _capacity(b) && return -1
        b_count += 1
        _set_node!(b, b_count, _poly_node(T; point = b_pt1))
        if intr_curr ≤ n_idx
            curr_idx = a_idx[intr_curr]
            curr_node = _get_node(a, curr_idx)
            prev_counter = b_count
            # add all intersection points on edge i - 1
            while curr_node.neighbor == i - 1
                new_intr = merge(curr_node, (; neighbor = curr_idx))
                if curr_node.fracs[2] == 0  # the intersection point is B's vertex
                    b_idx = prev_counter
                else
                    b_count == _capacity(b) && return -1
                    b_count += 1
                    b_idx = b_count
                end
                _set_node!(b, b_idx, new_intr)
                a.neighbor[a.k, curr_idx] = b_idx
                intr_curr += 1
                intr_curr > n_idx && break
                curr_idx = a_idx[intr_curr]
                curr_node = _get_node(a, curr_idx)
            end
        end
        b_pt1 = b_pt2
    end
    _sort_indices!(a_idx, n_idx)  # back to the order of the nodes in A's list
    return b_count
end
