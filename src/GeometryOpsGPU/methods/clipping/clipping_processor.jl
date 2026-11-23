#=
Foster-Hormann polygon clipping, as GeometryOps' `clipping_processor.jl`, for `ArrayPolygon`s
and with the node lists and results in fixed-size buffers (`ClipLists`, `PolyRegions`).
GeometryOps clips two polygons A and B in steps:
1. `_build_a_list`: a list of nodes for A: its vertices, with the points where B's edges
   intersect A's edges inserted in order along each edge.
2. `_build_b_list`: the same list for B, with each intersection node linked to the
   equivalent node in the other list (`neighbor`).
3. `_classify_crossing!`, `_flag_ent_exit!`, `_index_crossing_intrs!`: mark which
   intersections cross from the inside to the outside, and in which direction.
4. `_trace_polynodes`: walk the lists from crossing to crossing to form the regions.
5. `_find_non_cross_orientation`: if the boundaries don't cross, whether one polygon is
   inside the other.
6. `_remove_collinear_points!`: remove the vertices of the regions that are on a straight
   line between their neighbours.
The functions here have the names of those steps, with a `!` if they write into a buffer
instead of returning a new vector. `intersection!` (`intersection.jl`) runs them. Polygons
have no holes, so GeometryOps' steps for holes aren't needed.

GeometryOps keeps the lists in growing vectors of `PolyNode`s and the regions in a vector of
polygons; here the functions return `CLIP_TOO_SMALL` if the buffers are too small. Where
GeometryOps throws an error (a `TracingError`, or an index that isn't found), they return
`CLIP_FAILED` instead.
=#

#= What the clipping functions return instead of a count when the buffers are too small, or
where GeometryOps would throw an error. GeometryOps has no equivalent. =#
const CLIP_TOO_SMALL = -1
const CLIP_FAILED = -2

# Which side of the edges a point is on, as GeometryOps' `PointEdgeSide`
@enum PointEdgeSide left=1 right=2 unknown=3

#= Where an intersection node is in an overlapping chain of edges, as GeometryOps'
`EndPointType`, with the same values, stored in an Int8 to keep the node lists small =#
@enum EndPointType::Int8 start_chain=1 end_chain=2 not_endpoint=3

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
- `endpoint[k, m]`: `start_chain`, `end_chain`, or `not_endpoint`
- `fracs[k, m, :]`: for an intersection point, the fractions along the edges of both
  polygons where it is
`k` is the list that this value refers to; make one for another list with
`_PolyNodes(nodes, k)`. It's isbits on the device, so kernels can create and pass it.
=#
struct _PolyNodes{PA<:AbstractArray{<:Any, 3}, BA<:AbstractMatrix{Bool},
        IA<:AbstractMatrix{<:Integer}, EA<:AbstractMatrix{EndPointType}}
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
    zeros(Int32, nlists, capacity), fill(not_endpoint, nlists, capacity), 0,
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
    ent_exit = false, crossing = false, endpoint = not_endpoint,
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

# The point (x, y) converted to T, as GeometryOps' `_tuple_point(p, T)`
_tuple_point(p, ::Type{T}) where T = (T(p[1]), T(p[2]))

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
function _build_a_list!(a::_PolyNodes, a_idx::_IndexList, ::Type{T}, poly_a::ArrayPolygon,
        poly_b::ArrayPolygon) where T
    a_count = 0  # number of nodes in a
    n_idx = 0  # number of indices in a_idx
    n_b_intrs = 0
    prev_counter = 0  # index of the node of the current edge's first point
    for i in 1:(npoint(poly_a) - 1)
        a_pt1 = _tuple_point(getpoint(poly_a, i), T)
        a_pt2 = _tuple_point(getpoint(poly_a, i + 1), T)
        a_pt1 == a_pt2 && continue
        a_count == _capacity(a) && return -1, 0, 0
        a_count += 1
        _set_node!(a, a_count, _poly_node(T; point = a_pt1))
        prev_counter = a_count
        for j in 1:(npoint(poly_b) - 1)
            b_pt1 = _tuple_point(getpoint(poly_b, j), T)
            b_pt2 = _tuple_point(getpoint(poly_b, j + 1), T)
            b_pt1 == b_pt2 && continue
            line_orient, intr1, intr2 = _intersection_point(T, (a_pt1, a_pt2), (b_pt1, b_pt2))
            line_orient == line_out && continue
            if line_orient == line_cross  # intersection point that isn't a vertex
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
                if line_orient == line_over
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
        ::Type{T}, poly_b::ArrayPolygon) where T
    # order the intersection nodes by their place along B
    _sort_indices!(a_idx, n_idx, x -> a.neighbor[a.k, x] + a.fracs[a.k, x, 2])
    intr_curr = 1
    b_count = 0
    b_pt1 = (zero(T), zero(T))
    for i in 1:npoint(poly_b)
        b_pt2 = _tuple_point(getpoint(poly_b, i), T)
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

# Point m of the list
_node_point(l::_PolyNodes, m) = (l.point[l.k, m, 1], l.point[l.k, m, 2])

# Whether two nodes are the same, as GeometryOps' `equals` (`==`) of two `PolyNode`s
_same_node(n1, n2) = n1.point == n2.point && n1.inter == n2.inter && n1.fracs == n2.fracs

# Whether the node is a vertex of one of the polygons, as GeometryOps' `_is_vertex`
_is_vertex(node) = !node.inter || node.fracs[1] == 0 || node.fracs[1] == 1 ||
    node.fracs[2] == 0 || node.fracs[2] == 1

# Midpoint of nodes m1 and m2, as GeometryOps' `(p1 .+ p2) ./ 2`
function _node_midpoint(l::_PolyNodes, m1, m2)
    (x1, y1), (x2, y2) = _node_point(l, m1), _node_point(l, m2)
    return ((x1 + x2) / 2, (y1 + y2) / 2)
end

#= The last vertex node at or before node `m` of the list's `n` nodes, or else the last one
in the list, as GeometryOps' `findprev(_is_vertex, list, m)` and then
`findlast(_is_vertex, list)`. Every list has vertex nodes, so one is always found. =#
function _prev_vertex(l::_PolyNodes, m, n)
    for j in m:-1:1
        _is_vertex(_get_node(l, j)) && return j
    end
    for j in n:-1:1
        _is_vertex(_get_node(l, j)) && return j
    end
    return 1  # not reached
end

# The first vertex node at or after node `m`, or else the first one in the list, as
# `findnext(_is_vertex, list, m)` and then `findfirst(_is_vertex, list)`
function _next_vertex(l::_PolyNodes, m, n)
    for j in m:n
        _is_vertex(_get_node(l, j)) && return j
    end
    for j in 1:n
        _is_vertex(_get_node(l, j)) && return j
    end
    return 1  # not reached
end

# Which side of the line P1-P2-P3 the point Q is on, as `GO._get_side(Q, P1, P2, P3;
# exact = True())`
function _get_side(Q, P1, P2, P3)
    s1 = Predicates.orient(Q, P1, P2; exact = True())
    s2 = Predicates.orient(Q, P2, P3; exact = True())
    s3 = Predicates.orient(P1, P2, P3; exact = True())
    return if s3 ≥ 0
        (s1 < 0) || (s2 < 0) ? right : left
    else
        (s1 > 0) || (s2 > 0) ? left : right
    end
end

#= Which side of the line through a_prev, curr_pt and a_next the neighbours of curr_pt in B's
list are on, as `GO._get_sides(b_prev, b_next, a_prev, curr_pt, a_next, i, j, a_list, b_list;
exact = True())`. Intersection nodes that aren't vertices are replaced by the nearest vertex
along the polygon, because their rounded points may not be on the edge exactly. =#
function _get_sides(b_prev, b_next, a_prev, curr_pt, a_next, i, j, a::_PolyNodes, na,
        b::_PolyNodes, nb)
    b_prev_pt = _is_vertex(b_prev) ? b_prev.point : _node_point(b, _prev_vertex(b, j - 1, nb))
    b_next_pt = _is_vertex(b_next) ? b_next.point : _node_point(b, _next_vertex(b, j + 1, nb))
    a_prev_pt = _is_vertex(a_prev) ? a_prev.point : _node_point(a, _prev_vertex(a, i - 1, na))
    a_next_pt = _is_vertex(a_next) ? a_next.point : _node_point(a, _next_vertex(a, i + 1, na))
    b_prev_side = _get_side(b_prev_pt, a_prev_pt, curr_pt.point, a_next_pt)
    b_next_side = _get_side(b_next_pt, a_prev_pt, curr_pt.point, a_next_pt)
    return b_prev_side, b_next_side
end

#= Mark nodes `start_idx` and `end_idx` of A's list, and their neighbours in B's list, as the
start and end of an overlapping chain of edges, crossing or not, as the two places in
`GO._classify_crossing!` that close a chain. The end is written first, as there. =#
function _mark_chain!(a::_PolyNodes, b::_PolyNodes, start_idx, end_idx, crossing, same_winding)
    j = a.neighbor[a.k, end_idx]
    a.crossing[a.k, end_idx] = crossing
    a.endpoint[a.k, end_idx] = end_chain
    b.crossing[b.k, j] = crossing
    b.endpoint[b.k, j] = same_winding ? end_chain : start_chain
    j = a.neighbor[a.k, start_idx]
    a.crossing[a.k, start_idx] = crossing
    a.endpoint[a.k, start_idx] = start_chain
    b.crossing[b.k, j] = crossing
    b.endpoint[b.k, j] = same_winding ? start_chain : end_chain
    return
end

#=
Mark the intersection nodes where the boundaries cross rather than touch, as
`GO._classify_crossing!(alg, T, a_list, b_list; exact = True())`, for the `na` nodes of A's
list and the `nb` nodes of B's list. Where edges overlap, in a chain of edges, the start and
end of the chain are marked as such, and as crossing if the boundaries cross over the whole
chain.

Returns false where GeometryOps would throw an error (a chain whose start isn't found).
=#
function _classify_crossing!(a::_PolyNodes, na, b::_PolyNodes, nb)
    # start centred on the last node
    a_prev = _get_node(a, na - 1)
    curr_pt = _get_node(a, na)
    i = na
    # unmatched overlapping chains
    start_chain_edge, start_chain_idx = unknown, 0
    unmatched_end_chain_edge, unmatched_end_chain_idx = unknown, 0
    same_winding = true
    for next_idx in 1:na
        a_next = _get_node(a, next_idx)
        if curr_pt.inter && !curr_pt.crossing
            j = curr_pt.neighbor
            b_prev = _get_node(b, j == 1 ? nb : j - 1)
            b_next = _get_node(b, j == nb ? 1 : j + 1)
            # whether any edges are on top of each other
            a_prev_is_b_prev = a_prev.inter && _same_node(a_prev, b_prev)
            a_prev_is_b_next = a_prev.inter && _same_node(a_prev, b_next)
            a_next_is_b_prev = a_next.inter && _same_node(a_next, b_prev)
            a_next_is_b_next = a_next.inter && _same_node(a_next, b_next)
            b_prev_side, b_next_side = _get_sides(b_prev, b_next, a_prev, curr_pt, a_next,
                i, j, a, na, b, nb)
            if !a_prev_is_b_prev && !a_prev_is_b_next && !a_next_is_b_prev && !a_next_is_b_next
                # no edges overlap: the boundaries cross if B's edges are on different sides
                if b_prev_side != b_next_side
                    a.crossing[a.k, i] = true
                    b.crossing[b.k, j] = true
                end
            elseif !a_next_is_b_prev && !a_next_is_b_next  # end of an overlapping chain
                b_side = a_prev_is_b_prev ? b_next_side : b_prev_side
                if start_chain_edge == unknown  # started in the middle of the chain
                    unmatched_end_chain_edge = b_side
                    unmatched_end_chain_idx = i
                    same_winding = a_prev_is_b_prev
                else
                    _mark_chain!(a, b, start_chain_idx, i, b_side != start_chain_edge,
                        same_winding)
                end
            elseif !a_prev_is_b_prev && !a_prev_is_b_next  # start of an overlapping chain
                start_chain_edge = a_next_is_b_prev ? b_next_side : b_prev_side
                start_chain_idx = i
                same_winding = a_next_is_b_next
            end
        end
        a_prev = curr_pt
        curr_pt = a_next
        i = next_idx
    end
    # close the chain that the loop started in the middle of
    if unmatched_end_chain_edge != unknown
        start_chain_idx == 0 && return false
        _mark_chain!(a, b, start_chain_idx, unmatched_end_chain_idx,
            unmatched_end_chain_edge != start_chain_edge, same_winding)
    end
    return true
end

#= Where to start marking entries and exits, as `GO._pt_off_edge_status(alg, T, list, poly,
n; exact = True())`: the node after the first node that isn't an intersection point, or, if
all are, after the first node from which the edge leaves the other polygon's boundary.
Returns that node and whether it's an entry into the other polygon (that is, whether the
start is outside it), or (0, false) if there is no such node. =#
function _pt_off_edge_status(l::_PolyNodes, n, poly::ArrayPolygon)
    start_idx, is_non_intr_pt = 0, true
    for m in 1:n
        if !l.inter[l.k, m]
            start_idx = m
            break
        end
    end
    if start_idx == 0
        is_non_intr_pt = false
        for m in 1:n
            # GeometryOps' `_next_edge_off`
            endpoint = l.endpoint[l.k, m]
            if endpoint == end_chain || (l.crossing[l.k, m] && endpoint == not_endpoint)
                start_idx = m
                break
            end
        end
        start_idx == 0 && return 0, false
    end
    next_idx = start_idx < n ? (start_idx + 1) : 1
    start_pt = is_non_intr_pt ? _node_point(l, start_idx) : _node_midpoint(l, start_idx, next_idx)
    in_poly = _point_filled_curve_orientation(start_pt, poly; exact = True()) == point_in
    return next_idx, !in_poly
end

#=
Mark each crossing in the `n` nodes of the list as an entry into the other polygon `poly`
or an exit from it, as `GO._flag_ent_exit!(alg, T, GI.LinearRingTrait(), poly, list,
delay_cross_f, delay_bounce_f; exact = True())` with the functions of an intersection
(`_inter_delay_cross_f`, `_inter_delay_bounce_f`). Also decides which end of an overlapping
chain counts as the crossing for tracing.

Returns false where GeometryOps would throw an error (the end of a chain whose start
hasn't been found).
=#
function _flag_ent_exit!(l::_PolyNodes, n, poly::ArrayPolygon)
    next_idx, status = _pt_off_edge_status(l, n, poly)
    next_idx == 0 && return true
    start_chain_idx = 0
    # all nodes, starting from next_idx
    for s in 0:(n - 1)
        ii = next_idx + s > n ? next_idx + s - n : next_idx + s
        curr_pt = _get_node(l, ii)
        if curr_pt.endpoint == start_chain
            start_chain_idx = ii
        elseif curr_pt.crossing || curr_pt.endpoint == end_chain
            start_crossing, end_crossing = curr_pt.crossing, curr_pt.crossing
            if curr_pt.endpoint == end_chain  # end of an overlapping chain
                start_chain_idx == 0 && return false
                if curr_pt.crossing  # delayed crossing: `_inter_delay_cross_f(status)`
                    start_crossing, end_crossing = !status, status
                else  # delayed bouncing: `_inter_delay_bounce_f(pt_in_poly, _)`
                    next_m = ii < n ? (ii + 1) : 1
                    next_val = _node_midpoint(l, ii, next_m)
                    start_crossing =
                        _point_filled_curve_orientation(next_val, poly; exact = True()) == point_in
                    end_crossing = start_crossing
                end
                l.ent_exit[l.k, start_chain_idx] = status
                l.crossing[l.k, start_chain_idx] = start_crossing
                if !curr_pt.crossing
                    status = !status
                end
            end
            l.ent_exit[l.k, ii] = status
            l.crossing[l.k, ii] = end_crossing
            status = !status
        end
    end
    return true
end

#= Keep only the crossings among the `n_idx` indices of A's intersection nodes, and number
them in both lists (the nodes' `idx`), as `GO._index_crossing_intrs!(alg, a_list, b_list,
a_idx_list)`. Returns the number of crossings. =#
function _index_crossing_intrs!(a::_PolyNodes, b::_PolyNodes, a_idx::_IndexList, n_idx)
    n = 0
    for m in 1:n_idx
        x = a_idx[m]
        if a.crossing[a.k, x]
            n += 1
            a_idx[n] = x
        end
    end
    for m in 1:n
        x = a_idx[m]
        a.idx[a.k, x] = m
        b.idx[b.k, a.neighbor[a.k, x]] = m
    end
    return n
end

"""
    PolyRegions(T, nlists, max_polys, capacity)
    PolyRegions(regions, k)

Buffer for the polygons that clipping returns, as the vector of polygons that GeometryOps
returns: `nlists` lists of up to `max_polys` polygons, each with up to `capacity` points of
type `T`. The first form makes the arrays on the CPU (move them to a device with `adapt`);
the second refers to list `k` of the same arrays, e.g. the list of one pair of polygons in a
kernel. It's isbits on a device, so kernels can create and pass it.

Polygon `m` of list `k` is [`region(regions, m)`](@ref GeometryOpsGPU.region), an
`ArrayPolygon`. GeometryOps has no equivalent.
"""
struct PolyRegions{A<:AbstractArray{<:Any, 3}, IA<:AbstractVector{<:Integer}}
    coords::A
    n::IA
    max_polys::Int
    k::Int
end

Adapt.@adapt_structure PolyRegions

PolyRegions(::Type{T}, nlists, max_polys, capacity) where T = PolyRegions(
    zeros(T, nlists * max_polys, capacity, 2), zeros(Int32, nlists * max_polys),
    Int(max_polys), 0)

PolyRegions(r::PolyRegions, k) = PolyRegions(r.coords, r.n, r.max_polys, Int(k))

_capacity(r::PolyRegions) = size(r.coords, 2)
_region_row(r::PolyRegions, m) = (r.k - 1) * r.max_polys + m

"""
    region(regions::PolyRegions, m)

Polygon `m` of the list that `regions` refers to, as an `ArrayPolygon` (row
`(k - 1) * max_polys + m` of the arrays). GeometryOps has no equivalent.
"""
function region(r::PolyRegions, m)
    row = _region_row(r, m)
    return ArrayPolygon(r.coords, row, Int(r.n[row]))
end

# Point j of row `row` of the regions' arrays
@inline _region_point(r::PolyRegions, row, j) = (r.coords[row, j, 1], r.coords[row, j, 2])

# Set point j of polygon m of the list
@inline function _set_region_point!(r::PolyRegions, m, j, (x, y))
    row = _region_row(r, m)
    r.coords[row, j, 1] = x
    r.coords[row, j, 2] = y
    return
end

# Copy polygon m1 of the list to polygon m2
function _copy_region!(r::PolyRegions, m2, m1)
    row1, row2 = _region_row(r, m1), _region_row(r, m2)
    for j in 1:r.n[row1]
        r.coords[row2, j, 1] = r.coords[row1, j, 1]
        r.coords[row2, j, 2] = r.coords[row1, j, 2]
    end
    r.n[row2] = r.n[row1]
    return
end

#=
Walk the lists from crossing to crossing to form the overlap regions, as
`GO._trace_polynodes(alg, T, a_list, b_list, a_idx_list, _inter_step, poly_a, poly_b)`:
start at a crossing of A's list, and follow the list forwards from an entry and backwards
from an exit, switching to the other list at each crossing, until it's back at the start.
Uses up the `n_cross` indices of the crossings in `a_idx`.

Writes the regions to `out` and returns their number, `CLIP_TOO_SMALL` if `out` is too
small, or `CLIP_FAILED` where GeometryOps throws a `TracingError`.
=#
function _trace_polynodes!(out::PolyRegions, a::_PolyNodes, na, b::_PolyNodes, nb,
        a_idx::_IndexList, n_cross)
    total_pts = na + nb
    visited_pts = 0
    processed_pts = 0
    first_idx = 1
    n_polys = 0
    while processed_pts < n_cross
        on_a_list = true
        visited_pts += 1
        processed_pts += 1
        # the first crossing not processed yet, as `findnext(x -> x != 0, a_idx_list, ...)`
        while first_idx ≤ n_cross && a_idx[first_idx] == 0
            first_idx += 1
        end
        first_idx > n_cross && return CLIP_FAILED
        idx = a_idx[first_idx]
        a_idx[first_idx] = 0
        start_pt = _get_node(a, idx)
        start_neighbor = _get_node(b, start_pt.neighbor)
        curr = start_pt
        # start the region with the first point
        n_polys == out.max_polys && return CLIP_TOO_SMALL
        n_polys += 1
        n_pts = 1
        _set_region_point!(out, n_polys, n_pts, curr.point)
        curr_not_start = true
        while curr_not_start
            curr_list, curr_npoints = on_a_list ? (a, na) : (b, nb)
            # `_inter_step`: forwards from an entry, backwards from an exit
            step = curr.ent_exit ? 1 : -1
            same_status, prev_status = true, curr.ent_exit
            while same_status
                visited_pts >= total_pts && return CLIP_FAILED
                idx += step
                idx = (idx > curr_npoints) ? mod(idx, curr_npoints) : idx
                idx = (idx == 0) ? curr_npoints : idx
                curr = _get_node(curr_list, idx)
                n_pts == _capacity(out) && return CLIP_TOO_SMALL
                n_pts += 1
                _set_region_point!(out, n_polys, n_pts, curr.point)
                if curr.crossing || curr.endpoint != not_endpoint
                    same_status = curr.ent_exit == prev_status
                    curr_not_start = !_same_node(curr, start_pt) &&
                        !_same_node(curr, start_neighbor)
                    !curr_not_start && break
                    if (on_a_list && curr.crossing) ||
                            (!on_a_list && a.crossing[a.k, curr.neighbor])
                        processed_pts += 1
                        a_idx[curr.idx] = 0
                    end
                end
                visited_pts += 1
            end
            # switch to the equivalent node in the other list
            on_a_list = !on_a_list
            idx = curr.neighbor
            curr = _get_node(on_a_list ? a : b, idx)
        end
        out.n[_region_row(out, n_polys)] = n_pts
    end
    return n_polys
end

#= If the boundaries don't cross, whether A is inside B and whether B is inside A (possibly
sharing edges or points), as `GO._find_non_cross_orientation(a_list, b_list, a_poly, b_poly;
exact = True())`, from a node of each list that isn't an intersection point. =#
function _find_non_cross_orientation(a::_PolyNodes, na, b::_PolyNodes, nb,
        poly_a::ArrayPolygon, poly_b::ArrayPolygon)
    a_pt_orient = point_on  # if all nodes are on the other polygon's boundary
    for m in 1:na
        if !a.inter[a.k, m]
            a_pt_orient = _point_filled_curve_orientation(_node_point(a, m), poly_b; exact = True())
            break
        end
    end
    b_pt_orient = point_on
    for m in 1:nb
        if !b.inter[b.k, m]
            b_pt_orient = _point_filled_curve_orientation(_node_point(b, m), poly_a; exact = True())
            break
        end
    end
    a_in_b = a_pt_orient != point_out && b_pt_orient != point_in
    b_in_a = b_pt_orient != point_out && a_pt_orient != point_in
    return a_in_b, b_in_a
end

# Copy the polygon to polygon m of `out`, converted to T, as `GI.Polygon([tuples(ring)])`.
# Returns false if `out` is too small.
function _copy_poly_to_region!(out::PolyRegions, m, poly::ArrayPolygon, ::Type{T}) where T
    npoint(poly) > _capacity(out) && return false
    for j in 1:npoint(poly)
        _set_region_point!(out, m, j, _tuple_point(getpoint(poly, j), T))
    end
    out.n[_region_row(out, m)] = npoint(poly)
    return true
end

#=
Remove the points of the `n_polys` regions that are on the straight line between their
neighbours (within round-off), as `GO._remove_collinear_points!(alg, polys, remove_idx,
poly_a, poly_b)`, and remove the regions that have fewer than 3 points left. Whether a point
is removed depends on its neighbours before any are removed. Returns the number of regions
left.
=#
function _remove_collinear_points!(out::PolyRegions, n_polys)
    n_kept = 0
    for m in 1:n_polys
        row = _region_row(out, m)
        n = Int(out.n[row])
        # whether the first point (repeated as the last point) is needed
        first_removed = Predicates.orient(_region_point(out, row, n - 1),
            _region_point(out, row, 1), _region_point(out, row, 2); exact = False()) == 0
        # keep the points in place, overwriting the ones that are removed
        n_left = 0
        if !first_removed
            n_left += 1  # point 1 stays where it is
        end
        p1, p2 = _region_point(out, row, 1), _region_point(out, row, 2)
        for j in 3:n
            p3 = _region_point(out, row, j)
            if Predicates.orient(p1, p2, p3; exact = False()) != 0  # keep point j - 1
                n_left += 1
                _set_region_point!(out, m, n_left, p2)
            end
            p1, p2 = p2, p3
        end
        if !first_removed  # keep the last point
            n_left += 1
            _set_region_point!(out, m, n_left, _region_point(out, row, n))
        end
        # whether enough points are left to form a polygon
        n_left ≤ (first_removed ? 2 : 3) && continue
        if first_removed  # repeat the new first point as the last
            n_left += 1
            _set_region_point!(out, m, n_left, _region_point(out, row, 1))
        end
        out.n[row] = n_left
        n_kept += 1
        n_kept < m && _copy_region!(out, n_kept, m)
    end
    return n_kept
end

"""
    ClipLists(T, nlists, capacity)
    ClipLists(lists, k)

Buffer for the node lists that clipping two polygons needs, as GeometryOps' `a_list`,
`b_list` and `a_idx_list`: for `nlists` pairs of polygons, lists of up to `capacity` nodes
with points of type `T`. The first form makes the arrays on the CPU (move them to a device
with `adapt`); the second refers to pair `k` of the same arrays, e.g. in a kernel. It's
isbits on a device. GeometryOps has no equivalent.
"""
struct ClipLists{L<:_PolyNodes, I<:_IndexList}
    a::L
    b::L
    a_idx::I
end

Adapt.@adapt_structure ClipLists

ClipLists(::Type{T}, nlists, capacity) where T = ClipLists(_PolyNodes(T, nlists, capacity),
    _PolyNodes(T, nlists, capacity), _IndexList(nlists, capacity))
ClipLists(l::ClipLists, k) = ClipLists(_PolyNodes(l.a, k), _PolyNodes(l.b, k),
    _IndexList(l.a_idx, k))
