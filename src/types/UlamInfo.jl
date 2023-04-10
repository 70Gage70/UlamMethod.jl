"""
    UlamInfo{S, U}

A contained for some of the results of [`ulam_method`](@ref).
"""
struct UlamInfo{S<:AbstractString, U<:Integer}
    n_polys_requested::U
    n_polys_no_data::U
    n_polys_dis::U
    n_counts_total::U
    n_counts_removed_scc::U
    poly_type::S

    function UlamInfo(
        n_polys_requested::U,
        n_polys_no_data::U,
        n_polys_dis::U,
        n_counts_total::U,
        n_counts_removed_scc::U,
        poly_type::S) where {S<:AbstractString, U<:Integer}
    
        new{String, Int64}(
        n_polys_requested,
        n_polys_no_data,
        n_polys_dis,
        n_counts_total,
        n_counts_removed_scc,
        poly_type)
    end
end

function Base.show(io::IO, x::UlamInfo)
    print(io, " UlamInfo")
    println(io)
    print("  Polys requested: ")
    show(io, x.n_polys_requested)
    print("@" * x.poly_type)
    println(io)
    print("  Polys with no data: ")
    show(io, x.n_polys_no_data)
    println(io)
    print("  Polys disconnected: ")
    show(io, x.n_polys_dis)
    println(io)
    print("  Counts total: ")
    show(io, x.n_counts_total)
    println(io)
    print("  Counts removed by scc: ")
    show(io, x.n_counts_removed_scc)    
end