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
