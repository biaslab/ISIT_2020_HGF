macro symmetrical(orig::Expr)
    if orig.args[1].head == :where
        eval(orig)
        mirrored = deepcopy(orig)
        mirrored.args[1] = swap_arguments(orig.args[1])
    elseif orig.args[1].head == :call
        eval(orig)
        mirrored = swap_arguments(orig)
    else
        error("Invalid use of @symmetrical")
    end
    eval(mirrored)
end

function swap_arguments(orig::Expr)
    swap_arg_indexes = Int64[]
    for i=1:length(orig.args[1].args)
        (typeof(orig.args[1].args[i]) == Expr) || continue
        (orig.args[1].args[i].head == :(::)) || continue
        push!(swap_arg_indexes, i)
        (length(swap_arg_indexes) < 2) || break
    end

    mirrored = deepcopy(orig)
    mirrored.args[1].args[swap_arg_indexes] = orig.args[1].args[reverse(swap_arg_indexes)]
    return mirrored
end