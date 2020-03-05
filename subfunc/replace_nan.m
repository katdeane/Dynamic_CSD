function matrix = replace_nan(matrix)

    matrix(matrix >= 0) = nan;

end


