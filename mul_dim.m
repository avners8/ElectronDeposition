function A_mul = mul_dim(A,v,dim)
    
    % This function multiply (and sum) a matrix with a vector, according to
    % the specified dimension.
    dim_permute = 1:ndims(A); dim_permute(1) = dim; dim_permute(dim) = 1;
    A_permute = permute(A, dim_permute);
    size_A_permute = size(A_permute);
    v = repmat(v, [1, size_A_permute(2:end)]);
    A_mul = A_permute .* v;
    A_mul = permute(A_mul, dim_permute);
    A_mul = squeeze(sum(A_mul, dim));
    
end