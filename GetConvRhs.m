function [convRhs] = GetConvRhs(glob)
    convRhs = zeros(glob.nvars, glob.nx+2*glob.nguard, glob.ny+2*glob.nguard);
end