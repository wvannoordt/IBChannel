function [cons] = Prims2Cons(glob)
    cons = zeros(glob.nvars, glob.nx+2*glob.nguard, glob.ny+2*glob.nguard);
end