function [ghostswall, ghostscenter] = GetGhostValues(glob)
    ghostswall   = zeros(glob.nvars, glob.nx+2*glob.nguard, 2);
    ghostscenter = zeros(glob.nvars, glob.nx+2*glob.nguard, 2);
    
end

