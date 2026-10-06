function matrix = C_function_3D(phi,n,N)
% Compatibility wrapper for full symbolic reference assembly.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','common'));
matrix = ite_symbolic_matrix('C',3,phi,[],n,N);
end
