function opt = standardOptions()

% Plotting options
opt.plot_level = 0;
opt.plot_mesh = 0;
opt.plot_solution = 1;
opt.FS = 16;
opt.suffix = 'figures/test';
opt.plotGradient = 0;
opt.plotHessian = 0;

% Control parameters for adaptive refinement
opt.refinementType = 'global';
opt.number_of_refinements = 5; % Only relevant for ref_type 1 or 2
opt.refinement_ratio = 0.5; % Only relevant for ref_type 1 or 2

% Options for numerical solution algorithm
opt.compute_reduced = 1;
opt.upwinding = 1;
% opt.dual_mesh = 'voronoi';
opt.dual_mesh = 'centroid';

end % function opt = standardOptions()
