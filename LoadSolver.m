function X = LoadSolver(folder, solver_file, mat_file)

if exist(mat_file,'file') && ~isempty(mat_file)
    M_ = load(mat_file);
    X = T2_SolverParser(solver_file, M_); % TODO: save G_
else
    X = T2_SolverParser(solver_file);
end

if isfield(X,'snapshot_prefix') && exist('M_','var')
    [~,pr] = fileparts(X.snapshot_prefix);
    if isfield(X, 'iter')
        % saved path
        state = fullfile(folder, [pr,sprintf('_iter_%d.solverstate', X.iter)]);
        % tmp path
        X.state_file = [X.snapshot_prefix,sprintf('_iter_%d.solverstate', X.iter)];
        copyfile(state, X.state_file);
        % saved path
        modes = fullfile(folder, [pr,sprintf('_iter_%d.caffemodel', X.iter)]);
        % tmp path
        X.model_file = [X.snapshot_prefix, sprintf('_iter_%d.caffemodel', Solver.G.iter)];
        copyfile(modes, X.model_file);
    end
else
    % dump existing models if there is no solver mat info
    delete(fullfile(folder, [pr,'*']));
    X.state_file = [];
    fprintf('ALERT: please copy the pretrained caffemodel to %s.\n', X.model_file);
end

X = SR_init02(X, solver_file);

end