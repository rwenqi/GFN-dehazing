function Solver = modelconfig_test(solver_file, save_file)

if exist(save_file, 'file')
%     M_ = load(save_file);
%     fprintf('loading saved model.\n');
    Solver = SolverParser(solver_file, save_file);
else
    Solver = SolverParser(solver_file);
end

Solver.iter = 240000;

if isfield(Solver, 'snapshot_prefix')
    if isfield(Solver, 'iter')
        Solver.state_file = [Solver.snapshot_prefix,sprintf('_iter_%d.solverstate', Solver.iter)];
        Solver.model_file = [Solver.snapshot_prefix, sprintf('_iter_%d.caffemodel', Solver.iter)];
    else
        Solver.state_file = [];
        fprintf('ALERT: no pretrained snapshot found.\n');
    end
end

Solver = caffe_init(Solver, solver_file);
Solver.matfile = save_file;
% Solver.state_path = state_path;

end
