function X = caffe_init(X, solver_file)

if ~exist(solver_file,'file') || isempty(solver_file)
    error('You need a solver prototxt definition.')
end
X.Solver_ = caffe.Solver(solver_file); % to cpp
if isfield(X, 'state_file') && ~isempty(X.state_file) && exist(X.state_file, 'file')
    X.Solver_.restore(X.state_file);
    fprintf('resume from snapshot file %s...\n', X.state_file);
elseif isfield(X, 'model_file') && ~isempty(X.model_file)
    X.Solver_.net.copy_from(X.model_file);
    fprintf('resume from caffemodel file %s...\n', X.model_file);
end

% set caffe mode
if strcmp(X.solver_mode,'GPU')
    caffe.set_mode_gpu();
    caffe.set_device(X.device_id);
else
    caffe.set_mode_cpu();
end
fprintf('Done with init %s\n', solver_file);
end
