% parse all fields of a solver proto
% transfer all FC layers in the existing model to CONV layers (equvelent)
function [ Solver ] = SolverParser( solver_def_file, resume_file )
if ~exist(solver_def_file,'file')||isempty(solver_def_file)
    error('Solver definition file %s is not found.',solver_def_file);
end
if ~exist('resume_file','var')
    fprintf('Creating a face lmk solver structure based on %s.', solver_def_file);
    Solver = []; 
else
    load(resume_file);
end
fout = fopen(solver_def_file,'r');
tline = fgetl(fout);
while ischar(tline)
    disp(tline)
    pos = 1;
    while tline(pos)==' '
       pos = pos + 1;
    end 
    if tline(pos) == '#' 
            tline = fgetl(fout);
            continue;
    end 
    ind = find(tline == '"',1);
    if  ~isempty(ind)
        field = tline(ind + 1 : end - 1); 
        ind2 = find(tline == ':',1);
        name = tline(1:ind2-1);
    else
        ind2 = find(tline == ':',1);
        if isempty(ind2)
            error('incorrect format.')
        end 
        ctr = tline(ind2+2:end);
        if isempty(str2num(ctr))
            ctr(isspace(ctr)) = [];
            field = ctr;
        else
            field = str2double(ctr);
        end 
        name = tline(1:ind2-1);
    end 
    Solver = setfield(Solver, name, field);
    tline = fgetl(fout);
end
fclose(fout);
if ~isfield(Solver, 'solver_mode')
    Solver.solver_mode = 'GPU';
end
if ~isfield(Solver, 'device_id') && strcmp(Solver.solver_mode, 'GPU')
        Solver.device_id = 0;
end
if isfield(Solver,'model')
    lnum = length(Solver.model);
    for ind = 1:lnum
        if strcmp(Solver.model(ind).layer_names,sprintf('fc%d',ind))
            Solver.model(ind).layer_names = sprintf('conv%d',ind);
            weights = Solver.model(ind).weights{1};
            [s1,s2] = size(weights);
            [~,~,~,ch] = size(Solver.model(ind-1).weights{1});
            filtersize = sqrt(s1/ch);
            weights = reshape(weights,[filtersize,filtersize,ch,s2]);
            Solver.model(ind).weights{1} = weights;
        end 
    end 
end
end 