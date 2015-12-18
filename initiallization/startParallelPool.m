function startParallelPool(n_workers)

if ~nargin
    n_workers = 31;
end

%% Creating a parralel cluster, modifying the local profile
cluster_obj = parcluster('local');
cluster_obj.NumWorkers = n_workers;  % modifying NumWorkers property of cluster_obj, max is 12.
saveProfile(cluster_obj);    % 'local' profile now updated. 

%% Shutting down the currently open parallel pool
cprintf('*String', 'Shutting down the currently open parallel pool. \n');
% returns the current pool if one exists. 
% If no pool exists, the 'nocreate' option causes gcp not to create a pool
% regardless of your parallel preferences settings.
current_pool = gcp('nocreate');  
if isempty(current_pool)
    fprintf('No open parallel pool detected. \n');
else
    delete(current_pool);
    current_pool = gcp('nocreate');
    if isempty(current_pool)
        cprintf('*String', 'The already open parallel pool was shutdown successfuly. \n');
    end
end

%% Open the new pool with n_workers on cluster_obj
cprintf('*String', 'Opening a new pool with %g workers. \n', n_workers);
parpool(cluster_obj, n_workers);
% status check
current_pool = gcp('nocreate');
if isempty(current_pool)
    error('Error in creating the pool. \n');
elseif ~current_pool.Connected
    error('Pool created but unable to connect to. \n');
elseif current_pool.NumWorkers ~= n_workers
    error('Requested to connect to %g workers but currently connected to %g workers. \n\n', n_workers, cuurent_pool.NumWorkers);
end
 