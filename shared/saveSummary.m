function saveSummary(s)
%% saveSummary saves object s and Matlab's output (command window)
calss_copy_path = [s.imageStacksPathRoot, '\s.mat'];
save(calss_copy_path, 's');
summaryCopyPath = [s.imageStacksPathRoot '\summary.txt'];
diary(summaryCopyPath)
diary OFF