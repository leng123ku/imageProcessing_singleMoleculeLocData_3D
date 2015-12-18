for k = 1:num_setOfdate
    this_date = date(k,:);
    this_cellDescrip = cellDescrip(k,:);
    this_runNum_temp = runNum(k,:);
    this_runNum = this_runNum_temp (this_runNum_temp ~= 0);

    %% user's setting => creating s(setParam) object
    % setParam(date, cellDescrip, runNUm, compOrPlot)
    s = setParam(this_date, this_cellDescrip, this_runNum, 1);
end