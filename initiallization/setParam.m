classdef setParam
    properties
        rootPath = 'E:\'; %C:\Daily_STORM_Data\  'D:\'  D:\BCells\ 'F:\David\'
        date ='';
        cellDescrip = '';
        runNum = [];
        compOrPlot % 1: for a complete analysis, 2: plotting only
        imageStacksPathRoot
        imageStacksPathByRuns = {};
        % bead tracking calibration file
        %calFileBeadPath = 'D:\dSTORM_v2.0\calFiles_final_selection_ZCorrected_20140716\newtonCal\20140114\Coef.dat';
        calFileBeadPath = 'C:\Users\reza\Google Drive\dSTORM_v2.4\calFiles_final_selection_ZCorrected_20140716\newtonCal\20140114\Coef.dat';
        calFileBead
        % single molecule calibration file
        imagingDepthInfoPath
        imagingDepthInfo %imagingDepthInfo is a cell => {1}: run#, {2}: V0, {3}: V1
        imagingDepthInMicron
        %calFileSMPath_root ='D:\dSTORM_v2.0\calFiles_final_selection_ZCorrected_20140716';
        %calFileSMPath_root ='C:\Users\storm\Google Drive\dSTORM_v2.4\calFiles_final_selection_ZCorrected_20140716';
        calFileSMPath_root = 'C:\Users\reza\Google Drive\dSTORM_v2.4\calFiles_20150817';
        calFileSMPath = {};
        calFileSM_Rxy_Z = {};
        calFileSM_Z_Wxy = {};
        % ViSP path
        vispPath = '';
        % pixel sizes
        CCDpixely_SM = 106.1;%109.0;%106.3830; %(nm)   left -> right
        CCDpixelx_SM = 106.1;%107.2;%106.1;%102.5641; %(nm) up-> down (x for matlab coloum)
        CCDpixely_beadCam = 96.1538; %(nm) up-> down (x for matlab coloum)
        CCDpixelx_beadCam = 107.3826; %(nm)   left -> right
        % parallel pool workers
        n_workers = 31;
    end
    
    properties(Constant)
        conversionFac_v_to_um = 3.75;
    end
    
    
    
    methods
        function s = setParam(date, cellDescrip, runNum, compOrPlot)
            s.date = date;
            s.cellDescrip = cellDescrip;
            cprintf('*red', 'Cell''s info: %s, Imaged on: %s\n\n', s.cellDescrip, s.date);
            s.runNum = runNum;
            s.compOrPlot = compOrPlot;
            
            s.imageStacksPathRoot = fullfile(s.rootPath, s.date);
            s.check_exist(s.imageStacksPathRoot, 'dir');
            
            s.imagingDepthInfoPath = fullfile(s.imageStacksPathRoot, 'depth_info.txt');
            s.check_exist(s.imagingDepthInfoPath, 'file');
            
            s.calFileBead = load(s.calFileBeadPath);
            s.check_exist(s.calFileBeadPath, 'file');
            
            s.vispPath = fullfile(s.rootPath, date, 'VISP');
            create_vispFolder(s);
            
            fprintf('\n');
            check_exist_otherFiles(s)
            
        end
        
        function imageStacksPathByRuns = get.imageStacksPathByRuns(s)
            for i = 1:numel(s.runNum)
                imageStacksPathByRuns{i} = fullfile( s.imageStacksPathRoot, num2str(s.runNum(i)) );
            end
        end
        
        
        function imagingDepthInfo = get.imagingDepthInfo(s)
            file_id = fopen(s.imagingDepthInfoPath);
            if file_id == -1
                error('Unable to open depth_info.txt => check s.imagingDepthInfoPath');
            else
                % imagingDepthInfo is a cell => {1}: run#, {2}: V0, {3}: V1
                imagingDepthInfo = textscan(file_id, '%u %f %f', -1);
                fclose(file_id);
            end
        end
        
        function imagingDepthInMicron = get.imagingDepthInMicron(s)
            imagingDepthInMicron = zeros(numel(s.runNum),1);
            for i = 1:numel(s.runNum)
                imagingDepthInMicron(i,1) = round( (s.imagingDepthInfo{2}(s.runNum(i)) - s.imagingDepthInfo{3}(s.runNum(i))) / s.conversionFac_v_to_um);
            end
        end
        
        function calFileSMPath = get.calFileSMPath(s)
            %calFileSMPath = cell(s.runNum,1);
            for i = 1:numel(s.runNum)
                str = [num2str(s.imagingDepthInMicron(i)) 'um'];
                calFileSMPath{i} = fullfile(s.calFileSMPath_root, str);
            end
        end
        
        function calFileSM_Rxy_Z = get.calFileSM_Rxy_Z(s)
            %calFileSM_Rxy_Z = cell(s.runNum,1);
            for i = 1:numel(s.runNum) 
                %Rxy_Z_table = load(fullfile(s.calFileSMPath{i}, 'Rxy_Z_polyFit_wholeRange.dat'));
                Rxy_Z_table = load(fullfile(s.calFileSMPath{i}, 'Rxy_Z_curveFit_appRange.dat'));
                calFileSM_Rxy_Z{i} = Rxy_Z_table;
            end
        end
        
        function calFileSM_Z_Wxy = get.calFileSM_Z_Wxy(s)
            %calFileSM_Z_Wxy = cell(s.runNum,1);
            for i = 1:numel(s.runNum)
                ZWxy_table = load(fullfile(s.calFileSMPath{i}, '\ZWxy.dat'));
                calFileSM_Z_Wxy{i} = ZWxy_table;
            end
        end
        
        
        function check_exist_otherFiles(s)
            for i = 1:numel(s.runNum)

                PreSTORM_path = fullfile(s.imageStacksPathByRuns{i}, 'PreSTORM.tif');
                drift_path = fullfile(s.imageStacksPathByRuns{i}, 'drift.mat');
                newtonTime_path = fullfile(s.imageStacksPathByRuns{i}, 'newtonTime.mat');
                iXonTime_path = fullfile(s.imageStacksPathByRuns{i}, 'iXonTime.mat');
                
                s.check_exist(PreSTORM_path, 'file');
                s.check_exist(drift_path, 'file');
                s.check_exist(newtonTime_path, 'file');
                s.check_exist(iXonTime_path, 'file');
                fprintf('\n');

            end
        end
        
        function check_exist(s, path, fileOrDir)
            switch fileOrDir
                case 'file'
                    exist_id = 2;
                case 'dir'
                    exist_id = 7;
                otherwise
                    error('fileOrDir should be either file or dir');
            end
                
            if exist(path, 'file') ~= exist_id
                    error('%s file could not be found', path);
            else
                    cprintf('*green', 'EXISTS   ');
                    cprintf('*String', '%s\n', path); 
            end
        end
        
        
        function creat_copy_calibFiles(s)
            calFilesPath_copy = fullfile(s.imageStacksPathRoot, 'callibrationFiles');
            if exist(calFilesPath_copy, 'dir') ~= 7
                parentFolder = s.imageStacksPathRoot;
                mkdir(parentFolder, 'callibrationFiles');
                
                [status, message] = copyfile(s.calFileSMPath_root, calFilesPath_copy, 'f');
                if ~status
                    error(message);
                else
                    cprintf('*String', 'A copy of calibarion files was made in the destination path. \n');
                end
            else
                cprintf('*String', 'Did not make a copy of calibarion file in the destination path => A copy exists already. \n');
            end
            
        end
        
        
        function create_vispFolder(s)
            if exist(s.vispPath, 'dir') ~= 7
                parentFolder = fullfile(s.rootPath, s.date);
                mkdir(parentFolder, 'VISP');
            end
            
        end
        
    end
    
end