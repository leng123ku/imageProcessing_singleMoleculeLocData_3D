function n_beads = searchForBeads(directoryName)
%% searchFOrBeads determines number of beads (n_bead) used for active stage stabillization
        n_beads = 0;
        for n = 1:12
            filePath_bead = [directoryName,'\Bead',int2str(n),'.tif'];
            flag1 = exist(filePath_bead, 'file');
            if flag1 ~= 0
                I = (filePath_bead);
                info = imfinfo(I);
                if n == 1
                    n_beads = n_beads+1;
                    Nimage_bead1 = numel(info);
                end
                if n >= 2
                    Nimage_bead2 = numel(info);
                    if Nimage_bead2 == Nimage_bead1
                        n_beads = n_beads+1;
                    end
                end
            end
        end
end