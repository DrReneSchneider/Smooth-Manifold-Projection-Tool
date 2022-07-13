clear all
close all

% %% User interface: ask for folder location - MULTI %%
% Directory = uigetdir(pwd, 'Select a folder');
% Files = dir(fullfile(Directory, '*w2mCHERRY*'));
% % Display the names
% % Files.name;
% PathName = [Directory '/'];

% User interface: ask for file location - SINGLE %%
[FileName PathName] = uigetfile('*.TIF', 'Select file that contains "surface" information');
Files(1).name = FileName;
prompt = {'Direction of z-stack (IN or OUT):', 'Enter envelope stiffness [pixels]:', 'Enter final filter size [pixels]:', 'Enter number of ADDITIONAL stacks to be z-smoothed [e.g. 0, 1, 2]', 'Offset: N planes above (+) or below (-) blanket [pixels]', 'Depth: MIP for N pixels into blanket [pixels]'};
title = 'Parameter input';
dims = [1 35];
definput = {'OUT' '30', '30', '0', '0', '3'};
answer = inputdlg(prompt, title, dims, definput);
folder_save = [answer{1,1} '_' answer{2,1} 'px_' answer{3,1} 'px_' answer{4,1} ' additional channel(s)_z-offset ' answer{5,1} '_MIP-depth from SMP = ' answer{6,1} '/'];
mkdir([PathName folder_save]);

%% Process 1st file
for m = 1 : length(Files)
    FileName = Files(m).name;
    info = imfinfo([PathName '/' FileName]);
    num_images = numel(info);
    
    % Import file, i.e. all k frames of the stack %%
    Image1_raw = zeros(info(1).Height, info(1).Width);
    for k = 1:1:num_images
        I = imread([PathName FileName], k);
        
        %Image1_raw(:,:,k)= medfilt2(I, [3 3]);
        Image1_raw(:,:,k) = I;
    end
    [rows, columns, nSlices] = size(Image1_raw);
    I = round(im2double(I)*2^16);
    
    % Maximum Intensity Projection (MIP)
    MIP_Image = max(Image1_raw, [],3);
    imwrite(uint16(MIP_Image), [PathName folder_save num2str(FileName(1:end-4)) '_MIP.TIF']);

    % MIP z-map
    MIP_zmap = zeros(rows, columns, class(Image1_raw));
    for col = 1 : columns
        for row = 1 : rows
            z_vector = reshape(Image1_raw(row, col, :),[1, num_images]);
            locs_MIP = find(z_vector==max(z_vector));
            MIP_zmap(row, col) = locs_MIP(1);
        end
    end
%     %Export accurate MIP zmap
    imwrite(uint16(MIP_zmap), [PathName folder_save num2str(FileName(1:end-4)) '_MIP_zmap.TIF']);

    %Smooth Manifold Projection (SMP)
    SMP_Image = zeros(rows, columns, class(Image1_raw));
    SMP_zmap = zeros(rows, columns, class(Image1_raw));
    
%     %Smooth MIP zmap with 5-by-5 pixel area
%     MIP_zmap = medfilt2(MIP_zmap, [10 10]);
    
    %Place a smooth sheet over the zmap data
    [Env_1_up Env_1_low] = envelope(MIP_zmap, str2num(answer{2,1}), 'peak');
    [Env_2_up Env_2_low] = envelope(MIP_zmap', str2num(answer{2,1}), 'peak');

    Env_1_up = round(Env_1_up);
    Env_2_up = round(Env_2_up)';
    Env_1_low = round(Env_1_low);
    Env_2_low = round(Env_2_low)';
    
    %z-stack direction (into/out of) tissue
    if answer{1,1}(1,1) == 'O'
        Env_max = zeros(rows, columns, 2);
        Env_max(:, :, 1) = Env_1_up;
        Env_max(:, :, 2) = Env_2_up;
        Env_max = max(Env_max, [], 3);
    elseif answer{1,1}(1,1) == 'I'
        Env_max = zeros(rows, columns, 2);
        Env_max(:, :, 1) = Env_1_low;
        Env_max(:, :, 2) = Env_2_low;
        Env_max = min(Env_max, [], 3);
    end
    
    %Get rid of extreme values and replace them with lowest/highest slice
    outlier_high = [];
    outlier_high = find(Env_max > nSlices);
    Env_max(outlier_high) = nSlices;
    outlier_low = [];
    outlier_low = find(Env_max < 1);
    Env_max(outlier_low) = 1;
    Env_filt = ceil(medfilt2(Env_max, [str2num(answer{3,1}) str2num(answer{3,1})]));
        
    %Use smooth zmap to extract intensity from z-stack
    %Add planes above and below (=fuzziness) extracted z-position to allow for error
    %in finding ideal z-position
    zplanes = 1:1:nSlices;
    zmatrix = zeros(rows, columns, nSlices);
    offset = str2num(answer{5,1});
    fuzziness = 1;
    for i = 1:size(Image1_raw, 1)
        for j = 1:size(Image1_raw,2)
            ind_low = zplanes >= Env_filt(i,j)-fuzziness+offset;
            ind_high = zplanes <= Env_filt(i,j)+fuzziness+offset;
            ind = ind_low+ind_high-1;
            zmatrix(i, j, :) = ind;
        end
    end

%     %Save SMP to TIFF file
    Image1_SMP = Image1_raw.*zmatrix;
    SMP_Image = max(Image1_SMP, [], 3);
    if abs(offset) > 0
        imwrite(uint16(SMP_Image), [PathName folder_save num2str(FileName(1:end-4)) '_SMP_Offset(' num2str(offset) ').TIF']);
    else
        imwrite(uint16(SMP_Image), [PathName folder_save num2str(FileName(1:end-4)) '_SMP.TIF']);
    end
    
    %Find smooth z-map for z-stack and save to TIFF file
    SMP_zmap_opt = zeros(rows, columns, class(Image1_raw));
    for col = 1 : columns
        for row = 1 : rows
            z_vector = [];
            z_vector = reshape(Image1_SMP(row, col, :),[1, num_images]);
            z_vector = z_vector;
            if std(z_vector) == 0
                %use locs_MIP value from measurement before (neighboring
                %pixel)
                SMP_zmap_opt(row, col) = locs_MIP(1);
            else
                locs_MIP = [];
                locs_MIP = find(z_vector == max(z_vector));
            end
            SMP_zmap_opt(row, col) = locs_MIP(1);
        end
    end
    imwrite(uint16(SMP_zmap_opt), [PathName folder_save num2str(FileName(1:end-4)) '_SMP_zmap.TIF']);
    
    %Use optimized z-map to extract MIP-projection in a certain range of
    %z-levels
    MIPbased_value = str2double(answer{6,1});
    while abs(MIPbased_value) > 0
    depth = str2num(answer{6,1});
    ind2 = zeros(1, num_images);
    zmatrix2 = zeros(rows, columns, nSlices);
    for col = 1 : columns
        for row = 1 : rows
            ind2 = zeros(1, num_images);
            if answer{1,1}(1,1) == 'O'
                max_z = SMP_zmap_opt(row, col);
                min_z = max_z - depth;
                if min_z <= 1
                    min_z = 1;
                end
            elseif answer{1,1}(1,1) == 'I'   
                min_z = SMP_zmap_opt(row, col);
                max_z = min_z + depth;
                if max_z >= nSlices
                    max_z = nSlices-1;
                end
            end
            if MIPbased_value >= 0
                ind2(min_z:max_z) = 1;
            else
                ind2(max_z:min_z) = 1;
            end
            zmatrix2(row, col, :) = ind2;
        end
    end
    
    %Save SMP-based MIP to TIFF file
    Image1_SMPMIP = Image1_raw.*zmatrix2;
    SMP_MIP_Image1 = max(Image1_SMPMIP, [], 3);
    imwrite(uint16(SMP_MIP_Image1), [PathName folder_save num2str(FileName(1:end-4)) '_SMP-based (' answer{6,1} ') MIP.TIF']);

    %Find z-map for SMP-based MIP and save to TIFF file
    SMP_zmap_opt2 = zeros(rows, columns, class(Image1_raw));
    for col = 1 : columns
        for row = 1 : rows
            z_vector2 = [];
            z_vector2 = reshape(Image1_SMPMIP(row, col, :),[1, num_images]);
            z_vector2 = z_vector2;
            locs_MIP2 = [];
            locs_MIP2 = find(z_vector2==max(z_vector2));
            SMP_zmap_opt2(row, col) = locs_MIP2(1);
        end
    end
    imwrite(uint16(SMP_zmap_opt2), [PathName folder_save num2str(FileName(1:end-4)) '_SMP-based (' answer{6,1} ') MIP_zmap.TIF']);
    
    MIPbased_value = 0;
    end
    
%% (OPTIONAL) Process 2nd file    
    %Use SMP zmap to smooth another file
    if str2num(answer{4,1}(1,1)) > 0
        %% Use this block if your files are named MTs.. and PMs...
%         if mean(FileName(1:2) == 'MT') == 1
%             FileName2 = ['PMs' FileName(4:end)];
%         elseif mean(FileName(1:2) == 'PM') == 1
%             FileName2 = ['MTs' FileName(4:end)];
%         end
        
%         %% Use this block if you want this script to select the second file automatically
%         FileName2 = strrep(FileName, 'w1GFP', 'w2mCHERRY');
%         PathName2 = PathName;
        
        %% Use this block if you want to manually choose the second file

        [FileName2,PathName2] = uigetfile('*.tif','Select file to which SMP should be applied ...');
        info2 = imfinfo([PathName2 FileName2]);
        num_images2 = numel(info2);
        
        %Import file, i.e. all l frames of the stack %%
        Image2_raw = [];
        for l = 1:1:num_images
            I2 = imread([PathName2 FileName2], l);
            Image2_raw(:,:,l)=I2;
        end
        I2 = round(im2double(I2)*2^16);

        %Maximum Intensity Projection (MIP) for 2nd file
        MIP_Image2 = max(Image2_raw, [],3);
        imwrite(uint16(MIP_Image2), [PathName2 folder_save num2str(FileName2(1:end-4)) '_MIP.TIF']);
        
        %Apply SMP and save to TIFF file
        Image2_SMP = Image2_raw.*zmatrix;
        SMP_Image_2 = max(Image2_SMP, [], 3);
        if abs(offset) > 0
            imwrite(uint16(SMP_Image_2), [PathName2 folder_save num2str(FileName2(1:end-4)) '_SMP_Offset(' num2str(offset) ').TIF']);
        else
            imwrite(uint16(SMP_Image_2), [PathName2 folder_save num2str(FileName2(1:end-4)) '_SMP.TIF']);
        end
        
        %Find smooth z-map for 2nd file
        SMP_zmap_opt_2 = SMP_zmap_opt;
        imwrite(uint16(SMP_zmap_opt_2), [PathName2 folder_save num2str(FileName2(1:end-4)) '_SMP_zmap.TIF']);
        
        %Calculate MIP z-map of 2nd file and save to TIFF file
        MIP_zmap_2 = zeros(rows, columns, class(Image2_raw));
        for col = 1 : columns
            for row = 1 : rows
                z_vector = [];
                z_vector = reshape(Image2_raw(row, col, :),[1, num_images]);
                z_vector_double = z_vector;

                locs_MIP = [];
                locs_MIP = find(z_vector==max(z_vector));
                MIP_zmap_2(row, col) = locs_MIP(1);
            end
        end
        imwrite(uint16(MIP_zmap_2), [PathName2 folder_save num2str(FileName2(1:end-4)) '_MIP_zmap.TIF']);
        
        MIPbased_value = str2double(answer{6,1});
        while MIPbased_value > 0
        
        %Save SMP-based MIP to TIFF file
        Image2_SMPMIP = Image2_raw.*zmatrix2;
        SMP_MIP_Image2 = max(Image2_SMPMIP, [], 3);
        imwrite(uint16(SMP_MIP_Image2), [PathName folder_save num2str(FileName2(1:end-4)) '_SMP-based (' answer{6,1} ') MIP.TIF']);

        %Find z-map for SMP-based MIP and save to TIFF file
        SMP_zmap_opt2 = zeros(rows, columns, class(Image2_raw));
        for col = 1 : columns
            for row = 1 : rows
                z_vector2 = [];
                z_vector2 = reshape(Image2_SMPMIP(row, col, :),[1, num_images]);
                z_vector2 = z_vector2;
                locs_MIP2 = [];
                locs_MIP2 = find(z_vector2==max(z_vector2));
                SMP_zmap_opt2(row, col) = locs_MIP2(1);
            end
        end
        imwrite(uint16(SMP_zmap_opt2), [PathName folder_save num2str(FileName2(1:end-4)) '_SMP-based (' answer{6,1} ') MIP_zmap.TIF']);
        
        MIPbased_value = 0;
        end
    end
    
%% (OPTIONAL) Use SMP to z-smooth yet another file
    if str2num(answer{4,1}(1,1)) > 1
        [FileName3,PathName3] = uigetfile('*.tif','Select another file to which SMP should be applied ...');
        info3 = imfinfo([PathName3 '/' FileName3]);
        num_images3 = numel(info3);
        
        %Import file, i.e. all l frames of the stack %%
        Image3_raw = [];
        for l = 1:1:num_images
            I3 = imread([PathName3 FileName3], l);
            Image3_raw(:,:,l)=I3;
        end
        I3 = round(im2double(I3)*2^16);

        %Maximum Intensity Projection (MIP) of 3rd file
        MIP_Image3 = max(Image3_raw, [],3);
        imwrite(uint16(MIP_Image3), [PathName3 folder_save num2str(FileName3(1:end-4)) '_MIP.TIF']);
        
        %Apply SMP and save to TIFF file
        Image3_SMP = Image3_raw.*zmatrix;
        SMP_Image_3 = max(Image3_SMP, [], 3);
        imwrite(uint16(SMP_Image_3), [PathName3 folder_save num2str(FileName3(1:end-4)) '_SMP.TIF']);

        %Find smooth z-map for third channel
        SMP_zmap_opt_3 = SMP_zmap_opt;
        imwrite(uint16(SMP_zmap_opt_3), [PathName3 folder_save num2str(FileName3(1:end-4)) '_SMP_zmap.TIF']);
        
        %MIP z-map of second file
        MIP_zmap_3 = zeros(rows, columns, class(Image1_raw));
        for col = 1 : columns
            for row = 1 : rows
                z_vector = [];
                z_vector = reshape(Image3_raw(row, col, :),[1, num_images3]);
                z_vector_double = z_vector;

                locs_MIP = [];
                locs_MIP = find(z_vector==max(z_vector));
                MIP_zmap_3(row, col) = locs_MIP(1);
            end
        end
        imwrite(uint16(MIP_zmap_3), [PathName3 folder_save num2str(FileName3(1:end-4)) '_MIP_zmap.TIF']);
    end
end
%save([PathName folder_save 'workspace.mat']);