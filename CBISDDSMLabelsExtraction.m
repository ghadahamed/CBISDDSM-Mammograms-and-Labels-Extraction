% clear all; 
% close all; 
% clc;    % Clear the command window

% Define a starting folder.
% start_path = fullfile(CBISDDSMroot, 'E:\MY PHD inshaAllah\Datasets\Data\CBIS DDSM\All');
jpg_path = 'D:\CBIS\tiff2\';
% jpg_path = 'E:\CBIS DDSM\UpdatedVer\Masses_Training\NEW\PNG\';
mkdir(jpg_path, 'CC\');
mkdir(jpg_path, 'MLO\');
roi_path = 'E:\CBIS DDSM\UpdatedVer\Masses_Testing\ROI\';
mkdir(roi_path, 'CC\');
mkdir(roi_path, 'MLO\');
jpg_roi_path = 'E:\CBIS DDSM\UpdatedVer\Masses_Testing\JPGs-ROIs\';
mkdir(jpg_roi_path, 'CC\');
mkdir(jpg_roi_path, 'MLO\');
yolo_annotations_path = 'E:\CBIS DDSM\UpdatedVer\Masses_Testing\YOLO annotations\CBISBoundingBoxes_png.txt';

f_yolo_id = fopen(yolo_annotations_path,'w');
main_path = 'CBIS DDSM\All Mass Test\';
%google_drive_path = 'E:\CBIS DDSM\UpdatedVer\Masses_Testing\AllImgs\';
data_start_path = fullfile('E:\', main_path); %Temp path
d = dir(data_start_path);
isub = [d(:).isdir]; %# returns logical vector
mammosFolders = {d(isub).name}';
mammosFolders(ismember(mammosFolders,{'.','..'})) = [];
%display(mammosFolders);

ROI_cnt = 0;
% N = 28;
% figure;
% nS   = sqrt(N);
% nCol = ceil(nS);
% nRow = nCol - (nCol * nCol - N > nCol - 1);
no_of_folders = 0;
no_of_empty_cases = 0;
no_mammograms = 0;
for k = 1 : length(mammosFolders)
    no_of_folders = no_of_folders + 1;
    % Get list of all subfolders then access the last subfolder containing images.
    case_folder = strcat(main_path, mammosFolders{k});
    full_case_path = fullfile('E:\', case_folder);
    %display(full_case_path);    
    allSubFolders = genpath(full_case_path);
    %display(allSubFolders);
    imgs_folder_path = strsplit(allSubFolders,';'); %=number of subfolders + main folder + 1(since there is a semicolon at the end)
    %display(imgs_folder_path);
    folder_idx = 3;
    empty_case_flag = 0;
    separate_folders_flag = 0;
    if length(imgs_folder_path) < 3 
        display('h1');
        empty_case_flag = 1;
    end
    if length(imgs_folder_path) > 4 %in case there are 2 folders in the main one; one folder for the roi and the other for cropping area
        croppedimgORroi = strsplit(imgs_folder_path{3}, '\');
        croppedimgORroi = croppedimgORroi{6};
        %display(croppedimgORroi(1,3));
        if(croppedimgORroi(1,3) == 'c') %then this is the folder containing the cropping ROI
            folder_idx = 5;
        end
        separate_folders_flag = 1;
    end
    
    imgs_folder_path = (imgs_folder_path{folder_idx});
    imgs_folder_path = strcat(imgs_folder_path, '\');
    %display(imgs_folder_path);
    imgs = dir(imgs_folder_path);
    if any(size(dir([imgs_folder_path '*.dcm' ]),1)) == 0
         empty_case_flag = 1;
    end
    
    if empty_case_flag == 1
       no_of_empty_cases = no_of_empty_cases + 1;
       fprintf('xxxxxxxx %d- CASE PROBLEM %s.\n', no_of_folders, full_case_path);
       case_name = strsplit(full_case_path, '\');
       empty_cases{no_of_empty_cases} = case_name(4);
       display('=========================================================================');
       continue; 
    end
    
    %display(length(imgs));
    longer_w = 0;
    longer_h = 0;
    ROI_flag = 0;
    ROI_id = 0;
    ROI_file_name = '';
    for j = 1 : length(imgs)
        if strcmp(imgs(j).name, '.') || strcmp(imgs(j).name, '..')
            continue;
        end
        img_name = imgs(j).name;
        %display(img_name);
        %info = dicominfo([imgs_folder_path, img_name]);
        %display(info);
        [im] = dicomread([imgs_folder_path, img_name]);
%        disp(im(1001:1006,501:506));
        %im = uint8(255 * mat2gray(im));
        min1=min(min(im));
        max1=max(max(im));
%        disp(min1);
%        disp(max1);
%        im = uint8((255 .* (double(im)-double(min1)) ./ double(max1-min1)));
        im = uint16((1 .* (double(im)) ./ double(1)));
%        disp(im(1001:1006,501:506));
        [h,w] = size(im);
%         imshow(im, []);
        splitted_path = strsplit(imgs_folder_path, '\');
        fullORroi = strsplit(splitted_path{6}, '-');
        fullORroi = strsplit(fullORroi{2}, ' ');
        len = cellfun('length',splitted_path(4));
        view_type = cellfun(@(str) str(len:len), splitted_path(4), 'UniformOutput', false);
        %display(view_type);
        id = strcat(splitted_path(4), '=');
        id = strcat(id, fullORroi(1));
        %display(id);
        id = strcat(id{1}, '.png');
        if(ROI_flag == 0)
            ROI_id = id;
            bin_im = im;
            longer_h = h;
        end
        %fullFileName = fullfile(jpg_path, id);
        %display(fullFileName);
        %imwrite(im, fullFileName);
        if strcmp(fullORroi(1), 'ROI') || (strcmp(fullORroi(1), 'cropped') && folder_idx == 3)
            if(separate_folders_flag == 0 && size(dir([imgs_folder_path '*.dcm' ]),1) ~= 2)
                no_of_empty_cases = no_of_empty_cases + 1;
                fprintf('xxxxxxxx %d- CASE PROBLEM %s.\n', no_of_folders, id);
                case_name = strsplit(full_case_path, '\');
                empty_cases{no_of_empty_cases} = case_name(4);
                display('=========================================================================');
                break;
            end
            if h > longer_h % Then this is the ROI image extract then the mass/classification boundries
                bin_im = im;
                ROI_id = id;
                longer_h = h;
            end
            ROI_flag = 1;
        elseif (strcmp(fullORroi(1), 'cropped') == 0)
            if(strcmp(view_type, 'C') == 1)
                view_type = 'CC\';
            else
                view_type = 'MLO\';
            end
            idWithoutExt = strsplit(id, '.');
            idWithoutExt = idWithoutExt{1};
            idWithoutMass = strsplit(idWithoutExt, '-');
            idWithoutMass = idWithoutMass{2};
            display(idWithoutMass);
            fullFileName = fullfile(jpg_path, view_type, [idWithoutExt,'.tiff']);
            display(fullFileName);
            %imwrite(im, fullFileName, 'Mode','lossless');
            imwrite(im, fullFileName);
            %fullFileName = fullfile(google_drive_path, [idWithoutMass,'.jpg']);
            %imwrite(im, fullFileName);
            no_mammograms = no_mammograms + 1;
        end
    end
    if ROI_flag == 1
        %display(fullORroi(1));
        ROI_cnt = ROI_cnt + 1;
        %imshow(bin_im, []);
        bin_im1= bin_im;
        bin_im = bwconvhull(bin_im);
        [h,w]=size(bin_im);
        stats = regionprops(bin_im,'BoundingBox', 'MajorAxisLength','MinorAxisLength');
        % Show the result
        %subplot(nRow, nCol, ROI_cnt);
        bin_img_id = (strsplit(ROI_id, '='));
        original_im_id = bin_img_id(1);
        len = cellfun('length',original_im_id);
        original_im_id_1 = cellfun(@(str) str(1:len-2), original_im_id, 'UniformOutput', false);
        len = cellfun('length',original_im_id_1);
        view_type = cellfun(@(str) str(len:len), original_im_id_1, 'UniformOutput', false);
        %display(view_type);
        
        if(strcmp(view_type, 'C') == 1)
            view_type = 'CC\';
        else
            view_type = 'MLO\';
        end
        original_im_id = strcat(original_im_id_1, '=');
        original_im_id = strcat(original_im_id, 'full.png');
        original_im_path = fullfile(jpg_path, view_type, original_im_id);
        %display(original_im_path);
        
        original_im_id_1 = strcat(original_im_id_1, '.png');
        im_with_roi_path = fullfile(jpg_roi_path, view_type, original_im_id_1{1});
        %display(im_with_roi_path);
        if exist(im_with_roi_path, 'file') == 0  % This means that this is the second ROI for a specific mammogram
            %display('Not exist');        
            I = imread(original_im_path{1});
            %title(original_im_id);        
            %imshow(I);
            %hold on;
            %rectangle('Position', stats.BoundingBox,'EdgeColor','r');
        else
            %display('exist'); 
            I = imread(im_with_roi_path);
        end
        fullFileName = fullfile(jpg_roi_path, view_type, original_im_id_1{1});
        A = insertShape(I, 'rectangle', stats.BoundingBox,'LineWidth',30,'color','red');
        imwrite(A, fullFileName);

        ROI_cnt = ROI_cnt + 1;
        %subplot(nRow, nCol, ROI_cnt);
        %title(ROI_id);        
        %imshow(bin_im);
        %hold on;
        %rectangle('Position', stats.BoundingBox,'EdgeColor','r');
        fullFileName = fullfile(roi_path, view_type, ROI_id);
        A = insertShape(bin_im1, 'rectangle', stats.BoundingBox,'LineWidth',30,'color','red');
        imwrite(A, fullFileName);
            
            
        %display(stats);    %BoundingBox = xLeft, yTop, width, height]
        xMin = ceil(stats.BoundingBox(1));
        xMax = xMin + stats.BoundingBox(3) - 1;
        yMin = ceil(stats.BoundingBox(2));
        yMax = yMin + stats.BoundingBox(4) - 1;
        % to validate that the obtained coordinate are the correct ones,
        % draw on the image and show
        %rectangle('Position', [xMin, yMin, xMax - xMin, yMax - yMin], 'EdgeColor', 'y');
        x_yolo = (abs(xMax + xMin) / 2)/w;
        y_yolo = (abs(yMax + yMin) / 2)/h;
        width_yolo = abs(xMax - xMin) / w;
        height_yolo = abs(yMax - yMin) / h;

        fprintf(f_yolo_id, '%s ', ROI_id);
        fprintf(f_yolo_id, '%d %d %d %d %d %d', [xMin,yMin,xMax,yMax, h, w]);
        fprintf(f_yolo_id, '\n');
        %return;
    end
    fprintf('%d- PROCESSING DONE %s.\n', no_of_folders, full_case_path);
    display('=========================================================================');
end
fprintf('\n\n# of ALL images processed (=# of folders) = %d\n', no_of_folders);
fprintf('# of mammograms = %d\n', no_mammograms);
fprintf('# of ROI = %d\n', ROI_cnt/2);
fprintf('# of CASES with PROBLEMS = %d\n', no_of_empty_cases);
if no_of_empty_cases > 0
    display('Cases with PROBLEMS are:\n');
    for i = 1:no_of_empty_cases
        %fprintf('%d- %s', i, empty_cases{1:1});
        celldisp(empty_cases{i});
    end
end

