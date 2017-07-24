
pic_folder = 'YorkUrbanDB\'
%detected_hor_folder = 'hor_tard+EM+orth_validation_york\';
detected_hor_folder = 'result\'
true_hor_folder = 'YorkUrbanDB\';

mkdir(detected_hor_folder);
file = 'YorkUrbanDB\Manhattan_Image_DB_Names.mat';
names = load (file);

%fnames = dir(fullfile(true_hor_folder, '*.mat'));
NTrials = 5;

best_auc = 0;
best_Thr = 15; 
all_dist = [];

iThr = 5:5:30;
for Thr = best_Thr
  
    d_tard = [];
    for trial = 1:NTrials

        for idxPic = 1:size (names.Manhattan_Image_DB_Names, 1)
            s = char(names.Manhattan_Image_DB_Names(idxPic, 1));
            pict = sscanf(s, 'P%d\');
            pict0 = sprintf('P%d', pict);
            pict = sprintf ('P%d.jpg', pict);
    
            orthpoints_file = char(strcat('YorkUrbanDB\',names.Manhattan_Image_DB_Names(idxPic, 1), pict0, 'GroundTruthVP_Orthogonal_CamParams.mat'));
            pic_name = char(strcat('YorkUrbanDB\',names.Manhattan_Image_DB_Names(idxPic, 1), pict));
    
            %pic_name = [pic_folder fnames(idxPic).name(1:8) '.jpg'];
            pic = imread(pic_name);

            width = size(pic, 2);        
            height = size(pic, 1);
            imsize = max(size(pic));

            %load([true_hor_folder fnames(idxPic).name]);
            load ([orthpoints_file]);
            
            if (width > height)
                pp = [307, 251];
            else
                pp = [251, 307];
            end;

            vp_orthogonal(1, 3) = vp_orthogonal(1, 3)/vp_orthogonal(3, 3)*imsize + pp(1);
            vp_orthogonal(1, 2) = vp_orthogonal(1, 2)/vp_orthogonal(3, 2)*imsize + pp(1);
            vp_orthogonal(1, 1) = vp_orthogonal(1, 1)/vp_orthogonal(3, 1)*imsize + pp(1);

            vp_orthogonal(2, 3) = -vp_orthogonal(2, 3)/vp_orthogonal(3, 3)*imsize + pp(2);
            vp_orthogonal(2, 2) = -vp_orthogonal(2, 2)/vp_orthogonal(3, 2)*imsize + pp(2);
            vp_orthogonal(2, 1) = -vp_orthogonal(2, 1)/vp_orthogonal(3, 1)*imsize + pp(2);

            vp_orthogonal = vp_orthogonal';


            %detected_hor_file = sprintf('%sHOR%d.mat', detected_hor_folder, idxPic);
            %load(detected_hor_file);

            hor_tard = main_tardif(pic_name, idxPic, Thr, detected_hor_folder, 1);    
            hor_true = drawHorizont2(vp_orthogonal, width, 'Magenta');

            d_tard(end+1) = distanceBetweenLines( hor_tard, hor_true, width ) / height;    

            filename= sprintf( '%s\\%03dcomparison.jpg', detected_hor_folder, idxPic);
            print('-noui' ,'-r0', '-djpeg',  filename);

            %pause;

        end;
    end;

    distFilename= sprintf( '%s\\%03dTardifYorkerrorsThr%d.mat', detected_hor_folder, idxPic, Thr);
    save (distFilename, 'd_tard');
    mean_d_tard = mean(d_tard);
    std_d_tard = std(d_tard);

    [F, X] = ECDF(d_tard);
    auc = sum(0.5*(F(1:end-1) + F(2:end)) .* (X(2:end)-X(1:end-1))) + (1-max(X))*1;
    
    if (auc > best_auc)
       best_auc = auc;
       best_Thr = Thr;
    end;
    all_dist = [all_dist; d_tard]; 
   
end;

 answer = sprintf( '%s\\TheAnswer.mat', detected_hor_folder);
 save (answer, 'best_Thr', 'best_auc');
  


