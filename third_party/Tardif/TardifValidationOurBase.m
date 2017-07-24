
%pic_folder = 'EurasianCitiesBase\'
pic_folder = 'I:\11l_tre\VanishingPoints\baserenamed\';
%detected_hor_folder = 'tard_validation_eur\';
detected_hor_folder = 'I:\11l_tre\VanishingPoints\CODE_Tardif\validation';


%true_hor_folder = 'EurasianCitiesBase\';
true_hor_folder = 'I:\11l_tre\VanishingPoints\baserenamed\';
addpath ('../')

mkdir(detected_hor_folder);

%fnames = dir(fullfile(true_hor_folder, '*.mat'));
NTrials = 5;
best_auc = 0;
best_Thr = 5; 
all_dist = [];


BEST_THR = 25; %%%%%!!!!!!!! 4/06/2010 for eurasian !! the same
for Thr = BEST_THR %%%55:5:40
  
    d_tard = [];
    for trial = 1:NTrials

        for idxPic = 1:25
            pict = sprintf('%03d.jpg', idxPic);
            pic_name = [pic_folder pict];
    
            pic = imread(pic_name);
            width = size(pic, 2);        
            height = size(pic, 1);
            imsize = max(size(pic));

            %load([true_hor_folder fnames(idxPic).name]);
            
            %detected_hor_file = sprintf('%sHOR%d.mat', detected_hor_folder, idxPic);
            %load(detected_hor_file);

            hor_tard = main_tardif(pic_name, idxPic, Thr, detected_hor_folder, 0);    
            
            pict = sprintf('%03dhor.mat', idxPic);
            truehorfile = [pic_folder pict];
            st = load (truehorfile);
            %draw true horizon
            
            
            %hor_true = st.horizon; % for eurasian
            %l = st.horizon; % for eurasian
            hor_true = st.l;
            l = st.l;
            
            
            x = 1:width;
            l = l / l(2);
            y = -l(1)*x - l(3);
            plot( x, y, 'Color',  'Magenta', 'LineWidth', 3);
            %plot ([1 600], [700 700], 'Color', 'Magenta');

            %hor_true = drawHorizont2(st.VPs, width, 'Magenta');

            
            d_tard(end+1) = distanceBetweenLines( hor_tard, hor_true, width ) / height;    

            filename= sprintf( '%s\\%03dcomparison.jpg', detected_hor_folder, idxPic);
            print('-noui' ,'-r0', '-djpeg',  filename);

            %pause;

        end;
    end;

    distFilename= sprintf( '%s\\%03derrorsThr%d.mat', detected_hor_folder, idxPic, Thr);
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

 answer = sprintf( '%s\\TheAnswerTardif.mat', detected_hor_folder);
 save (answer, 'best_Thr', 'best_auc');
  


