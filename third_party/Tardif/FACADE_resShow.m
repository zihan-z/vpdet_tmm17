function FACADE_resShow(RES, testName, vIdx)
    
    if nargin < 2
        FCprintf('Provide a testName\n');
        return
    end
    
    if nargin >= 3
        %sub results
        RES.cpu = RES.cpu(vIdx,:);
        RES.f123 = RES.f123(vIdx,:);
        RES.f12 = RES.f12(vIdx,:);
        RES.algo = {RES.algo{vIdx}};
        
        
        tmp = reshape({RES.imgGtErr{vIdx,:}}, length(vIdx), size(RES.imgGtErr,2) );
        RES.imgGtErr = tmp;
    end
    
%get # of approach
    nbAlgo = length(RES.algo)
    %-------------------------------------------------------
    %focal length
    fprintf('----------------------Focal-----------------------\n')
    analyseF(RES, testName,  200, nbAlgo);

    fprintf('------------------classification------------------\n')
    analyse_classErr(RES, testName, 400, nbAlgo);
    %-------------------------------------------------------    
    %fitting error
    fprintf('-------------------Fitting error------------------\n')
    analyse_imgGtErr(RES,testName, 203, nbAlgo);

    fprintf('------------------computation time----------------\n')
    analyse_cpu(RES, testName, 300, nbAlgo);
  
return

%-------------------------------------------------------    
%-------------------------------------------------------    
%-------------------------------------------------------    
function vb = isValidF(vF)    
vb = vF>100 & vF<1500 & imag(vF)==0;
vb = imag(vF)==0;
%vb = vF>300 & vF<1000 & imag(vF)==0;
%vb = vF>500 & vF<800 & imag(vF)==0;
    return
%-------------------------------------------------------    
%-------------------------------------------------------    
%-------------------------------------------------------    
function analyse_cpu(RES, testName, fig, nbAlgo)
      
    
    clrMap = varycolor(nbAlgo); %FIXME
    clrMap = clrMap(2:end,:);
    vCpu = RES.cpu(:,end:-1:1)';
    vCpu = vCpu(:,2:end);
    
    Y = [min(vCpu); mean(vCpu); max(vCpu)];
    sfigure(fig)
    bh = bar(Y);
    for i=1:length(bh)
        set(bh(i),'facecolor',clrMap(i,:) )
    end
    
    legend();
    set(gca,'XTickLabel',{'Min.', 'Mean', 'Max.'})
    daspect([1,0.9,1]);
    filename = sprintf('fig/%s-computationTime.eps', testName);
    
    FACADE_printGraphics('',...
                         '', 'Compution time (sec.)',...
                         {RES.algo{2:end}},...
                         filename,...
                         'NorthWest')
    
%-------------------------------------------------------    
%-------------------------------------------------------    
%-------------------------------------------------------    
function analyse_classErr(RES, testName, fig, nbAlgo)

clrMap = hsv(nbAlgo); %FIXME
clrMap = varycolor(nbAlgo); %FIXME

% $$$ leg = {};
% $$$ hold on
% $$$ ii=1;
% $$$ for idx=1:nbAlgo
% $$$   vF123 =real(RES.f123(idx,:));
% $$$   vF12  =real(RES.f12(idx,:));
% $$$   
% $$$   vb123 = isValidF(vF123);
% $$$   vb12 = isValidF(vF12);
% $$$   vb = vb123 | vb12; %either one is good
% $$$ 
% $$$   vT = 1:length(vF123);
% $$$   
% $$$   %vErrF123(ii) = mean(abs(vF123(vb) - fGroundTruth));
% $$$   %vErrF12(ii) = mean(abs(vF12(vb)  - fGroundTruth));
% $$$   vFboth = [vF123;vF12];
% $$$   [vFer,vIdx] = min(abs(vFboth - fGroundTruth));
% $$$   vF = vFboth(1,:);
% $$$   vF(vIdx==2) = vFboth(2,vIdx==2);
% $$$   vFerr{ii}    = abs(vF(vb)-fGroundTruth);
% $$$   vFerr123{ii} = abs(vF123(vb123)-fGroundTruth);
% $$$   vErrF(ii) = mean(vFerr{ii});
% $$$   
% $$$   %assert(all(vFerr123{ii} >= vFerr{ii}));
% $$$   
% $$$   vNbF(ii)  = sum(vb)/length(vb);
% $$$   plot(vT(vb), vF(vb), '.-', 'Color', clrMap(ii,:) );
% $$$   ii=ii+1;
% $$$ end
% $$$ legend(RES.algo);
% $$$ hold off
% $$$ 
% $$$ disp([vErrF; vNbF])

if false
  v=0:300;
else
  v=0:0.05:1;
end

invfigure(fig+1); clf
hold on
for idx =1:nbAlgo
  vErr = real(RES.class_InOut(idx,:));
  vV = arrayfun(@(x) sum(vErr <= x)/length(vErr), v);
  plot(v, vV, '-', 'Color', clrMap(idx,:), 'LineStyle', linestyle(idx), 'LineWidth', 1.5 );
end
daspect([v(end)/2 1 1]);
hold off
filename = sprintf('fig/%s-cumulHistoInOut.eps', testName);
FACADE_printGraphics('', 'Classification error (%)', '', RES.algo, filename, 'SouthEast')


% $$$ invfigure(fig+2); clf
% $$$ hold on
% $$$ for idx =1:length(vFerr123)
% $$$   vErr = vFerr123{idx};
% $$$   vV = arrayfun(@(x) sum(vErr < x)/nbEdges, v);
% $$$   plot(v, vV, '-', 'Color', clrMap(idx,:), 'LineStyle', linestyle(idx), 'LineWidth', 1.5 );
% $$$ end
% $$$ daspect([v(end)/2 1 1]);
% $$$ hold off
% $$$ filename = sprintf('fig/%s-cumulHistoF123.eps', testName);
% $$$ FACADE_printGraphics('', 'Focal Length error (pixel)', '', RES.algo, filename, 'SouthEast')

return

%-------------------------------------------------------    
%-------------------------------------------------------    
%-------------------------------------------------------    
function analyse_imgGtErr(RES, testName, fig, nbAlgo)

    clrMap = hsv(nbAlgo); %FIXME
    clrMap = varycolor(nbAlgo); %FIXME
    
    
    sfigure(fig); clf;
    sfigure(fig+1); clf;
    sfigure(fig+2); clf;
    mx = 3*ones(3,1);
    leg = {};
    for i=1:nbAlgo
        sol = [RES.imgGtErr{i, :}] ; 

        mnErr = zeros(1,3);
        mxErr = zeros(1,3);
        vIdx = 1:size(sol,2);
        for j=1:3
            vb = sol(j,:)<inf;
            if sum(vb)==0, continue; end;
            mx(j) = max(mx(j), max(sol(j,vb)));
            sfigure(fig+j-1); 
            %subplot(3,1,j); 
            hold on
            plot(vIdx(vb), sol(j,vb), '.-', 'Color', clrMap(i,:) );
            hold off
        
            mnErr(j)  = mean(sol(j,vb));
            mxErr(j) = max(sol(j,vb));
        end
        fprintf('%s    \t  & %0.1f (%0.1f)\t& %0.1f (%0.1f)\t& %0.1f (%0.1f) \t \\\\ \n',...
                RES.algo{i},...
                mnErr(1), mxErr(1),...
                mnErr(2), mxErr(2),...
                mnErr(3), mxErr(3) );
        
    end 
    sfigure(fig); 
    for j=1:3
        subplot(3,1,j); hold on
        %legend(RES.algo);
        axis([0,length(sol),0, min(mx(j),10)  ]);
        hold off
    end

    %cumulative histogram    
    if true
        v=0:0.3:10;
        X=4
    else
        v=0:0.3:7;
        X = 3;
    end
    
    invfigure(fig+3); clf
    invfigure(fig+4); clf
    invfigure(fig+5); clf
    
    gt =  [RES.imgGtErr{1,:}];
    vb = gt < inf;
    
    for i=1:nbAlgo
        sol = [RES.imgGtErr{i,:}]; 
        
        vIdx = 1:size(sol,2);
        for j=1:3
            b = vb(j,:);
            vErr = sol(j,b);
            vV = arrayfun(@(x) sum(vErr <= x)/length(vErr), v);
            invfigure(fig+2+j);; %subplot(3,1,j); 
            hold on
            plot(v, vV, '-', 'Color', clrMap(i,:), 'LineStyle', linestyle(i), 'LineWidth', 1.5 );
            hold off
        end
    end
    
    lgLocation = 'SouthEast'; %'SouthEast'
    invfigure(fig+3);

    daspect([X 1 1]);
    filename = sprintf('fig/%s-cumulHistoVP-1.eps', testName);
    FACADE_printGraphics('', 'Consistency error (pixel)', '', RES.algo, filename, lgLocation)
    
    invfigure(fig+4);
    daspect([X 1 1]);        
    filename = sprintf('fig/%s-cumulHistoVP-2.eps', testName);
    FACADE_printGraphics('', 'Consistenty error (pixel)', '', {}, filename, lgLocation)
    
    invfigure(fig+5);
    daspect([X 1 1]);
    filename = sprintf('fig/%s-cumulHistoVP-3.eps', testName);
    FACADE_printGraphics('', 'Consistency error (pixel)', '', {}, filename, lgLocation)

    

    
    
    return
%-------------------------------------------------------    
%-------------------------------------------------------    
%-------------------------------------------------------    
function analyseF(RES, testName, fig, nbAlgo)

    fGroundTruth = RES.focal;
    nbEdges = size(RES.f123,2);
    sfigure(fig); clf;
    
    clrMap = hsv(nbAlgo); %FIXME
    clrMap = varycolor(nbAlgo); %FIXME
    
    leg = {};
    hold on
    ii=1;
    for idx=1:nbAlgo
        vF123 =real(RES.f123(idx,:));
        vF12  =real(RES.f12(idx,:));
        
        vb123 = isValidF(vF123);
        vb12 = isValidF(vF12);
        vb = vb123 | vb12; %either one is good

        vT = 1:length(vF123);
        
        %vErrF123(ii) = mean(abs(vF123(vb) - fGroundTruth));
        %vErrF12(ii) = mean(abs(vF12(vb)  - fGroundTruth));
        vFboth = [vF123;vF12];
        [vFer,vIdx] = min(abs(vFboth - fGroundTruth));
        vF = vFboth(1,:);
        vF(vIdx==2) = vFboth(2,vIdx==2);
        vFerr{ii}    = abs(vF(vb)-fGroundTruth);
        vFerr123{ii} = abs(vF123(vb123)-fGroundTruth);
        vErrF(ii) = mean(vFerr{ii});

        %assert(all(vFerr123{ii} >= vFerr{ii}));
        
        vNbF(ii)  = sum(vb)/length(vb);
        plot(vT(vb), vF(vb), '.-', 'Color', clrMap(ii,:) );
        ii=ii+1;
    end
    legend(RES.algo);
    hold off
    
    disp([vErrF; vNbF])

    if false
        v=0:300;
    else
        v=0:150;
    end
        
    invfigure(fig+1); clf
    hold on
    for idx =1:length(vFerr)
        vErr = vFerr{idx};
        vV = arrayfun(@(x) sum(vErr <= x)/nbEdges, v);
        plot(v, vV, '-', 'Color', clrMap(idx,:), 'LineStyle', linestyle(idx), 'LineWidth', 1.5 );
    end
    daspect([v(end)/2 1 1]);
    hold off
    filename = sprintf('fig/%s-cumulHistoF.eps', testName);
    FACADE_printGraphics('', 'Focal Length error (pixel)', '', RES.algo, filename, 'SouthEast')
    
    
    invfigure(fig+2); clf
    hold on
    for idx =1:length(vFerr123)
        vErr = vFerr123{idx};
        vV = arrayfun(@(x) sum(vErr <= x)/nbEdges, v);
        plot(v, vV, '-', 'Color', clrMap(idx,:), 'LineStyle', linestyle(idx), 'LineWidth', 1.5 );
    end
    daspect([v(end)/2 1 1]);
    hold off
    filename = sprintf('fig/%s-cumulHistoF123.eps', testName);
    FACADE_printGraphics('', 'Focal Length error (pixel)', '', RES.algo, filename, 'SouthEast')

    return

% $$$ %-------------------------------------------------------
% $$$ %-------------------------------------------------------
% $$$ %-------------------------------------------------------
% $$$ function b = strBeginWith(str,st)
% $$$     
% $$$     
% $$$     if length(str) < length(st)
% $$$         b = false; 
% $$$         return;
% $$$     end
% $$$     
% $$$     b = strcmp(str(1:length(st)), st);