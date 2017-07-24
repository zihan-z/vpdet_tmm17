function vp_onL = FACADE_closest_x_on_l(vL, vp)
    
    
%DEBUG
    if nargin==0
        vp1 = [1,2,1]';
        vp2 = [2,4,1]';
        vpout = [3,1,1]';
        
        vL = cross(vp1,vp2);
        
        vp_onL = FACADE_closest_x_on_l(vL, vpout);
        
        figure(1);clf(1)
        hold on
        plot(vp1(1),vp1(2),'.', 'Color', [1,0,0]);
        plot(vp2(1),vp2(2),'.', 'Color', [1,0,0]);
        plot(vpout(1), vpout(2), '.', 'Color', [0,1,0]);
        plot(vp_onL(1)/vp_onL(3), vp_onL(2)/vp_onL(3), 'x', 'Color', [0,0,1]);
        hline(vL, 'Color', [1,0,0] );
        hline(cross(vp_onL, vpout));
        hold off
        daspect([1,1,1]);
        return
    end
    %END
    
    a=vL(1);
    b=vL(2);
    c=vL(3);
    
    if min(size(vp)==1)
        %only a vector
        d = [b,-a]*vp(1:2) / vp(3);
        vp_onL = cross(vL, [-b,a,d])';
        
        assert(abs(vp_onL'*vL)<0.001);
    
        %DEBUG
        %verify that distance vp to vL is same as distance vp to vp_onL
        %norm([vp_onL(1)/vp_onL(3);    vp_onL(2)/vp_onL(3)] - [vp(1)/vp(3); vp(2)/vp(3)])   
        %abs(vp'*vL / norm(vL(1:2)))
    else
        %matrix
        nbPt = size(vp,2);
        if size(vp,1)==2
            vd = [b,-a]*vp(1:2,:);
        else
            vd = [b,-a]*vp(1:2,:) ./ vp(3,:);
        end
        msk = skew_mat(vL);
        v0 = zeros(1,nbPt);
        vLperp = [v0-b;v0+a;vd];  %perpendicular line passing by point
        vp_onL = msk*vLperp;              % intersection
        vp_onL(1,:) = vp_onL(1,:) ./ vp_onL(3,:);
        vp_onL(2,:) = vp_onL(2,:) ./ vp_onL(3,:);
        vp_onL = vp_onL(1:2,:);
        
% $$$         %DEBUG
% $$$         sfigure(100); clf(100);
% $$$         hold on
% $$$         plot(vp(1,:), vp(2,:), '.');
% $$$         plot(vp_onL(1,:), vp_onL(2,:), '.', 'Color', [1,0,0]);
% $$$         %hline(vL);
% $$$         hold off
% $$$         pause
    end
        
        
        
