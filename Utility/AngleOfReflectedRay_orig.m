function [interceptPoint, outRay, transmission ] = AngleOfReflectedRay(coarseIntPt, originPt, originRay, surface, approximateNormal, surfaceFitSize)
    %originPt and originRay should be in x,y,z coordinate vector
    %coarseIntPt is the surface voxel that is intersected.
    %polAngle is not used currently.
    %approximateNormal is used to specify the direction of ray propogation as normal direction may be flipped
    %surfaceFitSize describes how many points are incorperated into
    %surface. 
    %This needs tuning. Large value may smooth too much, small value may be too rough.
    
    debugPlots = 0; scler = 50;
    
    %get all points within some range
    surfaceFitSize = 100;  
    surfInds = find(sqrt((surface(:,1)-coarseIntPt(1)).^2 + (surface(:,2)-coarseIntPt(2)).^2 + (surface(:,3)-coarseIntPt(3)).^2) < surfaceFitSize);

    %use pca to align 3rd principle axis with z axis - should make it flat (enough)
    surfPts = [surface(surfInds,1)-coarseIntPt(1),surface(surfInds,2)-coarseIntPt(2),surface(surfInds,3)-coarseIntPt(3)];
    pcaVecs = pca(surfPts);

    rer = vectorRotationFromAtoB([0 0 1]',pcaVecs(:,3)');
    rer2 = vectorRotationFromAtoB(pcaVecs(:,3)',[0 0 1]');  
    newPts = surfPts*rer;
    approximateNormal = approximateNormal'*rer;
    
    warning('should shift intersect point to zero - normalToSurface');
    
    %fit surface to rotated points - use cftool to tests different fns (poly22 usually works well)
    fittedPoly = fit(newPts(:,1:2), newPts(:,3),'poly22');
    
    %test plots
    if debugPlots == 1;
        close all; figure;
        subplot(1,3,1);
        plot3(surface(surfInds,1),surface(surfInds,2),surface(surfInds,3), '.'); hold on
        scatter3(coarseIntPt(:,1), coarseIntPt(:,2), coarseIntPt(:,3),50,'g'); 
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,1)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,1)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,1)], 'color', 'b')
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,2)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,2)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,2)], 'color', 'g')
        line([coarseIntPt(1) coarseIntPt(1)+pcaVecs(1,3)], [coarseIntPt(2) coarseIntPt(2)+pcaVecs(2,3)], [coarseIntPt(3) coarseIntPt(3)+pcaVecs(3,3)], 'color', 'r')
        subplot(1,3,2);
        plot3(newPts(:,1), newPts(:,2), newPts(:,3),'.'); hold on
        plot(fittedPoly);
    end
    
    rotatedOriginPt = (originPt-coarseIntPt)*rer;
    rotOriginRay = (originRay)*rer;
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); 
        line([-originRay(1) 0]*scler+coarseIntPt(1), [-originRay(2) 0]*scler+coarseIntPt(2), [-originRay(3) 0]*scler+coarseIntPt(3),'color','g');
        subplot(1,3,2);
        line([-rotOriginRay(1) 0]*scler, [-rotOriginRay(2) 0]*scler, [-rotOriginRay(3) 0]*scler,'color','g');
        line([rotOriginRay(1) 0]*scler, [rotOriginRay(2) 0]*scler, [rotOriginRay(3) 0]*scler,'color','g');
    end
    
    %function to solve to calculate fine intersection of ray with surface
    function F = solveTheFn(in,poly,pt,ang)
        F = [in(3) - poly(in(1),in(2)) ;
            in(1) - pt(1) - in(4)*ang(1) ;
            in(2) - pt(2) - in(4)*ang(2) ;
            in(3) - pt(3) - in(4)*ang(3)];
    end
    options = optimoptions('fsolve','Display','none','TolFun',1e-6);
    f = @(x)solveTheFn(x,fittedPoly,rotatedOriginPt,rotOriginRay);
    
    %use optimizer to find fine intersection
    [fineInt] = fsolve(f,[0 0 0 sqrt((originPt(1)-coarseIntPt(1))^2+(originPt(2)-coarseIntPt(2))^2+(originPt(3)-coarseIntPt(3))^2)],options);

    warning('Should not need to calculate intersect - see normalToSurface')
    
    %rotate intercept point back
    interceptPoint = fineInt(1:3)*rer2+coarseIntPt;
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); 
        scatter3(interceptPoint(1),interceptPoint(2),interceptPoint(3),5,'r');
        subplot(1,3,2);
        scatter3(fineInt(1),fineInt(2),fineInt(3),5,'r');
        line([-rotOriginRay(1)*scler 0]+fineInt(1), [-rotOriginRay(2)*scler 0]+fineInt(2), [-rotOriginRay(3)*scler 0]+fineInt(3),'color','m');
        scatter3(0,0,0,50,'g');
    end
    
    %get normal of surface
    [dX, dY] = differentiate(fittedPoly, fineInt(1), fineInt(2));
    if isnan(dX)
        dX = 0;
    end
    if isnan(dY)
        dY = 0;
    end
    
    N = [dX, dY, -1]; N = N/sqrt(N(1)^2+N(2)^2+N(3)^2);
    
    %flip normal to same side as ray propogation if they are more than 90deg apart
    if sqrt((approximateNormal(1)-N(1))^2+(approximateNormal(2)-N(2))^2+(approximateNormal(3)-N(3))^2) > sqrt(2)
        N = N*-1;
    end
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,3);
        plot(fittedPoly); hold on
        scatter3(fineInt(1),fineInt(2),fineInt(3),5,'r'); xlabel('x'); ylabel('y'); zlabel('z');
        line([0 N(1)]*scler+fineInt(1), [0 N(2)]*scler+fineInt(2), [0 N(3)]*scler+fineInt(3),'color','b');
        line([-rotOriginRay(1) 0]*scler+fineInt(1), [-rotOriginRay(2) 0]*scler+fineInt(2), [-rotOriginRay(3) 0]*scler+fineInt(3),'color','g');
        line([-N(1) 0]*scler+fineInt(1), [-N(2) 0]*scler+fineInt(2), [-N(3) 0]*scler+fineInt(3),'color','k');
    end
    
    %check if angle within 90 deg of normal (should be solved in previous step...)
    if sqrt((-rotOriginRay(1)-N(1))^2+(-rotOriginRay(2)-N(2))^2+(-rotOriginRay(3)-N(3))^2) < sqrt(2)
        
%         outRay = rIn/rOut*cross(N,cross(-N,rotOriginRay)) -N*sqrt(1-((rIn/rOut)^2)*dot(cross(N,rotOriginRay),cross(N,rotOriginRay)));
        outRay = rotOriginRay-2*dot(N,rotOriginRay)*N;
        
        transmission = 1; %assume for now
    else
        %problems if we get here...
        outRay = [NaN NaN NaN];
        transmission = 0;
    end
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,3); 
        line([0 outRay(1)]*scler+fineInt(1), [0 outRay(2)]*scler+fineInt(2), [0 outRay(3)]*scler+fineInt(3),'color','r');
        line([0 -outRay(1)]*scler+fineInt(1), [0 -outRay(2)]*scler+fineInt(2), [0 -outRay(3)]*scler+fineInt(3),'color','r');
    end   
    
    %rotate the ray back
    outRay = outRay*rer2;
    
    %test plots
    if debugPlots == 1;
        subplot(1,3,1); 
        line([0 outRay(1)]*scler+interceptPoint(1), [0 outRay(2)]*scler+interceptPoint(2), [0 outRay(3)]*scler+interceptPoint(3),'color','r');
    end
end