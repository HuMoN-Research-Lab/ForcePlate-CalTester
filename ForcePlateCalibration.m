%FORCEPLATECALIBRATION Find the error between force plate reference frame 
%and motion capture reference frame
%
%Instructions for collecting calibration data
%
%IMPORTANT NOTE: Do not step on any of the active force plates while 
% recording
%Turn on force plates and cameras and calibrate cameras
%Assemble cal tester and start a qualisys recording set for around 6 
% minutes
%Place cal tester plate on first force plate
%Place the point at the end of the cal tester into the divot on the cal
% tester plate 
%
%Press down on the cal tester with a good amount of force for at least five
%seconds 
%   -There is no minimum amount of force required but give a decent push   
%   -Not required, but recommended to only have one continuous push and not
%   to push down with a lot of force then stop and push again
%Move the cal tester plate to a new spot on the same force plate and repeat
% the previous step
%Move the cal tester plate to a third location press down on the cal tester
%in that location
%
%Repeat the previous three steps on each force plate to have three
% locations where force was generated on each force plate
%
%End qualisys recording, label cal tester trajectories and export the data
%to mat file that will be loaded into this script
%
%
%Instructions for processing data
%
%Change the variable filename to the .mat file of the cal tester recordine
%If needed, change the plateData variable to the correct filepath
%
%This script takes the raw qualisys data of the trajectories of the cal
%tester and the center of pressure(cop) of the forces that were generated 
%in the recording to find the location of the forces that were generated 
%based on the location of the trajectory of the cal tester
%
%This is done by utilizing an optimization function, fminunc to optimize
%the error between the cop and the location where the cal tester
%trajectories are
%
%A function, calTestErrorFun.m, is called during the script and needs to be
%added to the path of this script in order for it to run 
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% 
close all
%closes all open graphs

clearvars
%deletes all variables

restoredefaultpath
addpath(genpath(cd))
%clears the path and establishes a new folder pathw

%% Read in the Struct/ set variables

%Place all data into a struct named data
filename = 'CalTest0001';
plateData = load(['Data/' filename '.mat']);
data = plateData.(strrep(filename,'-','_'));

%Set variables for certain values that will be used later on
numForcePlates = length(data.Force);
numFrames = data.Frames;
numForceFrames = data.Force.NrOfSamples;
mocapFramerate = data.FrameRate;
forceFramerate = data.Force(1).Frequency;
frameF = (1:numForceFrames);
framet = (1:numFrames);
calTesterRadius = 151.25; %%mm


%% Read in Force Plates

order = 4;
cutoff = 10;

for ff = 1:numForcePlates
    
    %Change any NaN values to zero
    copData = data.Force(ff).COP;
    copData(isnan(copData))=0;
    momentData = data.Force(ff).Moment;
    momentData(isnan(momentData))=0;
    forceData = data.Force(ff).Force;
    forceData(isnan(forceData))=0;
    
    %Filter cop, moment and force
    cop_fp_fr_dim(ff,:,:)= butterLowZero(order,cutoff,forceFramerate,copData);
    moment_fp_fr_dim(ff,:,:)= butterLowZero(order,cutoff,forceFramerate,momentData);
    force_fp_fr_dim(ff,:,:)= butterLowZero(order,cutoff,forceFramerate,forceData);
    
    cop_fp_fr_dim(ff,:,1) = -cop_fp_fr_dim(ff,:,1);
    
    %Get force plate location data
    fpLoc{ff} = data.Force(ff).ForcePlateLocation;
end

%% Restructure trajectories into XYZ columns
traj = data.Trajectories;
traject = traj.Labeled.Data;
traject_traj_dim_frame = data.Trajectories.Labeled.Data;
markerLabels = data.Trajectories.Labeled.Labels;
numMarkers = data.Trajectories.Labeled.Count;
topLeft = traject(1,:,:);
topRight = traject(2,:,:);
center = traject(3,:,:);
bottomLeft= traject(4,:,:);
bottomRight = traject(5,:,:);

topLeft = (squeeze(topLeft))';
topRight = (squeeze(topRight))';
center = (squeeze(center))';
bottomRight =(squeeze(bottomRight))';
bottomLeft = (squeeze(bottomLeft))';
%% Plot
figure(2)
subplot(2,1,1)
plot(framet, bottomLeft(:,3))
%legend('x','y','z')
title('Caltester - bottomLeft')
subplot(2,1,2)
plot(framet,bottomRight(:,3))
%legend('x','y','z')
title('Caltester bottomRight')
%% Calcs
%find the rod coordinats of the midpoint between the top and bottom parts
topMean = (topLeft + topRight)/2;
botMean = (bottomLeft + bottomRight)/2;
botMean = botMean(:,[ 1 2 3]);
topMean = topMean(:,[ 1 2 3]);
differenceIntTopAndBottomMean = botMean - topMean;

%% Motion capture data plot
debug = false;
if debug
    %must number figures for future reference
    f = figure(43782);
    %facilitates output of figures across multiple mediums
    f.Units = 'normalized';
    f.Position = [-0.0042 0.4306 0.9401 0.4741];
    
    %for frames 1:numFrames at interval of 10
    for fr = 1000:40:numFrames
        %clears the current figures to avoid plotting data over each other
        clf
        
        %determines how many rows and columns a figure will have
        numRows = 2;
        numCols = 3;
        
        
        %subplot(numRows, numCols, 1:3)
        
        %column 1(x), 2(y), 3(z)
        plot3(traject_traj_dim_frame(:, 1, fr),...
            traject_traj_dim_frame(:, 2, fr),...
            traject_traj_dim_frame(:, 3, fr),'k.','MarkerFaceColor','k')
        
        %facilitates plotting hipIDs on top of points
        hold on
        plot3(botMean(fr,1), botMean(fr,2), botMean(fr,3),'k.','MarkerFaceColor','k')
        plot3(topMean(fr,1), topMean(fr,2), topMean(fr,3),'k.','MarkerFaceColor','k')
        %Plot the force plates and get the cop relative to the location in the
        %room
        copLocationX =[];
        copLocationY =[];
        copLocationZ =[];
        for ff = 1:numForcePlates
            patch(fpLoc{ff}(:,1),fpLoc{ff}(:,2),fpLoc{ff}(:,3) )
            
            copLocationX(ff,:) = (cop_fp_fr_dim(ff,:,1)) + (mean(fpLoc{ff}(:,1)));
            copLocationY(ff,:) = (cop_fp_fr_dim(ff,:,2)) + (mean(fpLoc{ff}(:,2)));
            copLocationZ(ff,:) = (cop_fp_fr_dim(ff,:,3)) + (mean(fpLoc{ff}(:,3)));
        end
        
        for ff = 1:numForcePlates
            
            quiver3((copLocationX(ff,(fr*4))),...
                (copLocationY(ff,(fr*4))),...
                (copLocationZ(ff,(fr*4))),...
                (squeeze(force_fp_fr_dim(ff,(fr*4),1))),...
                (-squeeze(force_fp_fr_dim(ff,(fr*4),2))),...
                (squeeze(force_fp_fr_dim(ff,(fr*4),3))),'LineWidth',5,'Color','g')
            
        end
        title(['frame# ' num2str(fr)])
        
        %Plot Line through center of cal tester
        
        calTesterTrajLine = [topMean(fr,:) differenceIntTopAndBottomMean(fr,:)];
        fpPlane = [fpLoc{1}(1,:) fpLoc{1}(2,:) fpLoc{1}(3,:)];
        calTesterTrajectIntersect(fr,[1 2 3]) = intersectLinePlane(calTesterTrajLine, fpPlane);
        
        plot3(calTesterTrajectIntersect(fr,1), calTesterTrajectIntersect(fr,2), calTesterTrajectIntersect(fr,3),'k.','MarkerFaceColor','k','MarkerSize',25)
        drawLine3d([topMean(fr,:), differenceIntTopAndBottomMean(fr,:)],'LineWidth',3,'Color','r')
        
        title(['frame# ' num2str(fr)])
        axis equal
        
        grid on
        
        %rough lab x y z limits
        xlim([0 1.5e3])
        ylim([-1e3 4e3])
        zlim([-5 1e3])
        
        %Orientation of where you view the plot from
        az = -120.362;
        el =  19.417;
        
        view(az,el)
        
        drawnow
    end
end
%% Find the point that a line through the top and bottom means of cal tester pass through force plate
for fr = 1:numFrames
    calTesterTrajLine = [topMean(fr,:) differenceIntTopAndBottomMean(fr,:)];
    fpPlane = [fpLoc{1}(1,:) fpLoc{1}(2,:) fpLoc{1}(3,:)];
    calTesterTrajectIntersect(fr,[1 2 3]) = intersectLinePlane(calTesterTrajLine, fpPlane);
    
end

%% Find the three points on each force plate that the Cal Tester was placed on 
clear fp

for ff = 1:numForcePlates
    
    fp{ff}  = squeeze(force_fp_fr_dim(ff,:,:));
    fpz{ff} = fp{ff}(:,3); %just the Z value of Force plate 1
    
    fpThreshValue = 100; %set threshold at 100
    
    fpThresh{ff} = fpz{ff} > fpThreshValue;
    diffFpThresh = diff(fpThresh{ff});
    
    startFr = find(diffFpThresh == 1);
    endFr = find(diffFpThresh == -1);
    amountOfClips = length(startFr);
    %If there is more than three clips that are above the threshold value 
    %the smallest clips must be removed
    k=0;
    p =0;
    for sf = 2:amountOfClips;
        if (startFr(sf) - startFr(sf-1)) < 10000
            k = k+1;
            startFrameErr(k,:) = startFr(sf);
        end
        if (endFr(sf) - endFr(sf-1)) < 10000
            p = p+1;
            endFrameErr(p,:) = endFr(sf-1);
        end
        
    end
    for ii =1: length(startFrameErr)
        startFr = startFr(startFr~=startFrameErr(ii));
        endFr = endFr(endFr~= endFrameErr(ii));
    end
    startFr;
    endFr;

%     %Change the clips from force frame rate to camera frame rate
    clip1Long{ff} = (startFr(1):endFr(1));
    clip2Long{ff} = (startFr(2):endFr(2));
    clip3Long{ff} = (startFr(3):endFr(3));
    first1 = find(mod(clip1Long{ff},4) ==0,1,'first');
    first2 = find(mod(clip2Long{ff},4) ==0,1,'first');
    first3 = find(mod(clip3Long{ff},4) ==0,1,'first');
    clip1{ff} = clip1Long{ff}(first1:4:end)/4;
    clip2{ff} = clip2Long{ff}(first2:4:end)/4;
    clip3{ff} = clip3Long{ff}(first3:4:end)/4;
    
    %numTest = 3;
    %assert(numel(startFr) == numTest, 'Found the wrong number of start Frames')
    %assert(numel(endFr) == numTest, 'Found the wrong number of end Frames')
    
    %Put CalTestintersection into three cell array clips
    caltestIntersectXYZ{1}(:,[1:3]) = calTesterTrajectIntersect(clip1{ff},:);
    caltestIntersectXYZ{2}(:,[1:3]) = calTesterTrajectIntersect(clip2{ff},:);
    caltestIntersectXYZ{3}(:,[1:3]) = calTesterTrajectIntersect(clip3{ff},:);
      
    %Put cop into three cell array clips
    copXYZ{1}(:,[1:3]) = squeeze(cop_fp_fr_dim(ff,(clip1{ff}*4),:));
    copXYZ{2}(:,[1:3]) = squeeze(cop_fp_fr_dim(ff,(clip2{ff}*4),:));
    copXYZ{3}(:,[1:3]) = squeeze(cop_fp_fr_dim(ff,(clip3{ff}*4),:));
  
    %find mean of each caltestIntersection clip and put the mean into a cell array
    for ii = 1:3
        meanCalTester_fp_clip_XYZ{ff}{ii} = mean(caltestIntersectXYZ{ii});
    end
    %find mean of each cop clip and put the mean into a cell array
    for ii = 1:3
        meanCOPXYZ{ff}{ii} = mean(copXYZ{ii});
    end
    
    
    clear caltestIntersectXYZ
    clear copXYZ
     %clear startFrameErr
     %clear startFrameFilt
end

%% debug plot , plots the threshold values
figure(4372)
    clf
for ff = 1:numForcePlates    
    numCols = numForcePlates;
    subplot(2,3,ff)
    plot(fpz{ff})
    hold on
    plot(fpThresh{ff}*100,'r.-')
    plot(diff(fpThresh{ff})*100,'m-')
    title(['ForcePlate #' num2str(ff)])
end
%% debug plug #2 , plots the locations of each force clip
figure(27119)
clf
for ff = 1:numForcePlates
    
    subplot(2,3,ff)
    plot(meanCalTester_fp_clip_XYZ{ff}{1}(1), meanCalTester_fp_clip_XYZ{ff}{1}(2), 'kx')
    hold on;
    plot(meanCalTester_fp_clip_XYZ{ff}{2}(1), meanCalTester_fp_clip_XYZ{ff}{2}(2), 'gx')
    plot(meanCalTester_fp_clip_XYZ{ff}{3}(1), meanCalTester_fp_clip_XYZ{ff}{3}(2), 'bx')
    
    plot(meanCOPXYZ{ff}{1}(1), meanCOPXYZ{ff}{1}(2), 'ko')
    plot(meanCOPXYZ{ff}{2}(1), meanCOPXYZ{ff}{2}(2), 'go')
    plot(meanCOPXYZ{ff}{3}(1), meanCOPXYZ{ff}{3}(2), 'bo')
    
    title(['ForcePlate #' num2str(ff)])
    xlim([-500 2000])
    ylim([-500 3000])
end
%% caltest error
%meanPoint bool value for when you want to run it based on the mean of the
%three force clips



%allThreePoints bool value for when you want to run it for all three force 
%clips
allThreePoints = true;
if allThreePoints
    for  ff = 1:numForcePlates
        for ClipNum = 1:3
            
            initialFpLocGuess =  [0 0 0];
            meanCalTestXYZ = meanCalTester_fp_clip_XYZ{ff};
            meanCopXYZ = meanCOPXYZ{ff};
            
            
            %Call Function and run through optimizer
            calTestErrorAnonFun = @(fpLocGuess) calTestErrorFun(meanCalTestXYZ, meanCopXYZ, fpLocGuess,ClipNum);
            opts = optimset('Display', 'iter', 'MaxFunEvals',5000, 'PlotFcns',{@optimplotx, @optimplotfval,@optimplotfirstorderopt});
            [fpLocOptClips{ClipNum}(ff,:), fpLocErrClips{ClipNum}(ff)] = fminunc(calTestErrorAnonFun, initialFpLocGuess, opts);
            
            finalAdjCopXYZ{ff}{ClipNum}(1:3) = fpLocOptClips{ClipNum}(ff,:)+ meanCOPXYZ{ff}{ClipNum}(1:3);
        end
    end
end

%% Plots the original cop and cal tester data and the translated cop data


%Non-Calibrated force plate corner coordinates. (0,0) is the center of the
%non-calibrated force plate and they are 400x600 mm
initialNoCalibrateFpLoc = [-200 300 0; 200 300 0; 200 -300 0; -200 -300 0];  


figure(45678) 

clf
for ff = 1:numForcePlates
   
    %find the mean difference in uncalibrated cop and calibrated cop XYZ for each plate
    meanfpLocOptClips(ff,1:3) = mean([fpLocOptClips{1}(ff,:); fpLocOptClips{2}(ff,:); fpLocOptClips{3}(ff,:)]);
        
    %Add the mean to the uncalibrated force plates
    calibratedFpLoc{ff} = initialNoCalibrateFpLoc + meanfpLocOptClips(ff,:);

    subplot(2,3,ff)
    plot(meanCalTester_fp_clip_XYZ{ff}{1}(1), meanCalTester_fp_clip_XYZ{ff}{1}(2), 'kx')
    hold on;
    plot(meanCalTester_fp_clip_XYZ{ff}{2}(1), meanCalTester_fp_clip_XYZ{ff}{2}(2), 'gx')
    plot(meanCalTester_fp_clip_XYZ{ff}{3}(1), meanCalTester_fp_clip_XYZ{ff}{3}(2), 'bx')
    
    plot(meanCOPXYZ{ff}{1}(1), meanCOPXYZ{ff}{1}(2), 'ko')
    plot(meanCOPXYZ{ff}{2}(1), meanCOPXYZ{ff}{2}(2), 'go')
    plot(meanCOPXYZ{ff}{3}(1), meanCOPXYZ{ff}{3}(2), 'bo')
    
    plot(finalAdjCopXYZ{ff}{1}(1), finalAdjCopXYZ{ff}{1}(2), 'rp')
    plot(finalAdjCopXYZ{ff}{2}(1), finalAdjCopXYZ{ff}{2}(2),'yp')
    plot(finalAdjCopXYZ{ff}{3}(1), finalAdjCopXYZ{ff}{3}(2), 'cp')
   
    patch(initialNoCalibrateFpLoc(:,1),initialNoCalibrateFpLoc(:,2),initialNoCalibrateFpLoc(:,3),'FaceColor','none','EdgeColor','r')
    patch(calibratedFpLoc{ff}(:,1), calibratedFpLoc{ff}(:,2),calibratedFpLoc{ff}(:,3),'FaceColor','none')
    
    title(['ForcePlate #' num2str(ff)])
    xlim([-500 2000])
    ylim([-500 3000])
end

figure(23482)
for ff = 1:numForcePlates
    
    
    plot(meanCalTester_fp_clip_XYZ{ff}{1}(1), meanCalTester_fp_clip_XYZ{ff}{1}(2), 'kx')
    hold on;
    plot(meanCalTester_fp_clip_XYZ{ff}{2}(1), meanCalTester_fp_clip_XYZ{ff}{2}(2), 'gx')
    plot(meanCalTester_fp_clip_XYZ{ff}{3}(1), meanCalTester_fp_clip_XYZ{ff}{3}(2), 'bx')
    
    plot(meanCOPXYZ{ff}{1}(1), meanCOPXYZ{ff}{1}(2), 'ko')
    plot(meanCOPXYZ{ff}{2}(1), meanCOPXYZ{ff}{2}(2), 'go')
    plot(meanCOPXYZ{ff}{3}(1), meanCOPXYZ{ff}{3}(2), 'bo')
    
    plot(finalAdjCopXYZ{ff}{1}(1), finalAdjCopXYZ{ff}{1}(2), 'rp')
    plot(finalAdjCopXYZ{ff}{2}(1), finalAdjCopXYZ{ff}{2}(2),'yp')
    plot(finalAdjCopXYZ{ff}{3}(1), finalAdjCopXYZ{ff}{3}(2), 'cp')
    patch(initialNoCalibrateFpLoc(:,1),initialNoCalibrateFpLoc(:,2),initialNoCalibrateFpLoc(:,3),'FaceColor','none','EdgeColor','r')
    patch(calibratedFpLoc{ff}(:,1), calibratedFpLoc{ff}(:,2),calibratedFpLoc{ff}(:,3),'FaceColor','none')
    
    title(['All Force Plates'])
    xlim([-500 2000])
    ylim([-500 3000])    
end
%% Find the rotation  of the cop forceplate and the caltester force plate
for ff = 1:numForcePlates
    
    calTestCoordinates = [meanCalTester_fp_clip_XYZ{ff}{1}; meanCalTester_fp_clip_XYZ{ff}{2}; meanCalTester_fp_clip_XYZ{ff}{3}]';
    copCoordinates = [meanCOPXYZ{ff}{1}; meanCOPXYZ{ff}{2}; meanCOPXYZ{ff}{3}]';
    
    
    [regParams{ff}, Bfit{ff}, ErrorStats{ff}] = absor(copCoordinates, calTestCoordinates);
    %The absor function returns the struct regParams
    %RegParams gives the rotation/translation matrix, scale factor(in this
    %case the scale factor is 1), and a unit quarternion 
    rotationMatrix{ff}= regParams{ff}.R;
    rotationAngleXYZ(ff,1) = asind(rotationMatrix{ff}(2,3));
    rotationAngleXYZ(ff,2) = asind(rotationMatrix{ff}(1,3));
    rotationAngleXYZ(ff,3) = asind(rotationMatrix{ff}(2,1)); 
end
%% Testing for Difference Between Classic Cal Tester and This Code Cal Tester
for ff = 1:numForcePlates
    fpLocdiff{ff} = fpLoc{ff} - calibratedFpLoc{ff};
end
%% Create vectors for plotting

for ff = 1:numForcePlates
    
    timebeforeClip1 = NaN(3,clip1{ff}(1),'double');
    timeduringClip1 = repmat(finalAdjCopXYZ{ff}{1},length(clip1{ff}),1);
    idxAfterClip = numFrames - (length(timebeforeClip1) + length(timeduringClip1));
    timeafterClip1 = NaN(3,idxAfterClip,'double');  
    finalAdjCopWithTimeClip1{ff}= [timebeforeClip1'; timeduringClip1; timeafterClip1';];
    
    timebeforeClip2 = NaN(3,clip2{ff}(1),'double');
    timeduringClip2 = repmat(finalAdjCopXYZ{ff}{2},length(clip2{ff}),1);
    idxAfterClip = numFrames - (length(timebeforeClip2) + length(timeduringClip2));
    timeafterClip2 = NaN(3,idxAfterClip,'double');  
    finalAdjCopWithTimeClip2{ff}= [timebeforeClip2'; timeduringClip2; timeafterClip2';];

    timebeforeClip3 = NaN(3,clip3{ff}(1),'double');
    timeduringClip3 = repmat(finalAdjCopXYZ{ff}{3},length(clip3{ff}),1);
    idxAfterClip = numFrames - (length(timebeforeClip3) + length(timeduringClip3));
    timeafterClip3 = NaN(3,idxAfterClip,'double');  
    finalAdjCopWithTimeClip3{ff}= [timebeforeClip3'; timeduringClip3; timeafterClip3';];

end
  
%% Motion capture data plot
debug = true;


figure(43782);

%for frames 1:numFrames at interval of 10
v = VideoWriter('CalTestScriptComparison', 'MPEG-4');

open(v)
for fr = 1000:10:numFrames
    %clears the current figures to avoid plotting data over each other
    
    s = subplot(2,1,1);
    cla(s)
    %column 1(x), 2(y), 3(z)
    
    plot3(traject_traj_dim_frame(:, 1, fr),...
        traject_traj_dim_frame(:, 2, fr),...
        traject_traj_dim_frame(:, 3, fr),'k.','MarkerFaceColor','k')
    
    %facilitates plotting hipIDs on top of points
    hold on
    plot3(botMean(fr,1), botMean(fr,2), botMean(fr,3),'k.','MarkerFaceColor','k')
    plot3(topMean(fr,1), topMean(fr,2), topMean(fr,3),'k.','MarkerFaceColor','k')
   
    
    for ff = 1:numForcePlates
        patch(calibratedFpLoc{ff}(:,1),calibratedFpLoc{ff}(:,2),calibratedFpLoc{ff}(:,3))
        
        quiver3((finalAdjCopWithTimeClip1{ff}(fr,1)),...
            (finalAdjCopWithTimeClip1{ff}(fr,2)),...
            (finalAdjCopWithTimeClip1{ff}(fr,3)),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),1))),...
            (-squeeze(force_fp_fr_dim(ff,(fr*4),2))),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),3))),'LineWidth',3,'Color','g')
        
        quiver3((finalAdjCopWithTimeClip2{ff}(fr,1)),...
            (finalAdjCopWithTimeClip2{ff}(fr,2)),...
            (finalAdjCopWithTimeClip2{ff}(fr,3)),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),1))),...
            (-squeeze(force_fp_fr_dim(ff,(fr*4),2))),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),3))),'LineWidth',3,'Color','g')
        
        quiver3((finalAdjCopWithTimeClip3{ff}(fr,1)),...
            (finalAdjCopWithTimeClip3{ff}(fr,2)),...
            (finalAdjCopWithTimeClip3{ff}(fr,3)),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),1))),...
            (-squeeze(force_fp_fr_dim(ff,(fr*4),2))),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),3))),'LineWidth',3,'Color','g')
    end 
   %Plot Line through center of cal tester
    
    calTesterTrajLine = [topMean(fr,:) differenceIntTopAndBottomMean(fr,:)];
    fpPlane = [fpLoc{1}(1,:) fpLoc{1}(2,:) fpLoc{1}(3,:)];
    calTesterTrajectIntersect(fr,[1 2 3]) = intersectLinePlane(calTesterTrajLine, fpPlane);
    
    plot3(calTesterTrajectIntersect(fr,1), calTesterTrajectIntersect(fr,2), calTesterTrajectIntersect(fr,3),'k.','MarkerFaceColor','k','MarkerSize',5)
    drawLine3d([topMean(fr,:), differenceIntTopAndBottomMean(fr,:)],'LineWidth',1,'Color','r')
    
    title(['Data from Matlab Script '])
    axis equal
    
    grid on
    
    %rough lab x y z limits
    xlim([0 1.5e3])
    ylim([-1e3 4e3])
    zlim([-5 1e3])
    
    %Orientation of where you view the plot from
    az = -120.362;
    el =  19.417;
    
    view(az,el)
    
    drawnow
    
    
    o = subplot(2,1,2);
    cla(o)
    
    %column 1(x), 2(y), 3(z)
    plot3(traject_traj_dim_frame(:, 1, fr),...
        traject_traj_dim_frame(:, 2, fr),...
        traject_traj_dim_frame(:, 3, fr),'k.','MarkerFaceColor','k')
    
    %facilitates plotting hipIDs on top of points
    hold on
    plot3(botMean(fr,1), botMean(fr,2), botMean(fr,3),'k.','MarkerFaceColor','k')
    plot3(topMean(fr,1), topMean(fr,2), topMean(fr,3),'k.','MarkerFaceColor','k')
    %Plot the force plates and get the cop relative to the location in the
    %room
    copLocationX =[];
    copLocationY =[];
    copLocationZ =[];
    for ff = 1:numForcePlates
        patch(fpLoc{ff}(:,1),fpLoc{ff}(:,2),fpLoc{ff}(:,3) )
        
        copLocationX(ff,:) = (cop_fp_fr_dim(ff,:,1)) + (mean(fpLoc{ff}(:,1)));
        copLocationY(ff,:) = (cop_fp_fr_dim(ff,:,2)) + (mean(fpLoc{ff}(:,2)));
        copLocationZ(ff,:) = (cop_fp_fr_dim(ff,:,3)) + (mean(fpLoc{ff}(:,3)));
    end
    
    for ff = 1:numForcePlates
        
        quiver3((copLocationX(ff,(fr*4))),...
            (copLocationY(ff,(fr*4))),...
            (copLocationZ(ff,(fr*4))),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),1))),...
            (-squeeze(force_fp_fr_dim(ff,(fr*4),2))),...
            (squeeze(force_fp_fr_dim(ff,(fr*4),3))),'LineWidth',3,'Color','g')
        
    end
    title(['Data from Cal Tester '])
    
    %Plot Line through center of cal tester
    
    calTesterTrajLine = [topMean(fr,:) differenceIntTopAndBottomMean(fr,:)];
    fpPlane = [fpLoc{1}(1,:) fpLoc{1}(2,:) fpLoc{1}(3,:)];
    calTesterTrajectIntersect(fr,[1 2 3]) = intersectLinePlane(calTesterTrajLine, fpPlane);
    
    plot3(calTesterTrajectIntersect(fr,1), calTesterTrajectIntersect(fr,2), calTesterTrajectIntersect(fr,3),'k.','MarkerFaceColor','k','MarkerSize',5)
    drawLine3d([topMean(fr,:), differenceIntTopAndBottomMean(fr,:)],'LineWidth',1,'Color','r')
    
    axis equal
    grid on
    
    %rough lab x y z limits
    xlim([0 1.5e3])
    ylim([-1e3 4e3])
    zlim([-5 1e3])
    
    %Orientation of where you view the plot from
    az = -120.362;
    el =  19.417;
    
    view(az,el)
    
    drawnow
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)







