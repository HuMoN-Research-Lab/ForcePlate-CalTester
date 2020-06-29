function [error] = calTestErrorFun(meanCalTestXYZ,meanCopXYZ, fpLocGuess,ClipNum)
%calTestErrorFun - Error function for cal tester calibration procedure
%Function takes a guess of the difference between the cop and cal tester
%location and then takes the error of that guess. That error is then
%outputed to a the function fminunc which returns a new guess. This is
%reeated until the error is very low (~0).


%All three point bool values are to run the calibration script with values
%of all three clips
allThreePoints = true;
allThreePointsDebugPlot = true;


if allThreePoints
        
    %Find the new COP based on the guess of the fpLoc
    adjCOPxyz(1) = meanCopXYZ{ClipNum}(1)+fpLocGuess(1); 
    adjCOPxyz(2) = meanCopXYZ{ClipNum}(2)+fpLocGuess(2);
    adjCOPxyz(3) = meanCopXYZ{ClipNum}(3)+fpLocGuess(3);
    
    %Find the error of the caltester location and cop location
    error= sqrt(mean([meanCalTestXYZ{ClipNum}(1)-adjCOPxyz(1) meanCalTestXYZ{ClipNum}(2)-adjCOPxyz(2) meanCalTestXYZ{ClipNum}(3)-adjCOPxyz(3)].^2));
 
    % debug plot
   
    if allThreePointsDebugPlot
        figure(432773) 
        %Plot the three generated forces and the adjusted cop
        clf
        plot3(meanCalTestXYZ{1}(1), meanCalTestXYZ{1}(2), meanCalTestXYZ{1}(3), 'kx')
        hold on;
        plot3(meanCalTestXYZ{2}(1), meanCalTestXYZ{2}(2), meanCalTestXYZ{2}(3), 'gx')
        plot3(meanCalTestXYZ{3}(1), meanCalTestXYZ{3}(2), meanCalTestXYZ{3}(3), 'bx')
        
        plot3(meanCopXYZ{1}(1), meanCopXYZ{1}(2), meanCopXYZ{1}(3), 'ko')
        plot3(meanCopXYZ{2}(1), meanCopXYZ{2}(2), meanCopXYZ{2}(3), 'go')
        plot3(meanCopXYZ{3}(1), meanCopXYZ{3}(2), meanCopXYZ{3}(3), 'bo')
        
        plot3(adjCOPxyz(1), adjCOPxyz(2), adjCOPxyz(3), 'rp')
        
        axis equal
        drawnow
        pause(.2)
    end
end

