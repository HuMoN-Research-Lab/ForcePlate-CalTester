


%find angle between cal tester and ground through trig
diffMeans = (((topMean(:,3)) - (botMean(:,3))));
diffHighLowTraj = 600;                          %Difference in Z coordinate of top part and bottom part of tester when standing straight up
theta = abs(asind(diffMeans/diffHighLowTraj));  %Use inverse sin to find the angle between ground and caltester

%%Find angle on xy plane
diffLeftRight = topLeft - topRight;
deltaX = diffLeftRight(:,1)/2;                  %Divide by 2 because just using one side of the trajectory
deltaY = diffLeftRight(:,2)/2;
alpha = atand(deltaY./deltaX);

% Find difference from COP and trajectories
botDiff2cop = abs(botMean(:,3)./(tand(theta)));
XDiff2cop = cos(alpha).*botDiff2cop;
YDiff2cop = sin(alpha).*botDiff2cop;
Xtrajcop = botMean(:,1) + XDiff2cop;
Ytrajcop = botMean(:,2) + YDiff2cop;
