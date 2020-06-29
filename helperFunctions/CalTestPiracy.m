close all
clearvars
restoredefaultpath
addpath(genpath(cd))
%% load data

fileName = 'TopRightPlate4';

data  = load([fileName '.mat']);

data = data.(fileName);

numForcePlates = length(data.Force);

assert(numForcePlates == 5)


numFrames = data.Frames;
framerate = data.FrameRate;

order = 4
cutoff = 10;

%pull out force data for each force plate and butterworth filter
for ff = 1:numForcePlates

    cop_fp_dim_fr(ff,:,:)= butterLowZero(order,cutoff,framerate,data.Force(ff).COP);
    moment_fp_dim_fr(ff,:,:)= butterLowZero(order,cutoff,framerate,data.Force(ff).Moment);
    
    frLoc{ff} = data.Force(ff).ForcePlateLocation;           
end


%% plot some stuff

figure(3626)

for ff = 1:numForcePlates
subplot(2,3,ff)
plot(squeeze(cop_fp_dim_fr(ff,:,:)))
end

figure(764)

for ff = 1:numForcePlates
subplot(2,3,ff)
plot(squeeze(moment_fp_dim_fr(ff,:,:)))
end
