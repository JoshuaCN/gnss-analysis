%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you: 
prFileName = 'pseudoranges_log_2016_06_30_21_26_07.txt'; %with duty cycling, no carrier phase
% prFileName = 'pseudoranges_log_2016_08_22_14_45_50.txt'; %no duty cycling, with carrier phase
% as follows

% 1) copy everything from GitHub google/gps-measurement-tools/ to 
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
dirName = '~/Documents/MATLAB/gpstools/opensource/demoFiles';

dirName=['C:\Users\DELL\Desktop\GnssAnalysisFiles\demoStationaryNoEphemerisFiles\' ...
    'xiaomi_02_-1walk'];
FileFolder=fullfile(dirName);
DirOutput=dir(fullfile(FileFolder,'*.txt'));
prFileName=DirOutput.name;
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

%% data
%To add your own data:
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
gnssRaw.Svid(gnssRaw.ConstellationType==5) = gnssRaw.Svid(gnssRaw.ConstellationType==5) + 100;
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
[allGpsEph,allBdsEph] = GetNasaHourlyEphemeris(utcTime,dirName);
if isempty([allGpsEph,allBdsEph]), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw);
% gnssMeas = PseudorangeSmoother(gnssMeas);
%% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,allBdsEph,true);

poiLla = [[39.07214209,115.93713596,36.337];
[39.07236398,115.93627322,59.723];
[39.07201588,115.93641712,58.654];
[39.07159847,115.93663887,44.849];
[39.07267446,115.93601082,34.572]];

poiXyz = Lla2Xyz(poiLla);
minDist = zeros(size(gpsPvt.allBcMeters));
minId = zeros(size(gpsPvt.allBcMeters));
gnssMeas.outAntXyz = zeros(size(gpsPvt.allLlaDegDegM));
for iPos = 1:length(gpsPvt.allBcMeters)
    % [minDist(iPos),minId(iPos)] = min(vecnorm(Lla2Xyz(gpsPvt.allLlaDegDegM(iPos,:))-poiXyz,2,2));
    [minDist(iPos),minId(iPos)] = min(vecnorm(gpsPvt.allLlaDegDegM(iPos,1:2)-poiLla(:,1:2),2,2)*1.1e5);
    gnssMeas.outAntXyz(iPos,:) = poiXyz(minId(iPos),:);
end

% scatter3(gpsPvt.allLlaDegDegM(:,1),gpsPvt.allLlaDegDegM(:,2),gpsPvt.allLlaDegDegM(:,3)); hold on ;
% scatter3(apiLla(:,1),apiLla(:,2),apiLla(:,3)); hold on ;
% scatter3(poiLla(1:3,1),poiLla(1:3,2),poiLla(1:3,3),"LineWidth",5);

gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,allBdsEph,false);

wireLength = gpsPvt.allBcMeters - cumsum(gpsPvt.allBcDotMps,"omitnan");
plot(wireLength)

% z = [gpsPvt.allBcMeters,gpsPvt.allBcDotMps,wireLength];
% y = zeros(size(z));
% x_est = zeros(3,1);
% p_est = zeros(3);
% 
% for iz = 1:size(z,1)
%     [x_est,p_est,y(iz,:)] = Kalman_filter(1,1,1,1,x_est,p_est,z(iz,:)');
% end
% plot(y(:,3))

%% plot Pvt results
% h4 = figure;
% ts = 'Raw Pseudoranges, Weighted Least Squares solution';
% PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
% h5 = figure;
% PlotPvtStates(gpsPvt,prFileName);

%% Plot Accumulated Delta Range 
% if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
%     [gnssMeas]= ProcessAdr(gnssMeas);
%     h6 = figure;
%     PlotAdr(gnssMeas,prFileName,colors);
%     [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
%     h7 = figure;
%     PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
% end
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
