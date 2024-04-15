clear;
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
    '1.5h'];
ephdirName=['C:\Users\DELL\Desktop\GnssAnalysisFiles\demoStationaryNoEphemerisFiles\' ...
    'eph'];
FileFolder=fullfile(dirName);
DirOutput=dir(fullfile(FileFolder,'*.txt'));

for i=1:1%length(DirOutput)


prFileName=DirOutput(i).name;
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements
...................................................................................
% save data from GnssLogger App, and edit dirName and prFileName appropriately
%dirName = 'put the full path for your directory here';
%prFileName = 'put the pseuoranges log file name here';

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter;
[gnssRaw,gnssStatus,gnssFix,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
gnssRaw.Svid(gnssRaw.ConstellationType==5) = gnssRaw.Svid(gnssRaw.ConstellationType==5) + 100;
if ~isempty(gnssStatus)
    gnssStatus.Svid(gnssStatus.ConstellationType==5) = gnssStatus.Svid(gnssStatus.ConstellationType==5) + 100;
    [gnssStatus] = ProcessGnssStatus(gnssStatus);
end
if isempty(gnssRaw), return, end

%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end));
utcTime = Gps2Utc([],fctSeconds);
[allGpsEph,allBdsEph,iono] = GetNasaHourlyEphemeris(utcTime,ephdirName);
if isempty([allGpsEph,allBdsEph]), return, end

%% process raw measurements, compute pseudoranges:
[gnssMeas] = ProcessGnssMeas(gnssRaw,gnssStatus,gnssFix);
% gnssMeas = PseudorangeSmoother(gnssMeas);
%% plot pseudoranges and pseudorange rates
% h1 = figure;
% [colors] = PlotPseudoranges(gnssMeas,prFileName);
% h2 = figure;
% PlotPseudorangeRates(gnssMeas,prFileName,colors);
% h3 = figure;
% PlotCno(gnssMeas,prFileName,colors);

%% compute WLS position and velocity
gnssFix = [];
if(isempty(gnssFix))
    gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,allBdsEph,true);
else
    gpsPvt.allLlaDegDegM = [gnssFix.LatitudeDegrees,gnssFix.LongitudeDegrees,gnssFix.AltitudeMeters];
    % iono corr
    % if(~isempty(gnssStatus))
    %     D_iono = zeros(size(gnssMeas.tRxSeconds));
    %     minLen = min([length(gnssMeas.ClkDCount),length(gnssFix.AltitudeMeters),length(gnssStatus.UnixTimeMillis)]);
    %     for iCycle=1:minLen
    %         [~, D_iono(iCycle,:), ~, ~] = IonosphericDelay(gnssMeas.tRxSeconds(iCycle,1),gnssFix.LatitudeDegrees(iCycle),gnssFix.LongitudeDegrees(iCycle),gnssStatus.AzDeg(iCycle,:),gnssStatus.ElDeg(iCycle,:),iono.alpha,iono.beta);
    %     end
    %     D_iono(isnan(D_iono)) = 0;
    %     gnssMeas.PrM = gnssMeas.PrM - D_iono;
    % end
end

poiXyz = Lla2Xyz(POIPosition.poiLla);
minDist = zeros(size(gpsPvt.allLlaDegDegM,1),1);
minId = zeros(size(minDist));
gnssMeas.outAntXyz = zeros(size(gpsPvt.allLlaDegDegM));
for iPos = 1:size(gpsPvt.allLlaDegDegM,1)
    % [minDist(iPos),minId(iPos)] = min(vecnorm(Lla2Xyz(gpsPvt.allLlaDegDegM(iPos,:))-poiXyz,2,2));
    [minDist(iPos),minId(iPos)] = min(vecnorm(gpsPvt.allLlaDegDegM(iPos,1:2)-POIPosition.poiLla(:,1:2),2,2)*1.1e5*sqrt(2));
    gnssMeas.outAntXyz(iPos,:) = poiXyz(minId(iPos),:);
end

% scatter3(gpsPvt.allLlaDegDegM(:,1),gpsPvt.allLlaDegDegM(:,2),gpsPvt.allLlaDegDegM(:,3)); hold on ;
% scatter3(apiLla(:,1),apiLla(:,2),apiLla(:,3)); hold on ;
% scatter3(poiLla(1:3,1),poiLla(1:3,2),poiLla(1:3,3),"LineWidth",5);

% gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,allBdsEph,false); 
wireLength = gpsPvt.allBcMeters - cumsum(gpsPvt.allBcDotMps,"omitnan");

% % figure(i);
% % subplot(1,2,1)
% plot(gpsPvt.allBcMeters/GpsConstants.LIGHTSPEED);
% % subplot(1,2,2)
% sp = split(dirName,'\');
% sp = split(sp(end),'0');
% title(sp(1))
% plot(wireLength)
% hold on
% [b,wireLength_f] = stepDetect(wireLength,4,2,20);
% plot(wireLength_f)
% plot(b)

% subplot(2,2,1);
% plot(gpsPvt.allBcMeters);
% sp = split(dirName,'\');
% sp = split(sp(end),'0');
% title(sp(1)+"钟差")
% subplot(2,2,2);
% figure(2);
initclkbias = mean(gnssMeas.Fbns-gnssMeas.Fbns(1)-cumsum(gnssMeas.ClkDrift));
plot(wireLength-initclkbias);
% figure(2)
% plot(gpsPvt.allLlaDegDegM(:,3))


% hold on;
% [b,wireLength_f] = stepDetect(wireLength,4,2,20);
% % [b,wireLength_f] = extremeDetect(wireLength,4,2,3);
% 
% % plot(wireLength_f)
% % plot(b)
% title("线长")
% 
% %%%%%%%%%%%%%%% polyfit %%%%%%%%%%%%%%%%
% subplot(2,2,3);

% figure(2)
% winsize = 30;
% fitorder = 2;
% fitparam = zeros(length(gpsPvt.allBcMeters)-winsize+1,fitorder+1);
% fitresult = zeros(length(gpsPvt.allBcMeters)-winsize+1,1);
% for iBc = winsize:length(gpsPvt.allBcMeters)-1
%     bc = gpsPvt.allBcMeters(iBc-winsize+1:iBc);
%     fitparam(iBc-winsize+1,:) = polyfit(iBc-winsize+1:iBc,bc,fitorder);
%     pred = [(iBc+1)^2,iBc+1,1]*fitparam(iBc-winsize+1,:)';
%     fitresult(iBc-winsize+1) = gpsPvt.allBcMeters(iBc+1)-pred;
% end
% 
% plot(fitresult);
% title("window size = "+int2str(winsize))
% 
% load('bcMeters.mat'); % 导入钟差数据
% winsize = 15; % 窗口大小
% fitorder = 2; % 多项式拟合阶数
% fitresult = zeros(length(bcMeters)-winsize+1,1);
% for i = winsize:length(bcMeters)-1
%     bc = bcMeters(i-winsize+1:i); % 取等同于窗口大小的数据
% fitparam = polyfit(i-winsize+1:i,bc,fitorder); % 拟合
% pred = [(i+1)^2,i+1,1]*fitparam'; % 估计下一时刻钟差
%     fitresult(i-winsize+1) = bcMeters(i+1)-pred; % 测量值-预测值
% end
% 
% plot(fitresult);


    % % plot(minId);
% plot(gpsPvt.allBcDotMps);
% plot(movvar(gpsPvt.allBcDotMps,3));hold on;
% driftnanos = unique(gnssRaw.DriftNanosPerSecond,"stable")*GpsConstants.LIGHTSPEED*1e-9;
% plot(driftnanos)
% plot(movvar(driftnanos,3))
% title("解算钟漂滑动方差")
% 
% subplot(2,2,4);
% % plot(mean(gnssMeas.Cn0DbHz,2,"omitnan"))
% % plot(gpsPvt.allLlaDegDegM(:,3)) 
% % plot(gpsPvt.allBcDotMps)
% plot(vecnorm(gpsPvt.allVelMps,2,2));
% % plot(gpsPvt.allVelMps);
% title("速度模")


end
% ind = [[0,6];[6,12];[12,18]];
% xaxis = 1:length(gpsPvt.allBcDotMps);
% for i = 1:size(ind,1)
%     plot(xaxis(ind(i,1)*60+1:ind(i,2)*60),gpsPvt.allBcDotMps(ind(i,1)*60+1:ind(i,2)*60)); hold on;
% end
% for i = 18:29
%     plot(xaxis(i*60+1:(i+1)*60),gpsPvt.allBcDotMps(i*60+1:(i+1)*60)); hold on;
% end

% for i = 0:2:18
%     plot(xaxis(i*60+1:(i+2)*60),gpsPvt.allBcDotMps(i*60+1:(i+2)*60)); hold on;
% end

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
% 
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
