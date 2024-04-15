function gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,allBdsEph,bWls,bRaw)
%gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,bRaw)
%compute PVT from gnssMeas
% Input: gnssMeas, structure of pseudoranges, etc. from ProcessGnssMeas
%        allGpsEph, structure with all ephemeris
%        [bRaw], default true, true => use raw pr, false => use smoothed
%
% Output: 
% gpsPvt.FctSeconds    Nx1 time vector, same as gnssMeas.FctSeconds
%       .allLlaDegDegM Nx3 matrix, (i,:) = [lat (deg), lon (deg), alt (m)]
%       .sigmaLlaM     Nx3 standard deviation of [lat,lon,alt] (m)
%       .allBcMeters   Nx1 common bias computed with llaDegDegM
%       .allVelMps     Nx3 (i,:) = velocity in NED coords
%       .sigmaVelMps   Nx3 standard deviation of velocity (m/s)
%       .allBcDotMps   Nx1 common freq bias computed with velocity
%       .numSvs        Nx1 number of satellites used in corresponding llaDegDegM
%       .hdop          Nx1 hdop of corresponding fix
%
%Algorithm: Weighted Least Squares

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

if nargin<5
    bRaw = true;
else
    %check that smoothed pr fields exists in input
    if any(~isfield(gnssMeas,{'PrSmM','PrSmSigmaM'}))
       error('If bRaw is false, gnssMeas must have fields gnssMeas.PrSmM and gnssMeas.PrSmSigmaM')
    end
end

xo =zeros(10,1);%initial state: [center of the Earth, bc=0, velocities = 0]'

weekNum     = floor(gnssMeas.FctSeconds/GpsConstants.WEEKSEC);
%TBD check for week rollover here (it is checked in ProcessGnssMeas, but
%this function should stand alone, so we should check again, and adjust 
%tRxSeconds by +- a week if necessary)
%btw, Q. why not just carry around fct and not worry about the hassle of
%weeknumber, and the associated week rollover problems?
% A. because you cannot get better than 1000ns (1 microsecond) precsision
% when you put fct into a double. And that would cause errors of ~800m/s * 1us
% (satellite range rate * time error) ~ 1mm in the range residual computation
% So what? well, if you start processing with carrier phase, these errors
% could accumulate.

N = length(gnssMeas.FctSeconds);
gpsPvt.FctSeconds      = gnssMeas.FctSeconds;
gpsPvt.allLlaDegDegM   = zeros(N,3)+NaN; 
gpsPvt.sigmaLLaM       = zeros(N,3)+NaN;
gpsPvt.allBcMeters     = zeros(N,1)+NaN;
gpsPvt.allVelMps       = zeros(N,3)+NaN;
gpsPvt.sigmaVelMps     = zeros(N,3)+NaN;
gpsPvt.allBcDotMps     = zeros(N,1)+NaN;
gpsPvt.numSvs          = zeros(N,1);
gpsPvt.hdop            = zeros(N,1)+inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSvsize = 5;
zs = zeros(N,2*nSvsize);
zss = zeros(N,3*nSvsize);
x_est = zeros(nSvsize,3);
p_est = ones(nSvsize,3,3);
z_unfilt = zeros(N,nSvsize,3);
z_filt = zeros(N,nSvsize,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    iValid = find(isfinite(gnssMeas.PrM(i,:)) .* (gnssMeas.PrrSigmaMps(i,:)~=0) .*(gnssMeas.PrSigmaM(i,:)~=0)); %index into valid svid
    svid    = gnssMeas.Svid(iValid)';
    
    [gpsEph,iSv1] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(i,1));
    [bdsEph,iSv2] = ClosestGpsEph(allBdsEph,svid-100,gnssMeas.FctSeconds(i,2));
    if(isempty(bdsEph))
        keys = fieldnames(gpsEph);
        for key = keys'
            bdsEph.(key{1}) = [];
        end
    elseif(isempty(gpsEph))
        keys = fieldnames(bdsEph);
        for key = keys'
            gpsEph.(key{1}) = [];
        end
    else
        keys = fieldnames(bdsEph);
    end

    for key = keys'
        if(strcmp(key,"PRN"))
            allEph.(key{1}) = [gpsEph.(key{1}),[bdsEph.(key{1})]+100];
        else
        allEph.(key{1}) = [gpsEph.(key{1}),bdsEph.(key{1})];
        end
    end

    iSv = [iSv1,iSv2];
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid); %number of satellites this epoch
    gpsPvt.numSvs(i) = numSvs;

    
    %% WLS PVT -----------------------------------------------------------------
    %for those svIds with valid ephemeris, pack prs matrix for WlsNav
    prM     = gnssMeas.PrM(i,iValid(iSv))';
    prSigmaM= gnssMeas.PrSigmaM(i,iValid(iSv))';
    
    prrMps  = gnssMeas.PrrMps(i,iValid(iSv))';
    prrSigmaMps = gnssMeas.PrrSigmaMps(i,iValid(iSv))';
    
    tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];
    tRx(svid>100,1) = tRx(svid>100,1) - 1356;

    prs = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];
    
    xo(6:8) = zeros(3,1); %initialize speed to zero
    if(bWls)
        if numSvs<4
            continue;%skip to next epoch
        end
        [xHat,~,~,H,Wpr,Wrr] = WlsPvt(prs,allEph,xo);%compute WLS solution
        % [xHat,~,~,H,Wpr,Wrr] = WlsPvtAltHold(prs,allEph,xo);%compute WLS solution
    else
        xo(1:3) = Lla2Xyz([39.07267446,115.93601082,34.572]);
        xo(1:3) = gnssMeas.outAntXyz(i,:);
        [xHat,~,~,H,Wpr,Wrr,output{i}] = NaivePvt(prs,allEph,xo);%compute WLS solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % zss(i,:) = [zs(i,1:nSvsize),zs(i,nSvsize+1:2*nSvsize),zs(i,1:nSvsize)-sum(zs(1:i,nSvsize+1:2*nSvsize),1)];
        % z_unfilt(i,:,:) = reshape(zss(i,:),nSvsize,3);
        % z = reshape(zss(i,:),nSvsize,3);
        % for nSv=1:nSvsize
        %     [x_est(nSv,:),p_est(nSv,:,:),z(nSv,:)] = Kalman_filter(prSigmaM(nSv),prrSigmaMps(nSv),prSigmaM(nSv),1,x_est(nSv,:)',squeeze(p_est(nSv,:,:)),z(nSv,:)');
        % end
        % z_filt(i,:,:) = z;
    end
    xo = xo + xHat;
    
    %extract position states
    llaDegDegM = Xyz2Lla(xo(1:3)');
    gpsPvt.allLlaDegDegM(i,:) = llaDegDegM;
    gpsPvt.allBcMeters(i) = xo(4);
    gpsPvt.allBcMetersBtwSys(i) = xo(5);
    %extract velocity states
    RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
    %NOTE: in real-time code compute RE2N once until position changes
    vNed = RE2N*xo(6:8); %velocity in NED
    gpsPvt.allVelMps(i,:) = vNed;
    gpsPvt.allBcDotMps(i) = xo(9);
    gpsPvt.allBcDotMpsBtwSys(i) = xo(10);
    
    % if ~bWls
    %     gpsPvt.allBcMeters(i) = mean(z(:,1));
    %     gpsPvt.allBcDotMps(i) = mean(z(:,2));
    %     gpsPvt.allWireLength(i) = mean(z(:,3));
    % end

    % %compute HDOP
    % H = [H(:,1:3)*RE2N', ones(numSvs,1)]; %observation matrix in NED
    % P = inv(H'*H);%unweighted covariance
    % gpsPvt.hdop(i) = sqrt(P(1,1)+P(2,2));
    % 
    % %compute variance of llaDegDegM
    % %inside LsPvt the weights are used like this:
    % %  z = Hx, premultiply by W: Wz = WHx, and solve for x:
    % %  x = pinv(Wpr*H)*Wpr*zPr;
    % %  the point of the weights is to make sigma(Wz) = 1
    % %  therefore, the variances of x come from  diag(inv(H'Wpr'WprH))
    % P = inv(H'*(Wpr'*Wpr)*H); %weighted covariance
    % gpsPvt.sigmaLLaM(i,:) = sqrt(diag(P(1:3,1:3)));
    % 
    % %similarly, compute variance of velocity
    % P = inv(H'*(Wrr'*Wrr)*H); %weighted covariance
    % gpsPvt.sigmaVelMps(i,:) = sqrt(diag(P(1:3,1:3)));
    %%end WLS PVT --------------------------------------------------------------
end




end %end of function GpsWlsPvt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2016 Google Inc.
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

