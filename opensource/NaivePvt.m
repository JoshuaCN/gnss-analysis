function [xHat,z,svPos,H,Wpr,Wrr,output] = NaivePvt(prs,gpsEph,xo)
% [xHat,z,svPos,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo)
% calculate a weighted least squares PVT solution, xHat
% given pseudoranges, pr rates, and initial state
%
% Inputs: 
%  prs: matrix of raw pseudoranges, and pr rates, each row of the form:
%  [trxWeek,trxSeconds,sv,prMeters,prSigmaMeters,prrMps,prrSigmaMps] 
%   trxWeek, trxSeconds: Rx time of measurement 
%      where trxSeconds = seconds in the current GPS week
%   sv: satellite id number
%   prMeters, prSigmaMeters: pseudorange and standard deviation (meters)
%   prrMps, prrSigmaMps: pseudorange rate and standard deviation (m/s)
%   gpsEph: matching vector of GPS ephemeris struct, defined in ReadRinexNav
%   xo: initial (previous) state, [x,y,z,bc,xDot,yDot,xDot,bcDot]'
%       in ECEF coordinates(meters and m/s)
%
% Outputs: xHat: estimate of state update
%          z = [zPr; zPrr] a-posteriori residuals (measured-calculated)
%          svPos: matrix of calculated sv positions and sv clock error: 
%                 [sv prn, x,y,z (ecef m), dtsv (s),xDot,yDot,zDot, dtsvDot]
%          H: H observation matrix corresponding to svs in svPos(:,1)
%          Wpr,Wrr Weights used in WlsPvt = 1/diag(sigma measurements)
%                  use these matrices to compute variances of xHat
%
% The PVT solution = xo + xHat, in ECEF coordinates
% For unweighted solution, set all sigmas = 1

%Author: Frank van Diggelen
%Open Source code for processing Android GNSS Measurements

jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

[bOk,numVal] = checkInputs(prs, gpsEph, xo);
if ~bOk
    error('inputs not right size, or not properly aligned with each other')
end

xHat=[]; z=[]; H=[]; svPos=[];
xyz0 = xo(1:3);
bc = xo(4);
bcBtwSys = xo(5);

ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
%goes negative, and it is handled in GpsEph2Pvt, where we work with fct
ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock 
% this is accurate satellite time of tx, because we use actual pseudo-ranges 
% here, not corrected ranges
% write the equation for pseudorange to see the rx clock error exactly cancel
% to get precise GPS time: we subtract the satellite clock error from sv time, 
% as done next:
dtsv = GpsEph2Dtsv(gpsEph,ttxSeconds);
dtsv = dtsv(:); %make into a column for compatibility with other time vectors
ttx = ttxSeconds - dtsv; %subtract dtsv from sv time to get tr0000000000000000000001ue gps time

%calculate satellite position at ttx
[svXyzTtx,dtsv,svXyzDot,dtsvDot]=GpsEph2Pvt(gpsEph,[ttxWeek,ttx]);
svXyzTrx = svXyzTtx; %initialize svXyz at time of reception

%Compute weights ---------------------------------------------------
Wpr = diag(1./prs(:,jPrSig));
Wrr = diag(1./prs(:,jPrrSig));
Wpr_vec = 1./(prs(:,jPrSig).*prs(:,jPrSig));
Wrr_vec = 1./(prs(:,jPrrSig).*prs(:,jPrrSig));
%iterate on this next part tilL change in pos & line of sight vectors converge
xHat=zeros(5,1);
dx=xHat+inf;
%we expect the while loop to converge in < 10 iterations, even with initial
%position on other side of the Earth (see Stanford course AA272C "Intro to GPS")
for i=1:length([gpsEph.PRN])
    % calculate tflight from, bc and dtsv
    dtflight = (prs(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv(i);
    % Use of bc: bc>0 <=> pr too big <=> tflight too big.
    %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
    % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
    %   i.e ttx = ttxsv - dtsv
    svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
end

%calculate line of sight vectors and ranges from satellite to xo
v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
range = sqrt( sum(v.^2) );
v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo

svPos=[prs(:,3),svXyzTrx,dtsv(:)];

%calculate the a-priori range residual
prHat = range(:) + bc -GpsConstants.LIGHTSPEED*dtsv;
% Use of bc: bc>0 <=> pr too big <=> rangehat too big
% Use of dtsv: dtsv>0 <=> pr too small

zPr = prs(:,jPr)-prHat; 
zeroOneVec = zeros(numVal,1);
zeroOneVec(prs(:,jSv)>100) = 1;
H = [v', ones(numVal,1), zeroOneVec]; % H matrix = [unit vector,1]

%z = Hx, premultiply by W: Wz = WHx, and solve for x:
% dx = pinv(Wpr*H)*Wpr*zPr;
bcBtwSys = mean(zPr(prs(:,jSv)>100)) - mean(zPr(prs(:,jSv)<100));
if isnan(bcBtwSys) % single constellation
    bcBtwSys = 0;
end
zPr(prs(:,jSv)>100) = zPr(prs(:,jSv)>100) - bcBtwSys;
dx = [0,0,0,sum(zPr.*Wpr_vec)/sum(Wpr_vec),bcBtwSys]';

% update xo, xhat and bc
xHat=xHat+dx;

%Now calculate the a-posteriori range residual
zPr = zPr-H*dx;

% Compute velocities ---------------------------------------------------------
rrMps = zeros(numVal,1);
for i=1:numVal
    %range rate = [satellite velocity] dot product [los from xo to sv]
    rrMps(i) = -svXyzDot(i,:)*v(:,i);
end
prrHat = rrMps + xo(9) - GpsConstants.LIGHTSPEED*dtsvDot;
zPrr = prs(:,jPrr)-prrHat;
%z = Hx, premultiply by W: Wz = WHx, and solve for x:
% vHat = pinv(Wrr*H)*Wrr*zPrr;

vHat = [0,0,0,sum(zPrr.*Wrr_vec/sum(Wrr_vec)),0]';
xHat = [xHat;vHat]; 

z = [zPr+bc,zPrr+xo(9)];

output = [prs(:,3),prs(:,5),prs(:,7),z];
end %end of function WlsPvt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bOk,numVal] = checkInputs(prs, gpsEph, xo)
%utility function for WlsPvt
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

bOk=true;
%check inputs
numVal=size(prs,1);  
if (max(prs(:,jSec))-min(prs(:,jSec)))> eps
  return
elseif length(gpsEph)~=numVal
    return
elseif any(prs(:,jSv) ~= [gpsEph.PRN]')
    return
elseif  any(size(xo) ~= [8,1])
    return
elseif size(prs,2)~=7
    return
else
    bOk = true;
end

%We insist that gpsEph and prs are aligned first.
%ClosestGpsEph.m does this, and passes back indices for prs - this is the way to
%do it right, so we don't have nested searches for svId

end %end of function checkInputs
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

