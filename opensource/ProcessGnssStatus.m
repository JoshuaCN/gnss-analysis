function gnssStatusMeas = ProcessGnssStatus(gnssStatus)

% Filter valid values first, so that rollover checks, etc, are on valid data
gnssStatus = FilterValid(gnssStatus);
gnssStatusMeas.Svid       = unique(gnssStatus.Svid)'; %all the sv ids found in gnssRaw
gnssStatusMeas.UnixTimeMillis = unique(gnssStatus.UnixTimeMillis);
M = length(gnssStatusMeas.Svid);

N = length(gnssStatusMeas.UnixTimeMillis);
gnssStatusMeas.AzDeg      = zeros(N,M)+NaN;
gnssStatusMeas.ElDeg      = zeros(N,M)+NaN;

%Now pack these vectors into the NxM matrices
for i=1:N %i is index into gnssMeas.FctSeconds and matrix rows
    %get index of measurements within 1ms of this time tag
    J = find(abs(gnssStatusMeas.UnixTimeMillis(i) - gnssStatus.UnixTimeMillis)<1);
    for j=1:length(J) %J(j) is index into gnssRaw.*
        k = find(gnssStatusMeas.Svid==gnssStatus.Svid(J(j)));
        %k is the index into gnssMeas.Svid and matrix columns
        gnssStatusMeas.AzDeg(i,k) = gnssStatus.AzimuthDegrees(J(j));
        gnssStatusMeas.ElDeg(i,k) = gnssStatus.ElevationDegrees(J(j));
    end
end

end %of function ProcessGnssMeas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gnssRaw = FilterValid(gnssRaw)
%utility function for ProcessGnssMeas, 
%remove fields corresponding to measurements that are invalid
%NOTE: this makes it simpler to process data. But it removes data,
% so if you want to investigate *why* fields are invalid, then do so 
% before calling this function

%check ReceivedSvTimeUncertaintyNanos, PseudorangeRateUncertaintyMetersPerSecond
%for now keep only Svid with towUnc<0.5 microseconds and prrUnc < 10 mps
iBad = gnssRaw.ElevationDegrees < GnssThresholds.MINELEDEG;
if any(iBad)
    numBad = sum(iBad);
    %assert if we're about to remove everything:
    assert(numBad<length(iBad),'Removing all measurements in gnssStatus')
    
    names = fieldnames(gnssRaw);
    for i=1:length(names)
        ts = sprintf('gnssRaw.%s(iBad) = [];',names{i});
        eval(ts); %remove fields for invalid meas
    end
    %explain to user what happened:
    fprintf('\nRemoved %d bad meas inside ProcessGnssStatus>FilterValid because:\n',...
        sum(iBad))
    if any(iBad)
        fprintf('eleDeg < %.0f \n',GnssThresholds.MINELEDEG)
    end
end

end %end of function FilterValid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



