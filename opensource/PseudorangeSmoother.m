function gnssMeas = PseudorangeSmoother(gnssMeas)

rawPr = gnssMeas.PrM;
smPr = rawPr;
prSigma = gnssMeas.PrSigmaM;
prr = gnssMeas.PrrMps;
prrSigma = gnssMeas.PrrSigmaMps;
n = size(rawPr,2);
for nSv=1:n
    rawPr_nSv = rmmissing(rawPr(:,nSv));
    prSigma_nSv = rmmissing(prSigma(:,nSv));
    prr_nSv = rmmissing(prr(:,nSv));
    prrSigma_nSv = rmmissing(prrSigma(:,nSv));
    deltaT = diff(rmmissing(gnssMeas.tRxSeconds(:,nSv)));

    m = length(rawPr_nSv);
    I = eye(m);
    D = ([zeros(m-1,1),eye(m-1)]-[eye(m-1),zeros(m-1,1)])./deltaT;
    A = [I;D];

    y = [rawPr_nSv;prr_nSv(2:end)];
    W = diag(1./[prSigma_nSv;prrSigma_nSv(2:end)]);
    smPr(~isnan(smPr(:,nSv)),nSv) = pinv(W*A)*W*y; 
end
gnssMeas.PrM = smPr;

% function gnssMeas = PseudorangeSmoother(gnssMeas)
% 
% rawPr = gnssMeas.PrM;
% smPr = zeros(size(rawPr));
% prSigma = gnssMeas.PrSigmaM;
% prr = gnssMeas.PrrMps(2:end,:);
% prrSigma = gnssMeas.PrrSigmaMps(2:end,:);
% 
% [m,n] = size(rawPr);
% 
% I = eye(m);
% D = [zeros(m-1,1),eye(m-1)]-[eye(m-1),zeros(m-1,1)];
% A = [I;D];
% for nSv=1:n
%     y = [rawPr(:,nSv);prr(:,nSv)];
%     W = diag(1./[prSigma(:,nSv);prrSigma(:,nSv)]);
%     smPr(:,nSv) = pinv(W*A)*W*y; 
% end
% gnssMeas.PrM = smPr;