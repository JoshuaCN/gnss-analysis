% 0714步行
% time = [
%     1,32,1*60+14,1*60+51,2*60+20,2*60+41;
%     3*60+2,3*60+28,4*60+11,4*60+45,5*60+13,5*60+29;
%     5*60+48,6*60+14,6*60+57,7*60+33,8*60+4,8*60+18;
%     8*60+39,9*60+6,9*60+50,10*60+27,10*60+58,11*60+13;
%     11*60+34,12*60+1,12*60+47,13*60+25,13*60+55,14*60+11
%     ];
% time_end = 14*60+30;

% 0714驾车
% time = [
%     4,23,49,1*60+11,1*60+26,1*60+36;
%     1*60+47,2*60+2,2*60+25,2*60+46,3*60+3,3*60+14;
%     3*60+25,3*60+39,4*60+3,4*60+26,4*60+43,4*60+54;
%     5*60+6,5*60+22,5*60+49,6*60+14,6*60+33,6*60+46;
%     6*60+58,7*60+15,7*60+42,8*60+5,8*60+23,8*60+35;
%     8*60+48,8*60+57,9*60+11,9*60+24,9*60+34,9*60+41;
%     9*60+49,9*60+57,10*60+11,10*60+23,10*60+33,10*60+40;
%     10*60+47,10*60+59,11*60+14,11*60+27,11*60+36,11*60+45;
%     11*60+54,12*60,12*60+13,12*60+26,12*60+36,12*60+44;
%     12*60+51,12*60+58,13*60+12,13*60+24,13*60+33,13*60+42;
%     ];
% time_end = 13*60+51;

% % 0717步行
% time = [
%     7,34,1*60+13,1*60+45,2*60+13,2*60+29;
%     2*60+47,3*60+12,3*60+50,4*60+22,4*60+49,5*60+4;
%     5*60+23,5*60+47,6*60+27,6*60+59,7*60+27,7*60+42;
%     7*60+59,8*60+22,9*60,9*60+31,9*60+57,10*60+11;
%     10*60+28,10*60+51,11*60+29,12*60+2,12*60+30,12*60+46;
% ];
% time_end = 13*60+3;
% % 0717驾车
% time = [
%     2,14,30,45,53,61;
%     1*60+14,1*60+23,1*60+39,1*60+52,2*60+1,2*60+9;
%     2*60+18,2*60+26,2*60+40,2*60+54,3*60+2,3*60+10;
%     3*60+19,3*60+28,3*60+43,3*60+57,4*60+6,4*60+13;
%     4*60+22,4*60+30,4*60+45,5*60,5*60+9,5*60+17;
%     5*60+25,5*60+33,5*60+45,5*60+57,6*60+3,6*60+17;
%     6*60+26,6*60+34,6*60+46,6*60+58,7*60+5,7*60+13;
%     7*60+23,7*60+30,7*60+43,7*60+56,8*60+4,8*60+12;
%     8*60+20,8*60+28,8*60+41,8*60+54,9*60+2,9*60+10;
%     9*60+18,9*60+25,9*60+38,9*60+50,9*60+57,10*60+3;
% ];
% time_end = 10*60+20;
% 
% 
% time_vec = time(:);
% [numRound,numAnt] = size(time);
% 
% figure()
% hold on;

% wireLength = load("C:\Users\DELL\Desktop\nova5r.mat").wireLength;

function [b,wireLength_f]=stepDetect(wireLength,windowsize,k,thres)
    %% moving median filter
    if(windowsize<2)
        wireLength_f = wireLength;
    else
        wireLength_f(1:windowsize-1)=wireLength(1:windowsize-1);
        for i=windowsize:length(wireLength)
            wireLength_f(i) = median(wireLength(i-windowsize+1:i));
        end
    end
    
    %% bicusum
    mu = 0;
    wireLength_diff = diff(wireLength_f);
    cusum_up = zeros(size(wireLength_diff));
    cusum_down = cusum_up;
    b = cusum_up;
    for i=4:length(wireLength_diff)
        cusum_up(i) = max([0,cusum_up(i-1)-k+wireLength_diff(i-1)-mu]);
        cusum_down(i) = min([0,cusum_down(i-1)+k+wireLength_diff(i-1)+mu]);
        if cusum_up(i)>thres
            b(i) = cusum_up(i);
            cusum_up(i) = 0;
        elseif cusum_down(i)<-thres
            b(i) = cusum_down(i);
            cusum_down(i) = 0;
        else
            b(i) = 0;
        end
        
        if any(b(i-3:i-1))
            cusum_up(i) = 0;
            cusum_down(i) = 0;
            b(i) = 0;
        end
    end
end
%% plot
% plotTarget = diff(wireLength_f);
% plotTarget = gpsPvt.allBcDotMps;
% for i=1:numRound
%     plot(time(i,:),plotTarget(time(i,:)),'o',MarkerSize=10,LineWidth=2)
% end
% for i=1:numRound-1
%     plot(time(i,1):time(i+1,1)-1,plotTarget(time(i,1):time(i+1,1)-1))
% end
% i = numRound;
% plot(time(i,1):time_end-1,plotTarget(time(i,1):time_end-1))
% plot(time_end,plotTarget(time_end),'o',MarkerSize=10,LineWidth=2)