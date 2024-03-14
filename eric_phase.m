function [avgPhase,asymPhase ] = eric_phase(EruptionTimes)
N=size(EruptionTimes,1);
avgPhase = zeros(N,1);
asymPhase = zeros(N,1);
% EruptionTimes(:,1) = -EruptionTimes(:,1);
    for k=1:N
        % phase of current tooth relative to the previous and next tooth in the position to the left 
        indLeft = find(EruptionTimes(:,1) == EruptionTimes(k,1)-1); % returns indices of all teeth in the position to the left of tooth k
        indPrev = find(EruptionTimes(indLeft,2)<EruptionTimes(k,2));
        prev = max(EruptionTimes(indLeft(indPrev),2));
        indNext = find(EruptionTimes(indLeft,2)>=EruptionTimes(k,2));
        next = min(EruptionTimes(indLeft(indNext),2));
        leftPhase = (EruptionTimes(k,2)-prev) / (next-prev);
        % phase of current tooth relative to the previous and next tooth in the position to the right 
        indRight = find(EruptionTimes(:,1) == EruptionTimes(k,1)+1); % returns indices of all teeth in the position to the left of tooth k
        indPrev = find(EruptionTimes(indRight,2)<EruptionTimes(k,2));
        prev = max(EruptionTimes(indRight(indPrev),2));
        indNext = find(EruptionTimes(indRight,2)>=EruptionTimes(k,2));
        next = min(EruptionTimes(indRight(indNext),2));
        rightPhase = (EruptionTimes(k,2)-prev) / (next-prev);
        if isempty(leftPhase)|isempty(rightPhase)
            avgPhase(k) = NaN;
            asymPhase(k) = NaN;
        else
            avgPhase(k) = (rightPhase+leftPhase)/2;
            asymPhase(k) = rightPhase-leftPhase;
        end
    end

end