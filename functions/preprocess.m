function X_pre=preprocess(X)
% Preprocess the input data array. For a three-way data, it centers across
% the first (subject) mode, and scales within the third (metabolites) mode.
% For a two-way array, it centers across the first (subject) mode, and
% scales within the second (metabolites) mode. 
%
% The three-way array X has modes: subjects, time, metabolites
% The two-way array X has modes: subjects, metabolites

s=size(X);
if length(s)==2
    X_center=X-repmat(nanmean(X,1),s(1),1);
    for i=1:s(2)
        temp=X_center(:,i);
        X_pre(:,i)=temp/sqrt(nanmean(temp.^2));
    end
elseif length(s)==3
    
    X_center=X-repmat(nanmean(X,1),s(1),1);
    for i=1:s(3)
        temp=X_center(:,:,i);
        X_pre(:,:,i)=temp/sqrt(nanmean(temp.^2,'all'));
    end
end







