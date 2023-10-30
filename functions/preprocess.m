%%Preprocess the data: for a three-way data, centering across the first (subject) mode, and scaling within
%%the third (metabolites) mode;
%%For a two-way data: centering across the first(subject)mode, and scaling
%%within the second(metabolites) mode

%%the three-way data X has mode subjects*time*metabolites

%%the two-way data X has mode subjects*metabolites

function X_pre=preprocess(X)
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







