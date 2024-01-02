
%%%%%% This code is to remove a subset of subjects from the considered dataset
% outlier_pid: the index of the subjects to be removed
% input Y: the input dataset

function Y=removeoutlier(Y,outlier_pid)

[subject_notoutlier,index] = setdiff(str2num(Y.label{1}), outlier_pid);
Y=Y(index,:,:);

end