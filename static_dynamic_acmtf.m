% This script shows how to fit an ACMTF model to jointly analyze a matrix
% (static) and tensor (dynamic) metabolomics data
%% add auxilary functions to path
addpath(genpath('./functions'))

% add the necessary toolboxes: Tensor Toolbox, CMTF Toolbox, Poblano toolbox and dataset object
% see https://github.com/Lu-source/Joint-analysis-of-static-and-dynamic-metabolomics-data/blob/main/README.md

%%  load dataset
load('data.mat','X_T0c','X_T0') % X_T0c and X_T0 are dataset objects 
X_correct=X_T0c.data; %dynamic data
X_T0=X_T0.data; % static data


%% preprocess the data: centering and scaling
Z1=preprocess(X_correct);
Z1=permute(Z1, [1 3 2]);
Z1=tensor(Z1);
Z2=preprocess(X_T0);
Z2=tensor(Z2);


%% data preparation ---  ACMTF model ----

W1 = tensor(~isnan(Z1.data)); %indicator tensor/matrix: 0s for missing values
W2 = tensor(~isnan(Z2.data));
Z1(find(W1==0))=0;
Z2(find(W2==0))=0;
Z.object{1} = tensor(Z1);
Z.object{2} = tensor(Z2);
Z.miss{1} = W1;
Z.miss{2} = W2;
Z.modes     = {[1 2 3],[1 4]};
Z.size      = [size(Z1) size(Z2,2)];
for p=1:2
    norms(p)    = norm(Z.object{p});
    Z.object{p} = Z.object{p}/norms(p);
end


%% fit ACMTF model
options = ncg('defaults');
options.Display ='final';
options.MaxFuncEvals = 100000;
options.MaxIters     = 10000;
options.StopTol      = 1e-9;
options.RelFuncTol   = 1e-9;
R = 2;
beta     = [1e-3 1e-3]; % l1-penalty parameter for the higher-order tensors
nb_runs  = 10;
goodness_X = zeros(nb_runs,1); %Stores objective function value
for i=1:nb_runs
    if i==1
        [Fac{i}, ~ , out{i}]    = acmtf_opt(Z,R,'init','nvecs','alg_options',options, 'beta', beta, 'alg','ncg');
    else
        [Fac{i}, ~ , out{i}]    = acmtf_opt(Z,R,'init','random','alg_options',options, 'beta', beta, 'alg','ncg');
    end
    goodness_X(i) = out{i}.F; %objective function value
end



%% order the best factor 
[ff, index] = sort(goodness_X,'ascend');
Fac_sorted = Fac(index);
out_sorted = out(index);

%check the stopping condition
if (out{index(1)}.ExitFlag~=2) && (out{index(1)}.ExitFlag~=1)
    flag_stop = true;
else
    flag_stop = false;
end


%% collect all informatin in one variable 
P = 2;
Zhat = Fac_sorted{1};
l_rec = zeros(P, R);
for p = 1:P
    temp        = normalize(Zhat{p});
    l_rec(p,:)  = temp.lambda;
end

for p=1:P
    data.Zhat{p}      = normalize(Zhat{p});
    data.fit(p) = 100- (norm(Z.miss{p}.*tensor(full(data.Zhat{p})-Z.object{p}))^2/norm(Z.miss{p}.*Z.object{p})^2*100);
end
data.out        = out_sorted{1};
data.func_eval  = ff; %value of the objective function for all runs
data.lambda     = l_rec;
data.Z          = Z;




%% check uniqueness
figure
[Fac_f100, T1, T2] = show_spread(R, Fac_sorted, data.func_eval, 0);
set(gca,'fontsize', 20)


%% BMI index
index_Under=find(X_T0c.class{1,2}==4);
index_normal=find(X_T0c.class{1,2}==1);
index_obesity=find(X_T0c.class{1,2}==2);
index_over=find(X_T0c.class{1,2}==3);
index_Under_normal=[index_Under,index_normal];
index_over_obesity=[index_over,index_obesity];
sub_normal=index_Under_normal;
sub_abnormal=index_over_obesity;


%% compute the p-value in terms of BMI group difference
nm_comp=R;
for r=1:nm_comp
    r1=Zhat{1}.U{1}(sub_normal,r);
    r2=Zhat{1}.U{1}(sub_abnormal,r);
    [~, p(r)]=ttest2(r1, r2, 'alpha', 0.05, 'vartype','unequal');
end


%% plot the factors by ACMTF model, 
% each column represents one component of the ACMTF model
% each row represents one mode

nm_comp=R;
time_value=[0.25 0.5 1 1.5 2 2.5 4];


f=figure
set(gcf, 'Position', get(0, 'Screensize'));
k=0;
i=1; % subjects mode
for j=1:nm_comp
    k_temp=1+4*(j-1)+k;
    subplot(nm_comp,4,k_temp)
    xvalue=1:length(Zhat{1}.U{i}(:,1));
    plot(xvalue(sub_normal),Zhat{1}.U{i}(sub_normal,j),'rs','MarkerSize',10,...
        'MarkerEdgeColor','r','MarkerFaceColor','r');
    hold on
    plot(xvalue(sub_abnormal),Zhat{1}.U{i}(sub_abnormal,j),'bs','MarkerSize',10,...
        'MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot([0,length(Zhat{1}.U{i}(:,1))],[0,0],'--k')
    if j==nm_comp
        legend('Lower BMI','Higher BMI')
    end
    ylim([min(min(Zhat{1}.U{i})),max(max(Zhat{1}.U{i}))])
    ylabel(['a_', num2str(j) ],'FontSize', 30,'FontWeight','bold')
    xlabel('subjects')
    set(gca,'Fontsize',20)
end

k=k+1;
i=2; % metabolites (dynamic) mode
for j=1:nm_comp
    k_temp=1+4*(j-1)+k;
    subplot(nm_comp,4,k_temp)
    [x_values,index_meta,plot_color,EdgeColor,FaceColor,plot_shape,Marker_size]=plot_metab_set(X_T0c);
    for kk=1:length(index_meta)
        plot(x_values{kk},Zhat{1}.U{i}(index_meta{kk},j),plot_shape{kk},'MarkerSize',Marker_size{kk},...
            'MarkerEdgeColor',EdgeColor{kk},'MarkerFaceColor',FaceColor{kk})
        hold on
    end
    plot([0,length(Zhat{1}.U{i}(:,1))],[0,0],'--k')
    ylim([min(min(Zhat{1}.U{i})),max(max(Zhat{1}.U{i}))])
    xlim([0,length(Zhat{1}.U{2}(:,1))])
    ylabel(['b_', num2str(j) ],'FontSize', 30,'FontWeight','bold')
    xlabel('metabolites (dynamic)' )
    set(gca,'Fontsize',20)
end
k=k+1;
i=3; % time mode
for j=1:nm_comp
    k_temp=1+4*(j-1)+k;
    subplot(nm_comp,4,k_temp)
    plot(time_value,Zhat{1}.U{3}(:,j),'-o','LineWidth',2.4)
    hold on
    plot([0,length(Zhat{1}.U{i}(:,1))],[0,0],'--k')
    xticks(0:4)
    ylim([min(min(Zhat{1}.U{i})),max(max(Zhat{1}.U{i}))])
    xlim([0 4.5])
    ylabel(['c_', num2str(j) ],'FontSize', 30,'FontWeight','bold')
    xlabel('time(h)')
    set(gca,'Fontsize',20)
end
k=k+1;
i=2; % metabolites (static) mode
for j=1:nm_comp
    k_temp=1+4*(j-1)+k;
    subplot(nm_comp,4,k_temp)
    [x_values,index_meta,plot_color,EdgeColor,FaceColor,plot_shape,Marker_size]=plot_metab_set(X_T0c);
    for kk=1:length(index_meta)
        plot(x_values{kk},Zhat{2}.U{i}(index_meta{kk},j),plot_shape{kk},'MarkerSize',Marker_size{kk},...
            'MarkerEdgeColor',EdgeColor{kk},'MarkerFaceColor',FaceColor{kk})
        hold on
    end
    plot([0,length(Zhat{2}.U{2}(:,1))],[0,0],'--k')
    ylim([min(min(Zhat{2}.U{i})),max(max(Zhat{2}.U{2}))])
    xlim([0,length(Zhat{2}.U{i}(:,1))])
    ylabel(['d_', num2str(j) ],'FontSize', 30,'FontWeight','bold')
    xlabel('metabolites (static)' )
    set(gca,'Fontsize',20)
    if j==nm_comp
        legend('-HDL','XL-HDL','L-HDL','M-HDL','S-HDL',...
            '-IDL','-LDL','L-LDL','M-LDL','S-LDL',...
            'Rest','Total','-VLDL','XXL-VLDL','XL-VLDL','L-VLDL','M-VLDL','S-LDL','XS-VLDL')
    end
end
