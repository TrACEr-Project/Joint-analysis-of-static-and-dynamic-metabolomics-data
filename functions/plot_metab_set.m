


% This function reorder the metabolite indexes to put the same group of metabolites toegther in the plot
% Also, the same maker color was assgined to the same group of metabolites in the plot
% Also, the same marker size was assigned to the same size of metabolites


function [x_values,index_meta_perm,plot_color_perm,EdgeColor_perm,FaceColor_perm,plot_shape_perm Marker_size_perm]=plot_metab_set(NMR_remove)
%%give marker color and size according to metabolite size and group
k=0;
shapes={'o','diamond','hexagram','pentagram','square','+','>'};
colors={'r','b','g','k','c','m'};
Marker_size_set={13-5 16-5 14-5 12-5 18-5 10-5  20-5};
unique_meta_type=unique(NMR_remove.class{3,3});
for ii=1:length(unique_meta_type)
    index_meta_type=find(NMR_remove.class{3,3}==unique_meta_type(ii));
    NMR_temp=NMR_remove(:,:,index_meta_type);
    unique_meta_size=unique(NMR_temp.class{3,4});
    for j=1:length(unique_meta_size)
        k=k+1;
        index_meta{1,k}=intersect(index_meta_type,index_meta_type(NMR_temp.class{3,4}==unique_meta_size(j)));
        plot_color{1,k}=colors{ii};
        EdgeColor{1,k}=colors{ii};
        FaceColor{1,k}=colors{ii};
        plot_shape{1,k}=shapes{j};
        Marker_size{1,k}=Marker_size_set{j};
        
    end
    clear index_meta_type NMR_temp unique_meta_size
end

%%reorder matbolite indexes to put the same group of metabolites together in the plot
index_orig=1:length(index_meta);
index_perm=[[1 5 2 3 4 ],6:12, [13 19 17 14 15 16 18]];

for i=1:length(index_perm)
    index_meta_perm{1,i}=index_meta{1,index_perm(i)};
    plot_color_perm{1,i}=plot_color{1,index_perm(i)};
    EdgeColor_perm{1,i}=EdgeColor{1,index_perm(i)};
    FaceColor_perm{1,i}=FaceColor{1,index_perm(i)};
    plot_shape_perm{1,i}=plot_shape{1,index_perm(i)};
    Marker_size_perm{1,i}=Marker_size{1,index_perm(i)};
    x_value(i)=length(index_meta_perm{1,i});
end
x_value=cumsum(x_value);
x_value=[0, x_value];
for i=1:length(index_perm)
    x_values{1,i}= x_value(i)+1:x_value(i+1);
end

end


