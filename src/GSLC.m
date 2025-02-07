function [inits, seeds, candidates_all, matches] = GSLC(tar_feat, tar_desc, ref_feat, ref_desc, tar_size, ref_size)

% compute descriptor distances
D = vl_alldist(single(tar_desc),single(ref_desc));
[N, M] = size(D);

% initial correspondences
k =2; tau =1;

[scorek, min_indk] = mink(D, k, 2); % top-k for each row

NN=[(1:N)',min_indk(:,1)];
NN1=NN((scorek(:,1)./scorek(:,2))<tau,:);
[~, min_ind2] = min(D, [], 1); % top-1 for each column
NN2 = [min_ind2',(1:M)'];

inits = intersect(NN1, NN2,'rows'); % MNN when tau =1

% grid filter
num_init = size(inits,1);
n_c = max(min(floor(sqrt(num_init/10)),20),8); 
alpha = 1;

interval1 = [tar_size(2),tar_size(1)]./n_c;
interval2 = [ref_size(2),ref_size(1)]./n_c;
coord1_init = tar_feat(1:2,inits(:,1));
coord2_init = ref_feat(1:2,inits(:,2));

grid1_init = ceil(coord1_init'./repmat(interval1, num_init,1));
grid2_init = ceil(coord2_init'./repmat(interval2, num_init,1));
grid1_uni_init = unique(grid1_init,'rows'); % indexes of grids contained
grid2_uni_init = unique(grid2_init,'rows');

S_threshold = alpha*sqrt(num_init/(n_c*n_c));

while true
    ind = []; grid_pairs = [];seeds_gridpair = {};
    for i = 1:size(grid1_uni_init,1)
        for j = 1:size(grid2_uni_init,1)
            ind_grid = find(sum(abs([grid1_init, grid2_init] - [grid1_uni_init(i,:), grid2_uni_init(j,:)]),2)==0);
            S_grid = length(ind_grid);
            if S_grid>S_threshold
                ind = [ind; ind_grid]; % indexs of seeds on inits
                seeds_gridpair{end+1} = inits(ind_grid,:); % coordinates of seeds on each valid grid pair
                grid_pairs = [grid_pairs;[grid1_uni_init(i,:)],[grid2_uni_init(j,:)]]; % indexes of grid pairs
            end
        end
    end
    seeds = inits(ind,:);
    if isempty(seeds)
        S_threshold = S_threshold-1;
        if S_threshold<0
            % disp('tau too low.');
            candidates_all = [];matches = [];
            return;
        end
    else
        break;
    end
end

% indexes of features other than seeds on target image
rest_tar_ind = 1:size(tar_feat,2); 
rest_tar_ind(seeds(:,1)) = [];

miu = ceil((n_c)/5);

% divide grid pairs for blocks
grid_pair_vector = grid_pairs(:,3:4)-grid_pairs(:,1:2); % vector difference 
tag_gridpairs = [0;find(sum(abs(diff(grid_pair_vector))>2*miu,2));size(grid_pair_vector,1)];
if tag_gridpairs(end)==tag_gridpairs(end-1)
    tag_gridpairs(end)=[];
end
num_blocks = length(tag_gridpairs)-1;

min_boundary = zeros(num_blocks,4);
max_boundary = zeros(num_blocks,4);
% indexes of ceil pairs for features other than seeds on target image
rest_tar_cell = ceil(tar_feat(1:2,rest_tar_ind)'./repmat(interval1, length(rest_tar_ind),1));

matches = [];
candidates_all = [];
% find candidates for each block
for i=1:num_blocks
    gridpairs_block = grid_pairs(tag_gridpairs(i)+1:tag_gridpairs(i+1),:);

    if size(gridpairs_block,1)>1
        min_boundary(i,:) = max(min(gridpairs_block)-miu*ones(1,4),ones(1,4));
        max_boundary(i,:) = min(max(gridpairs_block)+miu*ones(1,4),n_c*ones(1,4));
    else
        min_boundary(i,:) = max(gridpairs_block-miu*ones(1,4),ones(1,4));
        max_boundary(i,:) = min(gridpairs_block+miu*ones(1,4),n_c*ones(1,4));
    end

    % find rest features in left region
    indexs_left_region = rest_tar_ind(logical(all(rest_tar_cell>=min_boundary(i,1:2),2).*all(rest_tar_cell<=max_boundary(i,1:2),2)));

    indexs_right_ind = min_indk(indexs_left_region,1);
    candidate_ref = ceil(ref_feat(1:2,indexs_right_ind)'./repmat(interval2,length(indexs_right_ind),1));
    indicator = find(all(candidate_ref>=min_boundary(i,3:4),2).*all(candidate_ref<=max_boundary(i,3:4),2));
    if ~isempty(indicator)
        indexs_right_ind = indexs_right_ind(indicator);
        % add indexes of candidate matches
        candidates = [indexs_left_region(indicator)',indexs_right_ind];
    end

    if size(candidates,1)<2
        matches = [matches;candidates];
        continue;
    end

    % find seeds in this block pairs
    seeds_block = [];
    for m=tag_gridpairs(i)+1:tag_gridpairs(i+1)
        seeds_block = [seeds_block;seeds_gridpair{m}];
    end

    tar_coordinates = tar_feat(1:2,candidates(:,1));
    ref_coordinates = ref_feat(1:2,candidates(:,2));
    tar_control = tar_feat(1:2,seeds_block(:,1));
    ref_control = ref_feat(1:2,seeds_block(:,2));
    candidates_all_block = [seeds_block;candidates];
    filtered = coherence_filter(candidates_all_block, tar_coordinates, ref_coordinates, tar_control, ref_control);
    matches = [matches;filtered];
    candidates_all = [candidates_all;candidates_all_block];

end

matches = unique(matches,'rows');
candidates_all = unique(candidates_all,'rows');

end