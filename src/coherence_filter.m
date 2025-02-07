function matches = coherence_filter(candidate_matches, tar_coordinates, ref_coordinates, tar_control, ref_control)
    
    P0 = [ones(size(tar_control,2),1);0.0001*ones(size(tar_coordinates,2),1)];
    X = [tar_control, tar_coordinates]';
    Y = [ref_control, ref_coordinates]';

    [Xn, Yn] = normr(X,Y);
    Ind_nan = ceil(find(isnan(Xn)|isnan(Yn))/2);
    if numel(Ind_nan)>0
        Xn(Ind_nan,:)=[];Yn(Ind_nan,:)=[];
        P0(Ind_nan)=[];
    end

    conf = SLC_init([]);
    [indx, ~, ~] = SLC(Xn, Yn, conf,   P0);

    if numel(Ind_nan)>0
        a = 1:size(X,1);
        a(Ind_nan)=[];indx = a(indx);
    end

    matches = candidate_matches(indx,:);

end