function [sampling_inds, sampling_coords, point_dists] = getsamplingpoints(vol_coords, X)

    [point_dists, sampling_inds] = sort(sqrt((vol_coords(:,1)-X(1)).^2+(vol_coords(:,2)-X(2)).^2+...
        (vol_coords(:,3)-X(3)).^2));
    
    point_dists = point_dists(1:128);

    sampling_inds = sampling_inds(1:128);

    point_dists = point_dists/max(point_dists);

    sampling_coords = vol_coords(sampling_inds,:); %nx3 vector of sampling point coordinates
end

