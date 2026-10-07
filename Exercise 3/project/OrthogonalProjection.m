function projected_point = OrthogonalProjection(point,restriction_boundaries,varsvect)
    %point & varsvect are column vectors
    %restriction boundaries are set horizontically for each variable:
    %x1-->rb(1,:), x2-->rb(2,:)
    
    restriction = [varsvect<=restriction_boundaries(:,1),varsvect>=restriction_boundaries(:,2),...
        ~(varsvect<=restriction_boundaries(:,1)|varsvect>=restriction_boundaries(:,2))];
    restriction_handle = matlabFunction(restriction,'Vars',varsvect);
    pointinparg = num2cell(point);
    projected_point = feval(restriction_handle,pointinparg{:});
    
    restriction_boundaries = horzcat(restriction_boundaries,point);
    restriction_boundaries = transpose(restriction_boundaries);
    projected_point = restriction_boundaries(transpose(projected_point));
end