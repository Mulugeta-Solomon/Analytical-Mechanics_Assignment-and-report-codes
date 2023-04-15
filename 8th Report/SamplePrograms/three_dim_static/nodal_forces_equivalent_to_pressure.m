function nodal_forces = nodal_forces_equivalent_to_pressure( pressure, faces, npoints, positions )
    nodal_forces = zeros(3, npoints);
    nfaces = size(faces,1);
    for p=1:nfaces
        i=faces(p,1); j=faces(p,2); k=faces(p,3);
        ri = positions(:,i);
        rj = positions(:,j);
        rk = positions(:,k);
        fequiv = (1/6)*pressure*cross(rj-ri,rk-ri);
        nodal_forces(:,i) = nodal_forces(:,i) + fequiv;
        nodal_forces(:,j) = nodal_forces(:,j) + fequiv;
        nodal_forces(:,k) = nodal_forces(:,k) + fequiv;
    end
end
