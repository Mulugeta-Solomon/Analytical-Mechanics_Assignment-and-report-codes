function draw_mass_and_springs (body, xc, connect, disps)
    arguments
        body, xc, connect;
        disps = zeros(2, body.numNodalPoints);
    end
    
    persistent npoints xn;
    
    if isempty(npoints)
        npoints = body.numNodalPoints;
        for k=1:npoints
            xn = [ xn, body.NodalPoints(k).Coordinates ];
        end
    end
    
    for i=1:8
        j = connect(i);
        rj = xn(:,j) + disps(:,j);
        plot( [ xc(1); rj(1) ], [ xc(2); rj(2) ], 'k' );
    end
    plot(xc(1),xc(2),'ok', 'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerSize', 9);
end
