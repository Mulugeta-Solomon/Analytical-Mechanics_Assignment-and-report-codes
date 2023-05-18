function [ coef_order_one, coef_order_two, c2i, c2j, coef_order_three, c3i, c3j, c3k ] ...
    = calcuclate_coefficient_matrices_for_Green_strain_forces (body)

np = body.numNodalPoints;
np1 = np + 1;
id2p = @(i,j) (np1*i+j);
id3p = @(i,j,k) (np1*np1*i+np1*j+k);

id2 = [];
for k=1:body.numTriangles
    tri = body.Triangles(k);
    vs = tri.Vertices;
    i = vs(1); j = vs(2); k = vs(3);
    id2 = [ id2 ; ...
        id2p(i,i); id2p(i,j); id2p(i,k);
        id2p(j,i); id2p(j,j); id2p(j,k);
        id2p(k,i); id2p(k,j); id2p(k,k); ];
    id2 = sort(unique(id2));
end

id2 = sort(unique(id2));
c2j = rem(id2,np1);
c2i = (id2-c2j)/np1;

id3 = [];
for k=1:body.numTriangles
    tri = body.Triangles(k);
    vs = tri.Vertices;
    i = vs(1); j = vs(2); k = vs(3);
    id3 = [ id3 ; ...
        id3p(i,i,i); id3p(i,i,j); id3p(i,i,k);
        id3p(i,j,i); id3p(i,j,j); id3p(i,j,k);
        id3p(i,k,i); id3p(i,k,j); id3p(i,k,k);
        id3p(j,i,i); id3p(j,i,j); id3p(j,i,k);
        id3p(j,j,i); id3p(j,j,j); id3p(j,j,k);
        id3p(j,k,i); id3p(j,k,j); id3p(j,k,k);
        id3p(k,i,i); id3p(k,i,j); id3p(k,i,k);
        id3p(k,j,i); id3p(k,j,j); id3p(k,j,k);
        id3p(k,k,i); id3p(k,k,j); id3p(k,k,k) ];
    id3 = sort(unique(id3));
end


id3 = sort(unique(id3));
c3k = rem(id3,np1);
c3ij = (id3-c3k)/np1;
c3j = rem(c3ij,np1);
c3i = (c3ij-c3j)/np1;

%%%

np2 = size(id2,1);
np3 = size(id3,1);

coef_order_one = zeros(2*np, 2*np);
for k=1:body.numTriangles
    tri = body.Triangles(k);
    vs = tri.Vertices;
    i = vs(1); j = vs(2); k = vs(3);
    a = tri.vector_a; b = tri.vector_b;
    mat = tri.lambda * tri.Area * tri.Thickness * [ a*a', a*b'; b*a', b*b' ] + ...
          tri.mu     * tri.Area * tri.Thickness * [ 2*a*a'+b*b', b*a'; a*b', 2*b*b'+a*a' ];
    ijk = [ i, j, k ];
    row = [ ijk, np+ijk ];
    col = [ ijk, np+ijk ];
    coef_order_one(row,col) = coef_order_one(row,col) + mat;
end

%coef_order_one;
%term_order_one = coef_order_one * [ c_u; c_v ];

odot = @(x,y) [ x(1)*y; x(2)*y; x(3)*y ];
odot3 = @(x,y,z) [ x(1)*y(1)*z; x(1)*y(2)*z; x(1)*y(3)*z; ...
                   x(2)*y(1)*z; x(2)*y(2)*z; x(2)*y(3)*z; ...
                   x(3)*y(1)*z; x(3)*y(2)*z; x(3)*y(3)*z ];

coef_order_two = zeros(2*np, 3*np2);
for k=1:body.numTriangles
    tri = body.Triangles(k);
    vs = tri.Vertices;
    i = vs(1); j = vs(2); k = vs(3);
    a = tri.vector_a; b = tri.vector_b;
    cpl211 = (3/2)*a*odot(a,a)' + (1/2)*a*odot(b,b)' + b*odot(a,b)';
    cpl212 = a*odot(a,b)' + b*odot(b,b)';
    cpl213 = (1/2)*a*odot(a,a)' + (1/2)*a*odot(b,b)';
    cpl221 = (1/2)*b*odot(b,b)' + (1/2)*b*odot(a,a)';
    cpl222 = b*odot(a,b)' + a*odot(a,a)';
    cpl223 = (3/2)*b*odot(b,b)' + (1/2)*b*odot(a,a)' + a*odot(b,a)';
    cpm211 = 3*a*odot(a,a)' + a*odot(b,b)' + 2*b*odot(a,b)';
    cpm212 = b*odot(a,a)' + a*odot(b,a)' + 2*b*odot(b,b)';
    cpm213 = a*odot(a,a)' + b*odot(a,b)';
    cpm221 = b*odot(b,b)' + a*odot(b,a)';
    cpm222 = a*odot(b,b)' + b*odot(b,a)' + 2*a*odot(a,a)';
    cpm223 = 3*b*odot(b,b)' + b*odot(a,a)' + 2*a*odot(b,a)';
    mat = tri.lambda * tri.Area * tri.Thickness * [ cpl211, cpl212, cpl213; cpl221, cpl222, cpl223 ] + ...
          tri.mu     * tri.Area * tri.Thickness * [ cpm211, cpm212, cpm213; cpm221, cpm222, cpm223 ];
    ijk = [ i, j, k ];
    row = [ ijk, np+ijk ];
    ijk2 = [ find(id2==id2p(i,i)), find(id2==id2p(i,j)), find(id2==id2p(i,k)), ...
             find(id2==id2p(j,i)), find(id2==id2p(j,j)), find(id2==id2p(j,k)), ...
             find(id2==id2p(k,i)), find(id2==id2p(k,j)), find(id2==id2p(k,k)) ];
    col = [ ijk2, np2+ijk2, 2*np2+ijk2 ];
    coef_order_two(row,col) = coef_order_two(row,col) + mat;
end

%coef_order_one;
%term_order_one = coef_order_one * [ c_u; c_v ];

odot = @(x,y) [ x(1)*y; x(2)*y; x(3)*y ];
odot3 = @(x,y,z) [ x(1)*y(1)*z; x(1)*y(2)*z; x(1)*y(3)*z; ...
                   x(2)*y(1)*z; x(2)*y(2)*z; x(2)*y(3)*z; ...
                   x(3)*y(1)*z; x(3)*y(2)*z; x(3)*y(3)*z ];

