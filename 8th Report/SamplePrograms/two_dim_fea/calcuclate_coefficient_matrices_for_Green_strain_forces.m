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