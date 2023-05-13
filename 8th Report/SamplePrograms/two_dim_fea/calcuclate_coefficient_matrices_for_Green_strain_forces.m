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