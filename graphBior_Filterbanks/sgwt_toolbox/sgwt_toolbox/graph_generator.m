function G = graph_generator(loc,thresh)
err = dist(loc);%err1 + err2;
G = err < thresh;
G = G- diag(diag(G));
G = G+G';
G = G>0;
