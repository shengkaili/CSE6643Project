for nodes = [2, 10, 100, 500]
   
    [K,F,ymax,id] = mkbeamproblem(nodes);
    
    alldisp = K\F;
    
    err = ymax - alldisp(id)

end