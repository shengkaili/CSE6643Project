clear all;
nodeArray=10:10:200;
for i=1:length(nodeArray)
    nodes=nodeArray(i)
    [K,F,ymax,id] = mkbeamproblem(nodes);
    condK(i)=cond(K,2);
%     alldisp = K\F;
%     err = ymax - alldisp(id);
    
%     KS=sparse(K); FS=sparse(F);
%     alldisp = KS\FS;
%     errS = ymax - alldisp(id);
end
figure,plot(nodeArray,condK,'o-');
xlabel('nodes'); ylabel('cond_2(K)');

figure,plot(log(nodeArray),log(condK),'o-');
xlabel('log(nodes)'); ylabel('log(cond_2(K))');