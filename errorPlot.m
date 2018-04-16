clear all;
nodeArray=10:10:300;
for i=1:length(nodeArray)
    nodes=nodeArray(i);
    [K,F,ymax,id] = mkbeamproblem(nodes);
%     condK(i)=cond(K,2);
    alldisp = K\F;
    errM(i) = ymax - alldisp(id);
    errCh(i)= ymax - sCholesky(nodes);
    
%     KS=sparse(K); FS=sparse(F);
%     alldisp = KS\FS;
%     errS = ymax - alldisp(id);
end
figure,plot(nodeArray,abs(errM),'o-');
hold on,plot(nodeArray,abs(errCh),'o-');
xlabel('nodes'); ylabel('|err|');
legend('Matlab in-built','sparse Cholesky');

figure,plot(log(nodeArray),log(abs(errM)),'o-');
hold on,plot(log(nodeArray),log(abs(errCh)),'o-');
xlabel('log(nodes)'); ylabel('log(|err|)');
legend('Matlab in-built','sparse Cholesky');