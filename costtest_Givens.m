clear all;
nodeArray=10:10:300;
tM = nan(length(nodeArray),1);
tGf = nan(length(nodeArray),1);
tGC = nan(length(nodeArray),1);
errM = nan(length(nodeArray),1);
errG_full = nan(length(nodeArray),1);
errG_COO = nan(length(nodeArray),1);
for i=1:length(nodeArray)
    nodes=nodeArray(i);
    [K,F,ymax,id] = mkbeamproblem(nodes);
%     condK(i)=cond(K,2);
    tic;alldisp = K\F;tM(i) = toc;
    errM(i) = ymax-alldisp(id);
    if i == 1
        t0 = tM(i);
    end
    tic;[xf,Af] = Givens_full(K,F,3);tGf(i) = toc;
    errG_full(i) = ymax-xf(id);
    K_COO = full2sparse(K,'COO');
    F_COO = full2sparse(F,'COO');
    tic;[x_COO,A_COO] = Givens_COO(K_COO,F_COO);tGC(i) = toc;
    xC = sparse2full(x_COO,'COO');
    errG_COO(i) = ymax-xC(id);
    
    
%     KS=sparse(K); FS=sparse(F);
%     alldisp = KS\FS;
%     errS = ymax - alldisp(id);
end
FS = 18;
figure(1)
subplot(2,2,1)
plot(nodeArray,abs(errM),'o-');
hold on
plot(nodeArray,abs(errG_full),'o-');
hold off
xlabel('nodes'); ylabel('|err|');
legend('Matlab in-built','Givens','location','NorthWest');
grid on
set(gca,'FontSize',FS)

subplot(2,2,2)
plot(log(nodeArray),log(abs(errM)),'o-');
hold on
plot(log(nodeArray),log(abs(errG_full)),'o-');
hold off
xlabel('log(nodes)'); ylabel('log(|err|)');
legend('Matlab in-built','Givens','location','NorthWest');
grid on
set(gca,'FontSize',FS)

subplot(2,2,3)
plot(nodeArray,tM/t0,'o-');
hold on
plot(nodeArray,tGf/t0,'o-');
plot(nodeArray,tGC/t0,'o-');
hold off
xlabel('nodes'); ylabel('relative time cost');
legend('Matlab in-built','Givens for full matrix','Givens for COO matrix','location','NorthWest');
grid on
set(gca,'FontSize',FS)

subplot(2,2,4)
plot(nodeArray,log(tM/t0),'o-');
hold on
plot(nodeArray,log(tGf/t0),'o-');
plot(nodeArray,log(tGC/t0),'o-');
hold off
xlabel('log(nodes)'); ylabel('log(relative time cost)');
legend('Matlab in-built','Givens for full matrix','Givens for COO matrix','location','NorthWest');
grid on
set(gca,'FontSize',FS)