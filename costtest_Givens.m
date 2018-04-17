clear all;
nodeArray=10:10:300;
tM = nan(length(nodeArray),1);
tGf = nan(length(nodeArray),1);
tGf_bw = nan(length(nodeArray),1);
errM = nan(length(nodeArray),1);
errG_full = nan(length(nodeArray),1);
errG_full_bw = nan(length(nodeArray),1);
for i=1:length(nodeArray)
    nodes=nodeArray(i);
    [K,F,ymax,id] = mkbeamproblem(nodes);
    tic;alldisp = K\F;tM(i) = toc;
    errM(i) = ymax-alldisp(id);
    if i == 1
        t0 = tM(i);
    end
    tic;[xf,Af] = Givens_full(K,F,0);tGf(i) = toc;
    errG_full(i) = ymax-xf(id);
%     K_COO = full2sparse(K,'COO');
%     F_COO = full2sparse(F,'COO');
%     tic;[x_COO,A_COO] = Givens_COO(K_COO,F_COO);tGC(i) = toc;
%     xC = sparse2full(x_COO,'COO');
%     errG_COO(i) = ymax-xC(id);
    tic;[xf_bw,Af_bw] = Givens_full(K,F,3);tGf_bw(i) = toc;
    errG_full_bw(i) = ymax-xf_bw(id);
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
xlim([min(nodeArray) max(nodeArray)])
set(gca,'FontSize',FS)

subplot(2,2,2)
loglog(nodeArray,abs(errM),'o-');
hold on
loglog(nodeArray,abs(errG_full),'o-');
hold off
xlabel('nodes'); ylabel('|err|');
legend('Matlab in-built','Givens','location','NorthWest');
grid on
xlim([min(nodeArray) max(nodeArray)])
set(gca,'FontSize',FS)

subplot(2,2,3)
plot(nodeArray,tM/t0,'o-');
hold on
plot(nodeArray,tGf/t0,'o-');
% plot(nodeArray,tGC/t0,'o-');
plot(nodeArray,tGf_bw/t0,'o-');
hold off
xlabel('nodes'); ylabel('relative time cost');
legend('Matlab in-built','Givens (no bandwith information)','Givens (known bandwith)','location','NorthWest');
grid on
xlim([min(nodeArray) max(nodeArray)])
set(gca,'FontSize',FS)

subplot(2,2,4)
loglog(nodeArray,tM/t0,'o-');
hold on
loglog(nodeArray,tGf/t0,'o-');
% plot(nodeArray,log(tGC/t0),'o-');
loglog(nodeArray,tGf_bw/t0,'o-');
hold off
xlabel('nodes'); ylabel('relative time cost');
legend('Matlab in-built','Givens (no bandwith information)','Givens (known bandwith)','location','NorthWest');
grid on
xlim([min(nodeArray) max(nodeArray)])
set(gca,'FontSize',FS)