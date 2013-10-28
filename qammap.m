function qammap(M,mode)
% Rectangular QAM Modulator Baseband ブロックにおける
% コンスタレーションプロットに対応するマッピングパターンを表示
%
% qammod(M,mode)
% M : M値QAM
% mode : コンスタレーションオーダー(bin or gray) 

x = [0:M-1];
y = qammod(x,M,[],mode);
scatterplot(y)
hold on;
for jj=1:length(y)
   text(real(y(jj)),imag(y(jj)),[' ' num2str(jj-1)]);
end
hold off