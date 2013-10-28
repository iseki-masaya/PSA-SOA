function qammap(M,mode)
% Rectangular QAM Modulator Baseband �u���b�N�ɂ�����
% �R���X�^���[�V�����v���b�g�ɑΉ�����}�b�s���O�p�^�[����\��
%
% qammod(M,mode)
% M : M�lQAM
% mode : �R���X�^���[�V�����I�[�_�[(bin or gray) 

x = [0:M-1];
y = qammod(x,M,[],mode);
scatterplot(y)
hold on;
for jj=1:length(y)
   text(real(y(jj)),imag(y(jj)),[' ' num2str(jj-1)]);
end
hold off