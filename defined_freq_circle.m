%% ���[�_�[�`���[�g��\������v���O����
%�C���p���X������ǂݍ���Ŏw��̎��g���ɑ΂���1�����̃`���[�g��`�悷��v���O����
%�C���p���X�����̃t�@�C���́@[��ނ��킩��悤�ȕ���]_[�p�x�̒l(���p)]_imp.wav�@�œ���@��ӏ��̃t�H���_�Ɋi�[���Ă���
%���s�̓C���p���X������ۑ����Ă���t�H���_��I��
%20170824 Sakaguchi�C��

%%
clear;
close all;
clc;
%% --���̃t�@�C����I�����鏈��

m = 360;
num = m-10;
f = 10000;%���o�����g����ݒ�
%for h = from:1000:to
for v = 0:10:num%0����350�܂ŉ�(10�����Ƃ̃f�[�^���擾)
    
    tmp = v/10;%10���Âf�[�^��������̂Ńf�[�^����10�Ŋ��������ɂȂ�
    s = num2str(v);%v�𕶎��z��^�ɕϊ�(fname�Ŏg������)
    fname =['S_' s '_imp.wav'];
    PathName = './';
    T = strcat(PathName, fname);
    [data,Fs] = audioread(T);%�f�[�^�̓ǂݍ���
    disp(v)%���ǂ̊p�x�����s���Ă���̂����R�}���h�ɕ\��
    data = data(:,1);%���m�����f�[�^���ق����̂�L�����f�[�^�����o��

    %%FFT

    len_y = length(data);
    %�����̒���
    if (len_y < 2*Fs)
        add1 = 2*Fs - len_y;
        z_add1 = zeros(add1,1);
        yn = vertcat(data,z_add1);
    else
        yn = data;
        
    end
    len_yn = length(yn);
    Y = fft(yn);
    %�O���t�`��̂��߂̏���

    %���g���̈�̏���
    len_Y = length(Y);%���g������\�̌v�Z
    len_Y = len_Y - rem(len_Y,2);
    time_Y = (1:len_Y)/Fs;
    Y = Y(1:len_Y);
    f_vin = Fs/len_Y;
    fre = 0:f_vin:Fs/2; %���g������\�̌v�Z�����܂�

    %�w��̎��g�������o��
    deltaF = ( 0:len_Y-1).*Fs./ len_Y; %����
    d = deltaF(3)-deltaF(2);%deltaF�̍������߂�(deltaF�͓�������)
    n = fix(( f + d ) / d);%n�ԖڂɎ��o�����g��������̂������߂Ă���,n�͐����Ȃ̂�fix�Ŋۂߍ���
    
    
    %deltaF = deltaF(n);
    avr = abs(Y);%���f���Œl���o��̂ł��̐�Βl�����߂�
    
    %�������g���̃f�[�^�����т���
        if(v == 1)
            avr_h = zeros(1,tmp);%0�z����������ď�����
            avr_h(1,tmp+1) = avr(n);%tmp��0����n�܂邪�A�z���1����n�܂�̂�tmp+1����X�^�[�g
        else
            avr_h(1,tmp+1) = avr(n);
        end
        
    avr_av = 20.*log10(avr_h./(max(avr_h)));
    
end
disp(f)%�ǂ̎��g�����v�Z�����̂��Y��Ȃ��悤�ɃR�}���h�ɕ\��
%dffr = dffr + 1;
%end
%% ���߂ă`���[�g�\������Ƃ��Ɏ��s
flg =1;%�t���O�w��
%% �`���[�g�\��
%�}�ɕ\������
if(flg == 1)
    figure(1)
    polarplot(avr_av,'LineWidth',1.0);
    %f_st =num2str(f);
    %fre =[f_st 'Hz'];
    %legend({fre})
    title('Rader Chart')
    hold on
    ax = gca;
    pax = gca;
    rticks(pax,[-20 -15 -10 -5 0])
    ax.ThetaTickLabel = {'0��','30��','60��','90��','120��','150��','180��','210��','240��','270��','300��','330��'};
    ax.RTickLabel = {'-20dB','','-10dB','','0dB'};
    ax.ThetaZeroLocation = 'top';
    ax.Color = 'w';
    ax.ThetaDir = 'clockwise';
    %ax.RDir = 'reverse';
    
    %������Ɍ����₷���悤�ɂ���ݒ�
    ax.GridAlpha = 0.7;
    ax.FontSize = 13;
    
    %�\���͈͂̐ݒ�
    ax.RLim = [-25 0];
    ax.ALim = [0 3];
    flg = 2;%�t���O��2�ɕύX
else
    avr_av1 = avr_av;
    polarplot(avr_av1,'-.','LineWidth',1.0);
    %f_st =num2str(f);
    %fre =[f_st 'Hz'];
    %legend({fre})
    hold on
end

%%
lgd = legend({'5kHz','6kHz','7kHz','8kHz','9kHz','10kHz'},'FontSize',12);
%title(lgd);