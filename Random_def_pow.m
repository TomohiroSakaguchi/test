%% ����̎��g���͈͂̊g�U���˓���(����)��\������v���O����
%�C���p���X������ǂݍ���Ŏw��̎��g���ɑ΂���1�����̃`���[�g��`�悷��v���O����
%����Ώ̐�p
%�p���ƂɃt�H���_�����D�����K�w�ɓ���Ă����D�p�͐��ʂ�0���Ƃ������ɃR�l�N�^�����˂��o��������̊p�x�C���΂ɓ˂��o��������̊p�x
%�C���p���X�����̃t�@�C���́@[��ނ��킩��悤�ȕ���]_[���ʊp(���p)]_[�p(���p)]_imp.wav�@�œ���@��ӏ��̃t�H���_�Ɋi�[���Ă���
%���s�̓C���p���X������ۑ����Ă���K�w��I��
%20180116 Sakaguchi�쐬

%%
clear;
close all;
clc;
%% --���̃t�@�C����I�����鏈��
flg = 0;
for freq = 5000:1000:15000
    disp(freq)
    m = 360;
    num = m-10;
    for deg = -90:10:90 %�p�̕ϐ�
    d = num2str(deg);
    if(abs(deg) <= 50)
        mad = 10;
    elseif(abs(deg) <= 70)
        mad = 30;
    elseif(abs(deg) == 80)
        mad = 90;
    elseif(abs(deg) == 90)
        mad = 10;
        num = 0;
    end
        for v = 0:mad:num%���ʊp�̕ϐ��@0����350�܂ŉ�(10�����Ƃ̃f�[�^���擾)
            s = num2str(v);%v�𕶎��z��^�ɕϊ�(fname�Ŏg������)
            fname =['U87_' s '_' d '_imp.wav'];
            PathName = ['./' d '/'];
            T = strcat(PathName, fname);
            [data,Fs] = audioread(T);%�f�[�^�̓ǂݍ���
            disp(fname)%���ǂ̊p�x�����s���Ă���̂����R�}���h�ɕ\��
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

            %�w��̎��g�������o��
            len_yn = length(yn);
            deltaF = ( 0:len_yn-1).*Fs./ len_yn; %����
            dm = deltaF(3)-deltaF(2);%deltaF�̍������߂�(deltaF�͓�������)
            n = fix(( freq + dm ) / dm);%n�ԖڂɎ��o�����g��������̂������߂Ă���,n�͐����Ȃ̂�fix�Ŋۂߍ�
            
            Y = abs(fft(yn));
            p = Y(n);%���f���Œl���o��̂ł��̐�Βl�����߂�
             %�O���t�`��̂��߂̏���

            %���g���̈�̏���
%             len_Y = length(Y);%���g������\�̌v�Z
%             len_Y = len_Y - rem(len_Y,2);
%             time_Y = (1:len_Y)/Fs;
%             Y = Y(1:len_Y);
%             f_vin = Fs/len_Y;
%             fre = 0:f_vin:Fs/2; %���g������\�̌v�Z�����܂�

             %���\�ʐς̌���
            r = 2.0;
            S = 2*pi*r*r*sin(2*deg*pi/180)*cos(deg/2*pi/180)/2.5847;
            S_all = 4*pi*r*r;
            %deltaF = deltaF(n);

            %�ʐςɂ��d�݂Â�
            if(deg <= 50)
                p_sq = p*p*((S/36)/S_all);
            elseif(deg <= 70)
                p_sq = p*p*((S/12)/S_all);
            elseif(deg == 80)
                p_sq = p*p*((S/4)/S_all);
            elseif(deg == 90)
                p_sq= p*p*(S/S_all);
            end
            %�������g���̃f�[�^�����Z���Ă���
            if(flg == 0)
                p_sum = 0;%�ꉞ������
                p_sum = p_sum + p_sq;
                flg = 1;
            else
                p_sum = p_sum + p_sq;
            end
        end
       pow = sqrt(p_sum);
       if(abs(deg) == 90)
           num = m-10;
       end
    end
    if(freq == 5000)
        tmp = 0;
        POW = zeros(1,11);
        POW(1,tmp+1) = pow;
        tmp = tmp+1;
    else
        POW(1,tmp+1) = pow;
        tmp = tmp+1;
    end
end

%% �X�s�[�J�̓����̈����Z
[ file_sp, path] = uigetfile('*.wav', '�X�s�[�J�̃C���p���X������I�����Ă�������');%�C���p���X������ǂݍ���
SP = strcat( path, file_sp);
[ sp_imp, fs_sp] = audioread(SP);
len_sp_imp = length(sp_imp);
%�����̒���
if (len_sp_imp < 2*fs_sp)
    add1 = 2*fs_sp - len_sp_imp;
    z_add1 = zeros(add1,1);
    sp = vertcat(sp_imp,z_add1);
else
    sp = sp_imp;

end
len_sp = length(sp);
deltaF = ( 0:len_sp-1).*fs_sp./ len_sp; %����
dm = deltaF(3)-deltaF(2);%deltaF�̍������߂�(deltaF�͓�������)

SP = abs(fft(sp));
for freq = 5000:1000:15000
    n = fix(( freq + dm ) / dm);%n�ԖڂɎ��o�����g��������̂������߂Ă���,n�͐����Ȃ̂�fix�Ŋۂߍ�
    disp(freq)
    if(freq == 5000)
        tmp = 0;
        SP_ant = zeros(1,11);
        SP_ant(1,tmp+1) = SP(n);
        disp(SP(n))
        tmp = tmp+1;
    else
        SP_ant(1,tmp+1) = SP(n);
        disp(SP(n))
        tmp = tmp+1;
    end
end
%%
POW_mic = POW./ SP_ant;
pow_ori = 20.*log10(POW);
pow_mic = 20.*log10(POW_mic/mean(POW_mic));
sp = 20.*log10(SP_ant);
len_pow = length(pow_ori);
dlt_pow = 0:len_pow-1;
figure;
%plot(dlt_pow,pow_ori,'Color','r');
grid on;
hold on;
plot(dlt_pow,pow_mic,'Color','r');
%plot(dlt_pow,sp,'Color','green');
xticks([1 6 11 16]);
xticklabels({'5k';'10k';'15k';'20k';});
axis([-inf inf -60 20]);
set(gca, 'XScale', 'Log');
title('Random Frequency Response');
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
%%
%lgd = legend({'ori','SP','Mic'},'FontSize','southeast',12);
lgd = legend({'Mic'},'FontSize',12,'Location','southeast');
%title(lgd);
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

