%TSP���g�p����IR����Ǝ��g�������̎Z�o
%adobe audition���ɂĘ^�����s���Ă�������
% 20171002 sakaguchi �C��
%% ---  �����ݒ�  ------------------------
clear;
close all;
% ----- �p�����[�^�ݒ� -----
fs = 48000;% �T���v�����O���g�� windows�ł��ݒ�ς��Ȃ��ƃo�O�肻��(�B��)
nai=fs/2;
%ori_len = 2^18*2;                   % �M���� (1����) 2^18*2����
ori_len = 2^16*2; %�����葁��
ori_half = ori_len/2;               % �M���̎������iori_len/2���������߁j
amp = 0.5;                          % �U��
n_time = 3;                         % ����񐔁i�K�v�ɉ����đ������������j
rec_num=ori_len*(n_time+1);         % �^���� 

frame=2^16;                         % fft frame length 1.3652 ms    
fr_half=frame/2;
vin=fs/frame;                       % ���g���𑜓x

tim_s=0.005*fs;                     % �G�����F���[�v�̊Ԉ���

%% ----- �M���쐬�iTSP�j -----
%  TSP �̎��g������
ori_f=zeros(ori_len,1);
kk = 0 : ori_half;           % kk�F���U���g���ԍ�

ori_f(kk+1) = exp(-1i * 2 *pi * ori_half * (kk / ori_len).^2 );
ori_f( ori_len/2+2: ori_len) = conj( ori_f( ori_len/2: -1: 2));   % frame/2+2 �ȏ�͋����̐܂�Ԃ�
ori_sub = real(ifft(ori_f));                         % ���Ԕg�`
    
    
ori_sub = [ori_sub(ori_half+ori_half/2+1: ori_len); ori_sub(1:ori_half+ori_half/2 )];  % �g�`��M�����ori_len�̒����ɂ���~��V�t�g
ori = ori_sub / (max(abs(ori_sub)) /amp);               % �U���ő��AA�ɂ���
out = ori;
for k=1:n_time
    out=[out;ori];                            % �Đ��M������
end

%% IR�̏����o��
fname =['TSP_' num2str(fs) 'Hz_t16_x3.wav'];
audiowrite(fname,out,fs);

%% �^��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% ����������Adobe Audition�Ƃ��g���Ę^�� %%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% �Ñ����̎�荞��
%��������
[Path, FN] = uigetfile('*.wav', '�Ñ�����I�����Ă�������');   %�Ñ����̑I��
BG = strcat(FN, Path);
bgn = audioread(BG);

bgn1 = bgn(1:ori_len,1);
w_hann = hann(ori_len);
bgn1 = bgn1.*w_hann;
%% �X�s�[�J�C���p���X�����̎�荞��
[ file, path] = uigetfile('*.wav', '�X�s�[�J�̃C���p���X������I�����Ă�������');%�C���p���X������ǂݍ���
SP = strcat( path, file);
[ sp_imp, fs_sp] = audioread(SP);
[ row_sp, clm_sp] = size(sp_imp);

%sp = conv(sp_imp,ori);
%sp = sp(1:ori_len,1);
%sp = sp.*w_hann;
%sp = sp ./ max( abs(sp));
len_sp = length(sp_imp);
%t = (0:len_sp-1)./fs_sp;
%figure;plot(t,sp)
%ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1 1];
%% �^���f�[�^�̓ǂݍ���
%������10�������o��(��������)
n = 36;

ang = 0:10:n*10;

[IR, I] = uigetfile('*.wav', '�^���f�[�^��I�����Ă�������');   %�^���f�[�^�̑I��
Y = strcat(I, IR);
y = audioread(Y);
rec = y;

k=1;
len = length(rec(ori_len*k+1:ori_len*(k+1))); t=(0:len-1)'./fs;
figure;plot(t,rec(ori_len*k+1:ori_len*(k+1),1))
ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1 1];
%% ----- �C���p���X���� �Z�o -----
sub_1=zeros(ori_len,1);
    
for k=1:n_time
    sub_1 = sub_1+rec(ori_len*k+1:ori_len*(k+1),:);        % �Q�����ڂ�؂�o���ĂR�����ȍ~�͑����Ă����i�{�]�u�j
end
rec2=sub_1./n_time;%�������񐔂Ŋ����ĕ��ω�
figure;plot(t,rec2);
axis([-inf inf -1 1]);
%%         
IMP = fft(rec2)./[ori_f,ori_f] ;          % �^���M����DFT���ċt�t�B���^(1/SS)����Z����
IMP = IMP *(max(abs(ori_sub)) /amp);      % �U�������̋t���Z
imp = real(ifft(IMP));                    % �tDFT�����ăC���p���X�����𓾂�(�g���͎̂��������̂�)

% ���M���̃C���p���X�����̎Z�o
%IMPP = fft(ori)./ori_f;
%IMPP = IMPP *(max(abs(ori_sub)) /amp);      % �U�������̋t���Z
%impp = real(ifft(IMPP));                    %�����ł��g���͎̂��������̂�

%% ----- ���ʂ̕\�� -----
tf= 1; te = t(end); %�K�v�ɉ�����x���͈̔͐ݒ肵�Ă�������

figure
len = length(impp); t=(0:len-1)'./fs;
subplot(3,1,1); plot(t,impp); grid on
ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1.1 1.1];
title('OriginalImpluse'); xlabel('time [s]');

len = length(imp); t=(0:len-1)'./fs;
subplot(3,1,2); plot(t,imp(:,1)); grid on
ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1.1 1.1];
title('Lch'); xlabel('time [s]');
subplot(3,1,3); plot(t,imp(:,2)); grid on
ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1.1 1.1];
title('Rch'); xlabel('time [s]');

%% ----- rec & imp �ۑ� -----
% rec
%file_name = './ir_rec1.wav';
%audiowrite(file_name,rec,fs); % wav file�@�ŕۑ�
% imp
file_name = './imp/S_0_imp_mono.wav';
audiowrite(file_name,imp(:,1),fs); % wav file�@��IR��ۑ�
      
%% ----- �����������g�����������߂� -----
%������10�����Ƃ肾������
clear

[ file, path] = uigetfile('*.wav', 'Select the .wav file');%�C���p���X������ǂݍ���
loc = strcat( path, file);
[ data, fs] = audioread( loc);
[ row, clm] = size( data);
len_data = length( data);
data = 0.9 .* data ./ max( abs(data));
%% �X�s�[�J�C���p���X�����̎�荞��
[ file, path] = uigetfile('*.wav', '�X�s�[�J�̃C���p���X������I�����Ă�������');%�C���p���X������ǂݍ���
SP = strcat( path, file);
[ sp_imp, fs_sp] = audioread(SP);
[ row_sp, clm_sp] = size(sp_imp);

%sp = conv(sp_imp,ori);
%sp = sp(1:ori_len,1);
%sp = sp.*w_hann;
%sp = sp ./ max( abs(sp));
len_sp = length(sp_imp);
%t = (0:len_sp-1)./fs_sp;
%figure;plot(t,sp)
%ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1 1];
%%
w_data = hann(len_data);
data = data.*w_data;
    
w_sp = hann(len_sp);
sp_imp = sp_imp.*w_sp;
%% -----
%ms = 1000.*( 0:len_data-1)./fs;
%norm_data = 0.9 .* data ./ max( abs(data(:)));


if len_data < 2*fs
    if clm == 1
        emp = zeros( 2*fs-len_data, 1);
    elseif clm == 2
        emp = zeros( 2*fs-len_data, 2);
    end
    
    data_comp = [ data; emp];
else
    data_comp = data;
end
len_comp = length( data_comp);
deltaF = ( 0:len_comp-1 ).*fs./ len_comp;

if len_sp < len_comp
    if clm_sp == 1
        emp_sp = zeros( len_comp-len_sp, 1);
    elseif clm_sp == 2
        emp_sp = zeros( len_comp-len_sp, 2);
    end
    
    sp = [ sp_imp; emp_sp];
else
    sp = sp_imp;
end
%50Hz��HPF��������
%n = 4096;
%a = 1;
%cutoff = 2 * 50 / fs;
%X = fir1(n,cutoff,'high');
%data_comp = filter(X,a,data_comp);
%HPF�I���


%% 
%%data_comp = data_comp.*w_comp;
%%data_comp_sp = data_comp_sp.*w_comp;

%data_comp = 0.9 .* data_comp ./ max( abs(data_comp));
%data_comp_sp = 0.9 .* data_comp_sp ./ max( abs(data_comp_sp));
%% 
FFT_data = fft( data_comp(:,1));
%FFT_sp = fft(sp);
FR = 20.*log10( abs( FFT_data));
%FR_sp = 20.*log10(abs(FFT_sp));
%PFC = unwrap( angle( FFT_data));

xtick = [20 100 1000 10000 20000 fs/2];
x_axis = {'20';'100';'1k';'10k';'20k';' '};
figure;
plot( deltaF, FR, 'Color', 'b'); grid on
hold on
%plot( deltaF, FR_sp, 'Color', 'b');
axis([20, fs/2 -60 20]);
set(gca, 'XScale', 'Log', 'XTick', xtick, 'XTickLabel', x_axis);
title('Frequency Response');
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
%%
FR_mic = FR - FR_sp;
figure;
plot( deltaF, FR_mic, 'Color', 'r');grid on
axis([20, fs/2 -inf inf]);

set(gca, 'XScale', 'Log', 'XTick', xtick, 'XTickLabel', x_axis);
title('Frequency Response');
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
