%TSPを使用したIR測定と周波数特性の算出
%adobe audition等にて録音を行ってください
% 20171002 sakaguchi 修正
%% ---  初期設定  ------------------------
clear;
close all;
% ----- パラメータ設定 -----
fs = 48000;% サンプリング周波数 windowsでも設定変えないとバグりそう(曖昧)
nai=fs/2;
%ori_len = 2^18*2;                   % 信号長 (1周期) 2^18*2ｵｽｽﾒ
ori_len = 2^16*2; %手っ取り早く
ori_half = ori_len/2;               % 信号の実効長（ori_len/2がおすすめ）
amp = 0.5;                          % 振幅
n_time = 3;                         % 測定回数（必要に応じて増減ください）
rec_num=ori_len*(n_time+1);         % 録音長 

frame=2^16;                         % fft frame length 1.3652 ms    
fr_half=frame/2;
vin=fs/frame;                       % 周波数解像度

tim_s=0.005*fs;                     % エンヴェロープの間引き

%% ----- 信号作成（TSP） -----
%  TSP の周波数成分
ori_f=zeros(ori_len,1);
kk = 0 : ori_half;           % kk：離散周波数番号

ori_f(kk+1) = exp(-1i * 2 *pi * ori_half * (kk / ori_len).^2 );
ori_f( ori_len/2+2: ori_len) = conj( ori_f( ori_len/2: -1: 2));   % frame/2+2 以上は共役の折り返し
ori_sub = real(ifft(ori_f));                         % 時間波形
    
    
ori_sub = [ori_sub(ori_half+ori_half/2+1: ori_len); ori_sub(1:ori_half+ori_half/2 )];  % 波形を信号区間ori_lenの中央にする円状シフト
ori = ori_sub / (max(abs(ori_sub)) /amp);               % 振幅最大をAAにする
out = ori;
for k=1:n_time
    out=[out;ori];                            % 再生信号合成
end

%% IRの書き出し
fname =['TSP_' num2str(fs) 'Hz_t16_x3.wav'];
audiowrite(fname,out,fs);

%% 録音

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%% つくった音はAdobe Auditionとか使って録音 %%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 暗騒音の取り込み
%書きかけ
[Path, FN] = uigetfile('*.wav', '暗騒音を選択してください');   %暗騒音の選択
BG = strcat(FN, Path);
bgn = audioread(BG);

bgn1 = bgn(1:ori_len,1);
w_hann = hann(ori_len);
bgn1 = bgn1.*w_hann;
%% スピーカインパルス応答の取り込み
[ file, path] = uigetfile('*.wav', 'スピーカのインパルス応答を選択してください');%インパルス応答を読み込む
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
%% 録音データの読み込み
%自動で10°ずつ取り出す(書きかけ)
n = 36;

ang = 0:10:n*10;

[IR, I] = uigetfile('*.wav', '録音データを選択してください');   %録音データの選択
Y = strcat(I, IR);
y = audioread(Y);
rec = y;

k=1;
len = length(rec(ori_len*k+1:ori_len*(k+1))); t=(0:len-1)'./fs;
figure;plot(t,rec(ori_len*k+1:ori_len*(k+1),1))
ax = gca; ax.XLim=[0 t(end)]; ax.YLim=[-1 1];
%% ----- インパルス応答 算出 -----
sub_1=zeros(ori_len,1);
    
for k=1:n_time
    sub_1 = sub_1+rec(ori_len*k+1:ori_len*(k+1),:);        % ２周期目を切り出して３周期以降は足していく（＋転置）
end
rec2=sub_1./n_time;%足した回数で割って平均化
figure;plot(t,rec2);
axis([-inf inf -1 1]);
%%         
IMP = fft(rec2)./[ori_f,ori_f] ;          % 録音信号をDFTして逆フィルタ(1/SS)を乗算する
IMP = IMP *(max(abs(ori_sub)) /amp);      % 振幅調整の逆演算
imp = real(ifft(IMP));                    % 逆DFTをしてインパルス応答を得る(使うのは実数部分のみ)

% 原信号のインパルス応答の算出
%IMPP = fft(ori)./ori_f;
%IMPP = IMPP *(max(abs(ori_sub)) /amp);      % 振幅調整の逆演算
%impp = real(ifft(IMPP));                    %ここでも使うのは実数部分のみ

%% ----- 結果の表示 -----
tf= 1; te = t(end); %必要に応じてx軸の範囲設定してください

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

%% ----- rec & imp 保存 -----
% rec
%file_name = './ir_rec1.wav';
%audiowrite(file_name,rec,fs); % wav file　で保存
% imp
file_name = './imp/S_0_imp_mono.wav';
audiowrite(file_name,imp(:,1),fs); % wav file　でIRを保存
      
%% ----- ここから先周波数特性を求める -----
%自動で10°ずつとりだしたい
clear

[ file, path] = uigetfile('*.wav', 'Select the .wav file');%インパルス応答を読み込む
loc = strcat( path, file);
[ data, fs] = audioread( loc);
[ row, clm] = size( data);
len_data = length( data);
data = 0.9 .* data ./ max( abs(data));
%% スピーカインパルス応答の取り込み
[ file, path] = uigetfile('*.wav', 'スピーカのインパルス応答を選択してください');%インパルス応答を読み込む
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
%50HzでHPFをかける
%n = 4096;
%a = 1;
%cutoff = 2 * 50 / fs;
%X = fir1(n,cutoff,'high');
%data_comp = filter(X,a,data_comp);
%HPF終わり


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
