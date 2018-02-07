%% 特定の周波数範囲の拡散入射特性(推定)を表示するプログラム
%インパルス応答を読み込んで指定の周波数に対して1周分のチャートを描画するプログラム
%軸非対称専用
%仰角ごとにフォルダ分け．同じ階層に入れておく．仰角は正面を0°とした時にコネクタ側が突き出る方が負の角度，反対に突き出る方が正の角度
%インパルス応答のファイルは　[種類がわかるような文字]_[方位角(半角)]_[仰角(半角)]_imp.wav　で統一　一箇所のフォルダに格納しておく
%実行はインパルス応答を保存してある階層を選択
%20180116 Sakaguchi作成

%%
clear;
close all;
clc;
%% --音のファイルを選択する処理
flg = 0;
for freq = 5000:1000:15000
    disp(freq)
    m = 360;
    num = m-10;
    for deg = -90:10:90 %仰角の変数
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
        for v = 0:mad:num%方位角の変数　0から350まで回す(10°ごとのデータを取得)
            s = num2str(v);%vを文字配列型に変換(fnameで使うため)
            fname =['U87_' s '_' d '_imp.wav'];
            PathName = ['./' d '/'];
            T = strcat(PathName, fname);
            [data,Fs] = audioread(T);%データの読み込み
            disp(fname)%今どの角度を実行しているのかをコマンドに表示
            data = data(:,1);%モノラルデータがほしいのでLだけデータを取り出す

            %%FFT

            len_y = length(data);
            %長さの調整
            if (len_y < 2*Fs)
                add1 = 2*Fs - len_y;
                z_add1 = zeros(add1,1);
                yn = vertcat(data,z_add1);
            else
                yn = data;

            end

            %指定の周波数を取り出す
            len_yn = length(yn);
            deltaF = ( 0:len_yn-1).*Fs./ len_yn; %横軸
            dm = deltaF(3)-deltaF(2);%deltaFの差を求める(deltaFは等差数列)
            n = fix(( freq + dm ) / dm);%n番目に取り出す周波数があるのかを求めている,nは整数なのでfixで丸め込
            
            Y = abs(fft(yn));
            p = Y(n);%複素数で値が出るのでその絶対値を求める
             %グラフ描画のための準備

            %周波数領域の準備
%             len_Y = length(Y);%周波数分解能の計算
%             len_Y = len_Y - rem(len_Y,2);
%             time_Y = (1:len_Y)/Fs;
%             Y = Y(1:len_Y);
%             f_vin = Fs/len_Y;
%             fre = 0:f_vin:Fs/2; %周波数分解能の計算ここまで

             %球表面積の公式
            r = 2.0;
            S = 2*pi*r*r*sin(2*deg*pi/180)*cos(deg/2*pi/180)/2.5847;
            S_all = 4*pi*r*r;
            %deltaF = deltaF(n);

            %面積による重みづけ
            if(deg <= 50)
                p_sq = p*p*((S/36)/S_all);
            elseif(deg <= 70)
                p_sq = p*p*((S/12)/S_all);
            elseif(deg == 80)
                p_sq = p*p*((S/4)/S_all);
            elseif(deg == 90)
                p_sq= p*p*(S/S_all);
            end
            %同じ周波数のデータを加算していく
            if(flg == 0)
                p_sum = 0;%一応初期化
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

%% スピーカの特性の引き算
[ file_sp, path] = uigetfile('*.wav', 'スピーカのインパルス応答を選択してください');%インパルス応答を読み込む
SP = strcat( path, file_sp);
[ sp_imp, fs_sp] = audioread(SP);
len_sp_imp = length(sp_imp);
%長さの調整
if (len_sp_imp < 2*fs_sp)
    add1 = 2*fs_sp - len_sp_imp;
    z_add1 = zeros(add1,1);
    sp = vertcat(sp_imp,z_add1);
else
    sp = sp_imp;

end
len_sp = length(sp);
deltaF = ( 0:len_sp-1).*fs_sp./ len_sp; %横軸
dm = deltaF(3)-deltaF(2);%deltaFの差を求める(deltaFは等差数列)

SP = abs(fft(sp));
for freq = 5000:1000:15000
    n = fix(( freq + dm ) / dm);%n番目に取り出す周波数があるのかを求めている,nは整数なのでfixで丸め込
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
%% 初めてチャート表示するときに実行
flg =1;%フラグ指定
%% チャート表示
%図に表示する
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
    ax.ThetaTickLabel = {'0°','30°','60°','90°','120°','150°','180°','210°','240°','270°','300°','330°'};
    ax.RTickLabel = {'-20dB','','-10dB','','0dB'};
    ax.ThetaZeroLocation = 'top';
    ax.Color = 'w';
    ax.ThetaDir = 'clockwise';
    %ax.RDir = 'reverse';
    
    %印刷時に見えやすいようにする設定
    ax.GridAlpha = 0.7;
    ax.FontSize = 13;
    
    %表示範囲の設定
    ax.RLim = [-25 0];
    ax.ALim = [0 3];
    flg = 2;%フラグを2に変更
else
    avr_av1 = avr_av;
    polarplot(avr_av1,'-.','LineWidth',1.0);
    %f_st =num2str(f);
    %fre =[f_st 'Hz'];
    %legend({fre})
    hold on
end

