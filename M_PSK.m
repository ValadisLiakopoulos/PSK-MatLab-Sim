format long;
M = 16; %M-PSK , change between 8 and 16
Fc = 1/4; %Carrier's Frequency. Reduced to 1/4
Tsymbol = 40; %Symbol's Period
square_pulse = sqrt(2/Tsymbol); %Square pulse
all_errors =[]; %For BER plots
all_theor_errors=[];
symbol_length = 34000; %Length of message to be sent




%create mapper matrix
if (M==16)
    mapping_matrix = ['0000';'0001';'0011';'0010';'0110';'0111';'0101';'0100';'1100';'1101';'1111';'1110';'1010';'1011';'1001';'1000'];
elseif (M==8)
    mapping_matrix = ['000';'001';'011';'010';'110';'111';'101';'100'];
else
    disp('Works only on M=8,16')
end

for SNR=-4:2:20
    fourier_signal=zeros(1,symbol_length);
    sended_symbols = []; 
    received_symbols =[];
    disp('Loop');
    disp(SNR);
    Es=1;
    No = Es/(10^(SNR/10)*log2(M)); % Noise deviation based on the SNR
    map = [];
    for i=1:M
        map = [map;sqrt(Es)*cos((2*pi*(i-1))/M),sqrt(Es)*sin((2*pi*(i-1)/M))];    
    end   
    for i=1:symbol_length
        %disp(i);
        m = mapper(mapping_matrix,M);
        sended_symbols = [sended_symbols,m];
        rt = [];
        for t=1:Tsymbol
            Sm = psk_carrier(m,M,Fc,t,square_pulse);
            %send message through AWGN channel
            channel = awgn(Sm,No);
            rt = [rt;channel];
             %= deconfigurator(rt,Fc,t,square_pulse);
        end
        fourier_signal(i) = Sm;
        r = deconfigurator(rt,Fc,Tsymbol,square_pulse);
        exp_symbol = wearer(r,map);
        received_symbols = [received_symbols, exp_symbol];
    end
    sended_message = '';
    received_message = '';
    
    
    for i=1:symbol_length
        rec = received_symbols(i)+1;
        sen = sended_symbols(i)+1;
        sended_message = [sended_message,mapping_matrix(sen,:)];
        received_message = [received_message,mapping_matrix(rec,:)];
    end
    errors = 0;

    %calculate all wrong bits(error)
    for i=1:length(sended_message)
        if(sended_message(i)~=received_message(i))
            errors = errors +1;
        end
    end
    
    %theoretical error
    expected_error = 2*qfunc(sqrt(2/No)*sin(pi/M));
    all_theor_errors = [all_theor_errors;expected_error];
    error_prsnt = (errors/(symbol_length*log2(M)))*100;
    all_errors = [all_errors;error_prsnt];
    disp('Loop Completed');

    %{
    disp('Message sent');
    disp(sended_message);
    disp(sended_symbols);
    disp('Message received');
    disp(received_message);
    disp(received_symbols);
    disp('Error Percentage');
    disp(error_prsnt);
    disp('Correct bits');
    disp(length(sended_message) - errors);
    disp(unique(received_symbols));
    disp('Expected theoretical error');
    disp(expected_error);
    %}
end


fft_sample_max_len = fix(length(fourier_signal)/2048);
power_spectrum2=zeros(2048,1);
for i=1:2048
    start_index = i*fft_sample_max_len;
    end_index = start_index+15;
    if(end_index>length(fourier_signal))
        end_index = length(fourier_signal);
    end
    fft_vector = fourier_signal(start_index:end_index);
    fft_magnitude = fft(fft_vector);
    fft_magnitude = abs(fft_magnitude).^2;
    power_spectrum2(i) = sum(fft_magnitude)/length(fft_magnitude);
end
%plot(10*log10(power_spectrum));


%{
for i=1:fft_frequency:symbol_length
    power_spectrum=[];
    end_index = i+15;
    if(end_index>symbol_length)
        end_index = symbol_length;
    end
    fft_vector = fourier_signal(i:end_index);
    fft_magnitude = fft(fft_vector);
    fft_magnitude = fft_magnitude.^2;

    for j=1:length(fft_magnitude)
        power = sum(fft_magnitude)/len(fft_magnitude);
    end
end
%}


indices = (-4:2:20);
figure;
hold on;
plot(indices, all_errors,'-o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
plot(indices,all_theor_errors,'-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
ylabel('Bit Error Rate in percentage');
xlabel('Signal/Noise Ratio');
ylim([0 70]);
title('BER Plot');
legend('Measured BER 8-PSK','Theoretical BER 8-PSK')



%mapper's function
function m = mapper(matrix,M)
    source_output = 0;
    str_output = '';
    % Generate a random binary sequence with probability 0.5 for each bit
    for i=1:log2(M)
        %simulate binary numbers with decimals
        source_output = (rand(1, 1) > 0.5);
        str_output = [str_output,num2str(source_output)];
    end
    index = 0;
    for i=1:M
        if strcmpi(str_output,matrix(i,:))
            index=i;
        end
    end
    m = index-1;
end


%carrier and square pulse multiplication stage
function Sm_t = psk_carrier(m, M, Fc, t,pulse)
    % Carrier multiplication + pulse (pulse is steady for every different t)
    Sm_t_cos = cos((2*pi*m)/M) * cos(2*pi*Fc*t) * pulse;
    Sm_t_sin = sin((2*pi*m)/M) * sin(2*pi*Fc*t) * pulse;
    Sm_t = Sm_t_sin + Sm_t_cos;
end

%AWGN Channel Function
function rt = awgn(Sm,No)
    noise = sqrt(No/2) * randn; %1 bit at a time
    rt = Sm + noise;
end

%Deconfigurator's Function
function r = deconfigurator(rt,Fc,Ts,pulse)
    for t = 1:Ts
        rt_cos = rt(t)*pulse*cos(2*pi*Fc*t); %multiply with carrier's frequency

        rt_sin = rt(t)*pulse*sin(2*pi*Fc*t);
    end
    r1 = sum(rt_cos); %get the sum from all the samples for all symbol's period
    r2 = sum(rt_sin);
    r = [r1;r2]; 
end

function symbol = wearer(r,map)
    diff_vec = [];
    for i=1:length(map)
        %diff_cos = abs((r(1) - map(i,1)))^2;
        %diff_sin = abs((r(2) - map(i,2)))^2;
        %diff_sum = sqrt(diff_cos+diff_sin);
        diff_sum = norm(r - map(i,:));
        diff_vec = [diff_vec;diff_sum];
    end
        min_diff = min(diff_vec);
        symbol = find(diff_vec == min_diff);
        symbol = symbol - 1;
end
