%% EE676A: Simulation-Based Design of 5G New Radio (NR) Wireless Standard
% MATLAB Assignment-4: 5G Transceiver Implementation -II
% NAME: S. Srikanth Reddy; Roll No: 22104092

clear all;
clc;

TB_length = 20496; %Transport Block size = 20496

fprintf('\nA = D^24 + D^23 + D^18 + D^17 + D^14 + D^11 + D^10 + D^7 + D^6 + D^5 + D^4 + D^3 + D + 1\n');
fprintf('B = D^24 + D^23 + D^6 + D^5 + D + 1\n');
fprintf('C = D^24 + D^23 + D^21 + D^20 + D^17 + D^15 + D^13 + D^12 + D^8 + D^4 + D^2 + D + 1\n');

S = input("\nFor a CRC of length 24, select one of the above standard CRC polynomials...: \n","s"); % select one of the standard CRC polynomial 

poly_CRC24A = [24,23,18,17,14,11,10,7,6,5,4,3,1,0];
poly_CRC24B = [24,23,6,5,1,0];
poly_CRC24C = [24,23,21,20,17,15,13,12,8,4,2,1,0];

if strcmpi(S,'A')
    poly_CRC = poly_CRC24A;
elseif strcmpi(S,'B')
    poly_CRC = poly_CRC24B;
elseif strcmpi(S,'C')
    poly_CRC = poly_CRC24C;
end

random_input = randi([0,1],1,TB_length); % generating random input

TransmitTB = add_crc(random_input,poly_CRC); %function call to crc generation and appending


L = 24
B = length(TransmitTB) % B = 20496 + 24
C = ceil(B/(8448-L)) % Kcb = 8448 for base graph 1 
B_dash = B + (C*24) % since we already know TB size, we know that C>1
K_dash = (B_dash)/C %segmented code block size with CRC
K = K_dash + 176 % from base graph 1 we get Zc as 320 and 
                 % 22*320 - K_dash(=6864) = 176 filler bits
G = 100*162*2;
E = G/C;


cb = zeros(C,K_dash-L); % each row for each code blocks

s = 1;
il_cb = [];
for r = 1:C

    for k = 1 : K_dash-L % 6840 bits
        cb(r,k) = TransmitTB(s);
        s = s + 1;
    end

    temp = cb(r,1:end);
    tempcb = add_crc(temp,poly_CRC); %generated and CRC appended for CB-block
    tempcb = tempcb(K_dash-L+1:end); %capturing the CRC bits
    p_idx = 1;
    for k = K_dash-L+1:K_dash 
        cb(r,k)=tempcb(p_idx); %appending the CRC bits
        p_idx = p_idx + 1;
    end
    
    for k = K_dash+1 : K % appending 0's as null bits
        cb(r,k) = 0; 
    end

    nbg = 1; % Base graph 1
    nldpcdecits = 25; % Decode with maximum no of iteration
  
    
  
    ldpc_coded_bits = double(LDPCEncode(cb(r,1:end)',nbg)); %LDPC encoding
    N = length(ldpc_coded_bits);

    ratematched_bits = ldpc_coded_bits(1:E); % ratematching
    
    interleaved_bits = interleave(ratematched_bits,4,E/4); %interleaving - function defined @bottom

    il_cb = [il_cb interleaved_bits]; %CB concat
end

    cb_scr = [];
    %scrambling-start
    x1 = dec2bin(1,31) - '0';
    x2 = dec2bin(255,31) - '0';
    gs = [];    % gs - gold sequence
    for n = 1:1600+G
        gs = [gs xor(x1(31),x2(31))];
        tempx1 = xor(x1(28),x1(31));
        tempx2 = xor(xor(xor(x2(28),x2(29)),x2(30)),x2(31));
        x1 = circshift(x1,1);
        x2 = circshift(x2,1);
        x1(1) = tempx1;
        x2(1) = tempx2;
    end
    gs_dash = gs(1601:end); %first 1600 samples are rejected
    cb_scr = xor(il_cb,gs_dash);
    %scrambling-end

    
    cb_scr_NRZ  = 2*(cb_scr-0.5); % Non-Return to Zero form for QPSK modulation
    rs_cb_scr_NRZ = reshape(cb_scr_NRZ,2,G/2); % separating inphase and quadrature streams
    rx = [];
    noise_power = (10^-5);
    noise = sqrt(noise_power);
    for t = 1:G/2
        tx1  = rs_cb_scr_NRZ(1,t); %inphase data
        tx2 = rs_cb_scr_NRZ(2,t);  %quadrature data
        tx = tx1+j*tx2;
        rx_sig = tx + noise;   %adding noise
        rx1 = real(rx_sig);    %demod inphase 
        rx2 = imag(rx_sig);    %demod quadrature 
        rx = [rx rx1 rx2];     %demodulated output
    end
    
    
    llr0 =  abs(-1 + rx);   
    llr1 =  abs(1 + rx);    
    
    llr = log(llr0./llr1);      % ldpc decoder requires log(p(r/0)/p(r/1))
    demod_output = llr;

    %descrambling start
    descrambled_bits = [];
    ds_bit = 0;
    for q = 1:G
        if gs_dash(q) == 1
            ds_bit = -demod_output(q);
        else
            ds_bit = demod_output(q);
        end
        descrambled_bits = [descrambled_bits ds_bit];
    end
    %descrambling end

    rc_cb = zeros(C,N);
    neut_inf = zeros(1,N-E);     
    for r=1:C
        deinterleaved_bits = deinterleave(descrambled_bits((r-1)*E+1:r*E),4,E/4); %de-interleaving - function defined @bottom
        rc_cb(r,:) = [deinterleaved_bits neut_inf];       %CB segmentation and rate recovery    
        outputbits = double(LDPCDecode(rc_cb(r,:)',nbg,nldpcdecits));  %LDPC decoding
        errors = find(outputbits - cb(r,1:end)')
        output_no_filler = outputbits(1:K_dash); %removing the filler bits
        valid = validate_crc(output_no_filler,poly_CRC) % validating CRC of each code block
        rcv_cb = [];
        if valid
            x = output_no_filler(1:K_dash-L); % removing the CRC bits after validation
            rcv_cb = [rcv_cb x]; % concatinating the code blocks after CRC removal
        else
        %request for retransmit
        end
    end
    
valid = validate_crc(rcv_cb,poly_CRC); % crc validation of concatenated blocks

fprintf('TB-CRC validation returned...: %d\n',valid); 
% valid=1 => final CRC validated and transmit and receive code block match 
% valid=0 => final CRC not validated and transmit and receive code block do
% not match

% CRC addition and validation functions are as follows:
function TransmitTB = add_crc(random_input,poly_CRC)

TB_bits = random_input;

gen_poly = zeros(1,length(poly_CRC)+1);
for i=1:length(poly_CRC)
    gen_poly(1,poly_CRC(1,i)+1) = 1;
end
gen_poly = flip(gen_poly);

append_zeros = zeros(1,length(gen_poly)-1); %generation zeros equal to the highest degree of CRC polynomial
TB_bits = [random_input append_zeros]; %adding those zeros to the input bits

%mod2 division 
p = (length(TB_bits)-length(gen_poly));
for j = 1:p
    if(length(TB_bits)<length(gen_poly))
        break
    end
    T_r=xor(TB_bits(1,:),[gen_poly zeros(1,(length(TB_bits)-length(gen_poly)))]);
    TB_bits = T_r(find(T_r,1,'first'):length(T_r));
end

TransmitTB = [random_input zeros(1,poly_CRC(1,1)-length(TB_bits)) TB_bits]; %appending the remainder
end

function valid = validate_crc(RX_bits,poly_CRC)

gen_poly = zeros(1,length(poly_CRC)+1);
for i=1:length(poly_CRC)
    gen_poly(1,poly_CRC(1,i)+1)=1;    
end
gen_poly=flip(gen_poly);
%mod2 division for checksum
x=(length(RX_bits)-length(gen_poly));
for j=1:x
    if (length(RX_bits)<length(gen_poly)) 
        break 
    end
    R_x=xor(RX_bits(1,:),[gen_poly zeros(1,(length(RX_bits)-length(gen_poly)))]);
    RX_bits=R_x(find(R_x,1,'first'):length(R_x));
end

valid=isempty(RX_bits); % 1=>crc validated and 0=>crc not validated.
end

function interleaved_bits = interleave(x,M,Col)
rows = zeros(M,Col);
w = [];
for m = 1:M
    rows(m,:) = x((m-1)*Col+1:m*Col);
end
for j = 1:Col
    for i=1:M
        w = [w rows(i,j)];
    end
end
interleaved_bits = w;
end

function deinterleaved_bits = deinterleave(y,M,Col)
columns = zeros(M,Col);
z = [];
for col = 1:Col
    columns(:,col) = (y((col-1)*M+1:col*M))';
end
for i = 1:M
    z = [z columns(i,:)];
end
deinterleaved_bits = z;
end

