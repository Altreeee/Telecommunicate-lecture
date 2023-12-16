function BER=analytical_ber(snrdB)

snr=10.^(snrdB/10);
L=length(snr);

for m=1:L
    BER(m)=erfc(sqrt(snr(m)))/2;
end

end