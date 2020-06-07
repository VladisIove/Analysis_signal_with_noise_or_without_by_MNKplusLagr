
 function koef = find_koef(y,Y_noise_vector, signal_to_noise, delta)
    koef = 0.00000001;step=0.00001;
    while 1
        sign_to_noise_ratio = snr(y, koef*Y_noise_vector)
        if sign_to_noise_ratio>=signal_to_noise-delta && sign_to_noise_ratio<=signal_to_noise+delta
            break
        end

            
        if sign_to_noise_ratio<signal_to_noise
            koef = koef+step;
        else
            koef = koef-step;
        end
    end
end
