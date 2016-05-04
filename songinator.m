    
    srate = 44100;
    frequency = 1;
    files = dir('./songs/*.mp3');
    sLength = 0;
    longest = 0;
    
    %find the length of the longest song
    for k = 1:length(files)

        if length(audioread(['./songs/', files(k).name])) > sLength
            sLength = length(audioread(['./songs/', files(k).name]));
            longest = k;
        end
        
    end
    
    %stretch the songs to equal length
    for i = 1:length(files)

        x = audioread(['./songs/', files(i).name]);
        
        rsFactor = round(srate*length(audioread(['./songs/', files(longest).name]))/length(x));
        
        x = resample(x, round(rsFactor/10), round(srate/10));
        
        l = length(x);
        
        f = zeros(l, 2);
        
        x = x * .5;
        
        audiowrite(strcat(files(i).name, '.wav'), x, srate);
        
    end
    
    %create a container for the output files
    %output = zeros(length(files), length(audioread(['./songs/', files(longest).name])), 2);
    
    %create the amplitude envelopes
    envelope1 = zeros(l, 2);
    
    %fill the left amplitude envelope with the highest average RMS song
    maxRMS = rms(envelope1, 1);
    newFiles = dir;
    loudest = 0;
    
    for i = 1:length(newFiles)
        
        if strfind(newFiles(i).name, '.wav')
            z = audioread(newFiles(i).name);
            if rms(z, 1) > maxRMS
                maxRMS = rms(z, 1);
                loudest = i;
            end
        end
        
    end
    
    q = audioread(newFiles(loudest).name);
    
    if length(q) < length(envelope1)
        q = vertcat(q, zeros(abs(length(q)-length(envelope1)), 2));
    elseif length(envelope1) < length(q)
        envelope1 = vertcat(envelope1, zeros(abs(length(q)-length(envelope1)), 2));
    end
    
    envelope1(:, 1) = q(:, 1);
    
    %fill the right amplitude envelope with the polarity-inverted left side
    envelope1(:, 2) = envelope1(:, 1) * -1;
    
    %create alternate envelope
    envelope2 = envelope1 * -1;
    
    %now we do some fft shit
    
    %find the song with the most low frequency energy
    low = zeros(round((length(envelope1)*(120/srate)+1)-(length(envelope1)*(40/srate)+1)), 2);
    lowest = 0;
    energy = 0;
    
    for i = 1:length(newFiles)
        
        if strfind(newFiles(i).name, '.wav')
            w = audioread(newFiles(i).name);
            wFFT = fft(w((length(envelope1)*(40/srate)+1):(length(envelope1)*(120/srate)+1), 2));
            
            for j = 1:length(wFFT)
                energy = energy + wFFT(j);
            end
            
            energy = energy / length(wFFT);
        
            if energy > low
                low = energy;
                lowest = i;
            end 
        end
        
    end
    
    %filter the lowest song above 500Hz
    pF = audioread(newFiles(lowest).name);
    [B, A] = butter(4, 500/22050, 'low');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
    audiowrite((newFiles(lowest).name), pF, srate);
    
    %filter everything else below 120 Hz
    [B, A] = butter(2, 120/22050, 'high');
    
    for i = 1:length(newFiles)
        
        if strfind(newFiles(i).name, '.wav')
            if i~=lowest
                pF = audioread(newFiles(i).name);
                pF(:, 1) = filter(B, A, pF(:, 1));
                pF(:, 2) = filter(B, A, pF(:, 2));
                audiowrite((newFiles(i).name), pF, srate);
            end
        end
        
    end
    
    %find the song with the most low-mid frequency energy
    lowmid = zeros(round((length(envelope1)*(750/srate)+1)-(length(envelope1)*(120/srate)+1)), 2);
    lowmidest = 0;
    energy = 0;
    
    for i = 1:length(newFiles)
        if i~=lowest
            if strfind(newFiles(i).name, '.wav')
                w = audioread(newFiles(i).name);
                wFFT = fft(w((length(envelope1)*(120/srate)+1):(length(envelope1)*(750/srate)+1), 2));
            
                for j = 1:length(wFFT)
                    energy = energy + wFFT(j);
                end
            
                energy = energy / length(wFFT);
        
                if energy > lowmid
                    lowmid = energy;
                    lowmidest = i;
                end 
            end
        end
        
    end
    
    %filter the low-midest song above 1kHz and above 250 Hz
    pF = audioread(newFiles(lowmidest).name);
    
    [B, A] = butter(4, 1000/22050, 'low');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
    
    [B, A] = butter(4, 250/22050, 'high');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
    
    audiowrite((newFiles(lowmidest).name), pF, srate);
    
    
    %find the song with the most high-mid frequency energy
    highmid = zeros(round((length(envelope1)*(3000/srate)+1)-(length(envelope1)*(750/srate)+1)), 2);
    highmidest = 0;
    energy = 0;
    
    for i = 1:length(newFiles)
        if i~=lowest
            if i~=lowmidest
                if strfind(newFiles(i).name, '.wav')
                    w = audioread(newFiles(i).name);
                    wFFT = fft(w((length(envelope1)*(750/srate)+1):(length(envelope1)*(3000/srate)+1), 2));
            
                    for j = 1:length(wFFT)
                        energy = energy + wFFT(j);
                    end
            
                    energy = energy / length(wFFT);
        
                    if energy > highmid
                        highmid = energy;
                        highmidest = i;
                    end 
                end
            end
        end
        
    end
    
    %filter the high-midest song below 750 Hz and above 8 kHz
    pF = audioread(newFiles(highmidest).name);
    
    [B, A] = butter(2, 750/22050, 'high');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
    
    [B, A] = butter(2, 8000/22050, 'low');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
        
    audiowrite((newFiles(highmidest).name), pF, srate);
    
    %find the song with the most high frequency energy
    high = zeros(round((length(envelope1)*(16000/srate)+1)-(length(envelope1)*(3000/srate)+1)), 2);
    highest = 0;
    energy = 0;
    
    for i = 1:length(newFiles)
        if i~=lowest
            if i~=lowmidest
                if i~=highmidest
                    if strfind(newFiles(i).name, '.wav')
                        w = audioread(newFiles(i).name);
                        wFFT = fft(w((length(envelope1)*(3000/srate)+1):(length(envelope1)*(16000/srate)+1), 2));
            
                        for j = 1:length(wFFT)
                            energy = energy + wFFT(j);
                        end
            
                        energy = energy / length(wFFT);
        
                        if energy > high
                            high = energy;
                            highest = i;
                        end 
                    end
                end
            end
        end
        
    end
    
    %filter the highest song below 2kHz
    pF = audioread(newFiles(highest).name);
    [B, A] = butter(2, 2000/22050, 'high');
    pF(:, 1) = filter(B, A, pF(:, 1));
    pF(:, 2) = filter(B, A, pF(:, 2));
    audiowrite((newFiles(highest).name), pF, srate);
    
    %apply an envelope to each song and attenuate
    for i = 1:length(newFiles)
        
        if strfind(newFiles(i).name, '.wav')
            
            temp = audioread(newFiles(i).name);
            temp = temp * 0.3;
            
            if length(temp)<length(envelope1)
                temp = vertcat(temp, zeros(abs(length(temp)-length(envelope1)), 2));
            elseif length(temp)>length(envelope1)
                envelope1 = vertcat(envelope1, zeros(abs(length(temp)-length(envelope1)), 2));
                envelope2 = vertcat(envelope2, zeros(abs(length(temp)-length(envelope2)), 2));
            end
            
            if round(rand(1)) ~= 0              
                temp(:,1) = temp(:,1) .* envelope1(:,1);
                temp(:,2) = temp(:,2) .* envelope1(:,2);
            else
                temp(:,1) = temp(:,1) .* envelope2(:,1);
                temp(:,2) = temp(:,2) .* envelope2(:,2);
            end
           
        end
        
    end
    
    %do the mix down
    mix = zeros(length(envelope1), 2);
    
    for i = 1:length(newFiles)
        
        if strfind(newFiles(i).name, '.wav')
            
            p = audioread(newFiles(i).name);
            
            %chance to reverse
            if round(rand(1))~=0
                if round(rand(1))~=0
                    p(:, 1) = flipud(p(:, 1));
                    for n = 1:length(p)
                        p(n, 1) = p(n, 1) * (n/length(p));
                    end
                end
            end
            
            if round(rand(1))~=0
                if round(rand(1))~=0
                    p(:, 2) = flipud(p(:, 2));
                    for n = 1:length(p)
                        p(n, 2) = p(n, 2) * (n/length(p));
                    end
                end
            end
            
            %chance to crescendo
            if round(rand(1))~=0
                    for m = 1:length(p)
                        p(m, 1) = p(m, 1) * (m/length(p));
                    end
            end

            if round(rand(1))~=0
                    for m = 1:length(p)
                        p(m, 2) = p(m, 2) * (m/length(p));
                    end
            end
            
            %chance to decrescendo
            if round(rand(1))~=0
                if round(rand(1))~=0
                    for m = length(p):-1:1
                        p(m, 1) = p(m, 1) * (m/length(p));
                    end
                end
            end

            if round(rand(1))~=0
                if round(rand(1))~=0
                    for m = length(p):-1:1
                        p(m, 2) = p(m, 2) * (m/length(p));
                    end
                end
            end
                      
            %add the song
            if length(mix) > length(p)
                p = vertcat(p, zeros(abs(length(mix)-length(p)), 2));
                mix = mix + p;
            elseif length(mix) < length(p)
                mix = vertcat(mix, zeros(abs(length(mix)-length(p)), 2));
                mix = mix + p;
            else
                mix = mix + p;
            end
            
        end
        
    end
    
    audiowrite('mix.wav', mix, srate);
    
 
    
    