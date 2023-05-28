function env = envGen(fadeinTime,dur,fadeoutTime,fs,nCh,envType)
    nfadein = fadeinTime*fs;
    nfadeout = fadeoutTime*fs;
    nsteady = round((dur-(fadeinTime+fadeoutTime))*fs);

    % Linear fade
    if envType == "linear"
        envFadeIn = linspace(0,1,nfadein);
        envSteady = ones(1,nsteady);
        envFadeOut = linspace(1,0,nfadeout);
    end

    % Sine-squared fade
    if envType == "sinsq"
        alphain = (pi/2)*linspace(0,1,nfadein);
        alphaout = (pi/2)*linspace(0,1,nfadeout);
        envFadeIn = sin(alphain).^2;
        envFadeOut = fliplr(sin(alphaout).^2);
        envSteady = ones(1,nsteady);
    end
    
    env = [envFadeIn envSteady envFadeOut];
    env = repmat(env',1,nCh)';

end